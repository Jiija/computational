using ModalIntervalArithmetic
using Plots
using DelimitedFiles
using Random

include("sub_diff.jl")

function solve_first_problem()
    A = [ModalInterval(0.9, 1.1) ModalInterval(0.9, 1.1) ModalInterval(0.9, 1.1)
         ModalInterval(0.9, 1.1) ModalInterval(0.9, 1.1) ModalInterval(0.0, 0.0)
         ModalInterval(0.9, 1.1) ModalInterval(0.0, 0.0) ModalInterval(0.0, 0.0)]
    b = [ModalInterval(2.8, 3.2), ModalInterval(2.2, 1.8), ModalInterval(0.9, 1.1)]

    x0 = init_point(A, b)
    println("x0    = ", x0)
    (x, count) = sub_diff(A, b, x0, 1.e-4)
    println("x     = ", x)
    println("count = ", count)

    n = 8
    iters = []
    for k = 1:n
        prec = 10.0^-k
        (x, count) = sub_diff(A, b, x0, prec)
        push!(iters, count)
    end
    plot(1:n, iters,
        xlabel = "-lg(ε)",
        ylabel = "N",
        xticks = [0:1:n;],
        yticks = [0:2:max(count);],
        dpi=300,
        legend = false)
    savefig("plot.png")
end

function get_good_system(A, 𝐛, n_attempts)
    counter = 0
    C = []
    𝐝 = []
    while true
        counter += 1
        unused = collect(1:128)
        columns = []
        𝐝 = []
        for k = 1:18
            val = rand(unused)
            key = findfirst(isequal(val), unused)
            splice!(unused, key)
            push!(columns, transpose(A[val, :]))
            push!(𝐝, 𝐛[val])
        end
        C = vcat(columns...)
        if abs(det(C)) > 1e-9
            open("C.txt", "w") do io
                show(io, "text/plain", C)
            end
            break
        end
        if counter > n_attempts
            break
        end
    end
    return ModalInterval.(C), 𝐝
end

function heuristical_solver(A, 𝐛, n_attempts)
    a = fill(-Inf, 18)
    b = fill(+Inf, 18)
    for k = 1:35
        𝐂, 𝐝 = get_good_system(A, 𝐛, n_attempts)
        x0 = init_point(𝐂, 𝐝)
        x, iter = sub_diff(𝐂, 𝐝, x0, 0.001)
        for j = 1:18
            a[j] = max(a[j], x[j].inf)
            b[j] = min(b[j], x[j].sup)
        end
    end
    return ModalInterval.(a, b)
end

function solution_checker(A, 𝐛, x)
    prod = A * x
    n_correct = 0
    for k = 1:256
        correct = true
        correct &= (𝐛[k].inf < prod[k].inf)
        correct &= (prod[k].sup < 𝐛[k].sup)
        if correct
            n_correct += 1
        end
    end
    return n_correct
end


function plot_bounds(𝐱, s)
    inf_values = []
    sup_values = []
    mid_values = []
    index = collect(1:36)
    for i in index
        push!(inf_values, 𝐱[i].inf)
        push!(sup_values, 𝐱[i].sup)
        push!(mid_values, (𝐱[i].inf + 𝐱[i].sup) / 2)
    end
    plot(index, inf_values, label = "inf", dpi = 300)
    plot!(index, sup_values, label = "sup")
    plot!(index, s, label = "solution")
    plot!(index, mid_values, label = "mid")
    savefig("bounds_plot.png")
end


function solve_second_problem()
    ε = 0.4
    A = readdlm("data\\A.txt")
    b = readdlm("data\\bnew.txt")
    s = readdlm("data\\snew.txt")
    𝐛 = ModalInterval.(b .- ε, b .+ ε)

    U = A[1:128, 19:36]
    D = A[129:256, 1:18]

    x1 = heuristical_solver(U, 𝐛[1:128], 10000)
    x2 = heuristical_solver(D, 𝐛[129:256], 10000)
    x = vcat(x1, x2)
    open("x.txt", "w") do io
        show(io, "text/plain", x)
    end;
    println("Solved: ", solution_checker(A, 𝐛, x), " out of: ", 256)
    plot_bounds(x, s)
end

solve_first_problem()
solve_second_problem()
