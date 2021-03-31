using IntervalArithmetic
using LinearAlgebra
using IterTools
using Printf

function gen(n::Int)
    return Base.Iterators.product([0:1 for i in 1:n]...)
end

function determine_degeneracy_problem_1(ε)
    A = [interval(1-ε, 1+ε) interval(1-ε, 1+ε);
        interval(1.1-ε, 1.1+ε) interval(1-ε, 1+ε)]

    det_A = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    return 0 ∈ det_A
end


function determine_degeneracy_Baumans(A, n)
    list = gen(n^2)
    A1   = zeros((n, n))
    A2   = zeros((n, n))
    res  = false
    counter = 0
    next1 = iterate(list)
    while next1 !== nothing
        (elem1, state1) = next1
        for i = 1:n, j = 1:n
            A1[i, j] = (elem1[(i - 1) * n + j] == 0) ? A[i, j].hi : A[i, j].lo
        end
        next2 = iterate(list, state1)
        det_A1 = det(A1)
        while next2 !== nothing
            (elem2, state2) = next2
            for i = 1:n, j = 1:n
                A2[i, j] = (elem2[(i - 1) * n + j] == 0) ? A[i, j].hi : A[i, j].lo
            end
            counter = counter + 1
            if (det_A1 * det(A2) <= 0)
                res = true
                break
            end
            next2 = iterate(list, state2)
        end

        if res
            return res
        end
        next1 = iterate(list, state1)
    end
    return res
end

function determine_degeneracy_Beck(A, n)
    M = [mid(A[i, j]) for i=1:n, j=1:n]
    R = [radius(A[i, j]) for i=1:n, j=1:n]
    applicable = det(M) != 0
    if !applicable
        return (applicable, false)
    end

    M1 = [abs(inv(M)[i, j]) for i=1:n, j=1:n]
    F = M1 * R

    λ = eigvals(F)
    m = size(λ, 1)
    ρ = abs(λ[m])
    if ρ < 1
        return (applicable, false)
    end
    # Beck 2
    F = R * M1
    diag_max = F[1, 1]
    for i = 2:n
        if (diag_max < F[i, i])
            diag_max = F[i, i]
        end
    end
    if (diag_max >= 1)
        return (applicable, true)
    end

    return (false, false)
end

function determine_degeneracy_problem_2(n, ε)
    A = Array{Interval{Float64}, 2}(undef, n, n)
    for i = 1:n, j = 1:n
        A[i, j] = i == j ? 1 : interval(0.0, ε)
    end

    beck_applicable, degenerate = determine_degeneracy_Beck(A, n)
    if beck_applicable
        return degenerate
    end
    return determine_degeneracy_Baumans(A, n)
end
ε_list_1 = Any[]
ε_list_1 = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,]  # Reverse order for problem 2
ε_list_2 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
ε_list_3 = [0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029]
ε_list_4 = [0.024, 0.0241, 0.0242, 0.0243, 0.0244, 0.0245, 0.0246, 0.0247, 0.0248, 0.0249]
list_of_ε_lists = Any[]
push!(list_of_ε_lists, ε_list_1)

problem_num = 1
if problem_num == 1
    reverse!(ε_list_1)
    push!(list_of_ε_lists, ε_list_2)
    push!(list_of_ε_lists, ε_list_3)
    push!(list_of_ε_lists, ε_list_4)
    for ε_list in list_of_ε_lists
        for ε in ε_list
            @printf("ε: %f, degenerate: %s\n", ε, determine_degeneracy_problem_1(ε))
        end
        println();
    end
    println("\nProblem 1")
elseif problem_num == 2
    n = 4
    δ_max = 0.00001
    δ = list_of_ε_lists[1][1] - list_of_ε_lists[1][2]
    while (δ > δ_max)
        ε_list = pop!(list_of_ε_lists)
        global δ = ε_list[1] - ε_list[2]
        println(ε_list)
        println("δ: ", δ)
        index = 1
        while(index <= sizeof(ε_list))
            global ε = ε_list[index]
            degenerate = determine_degeneracy_problem_2(n, ε)
            @printf("ε: %f, degenerate: %s\n", ε, degenerate)
            if !degenerate
                global ε_ans = ε_list[index - 1]
                #new_list = 
                #println(new_list)
                push!(list_of_ε_lists, collect(ε_list[index - 1]:-(δ / 10):ε))
                break;
            end
            index += 1
        end
        println();
    end
    println("\nProblem 2")
    println("ε = ", ε_ans, ", n = ", n, ", with δ = ", δ_max)
end
