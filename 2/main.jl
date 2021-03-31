using IntervalArithmetic
using Plots

n = []
widthX = []
widthF = []

rosenbrock(X) = 100 * (X[1]^2 - X[2])^2 + (X[1] - 1)^2

# rasstrigin(X) = X[1]^2 + X[2]^2 - cos(18 * X[1]) - cos(18 * X[2])

shaffer(X) = 0.5 + (sin(X[1]^2 + X[2]^2)^2 - 0.5) / (1 + 0.001  * (X[1]^2 + X[2]^2))^2

function globOpt(X, func, iter_num)
    Y = func(X)
    iteration = 1
    lead_est = Y.lo
    R = copy(X)
    work_list = [Dict("Box" => X, "Est" => Y.lo)]

    while iteration <= iter_num
        list_lenght = size(work_list)[1]
        lead_est = work_list[1]["Est"]
        lead_ind = 1
        for k = 1:list_lenght
            est = work_list[k]["Est"]
            if (est < lead_est)
                lead_est = est
                lead_ind = k
            end
        end
        R = work_list[lead_ind]["Box"]

        D1 = copy(R)
        D2 = copy(D1)
        (rad_max, ind_max) = findmax(radius.(D1))

        s = D1[ind_max]
        sl = s.lo
        sh = s.hi
        sm = mid(s)
        D1[ind_max] = @interval(sl, sm)
        D2[ind_max] = @interval(sm, sh)

        Y1 = func(D1)
        Y2 = func(D2)

        rec1 = Dict("Box" => D1, "Est" => Y1.lo)
        rec2 = Dict("Box" => D2, "Est" => Y2.lo)

        println("iteration: ", iteration, ", F: ", func(R))

        push!(n, iteration)
        push!(widthX, 2 * sum(radius.(R)))
        push!(widthF, 2 * radius(func(R)))

        deleteat!(work_list, lead_ind)

        push!(work_list, rec1)
        push!(work_list, rec2)
        iteration += 1
    end

    for rec in work_list
        print("[ ")
        for cmp in rec["Box"]
            print(cmp, " ")
        end
        println("], ", rec["Est"])
    end

    return (lead_est, R)
end

X = [@interval(-30, 30), @interval(-30, 30)]
res = globOpt(X, rosenbrock, 20)
println("-----------------------------------------------")
println(res)
println("-----------------------------------------------")
# Plot leading interval widths
plot(n, widthX,
    xticks = ([0:1:20;]),
    yticks = ([0:15:120;]),
    xlabel = "n",
    ylabel = "WidthX",
    legend = false)
png("plot1")

X = [@interval(-15, 15), @interval(-15, 15)]


# Plot function value widths for rosenbrock
global widthF = []
global n = []
globOpt(X, rosenbrock, 120)
plot(n, widthF,
    xticks = ([0:10:120;]),
    xlabel = "n",
    ylabel = "WidthF",
    yscale = :log10,
    legend = false)
png("plot2")

#Plot function value widths for shaffer
global widthF = []
global n = []
globOpt(X, shaffer, 120)
plot(n, widthF,
    xticks = ([0:10:120;]),
    xlabel = "n",
    ylabel = "WidthF",
    yscale = :log10,
    legend = false)
png("plot3")
