function [min_cond] = get_min_cond(inf_A, sup_A)
    m = size(inf_A, 1);
    n = size(inf_A, 2);
    
    min_cond = Inf;
    num_of_tests = 20;
    for k = 1:num_of_tests
        M1 = ones(m, n);
        M2 = ones(m, n);
        epm = randi([0,1], m, n);
        for i = 1:m
            for j = 1:n
                if epm(i, j) == 0
                    M1(i, j) = inf_A(i, j);
                    M2(i, j) = sup_A(i, j);
                else
                    M1(i, j) = sup_A(i, j);
                    M2(i, j) = inf_A(i, j);
                end
            end
        end
        c1 = cond(M1, 2);
        c2 = cond(M2, 2);
        if min_cond > c1
            min_cond = c1;
        end
        if min_cond > c2
            min_cond = c2;
        end
    end
end