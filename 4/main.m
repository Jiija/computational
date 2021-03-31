A = [2, 1, 4;
     2.1, 1, 4;
     8, 0, 0];
inf_b = [0; 6;-2];
sup_b = [1; 7; 0];

[tolmax, argmax, envs, ccode] = tolsolvty(A, A, inf_b, sup_b);

display(tolmax)

diag_rad_b = diag(0.5 * (sup_b - inf_b));

C = [-A, -diag_rad_b;
     A,  -diag_rad_b];
 
 mid_b = 0.5 * (inf_b + sup_b);
 r = cat(1, -mid_b, mid_b);
 r = squeeze(r);
 
 min_vals = [-Inf, -Inf, -Inf, 0, 0, 0];
 max_vals = [Inf, Inf, Inf, Inf, Inf, Inf];
 options = optimoptions('linprog','Algorithm','interior-point');
 x = linprog([0, 0, 0, 1, 1, 1], C, r, [], [], min_vals, max_vals, options);
 display(x)
 display(det(C))