[A, b, c, x_0] = lp_generate(100, 500);
[x_opt, opt_val] = lp_solve_feasible(A, b, c, x_0);
x_opt
opt_val
