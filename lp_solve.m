% Solve the LP
%   minimize( c' * x )
%   subject to
%       A * x == b
%       x >= 0
%
% [x_opt, history] = lp_solve(A, b, c)
function [x_opt, opt_val] = lp_solve(A, b, c)
    % Phase I method
    n = length(c);
    x_0 = A \ b;
    t_0 = 2 - min(x_0);
    z_0 = [t_0; x_0 + (t_0 - 1)];
    [z, ~] = lp_solve_feasible(...
            A * [-ones(n, 1), eye(n)], ...
            b - sum(A, 2), ...
            [1; zeros(n, 1)], ...
            z_0);

    t_opt = z(1);
    x_0 = z(2:end) + (1 - t_opt);

    % Problem is strictly feasible iff t* < 1
    if t_opt < 1
        % Solve problem using feasible x_0
        [x_opt, ~] = lp_solve_feasible(A, b, c, x_0);
        opt_val = c' * x_opt;
    else
        % Report infeasibility
        x_opt = NaN;
        opt_val = +Inf;
    end
end
