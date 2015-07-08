% Solve the LP
%   minimize( c' * x )
%   subject to
%       A * x == b
%       x >= 0
% given a feasible starting point x_0
%
% [x_opt, history] = lp_solve_feasible(A, b, c, x_0)
function [x, history] = lp_solve_feasible(A, b, c, x, mu)
    if nargin < 5
        mu = 50;
    end

    t = 1;
    n = length(x);
    c_steps = 1;
    while true
        % Centering step
        [x, ~, n_steps, ~] = lp_center(A, b, t * c, x);
        history(1, c_steps) = n_steps;
        history(2, c_steps) = n / t;

        % Check stopping criterion
        if n / t < 1e-3
            return
        end

        % Update
        t = mu * t;
        c_steps = c_steps + 1;
    end

end

