% Solve the centering problem
%   minimize( c' * x - sum(log(x)))
%   subject to
%       A * x == b
%
% [x_opt, v_opt, n_steps] = lp_center(A, b, c, x_0)
function [x, v, n_steps, opt_gap] = lp_center(A, b, c, x, alpha, beta)
    if nargin < 6
        alpha = 0.25;
        beta = 0.9;
    end

    f = @(z) c' * z - sum(log(z));

    n_steps = 0;
    while true
        % Compute derivatives
        f_grad = c - (1 ./ x);
        f_hess_inv = diag(x.^2);

        % Compute Newton step using reduced equation
        S = A * f_hess_inv * A';
        b = A * (f_hess_inv * -f_grad);
        v = S \ b;
        x_step = f_hess_inv * (-f_grad - (A' * v));

        % Check stopping criterion
        sq_decr = -x_step' * f_grad;
        opt_gap(n_steps + 1) = sq_decr / 2;
        if sq_decr / 2 < 1e-6
            break
        end

        % Backtracking line search to choose t
        t = 1;
        fx = f(x);
        search_slope = alpha * f_grad' * x_step;
        while true
            if t < 1e-10
                fprintf('ERROR: t too small, giving up...\n');
                return
            end

            x_next = x + t * x_step;
            if all(x_next > 0) && f(x_next) < fx + t * search_slope
                break
            end
            t = beta * t;
        end

        % Update
        x = x_next;
        n_steps = n_steps + 1;
    end
end
