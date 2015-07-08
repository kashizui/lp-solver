function [A, b, c, x_0] = lp_generate(m, n)
    % Generate random problem data
    while true
        % Ensure A is full rank and sublevel sets are bounded
        A = [randn(m - 1, n); rand(1, n)];
        if rank(A) == m
            break
        end
    end
    x_0 = rand(n, 1);
    b = A * x_0;
    c = randn(n, 1);
end
