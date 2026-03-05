function [z,iter] = QPAS(H, f)
n = size(H, 1);
tol = 1e-6;
max_iter = 100;

active_set = false(n, 1);
z = ones(n, 1) / n; 

for iter = 1:max_iter
    F = find(~active_set); 
    H_F = 2 * H(F, F);
    f_F = f(F);
    A = [H_F, -ones(length(F), 1); ones(1, length(F)), 0]; 
    b = [-f_F; 1];

    x_lambda = A \ b;
    z_F = x_lambda(1:end-1);
    lambda = x_lambda(end);

    z_star = zeros(n, 1);
    z_star(F) = z_F;

    violation = z_star(F) < -tol;
    if ~any(violation)
        mu = 2 * H * z_star + f - lambda * ones(n, 1);
        mu_A = mu(active_set);

        if all(mu_A >= -tol) 
            z=z_star;
            return;
        else
            [min_mu, idx] = min(mu_A);
            if min_mu < -tol
                active_set(idx) = false;    
            end
        end
    else
        d = z_star - z;
        alpha = 1;
        blocking_idx = [];

        for i = 1:length(F)
            j = F(i);
            if d(j) < -tol  
                alpha_i = (0 - z(j)) / d(j);
                if alpha_i < alpha
                    alpha = alpha_i;
                    blocking_idx = j;
                end
            end
        end

        alpha = min(alpha, 1);
        z = z + alpha * d;
        if ~isempty(blocking_idx)
            active_set(blocking_idx) = true;
        end
    end
end 
z=z';
end