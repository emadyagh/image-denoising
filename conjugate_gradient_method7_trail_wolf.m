function [x, iterations] = conjugate_gradient_method7_trail_wolf(x0, delta, sigma, noisy_set, noisy, class_mat, ij_to_p)
    % Step 0
    x = x0';  % x0 is a vector [x1, x2]
    w = 0.01;
    t_bar = 0.3;
    k = 0;
    g0 = double(numerical_gradient_new_trial(x',noisy_set, noisy, class_mat, ij_to_p));  % input is horizontal [1 1 1 1 1 1] 
    d = -g0;
    
    % Initialize x_prev and g_prev here
    x_prev = x;
    g = g0;
    
    % Initialize t here
    t = t_bar;
    iterations  = 0 
    w_con = 0.01
    f_val = 10^9;       %initialization
    f_val_pre = -10^9; % initialization
    while iterations < 300  %&& (abs(f_val - f_val_pre)/abs(f_val)) > 10^-4
        iterations = iterations + 1 
        % Step 2
        if k > 0
             % Make sure d and y are column vectors
            w = max(max(w_con * norm(d_prev) * norm(y_prev), (d_prev' * y_prev)), -d_prev' * g_prev);
            
            beta = ((g' * y_prev) / w) - (norm(y_prev)^2*g'*d_prev/w^2);

            t = min(t_bar, max(0, (y_prev' * (y_prev - s_prev)) / (norm(y_prev)^2)));

            gamma = (t * (g' * d_prev)) / w;
            d = -g + beta*d_prev + gamma*y_prev;
        end
        
        % Step 3: Wolfe Line Search
        f_val_pre = f_val;
        %[alpha, f_val] = wolfe_line_search(x, d, delta, sigma, noisy_set, noisy, class_mat, ij_to_p)
        alpha = 1;
        % Step 4
        x = double(x) + double(alpha * d);
        g_new = double(numerical_gradient_new_trial(x',noisy_set, noisy, class_mat, ij_to_p));
        
        % Step 5
        s = double(x) - double(x_prev);
        y_prev = g_new - g;
        
        % Step 6
        k = k + 1;
        
        % Store previous values for the next iteration
        s_prev = s;
        x_prev = double(x);
        g_prev = g;
        g = g_new;
        d_prev = d;
        %fval = double(numerical_summation(x', noisy_set, noisy))
    end
    
    %fval = double(numerical_summation(x', noisy_set, noisy))
end


function [alpha, f0] = wolfe_line_search(x, d, delta, sigma, noisy_set, noisy, class_mat, ij_to_p)
    alpha = 10.0;
    max_iterations = 20;   % old value is 100
    %c1 = 1e-4;
    %c2 = 0.9;
    
    f0 = double(numerical_summation_new_trial(x', noisy_set, class_mat, noisy, ij_to_p));
    g0 = double(numerical_gradient_new_trial(x',noisy_set, noisy, class_mat, ij_to_p));
    
    for i = 1:max_iterations
        alpha;
        
        x_new = double(x) + double(alpha * d);
        f_new = double(numerical_summation_new_trial(x_new', noisy_set, class_mat, noisy, ij_to_p));
        %g_new = double(numerical_gradient_new_trial(x_new', noisy_set, noisy, class_mat, ij_to_p));
        
        if f_new <= f0 + delta * alpha * g0' * d %&& (g_new' * d >= sigma * g0' * d)
            break;
        end
        
        alpha = alpha * 0.9;  % Backtracking   old value is 0.5
     end
end

