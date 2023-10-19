% 4. find the summation function 
% inputs to the function: input_vector, noisy_set, noisy_image
function g = numerical_gradient_new_trial(input_vector, noisy_set, noisy, class_mat, ij_to_p)
  
    g = zeros(size(input_vector));

    [m, n] = size(noisy);
    L = size(noisy_set, 1);    % to be fixed   # m is the length of th noise set
    alfa = @(a) sqrt(double(a.^2 + 1));
    noisy_image = noisy;
   
    for p = 1:L                % for each point in the noise set 
        g1 = 0;
        g2 = 0;
        
        i = noisy_set(p,1);
        j = noisy_set(p,2);
        % finding the left summation f1 
        isInMatrix = class_mat(p,1);
        if ~isInMatrix    % first nieghboring point 
            %f1 = f1 + alfa(input_vector(p) - noisy_image(min(i+1, m), j));  
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(min(i+1, m), j)))/(2*((double(input_vector(p)) - double(noisy_image(min(i+1, m), j))).^2 + 1)^(1/2));
        else
            p2 = ij_to_p(min(i+1, m), j);
            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 

        end
    
        isInMatrix = class_mat(p,2);
        if ~isInMatrix   % second nieghboring point 
            %f1 = f1 + alfa(input_vector(p) - noisy_image(min(i+1, m),min(j+1, n))); 
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(min(i+1, m), min(j+1, n))))/(2*((double(input_vector(p)) - double(noisy_image(min(i+1, m), min(j+1, n))))^2 + 1)^(1/2));
        else
            p2 = ij_to_p(min(i+1, m), min(j+1, n));

            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 

        end
        isInMatrix = class_mat(p,3);
        if ~isInMatrix    %  nieghboring point 
           % f1 = f1 + alfa(input_vector(p) - noisy_image(max(i-1, 1), j)); 
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(max(i-1, 1), j)))/(2*((double(input_vector(p)) - double(noisy_image(max(i-1, 1), j)))^2 + 1)^(1/2));
        else
            p2 = ij_to_p(max(i-1, 1), j);
            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 

        end
    
        isInMatrix = class_mat(p,4);
        if ~isInMatrix   %  nieghboring point 
           % f1 = f1 + alfa(input_vector(p) - noisy_image(i, min(j+1, n)));
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(i, min(j+1, n))))/(2*((double(input_vector(p)) - double(noisy_image(i, min(j+1, n))))^2 + 1)^(1/2));
        else
            p2 = ij_to_p(i, min(j+1, n));
            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 

        end
    
        isInMatrix = class_mat(p,5);
        if ~isInMatrix    %  nieghboring point 
          %  f1 = f1 + alfa(input_vector(p) - noisy_image(i, max(j-1, 1)));
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(i, max(j-1, 1))))/(2*((double(input_vector(p)) - double(noisy_image(i, max(j-1, 1))))^2 + 1)^(1/2));
        else
            p2 = ij_to_p(i, max(j-1, 1));
            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 

        end
    
        isInMatrix = class_mat(p,6);
        if ~isInMatrix   %  nieghboring point 
            %f1 = f1 + alfa(input_vector(p) - noisy_image(max(i-1, 1), min(j+1, n)));   
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(max(i-1, 1), min(j+1, n))))/(2*((double(input_vector(p)) - double(noisy_image(max(i-1, 1), min(j+1, n))))^2 + 1)^(1/2));
        else
            p2 = ij_to_p(max(i-1, 1), min(j+1, n));
            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 

        end
    
        isInMatrix = class_mat(p,7);
        if ~isInMatrix    %  nieghboring point 
            %f1 = f1 + alfa(input_vector(p) - noisy_image(max(i-1, 1), max(j-1, 1)));  
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(max(i-1, 1), max(j-1, 1))))/(2*((double(input_vector(p)) - double(noisy_image(max(i-1, 1), max(j-1, 1))))^2 + 1)^(1/2));
        else
            p2 = ij_to_p(max(i-1, 1), max(j-1, 1));
            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 

        end
    
        isInMatrix = class_mat(p,8);
        if ~isInMatrix   % nieghboring point 
            %f1 = f1 + alfa(input_vector(p) - noisy_image(min(i+1,m), max(j-1, 1)));  
            g1 = g1 + (2*double(input_vector(p)) - 2*double(noisy_image(min(i+1,m), max(j-1, 1))))/(2*((double(input_vector(p)) - double(noisy_image(min(i+1,m), max(j-1, 1))))^2 + 1)^(1/2));
        else
            p2 = ij_to_p(min(i+1,m), max(j-1, 1));
            g2 = g2 + (2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2));
            g(p2) =g(p2) - 0.5*(2*double(input_vector(p)) - 2*double(input_vector(p2)))/(2*((double(input_vector(p)) - double(input_vector(p2)))^2 + 1)^(1/2)); 
        end
    
    
    
    
    
    g(p) = g(p) + g1 + 0.5*g2;
    
    end
  g = g';

end
