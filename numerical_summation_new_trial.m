% 4. find the summation function 
% inputs to the function: input_vector, noisy_set, noisy_image
function f = numerical_summation_new(input_vector, noisy_set, class_mat, noisy, ij_to_p)
  
    [m, n] = size(noisy);
    L = size(noisy_set, 1);    % to be fixed   # m is the length of th noise set
    alfa = @(a) sqrt(double(a.^2 + 1));
    noisy_image = noisy;
    f1 = 0 ;
    f2 = 0 ;
    for p = 1:L                % for each point in the noise set 
        i = noisy_set(p,1);
        j = noisy_set(p,2);
        % finding the left summation f1 
        isInMatrix = class_mat(p,1);
        if ~isInMatrix    % first nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(min(i+1, m), j));
        else
            p2 = ij_to_p(min(i+1, m), j);
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 
            
        end
    
        isInMatrix = class_mat(p,2);
        if ~isInMatrix   % second nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(min(i+1, m),min(j+1, n)));
        else
            p2 = ij_to_p(min(i+1, m),min(j+1, n));
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 

            
        end
        isInMatrix = class_mat(p,3);
        if ~isInMatrix    %  nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(max(i-1, 1), j)); 
        else
            p2 = ij_to_p(max(i-1, 1), j);
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 

        end
    
        isInMatrix = class_mat(p,4);
        if ~isInMatrix   %  nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(i, min(j+1, n))); 
        else
            p2 = ij_to_p(i, min(j+1, n));
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 


        end
    
        isInMatrix = class_mat(p,5);
        if ~isInMatrix    %  nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(i, max(j-1, 1))); 
        else
            p2 = ij_to_p(i, max(j-1, 1));
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 
  

        end
    
        isInMatrix = class_mat(p,6);
        if ~isInMatrix   %  nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(max(i-1, 1), min(j+1, n))); 
        else
            p2 = ij_to_p(max(i-1, 1), min(j+1, n));
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 
  

        end
    
        isInMatrix = class_mat(p,7);
        if ~isInMatrix    %  nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(max(i-1, 1), max(j-1, 1)));
        else
            
            p2 = ij_to_p(max(i-1, 1), max(j-1, 1));
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 
  

        end
    
        isInMatrix = class_mat(p,8);
        if ~isInMatrix   % nieghboring point 
            f1 = f1 + alfa(input_vector(p) - noisy_image(min(i+1,m), max(j-1, 1)));
        else
            p2 = ij_to_p(min(i+1,m), max(j-1, 1));
            f2 = f2 + alfa(input_vector(p) - input_vector(p2)); 
  

        end
    
    
    end
    % find eval_function 
    f = f1 + 0.5*f2;
end
