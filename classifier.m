% 4. find the summation function 
% inputs to the function: input_vector, noisy_set, noisy_image
function class_mat = classifier(noisy_set, noisy)
  
    [m, n] = size(noisy);
    L = size(noisy_set, 1);    % to be fixed   # m is the length of th noise set
    alfa = @(a) sqrt(double(a.^2 + 1));
    f1 = 0 ;
    f2 = 0 ;
    class_mat = []; 
    for p = 1:L                % for each point in the noise set 
        i = noisy_set(p,1);
        j = noisy_set(p,2);
        % finding the left summation f1 
        isInMatrix = any(all(noisy_set == [min(i+1, m), j], 2));
        if ~isInMatrix    % first nieghboring point 
            class_mat(p, 1) = 0;   % 1 for noise and 0 for not noise
        else
            class_mat(p, 1) = 1;
            
        end
    
        isInMatrix = any(all(noisy_set == [min(i+1, m), min(j+1, n)], 2));
        if ~isInMatrix   % second nieghboring point 
            class_mat(p, 2) = 0;
        else
            class_mat(p, 2) = 1;
            
        end
        isInMatrix = any(all(noisy_set == [max(i-1, 1), j], 2));
        if ~isInMatrix    %  nieghboring point 
            class_mat(p, 3) = 0;
        else
            class_mat(p, 3) = 1;
        end
    
        isInMatrix = any(all(noisy_set == [i, min(j+1, n)], 2));
        if ~isInMatrix   %  nieghboring point 
            class_mat(p, 4) = 0;
        else
            class_mat(p, 4) = 1;
        end
    
        isInMatrix = any(all(noisy_set == [i, max(j-1, 1)], 2));
        if ~isInMatrix    %  nieghboring point 
            class_mat(p, 5) = 0;
        else
            class_mat(p, 5) = 1;
        end
    
        isInMatrix = any(all(noisy_set == [max(i-1, 1), min(j+1, n)], 2));
        if ~isInMatrix   %  nieghboring point 
            class_mat(p, 6) = 0;
        else
            class_mat(p, 6) = 1;
        end
    
        isInMatrix = any(all(noisy_set == [max(i-1, 1), max(j-1, 1)], 2));
        if ~isInMatrix    %  nieghboring point 
            class_mat(p, 7) = 0;
        else
            class_mat(p, 7) = 1;
        end
    
        isInMatrix = any(all(noisy_set == [min(i+1,m), max(j-1, 1)], 2));
        if ~isInMatrix   % nieghboring point 
            class_mat(p, 8) = 0; 
        else
            class_mat(p, 8) = 1;
        end
 
    
    
    end
end
