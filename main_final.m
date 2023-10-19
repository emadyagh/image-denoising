close all;
clear all;
clc;
syms x


I = imread('/MATLAB Drive/denoising/20230519164757_00000.png');
%I = imread('/MATLAB Drive/denoising/test1.png');

%desiredSize = [50, 50]; % Replace newHeight and newWidth with the desired dimensions

% Resize the image to the desired size using imresize
%resizedImage = imresize(I, desiredSize);
% crop the image: 
% Define the cropping rectangle [x, y, width, height]
cropRect = [0, 0 , 300, 300]; % Adjust these values

% Crop the image
croppedImage = imcrop(I, cropRect);


original_image = rgb2gray(croppedImage);
Full_noisy = imnoise(original_image, 'salt & pepper', 0.3);

un_denoised = Full_noisy
% Define the size of the sliding window (e.g., 8x8)
window_size = 100;
step = 100;


% Get the size of the image
[rows, cols] = size(Full_noisy);
M = rows
N = cols
% Start the timer
tic;

% Define a step size for the sliding window (e.g., 1)
for loop = 1:1
    loop
       if loop == 1
            Smin = 0
            Smax = 255 
            n_threshold = 0 
       end
    % Iterate over rows and columns with the specified step
     for r = 1:step:rows
        
        for c = 1:step:cols
            
            % 1. Extract the current sliding window
            %noisy = Full_noisy(r:r+window_size-1, c:c+window_size-1);
            noisy = Full_noisy(r:min(r+window_size-1, rows), c:min(c+window_size-1, cols));
            
            
            % 2.  finding the adoptive median filter: 
            [m, n] = size(noisy);
    
            median_image = zeros(m, n);
            median_image = uint8(median_image);
            
            % Define initial and maximum window sizes
            initialWindowSize = 36; % Initial window size (3x3)
            maxWindowSize = 43;     % Maximum window size (7x7)
            
            for i = 1:m
                for j = 1:n
                    windowSize = initialWindowSize; % Start with the initial window size
            
                    while windowSize <= maxWindowSize
                        % Define the boundaries for the window around each pixel
                        xmin = max(1, i - floor(windowSize/2));
                        xmax = min(m, i + floor(windowSize/2));
                        ymin = max(1, j - floor(windowSize/2));
                        ymax = min(n, j + floor(windowSize/2));
            
                        % Extract the neighborhood matrix
                        temp = noisy(xmin:xmax, ymin:ymax);
            
                        % Calculate the median of this matrix
                        medianValue = median(temp(:));
            
                        % Check the conditions for adaptive median filter
                        if noisy(i, j) == 0 || noisy(i, j) == 255 || (noisy(i, j) < medianValue && noisy(i, j) > 0) || (noisy(i, j) > medianValue && noisy(i, j) < 255)
                            median_image(i, j) = medianValue;
                            break; % Stop filtering
                        else
                            windowSize = windowSize + 2; % Increase the window size
                        end
                    end
    
                    % If the while loop completes without satisfying conditions, keep the original pixel value
                    if windowSize > maxWindowSize
                        median_image(i, j) = noisy(i, j);
                    end
                end
            end
    
            % 3. identify the noisy pixels 
            noisy_set = [];
            normal_set = []; % the non noisy set
            
            for i = 1:m
                for j = 1:n
                    % Check the conditions for noisy pixels
                   
                    
                
                        if (noisy(i, j) <= Smin || noisy(i, j) >= Smax) && (noisy(i, j) ~= median_image(i, j))
                       % if (noisy(i, j) <= Smin || noisy(i, j) >= Smax) && (abs(noisy(i, j) - median_image(i, j))<=n_threshold)   
                            noisy_set = [noisy_set; i, j];
                         
                        end
                    
                        
                   

                end
            end

            
            class_mat = classifier(noisy_set, noisy);


                        % finding input_vector 
            input_vector = [];
            ij_to_p = [];
            L = size(noisy_set, 1);    % to be fixed   # m is the length of th noise set

            for p = 1:L                % for each point in the noise set 
                i = noisy_set(p,1);
                j = noisy_set(p,2);
                ij_to_p(i, j) = p; 
                %input_vector = [input_vector, noisy(i, j)];
                input_vector = [input_vector, 125];
            end

    
            % 4. find the summation function 
            %f = numerical_summation(input_vector, noisy_set, noisy)
            %g = numerical_gradient_2(input_vector, noisy_set, noisy)

            
            % 5. call the optimization algorithm
            x0 = [input_vector] % Your initial point here
            delta = 1e-4 %0.001%1e-4; % You can adjust these values as needed   1e-4
            sigma = 0.9 %0.009%0.9;  %0.9
            
            if ~isempty(x0)
                % Call the optimization function
                [x, iterations] = conjugate_gradient_method7_trail_wolf(x0, delta, sigma, noisy_set, noisy, class_mat, ij_to_p);
            
           
                 %6. update noisy image with the optimized solution:
                for p = 1:L
                    Full_noisy(noisy_set(p,1) + r - 1, noisy_set(p,2) + c - 1) = x(p);
                end
            end            
            
            
            
        end
    end

end

imshow(un_denoised)
figure
imshow(Full_noisy)
figure
imshow(original_image)
% Stop the timer
elapsedTime = toc;

% Display the processing time in seconds
disp(['Processing time: ', num2str(elapsedTime), ' seconds']);
iterations
% Calculate the Frobenius norm of the difference between X1 and X2
frobenius_norm = norm(double(Full_noisy) - double(original_image), 'fro');
PSNR = 10 * log10(255^2 / ((1 / (M * N)) * frobenius_norm^2))

%frobenius_norm2 = norm(double(un_denoised) - double(original_image), 'fro');
%PSNR2 = 10 * log10(255^2 / ((1 / (M * N)) * frobenius_norm2^2))


% Display or use the Frobenius norm value as needed
disp(['PSNR: ', num2str(PSNR)]);

% relative error 

% Calculate Mean Squared Errors (MSE)
%MSE_noisy = sum((double(original_image(:)) - double(un_denoised(:))).^2) / numel(original_image);
%MSE_denoised = sum((double(original_image(:)) - double(Full_noisy(:))).^2) / numel(original_image);

% Calculate the Relative Error (RE)
%RE =  MSE_noisy / MSE_denoised

% Calculate the relative error
%relativeError = 100*norm(double(Full_noisy) - double(original_image), 'fro') / norm(double(original_image), 'fro')
%...............

% Calculate the sum of squared pixel differences
sum_squared_pixel_diff = sum((double(Full_noisy(:)) - double(un_denoised(:))).^2);

% Calculate the sum of squared original pixel values
sum_squared_noisy_pixels = sum(double(un_denoised(:)).^2);

% Calculate the relative error
relativeError =   sum_squared_pixel_diff*5/sum_squared_noisy_pixels

%  ....................... MSE ..........................
MSE = (1 / (M * N)) * sum((double(Full_noisy(:)) - double(un_denoised(:))).^2)

ssimValue = ssim(Full_noisy, original_image)

