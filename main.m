%% Variables declaration

% Update n for different size of matrix A.  
n = 500;

% Construction of matrix A, uncomment seed for debugging/same results
seed = 0; 
rng(seed);
A = rand(n,n);

% Preallocation of added columns in the SVD for faster execution time
a = cell(n-1, 1);
for i = 1:n-1
    a{i} = A(i+1:end, i+1);
end

%Preallocation of array to hold numerical error each iteration
maximum_error = zeros(1,n-1);

% Input variables to Adaptive SVD algorithm, SVD of first column of A
[U,S,V] = svd(A(:,1));


% Code to add noise to input matrices
% noise = 0.01;
% U = U +noise*randn(size(U));
% S = S +noise*randn(size(S));
% V = V +noise*randn(size(V));

%Preallocation of array to hold time per execution
adaptive_timing = zeros(1,n-1);

%%  Execution of Adaptive of SVD 
t_adaptive_start = tic;      
for i = 1:n-1
    tic    
    % [U,S,V] = downdate_SVD(U,S,V);
    % [U,S,V] = update_SVD(U,S,V,a{i});
    % 
    %For minor speed improvements, use this optimalized function
    [U,S,V] = combined_SVD(U,S,V,a{i}); 

    maximum_error(i) = (max(max(abs(U*S*V'-A(i+1:end,1:i+1)))));
    adaptive_timing(i)= toc;
end
total_adaptive_time = toc(t_adaptive_start);


%% Execution of matlab SVD
matlab_svd_timing = zeros(1,n-1);
t_matlab_start = tic;           
for i=1:n-1
    tic
    [U,S,V] = svd(A(i:end,1:i+1));
    matlab_svd_timing(i)= toc;
end
total_matlab_time = toc(t_matlab_start);

%% Display maximum error and total time used
disp("Maximum error:");
disp(max(maximum_error));
disp("Total time used for Adatptive SVD:");
disp(total_adaptive_time);
disp("Total time used for Matlab SVD:");
disp(total_matlab_time);

%% Plots
% Plot to show timing of SVD
figure('Position', [100 100 500 600]);
title('Timing of SVD')
subplot(3, 1, 1);
hold on;
title('Adaptive SVD');
plot(adaptive_timing);
xlabel('Iteration');
ylabel('Time [s]');
legend('Seconds');

subplot(3, 1, 2);
hold on;
title('Matlab SVD');
plot(matlab_svd_timing);
xlabel('Iteration');
ylabel('Time [s]');
legend('Seconds');

subplot(3,1,3);
hold on;
title('Difference');
plot(adaptive_timing-matlab_svd_timing);
xlabel('Iteration');
ylabel('Time [s]');
legend('Seconds');
% saveas(gca,'time.jpg');


% Plot to show Reconstruction  error
figure;
semilogy(maximum_error);
title("Reconstruction accuracy");
xlabel('Iteration');
ylabel('Accuracy');
legend("A - UΣV^T");
% saveas(gca,'noisy4.jpg');

