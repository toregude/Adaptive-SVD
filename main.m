close all;
% %%Rows
% A = rand(200,200);
% [U,S,V] = svd(A);

% %%Columns
% A_full = rand(100,200);
% [U,S,V] = svd(A_full(:,1));

% % Both
n = 1000;
% seed = 0;  % Replace 123 with your desired seed value
rng(seed);
A = rand(n,n);
[U,S,V] = svd(A(:,1));

% noise = 0.001; % Adjust the noise level as desired
% A_noisy = A + noise * randn(size(A));


U_dash = U;
S_dash = S;
V_dash = V;


%Code for checking adding of columns
% m = zeros(1,199);
% for i = 2:200
%     a = A_full(:,i);
%     [U_dash,S_dash,V_dash] = update_SVD(U_dash,S_dash,V_dash,a);
% 
%     m(i-1) = (max(max(U_dash*S_dash*V_dash'-A_full(:,1:i))));
% end
% disp(max(m));
% disp(sum(isnan(m(:))));


% %Code for checking deletion of rows
% m = zeros(1,199);
% for i = 1:199
%     [U_dash,S_dash,V_dash] = downdate_SVD(U_dash,S_dash,V_dash);
% 
%     m(i) = (max(max(U_dash*S_dash*V_dash'-A(i+1:end,:))));
% end
% disp(max(m));
% disp(sum(isnan(m(:))));

% %Code for checking both

% T2 = zeros(1,n-1);
% tStart2 = tic;           % pair 2: tic
% for i=1:n-1
%     tic
%     [U,S,V] = svd(A(i:end,1:i+1));
%     T2(i)= toc;
% end
% 
% tEnd2 = toc(tStart2);
% disp(tEnd2);



m = zeros(1,n-1);

a = cell(n-1, 1);
for i = 1:n-1
    a{i} = A(i+1:end, i+1);
end

noise = 0.1;

[U,S,V] = svd(A(:,1));

U = U +noise*randn(size(U));
S = S +noise*randn(size(S));
V = V +noise*randn(size(V));


T = zeros(1,n-1);
tStart = tic;           % pair 2: tic
for i = 1:n-1
    tic    
    [U,S,V] = downdate_SVD(U,S,V);
    [U,S,V] = update_SVD(U,S,V,a{i});
    
    % [U,S,V] = combined_SVD(U,S,V,a{i});
    m(i) = (max(max(abs(U*S*V'-A(i+1:end,1:i+1)))));
    T(i)= toc;
end
tEnd = toc(tStart);
disp(max(abs(m)));
disp(sum(isnan(m(:))));
disp(tEnd);





% figure;
% hold on;
% set(0,'defaultfigurecolor',[1 1 1]) % white background
% set(0,'defaultaxesfontname','Ariel') % beautify the axes a bit
% title('Adaptive SVD');
% plot(T);
% xlabel('Iteration');
% ylabel('Time [s]');
% legend("seconds");

% figure;
% hold on;
% title('Matlab SVD');
% plot(T2);
% xlabel('seconds');
% ylabel('Time [s]');
% legend("A - UΣV^T");

% figure;
% hold on;
% title('Difference');
% plot(T-T2);

figure;
semilogy(m);
title("Reconstruction accuracy polluted by noise");
xlabel('Iteration');
ylabel('Accuracy');
legend("A - UΣV^T");
saveas(gca,'noisy4.jpg');

% 
% figure;
% hold on;
% title('North acceleration');
% plot(time, -diff8);
% legend('A priori - A posteriori');
% xlabel('Time [s]');
% ylabel('Acceleration [m/s^2]');
% % ylim([-200 20])
% saveas(gca,'North acceleration Predicted.jpg');
% 
% % Create a new figure and subplot
% figure('Position', [100 100 500 600]);
% title('Timing of SVD')
% subplot(3, 1, 1);
% hold on;
% set(gcf, 'Color', [1 1 1]); % white background
% set(gca, 'FontName', 'Arial'); % beautify the axes a bit
% title('Adaptive SVD');
% plot(T);
% xlabel('Iteration');
% ylabel('Time [s]');
% legend('Seconds');
% 
% % Create the second subplot and adjust its height
% subplot(3, 1, 2);
% hold on;
% title('Matlab SVD');
% plot(T2);
% xlabel('Iteration');
% ylabel('Time [s]');
% legend('Seconds');
% 
% subplot(3,1,3);
% hold on;
% title('Difference');
% plot(T-T2);
% xlabel('Iteration');
% ylabel('Time [s]');
% legend('Seconds');
% 
% saveas(gca,'times500_4.jpg');
% 
% 
