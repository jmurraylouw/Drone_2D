%% Implentation of Hankel Alternative View Of Koopman for 2D Drone
% close all;

% simulation_data_file = 'No_payload_data_5';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Extract data
u_data  = out.F_r.Data';
x_data  = out.x.Data';
y_data  = x_data([1,2,3],:); % Measurement data (x, z, theta)
t       = out.tout'; % Time

% Adjust for constant disturbance / mean control values
% u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; 4.5*9.81];
u_data  = u_data - u_bar; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
% N_test = 2000; % Num of data samples for testing
x_test = x_data(:,end-N_test+1:end);
y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Data dimentions
n = size(x_data,1); % number of states
m = size(y_data,1); % number of measurements
l = size(u_data,1); % number of inputs
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

% Add noise
rng('default');
rng(1); % Repeatable random numbers
% sigma = 0.001; % Noise standard deviation
y_data_noise = y_data + sigma*randn(size(y_data));

% Training data - Last sample of training is first sample of testing
N_train = 8000; % Number of sampels in training data
y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
t_train = t(:,end-N_test-N_train+2:end-N_test+1);

% Read previous results
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
results_file = ['Data/havok_results_', simulation_data_file, '_sig=', sig_str, '.mat'];

try
    load(results_file);
    results(~results.q,:) = []; % remove empty rows
catch
    disp('No saved results file')  
end

% Parameters
best_row = find(results.MAE_mean == min(results.MAE_mean));
best_results = results(best_row,:)
q = double(best_results.q);
p = double(best_results.p);

% q = 300; % Override
% p = 120; % Override

w = N_train - q + 1; % num columns of Hankel matrix
D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)

% Create Hankel matrix with measurements
% Y = zeros(q*m,w); % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
% for row = 0:q-1 % Add delay coordinates
%     Y(row*m+1:(row+1)*m, :) = y_train(:, row + (0:w-1) + 1);
% end

Y = zeros(q*m,w); % Augmented state with delay coordinates [..., Y(k-2), Y(k-1), Y(k)]
for row = 0:q-1 % Add delay coordinates
    Y((end - m*(row+1) + 1):(end - m*row), :) = y_train(:, row + (1:w));
end

Upsilon = u_train(:, q:end); % Leave out last time step to match V_til_1
YU_bar = [Y; Upsilon];

% SVD of the Hankel matrix
[U1,S1,V1] = svd(YU_bar, 'econ');
figure, semilogy(diag(S1), 'x'), hold on;
title('Singular values of Omega, showing p truncation')
plot(p,S1(p,p), 'ro'), hold off;

% Truncate SVD matrixes
U_tilde = U1(:, 1:p); 
S_tilde = S1(1:p, 1:p);
V_tilde = V1(:, 1:p);

% Setup V2 one timestep into future from V1
V_til_2 = V_tilde(2:end  , :)'; % Turnd on side (wide short matrix)
V_til_1 = V_tilde(1:end-1, :)';

% DMD on V
AB_tilde = V_til_2*pinv(V_til_1); % combined A and B matrix, side by side
AB_tilde = stabilise(AB_tilde,3);

% convert to x coordinates
AB_bar = (U_tilde*S_tilde)*AB_tilde*pinv(U_tilde*S_tilde);
A_bar = AB_bar(1:q*m, 1:q*m);
B_bar = AB_bar(1:q*m, q*m+1:end);
% A_bar = stabilise(A_bar,10);

% DMD of Y
Y2 = Y(:, 2:end  );
Y1 = Y(:, 1:end-1);

YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y and U, above and below
AB = Y2*pinv(YU); % combined A and B matrix, side by side
A  = AB(:,1:q*m); % Extract A matrix
B  = AB(:,(q*m+1):end);

% A = stabilise(A,10);

% Compare to testing data

%% Run with A_bar, B_bar and x
figure;
plot(U1(:,1:5))
title('First 5 modes of SVD')

% Initial condition (last entries of training data)
% y_hat_0 = zeros(q*m,1);
% for row = 0:q-1 % First column of spaced Hankel matrix
%     y_hat_0(row*m+1:(row+1)*m, 1) = y_train(:, end - ((q-1)+1) + row + 1);
% end

y_hat_0 = zeros(q*m,1);
for row = 0:q-1 % Add delay coordinates
    y_hat_0((end - m*(row+1) + 1):(end - m*row), 1) = y_train(:, end - ((q-1)+1) + row + 1);
end

% Run model
Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
Y_hat(:,1) = y_hat_0; % Initial condition
for k = 1:N_test-1
    Y_hat(:,k+1) = A_bar*Y_hat(:,k) + B_bar*u_test(:,k);
end

% y_hat_bar = Y_hat(end-m+1:end, :); % Extract only non-delay time series (last m rows)
y_hat_bar = Y_hat(1:m, :); % Extract only non-delay time series (first m rows)

% Vector of Mean Absolute Error on testing data
MAE_bar = sum(abs(y_hat_bar - y_test), 2)./N_test % For each measured state


%% Run with A and x

% Initial condition
% y_hat_0 = zeros(q*m,1);
% for row = 0:q-1 % First column of spaced Hankel matrix
%     y_hat_0(row*m+1:(row+1)*m, 1) = y_train(:, end - ((q-1)+1) + row + 1);
% end

y_hat_0 = zeros(q*m,1);
for row = 0:q-1 % Add delay coordinates
    y_hat_0((end - m*(row+1) + 1):(end - m*row), 1) = y_train(:, end - ((q-1)+1) + row + 1);
end
            
% Run model
Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
Y_hat(:,1) = y_hat_0; % Initial condition
for k = 1:N_test-1
    Y_hat(:,k+1) = A*Y_hat(:,k) + B*u_test(:,k);
end

% y_hat = Y_hat(end-m+1:end, :); % Extract only non-delay time series (last m rows)
y_hat = Y_hat(1:m, :); % Extract only non-delay time series (first m rows)

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state

% Compare MAE and MAE_til
MAE_error_percent = (MAE - MAE_bar)./MAE_bar*100

%% Plot data vs model
figure;
plot(t_train, y_train);
hold on;
plot(t_test, y_test);

plot(t_test, y_hat, 'k--', 'LineWidth', 1); % Plot only non-delay coordinate
plot(t_test, y_hat_bar, 'r--', 'LineWidth', 1); % Plot only non-delay coordinate
plot((D + t(N-N_test-N_train)).*[1, 1], ylim, 'r');
plot(t(N-N_test-N_train).*[1,1], ylim, 'k');
plot(t(N-N_test).*[1,1], ylim, 'k');
title('Training and Testing data vs Model');
% legend('','','','','','', 'x hat','z hat','theta hat', 'x hat bar','z hat bar','theta hat bar');
hold off;

function A = stabilise(A_unstable,max_iterations)
    % If some eigenvalues are unstable due to machine tolerance,
    % Scale them to be stable
    A = A_unstable;
    count = 0;
    while (sum(abs(eig(A)) > 1) ~= 0)       
        [Ve,De] = eig(A);
        unstable = abs(De)>1; % indexes of unstable eigenvalues
        De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-14 + count*2); % Normalize all unstable eigenvalues (set abs(eig) = 1)
        A = Ve*De/(Ve); % New A with margininally stable eigenvalues
        A = real(A);
        count = count+1;
        if(count > max_iterations)
            break
        end
    end

end