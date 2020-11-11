%% Implentation of DMD
% close all;

% % Extract data
% simulation_data_file = 'With_payload_and_noise_data_1';
% load(['Data/', simulation_data_file, '.mat']) % Load simulation data
% 
% Ts = 0.03;     % Desired sample time
% Ts_dmd = Ts;
% y_rows = 1:4;
% 
% % Adjust for constant disturbance / mean control values
% % u_bar = mean(out.u.Data,1); % Input needed to keep at a fixed point
% u_bar = [0, (m + M)*g];
% out.u.Data  = out.u.Data - u_bar; % Adjust for unmeasured input

% Training data
% % train_time = 0:Ts:200;
% x_train = resample(out.x, train_time );% Resample time series to desired sample time and training period  
% u_train = resample(out.u, train_time );  
% t_train = x_train.Time';
% N_train = length(t_train);
% 
% x_train = x_train.Data';
% y_train = x_train(y_rows,:);
% u_train = u_train.Data';
% 
% % Testing data
% % test_time = 400:Ts:500;
% x_test = resample(out.x, test_time );  
% u_test = resample(out.u, test_time );  
% t_test = x_test.Time';
% N_test = length(t_test); % Num of data samples for testing
% 
% x_test = x_test.Data';
% y_test = x_test(y_rows,:); % One sample of testing data overlaps for initial condition
% u_test = u_test.Data';
% 
% % Data dimentions
% nx = size(x_train,1); % number of states
% ny = size(y_train,1); % number of measurements
% nu = size(u_train,1); % number of inputs
% 
% % % Add noise
% rng('default');
% rng(1); % Repeatable random numbers
% sigma = 0; % Noise standard deviation
% y_data_noise = y_data + sigma*randn(size(y_data));

% comment = ''; % Extra comment to differentiate this run
% 
% % Read previous results
% sigma = 0;
% sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
% results_file = ['Data/havok_results_', comment, simulation_data_file, '_sig=', sig_str, '.mat'];

try
    load(results_file);
    results(~results.q,:) = []; % remove empty rows
    
    % Parameters
    best_row = find(results.MAE_mean == min(results.MAE_mean));
    best_results = results(best_row,:);
    q = double(best_results.q);
    p = double(best_results.p);
    
    only_q_Ts = 0; % Try best result for specific q
    if only_q_Ts
        q = 20;
        q_results = results((results.q == q & results.Ts == Ts),:);
        best_row = find(q_results.MAE_mean == min(q_results.MAE_mean));
        best_results = q_results(best_row,:)
        p = double(best_results.p);
    end
    
catch
    disp('No saved results file')  
end

% % Override parameters:
% q = 80
% p = 40

q
p

w = N_train - q + 1; % num columns of Hankel matrix
D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)

% Hankel matrix with delay measurements
if q == 1 % Special case if no delay coordinates
    Upsilon = u_train(:, q:end);
else
    Upsilon = zeros((q-1)*ny,w); % Augmented state Y[k] at top
    for row = 0:q-2 % Add delay coordinates
        Upsilon((end - ny*(row+1) + 1):(end - ny*row), :) = y_train(:, row + (1:w));
    end

    Upsilon = [Upsilon; u_train(:, q:end)]; % Leave out last time step to match V_til_1
end

% Matrix with time series of states
Y = y_train(:, q-1 + (1:w));

% DMD of Y
Y2 = Y(:, 2:end  );
Y1 = Y(:, 1:end-1);

YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y above and U below

% SVD of the Hankel matrix
[U1,S1,V1] = svd(YU, 'econ');

% Truncate SVD matrixes
U_tilde = U1(:, 1:p); 
S_tilde = S1(1:p, 1:p);
V_tilde = V1(:, 1:p);

% YU = \approx U_tilde*S_tilde*V_tilde'
AB = Y2*pinv(U_tilde*S_tilde*V_tilde'); % combined A and B matrix, side by side
% AB = Y2*pinv(YU); % combined A and B matrix, side by side

% System matrixes from DMD
A_dmd  = AB(:,1:ny); % Extract A matrix
B_dmd  = AB(:,(ny+1):end);
% A = stabilise(A,1);

% Initial condition
y_hat_0 = y_test(:,q);

% Initial delay coordinates
y_delays = zeros((q-1)*ny,1);
k = q; % index of y_data
for i = 1:ny:ny*(q-1) % index of y_delays
    k = k - 1; % previosu index of y_data
    y_delays(i:(i+ny-1)) = y_test(:,k);
end

% Run model
y_hat = zeros(ny,N_test); % Empty estimated Y
y_hat(:,1) = y_hat_0; % Initial condition
for k = 1:N_test-1
    upsilon = [y_delays; u_test(:,k)]; % Concat delays and control for use with B
    y_hat(:,k+1) = A_dmd*y_hat(:,k) + B_dmd*upsilon;
    if q ~= 1
        y_delays = [y_hat(:,k); y_delays(1:(end-ny),:)]; % Add y(k) to y_delay for next step [y(k); y(k-1); ...]
    end
end

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state

%% Plot data vs model
figure;
plot(t_train, y_train);
hold on;
plot(t_test, y_test);

% plot(t_test, y_hat, 'k--', 'LineWidth', 1); % Plot only non-delay coordinate
plot(t_test, y_hat, 'r--', 'LineWidth', 1); % Plot only non-delay coordinate
title('Training and Testing data vs Model (red = HAVOK, black = DMD)');
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