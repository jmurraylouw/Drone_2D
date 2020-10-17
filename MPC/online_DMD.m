%% Implentation of DMD for 2D Drone
% Simulate calculating a new DMD model at every timestep

close all;

simulation_data_file = 'No_payload_data_6';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Extract data
u_data  = out.F_r.Data';
x_data  = out.x.Data';
y_data  = x_data([1,2,3],:); % Measurement data (x, z, theta)
t       = out.x.Time'; % Time

% Adjust for constant disturbance / mean control values
u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; -4.5*9.81];
u_data  = u_data - u_bar; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 1000; % Num of data samples for testing

% Data dimentions
nx = size(x_data,1); % number of states
ny = size(y_data,1); % number of measurements
nu = size(u_data,1); % number of inputs
Ts = t(2)-t(1)     % Sample time of data
N  = length(t);     % Number of data samples

% Add noise
rng('default');
rng(1); % Repeatable random numbers
sigma = 0.001; % Noise standard deviation
y_data_noise = y_data + sigma*randn(size(y_data));

% Read previous results
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
results_file = ['Data/dmd_results_', simulation_data_file, '_sig=', sig_str, '.mat'];

try
    load(results_file);
    results(~results.q,:) = []; % remove empty rows
catch
    disp('No saved results file')  
end

% Parameters
% best_row = find(results.MAE_mean == min(results.MAE_mean));
% best_results = results(best_row,:)
% q = double(best_results.q);
% p = double(best_results.p);

q = 20; % Override

w = N_train - q + 1; % num columns of Hankel matrix

% Training data with only initial conditions
y_train = y_data(:,1) + zeros(ny,N_train); % Assume all data before simulation is at init conditions
u_train = u_bar + zeros(nu,N_train);

% Hankel matrix with delay measurements
Delays = zeros((q-1)*ny,w); % Augmented state with delay coordinates [...; Y(k-2); Y(k-1); Y(k)]
for row = 0:q-2 % Add delay coordinates
    Delays((end - ny*(row+1) + 1):(end - ny*row), :) = y_train(:, row + (1:w));
end
Upsilon = [Delays; u_train(:, q:end)]; % Leave out last time step to match V_til_1

% Matrix with time series of states
Y = y_train(:, q:end);

for k = 1:N
    
    % Current timestep data
    y_k = y_data_noise(:,k); % Sample from noisy data
    u_k = u_data(:,k);
    
    % Update Upsilon
    new_row = [Y(:,end); Upsilon(1:(end - ny - nu), end); u_k];
    Upsilon = [Upsilon(:,2:end), new_row]; % Add new row of updated values
    
    % Update Y
    Y = [Y(:,2:end), y_k]; % Forget oldest entry, add y_k
    
    % DMD of Y
    Y2 = Y(:, 2:end  ); 
    Y1 = Y(:, 1:end-1); % One time-step behin of Y2

    YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y and U, above and below
    AB = Y2*pinv(YU); % combined A and B matrix, side by side

    % System matrixes from DMD
    A  = AB(:,1:ny); % Extract A matrix
    B  = AB(:,(ny+1):end);

    % Adjust slightly unstable eigenvalues
    A = stabilise(A,3);

    %% Run with A and x

    k_start = 5000; % k from where to start applying model
    N_test = 1000; % Number of data points to run and test for

    % Start at end of initial condition k
    y_run = y_data(:, k_start + (1:N_test));
    u_run = u_data(:, k_start + (1:N_test));
    t_run = t(:, k_start + (1:N_test));
    N_run = length(y_run);

    % Initial condition
    y_hat_0 = y_data(:,k_start);

    % Initial delay coordinates
    y_delays = zeros((q-1)*ny,1);
    k = k_start; % index of y_data
    for i = 1:ny:ny*(q-1) % index of y_delays
        k = k - 1; % previosu index of y_data
        y_delays(i:(i+ny-1)) = y_data(:,k);
    end

    % Run model
    y_hat = zeros(ny,N_run); % Empty estimated Y
    y_hat(:,1) = y_hat_0; % Initial condition
    for k = 1:N_run-1
        upsilon = [y_delays; u_run(:,k)]; % Concat delays and control for use with B
        y_hat(:,k+1) = A*y_hat(:,k) + B*upsilon;
        y_delays = [y_hat(:,k); y_delays(1:(end-ny),:)]; % Add y(k) to y_delay for next step [y(k); y(k-1); ...]
    end

    % Vector of Mean Absolute Error on testing data
    MAE = sum(abs(y_hat - y_run), 2)./N_run % For each measured state

end % run through time steps
%% Plot data vs model
figure;
plot(t_run, y_run);
hold on;

plot(t_run, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
title('Training and Testing data vs Model');
legend()
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