%% Implentation of DMD for 2D Drone
close all;

% Load simulation data
simulation_data_file = 'No_payload_data_6';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Resample time series to desired sample time
Ts = 0.02;     % Desired sample time of data
x_resamp = resample(out.x,   0:Ts:out.x.Time(end));  
u_resamp = resample(out.F_r, 0:Ts:out.x.Time(end));  

% Extract data
u_data  = u_resamp.Data';
x_data  = x_resamp.Data';
y_data  = x_data(1:end,:); % Measurement data
t       = x_resamp.Time'; % Time

u_data_org = u_data;
% Adjust for constant disturbance / mean control values
u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; -4.5*9.81];
u_data  = u_data - u_bar; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 1000; % Num of data samples for testing
x_test = x_data(:,end-N_test+1:end);
y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Data dimentions
nx = size(x_data,1); % number of states
ny = size(y_data,1); % number of measurements
nu = size(u_data,1); % number of inputs
N  = length(t);     % Number of data samples

% Add noise
rng('default');
rng(1); % Repeatable random numbers
sigma = 0.001; % Noise standard deviation
y_data_noise = y_data + sigma*randn(size(y_data));

% Training data - Last sample of training is first sample of testing
N_train = 30/Ts % Number of sampels in training data
y_train = y_data_noise(:, end-N_test-N_train+2:end-N_test+1); % Use noisy data
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
t_train = t(:,end-N_test-N_train+2:end-N_test+1);

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

q = 10; % Override

w = N_train - q + 1; % num columns of Hankel matrix

if q == 1
    Upsilon = u_train(:, q:end);
else
    % Hankel matrix with delay measurements
    Upsilon = zeros((q-1)*ny,w); % Augmented state with delay coordinates [...; Y(k-2); Y(k-1); Y(k)]
    for row = 0:q-2 % Add delay coordinates
        Upsilon((end - ny*(row+1) + 1):(end - ny*row), :) = y_train(:, row + (1:w));
    end

    % Matrix with time series of states
    Y = y_train(:, q-1 + (1:w));

    Upsilon = [Upsilon; u_train(:, q:end)]; % Leave out last time step to match V_til_1
end

YU_bar = [Y; Upsilon];

% DMD of Y
Y2 = Y(:, 2:end  );
Y1 = Y(:, 1:end-1);

YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y and U, above and below
AB = Y2*pinv(YU); % combined A and B matrix, side by side

% System matrixes from DMD
A  = AB(:,1:ny); % Extract A matrix
B  = AB(:,(ny+1):end);

% A = stabilise(A,10);

%% Run with A and x

k_start = 25/Ts; % k from where to start applying model
N_test = 10/Ts; % Number of data points to run and test for

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
    if q ~= 1
        y_delays = [y_hat(:,k); y_delays(1:(end-ny),:)]; % Add y(k) to y_delay for next step [y(k); y(k-1); ...]
    end
end

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat - y_run), 2)./N_run % For each measured state

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