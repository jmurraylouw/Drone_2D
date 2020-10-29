function [A, B] = model_DMD(out, start_time, end_time, Ts, q, y_rows, sigma, plot_prediction, varargin)
%% Generate a DMD model from simulation data with specific Ts and q
% A = discrete state-space system matrix
% B = input matrix

% out        = simulation data from simulink To Workspace blocks
% Ts         = desireed sample time of model
% q          = delay coordinates present in model, including current timestep. 
%   q = 1 x(k) used in model
%   q = 2 x(k), x(k-1) used in model, etc
% y_rows     = vector of rows to measure from state, x, to produce y
% start_time = time when training starts
% end_time   = time when training ends
% sigma      = standard deviation of noise to add to data
% varargin   = {u_bar} if you want to normalise input by specific u_bar
%   else leave out param to use mean of input data as normalising u_bar

% Resample time series to desired sample time and training period
x_resamp = resample(out.x, start_time:Ts:end_time);  
u_resamp = resample(out.u, start_time:Ts:end_time);  

% Extract data
u_train  = u_resamp.Data';
x_train  = x_resamp.Data';
y_train  = x_train(y_rows,:); % Measurement data
N_train  = size(y_train,2); % Number of sampels in training data
t        = x_resamp.Time;

% Add noise
rng('default');
rng(1); % Repeatable random numbers
y_train = y_train + sigma*randn(size(y_train));

% Normalise input data
if isempty(varargin)
    u_bar = mean(u_train,2); % Input needed to keep at a fixed point
elseif length(varargin) == 1
    u_bar = varargin{1};
else
    error("Too many input arguments")
end

u_train  = u_train - u_bar; % Adjust for unmeasured input

% Data dimentions
ny = size(y_train,1); % number of measurements
w = N_train - q + 1; % num columns of Hankel matrix

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
AB = Y2*pinv(YU); % combined A and B matrix, side by side

% System matrixes from DMD
A  = AB(:,1:ny); % Extract A matrix
B  = AB(:,(ny+1):end);
% A = stabilise(A,1);

%% Run with A and x
if plot_prediction
    k_start = q+1; % k from where to start applying model
    N_test = N_train-k_start-1; % Number of data points to run and test for

    % Start at end of initial condition k
    y_run = y_train(:, k_start + (1:N_test));
    u_run = u_train(:, k_start + (1:N_test));
    t_run = t(k_start + (1:N_test));
    N_run = length(y_run);

    % Initial condition
    y_hat_0 = y_run(:,1);

    % Initial delay coordinates
    y_delays = zeros((q-1)*ny,1);
    k = k_start; % index of y_data
    for i = 1:ny:ny*(q-1) % index of y_delays
        k = k - 1; % previosu index of y_data
        y_delays(i:(i+ny-1)) = y_train(:,k);
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
    title('Training data vs Model');
    legend()
    hold off;

end

function A = stabilise(A_unstable,max_iterations)
    % If some eigenvalues are unstable due to machine tolerance,
    % Scale them to be stable
    A = A_unstable;
    count = 0;
    while (sum(abs(eig(A)) > 1) ~= 0)       
        [Ve,De] = eig(A);
        unstable = abs(De)>1; % indexes of unstable eigenvalues
        De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-16 + count*2); % Normalize all unstable eigenvalues (set abs(eig) = 1)
        A = Ve*De/(Ve); % New A with margininally stable eigenvalues
        A = real(A);
        count = count+1;
        if(count > max_iterations)
            break
        end
    end

end

end