%% Implentation of DMD for 2D Drone
% Simulate calculating a new DMD model at every timestep

% close all;

% Load simulation data
simulation_data_file = 'No_payload_data_5';
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

% Data dimentions
nx = size(x_data,1); % number of states
ny = size(y_data,1); % number of measurements
nu = size(u_data,1); % number of inputs
Ts = t(2)-t(1)     % Sample time of data
N  = length(t);     % Number of data samples

% Parameters
N_train = 30/Ts; % Num of data samples for training
N_test = 20; % Num of data samples for testing
sigma = 0.001; % Noise standard deviation
q = 10; % Override
model_intervals = 10; % Only dod DMD every so many time-steps

% Add noise
rng('default');
rng(1); % Repeatable random numbers
y_data_noise = y_data + sigma*randn(size(y_data));

% Training data with only initial conditions
y_train = y_data(:,1) + zeros(ny,N_train); % Assume all data before simulation is at init conditions
u_train = u_bar + zeros(nu,N_train);

% Hankel matrix with delay measurements
w = N_train - q + 1; % num columns of Hankel matrix
Delays = zeros((q-1)*ny,w); % Augmented state with delay coordinates [...; Y(k-2); Y(k-1); Y(k)]
for row = 0:q-2 % Add delay coordinates
    Delays((end - ny*(row+1) + 1):(end - ny*row), :) = y_train(:, row + (1:w));
end
Upsilon = [Delays; u_train(:, q:end)]; % Leave out last time step to match V_til_1

% Matrix with time series of states
Y = y_train(:, q:end);

% Create empty results table
VariableTypes = {'int16', 'int16', 'double'}; % id, q, MAE_mean
VariableNames = {'k',     'q',     'MAE_mean'};
for i = 1:ny % Mae column for each measured state
    VariableNames = [VariableNames, strcat('MAE_', num2str(i))];
    VariableTypes = [VariableTypes, 'double'];
end
Size = [N, length(VariableTypes)];

results = table('Size',Size,'VariableTypes',VariableTypes,'VariableNames',VariableNames);
emptry_row = 1; % Keep track of next empty row to insert results 
    
% for k = N_test:N - N_test
for k = N_test:20/Ts
    k
    tic;
    
    % Current timestep data
    y = y_data_noise(:,k); % Sample from noisy data
    u = u_data(:,k);
    
    % Update Upsilon
    new_row = [Y(:,end); Upsilon(1:(end - ny - nu), end); u];
    Upsilon = [Upsilon(:,2:end), new_row]; % Add new row of updated values
    
    % Update Y
    Y = [Y(:,2:end), y]; % Forget oldest entry, add y_k
    
    if (mod(k, model_intervals) == 0)
        
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


        if (k > q) % Need initial conditions for delays        
        %% Run with A and x
            % Start at end of initial condition k

            y_run = y_data(:, k + (1:N_test));
            u_run = u_data(:, k + (1:N_test));
            t_run = t(:, k + (1:N_test));
            N_run = length(y_run);

            % Initial condition
            y_hat_0 = y_data(:,k);

            % Initial delay coordinates
            y_delays = zeros((q-1)*ny,1);
            j = k; % index of y_data
            for i = 1:ny:ny*(q-1) % index of y_delays
                j = j - 1; % previosu index of y_data
                y_delays(i:(i+ny-1)) = y_data(:,j);
            end

            % Run model
            y_hat = zeros(ny,N_run); % Empty estimated Y
            y_hat(:,1) = y_hat_0; % Initial condition
            for j = 1:N_run-1
                upsilon = [y_delays; u_run(:,j)]; % Concat delays and control for use with B
                y_hat(:,j+1) = A*y_hat(:,j) + B*upsilon;
                y_delays = [y_hat(:,j); y_delays(1:(end-ny),:)]; % Add y(k) to y_delay for next step [y(k); y(k-1); ...]
            end

            % Vector of Mean Absolute Error on testing data
            MAE = sum(abs(y_hat - y_run), 2)./N_run; % For each measured state

            % Save results
            results(emptry_row,:) = [{k, q, mean(MAE)}, num2cell(MAE')]; % add to table of results
            emptry_row = emptry_row + 1; 

            % Plot and pause
            plot_and_pause = 0;
            if plot_and_pause
                clf('reset')
    %             plot(t_run, y_run)
                hold on;
                plot(t_run, y_hat, '--')
                hold off;
                disp('Pausing... Press enter to continue')
                pause
                disp('Continuing...')
            end

        end % Run model
    end % Every k_intervals
    toc; 
end % run through time steps


%% Plot spread of results
results(~results.q,:) = []; % remove empty rows

plot_results = 1;
if plot_results
    figure
    semilogy(t(results.k), results.MAE_mean)
    title("Model error over time")
    
    figure
    plot(t(results.k), results.MAE_mean*100)
    hold on
    title("Model error vs states")
    plot(t, y_data)
    legend('MAE_mean', 'x', 'z', 'theta')
    hold off
    
end

%% MAE over all time steps
valid_rows = find(results.k < N_train); % Only rows where full N_train sample were available for model
valid_results = results(valid_rows,:);
prev_MAE = total_MAE_mean
total_MAE_mean = mean(valid_results.MAE_mean)

%%
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