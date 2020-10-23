% Load simulation data
simulation_data_file = 'No_payload_data_5';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Desired sampling time
Ts = Ts_mpc;
N_test = 100;
pause_and_plot = 1;
plot_results = 1;
model_type = 'delay_A';
% A = A_dmd;
% B = B_dmd;

% Resample time series to desired sample time
x_resamp = resample(out.x, 0:Ts:out.x.Time(end));  
u_resamp = resample(out.F_r, 0:Ts:out.x.Time(end));  

% Extract data
u_data  = u_resamp.Data';
x_data  = x_resamp.Data';
y_data  = x_data(1:end,:); % Measurement data (x, z, theta)
t       = x_resamp.Time'; % Time

% Adjust for constant disturbance / mean control values
u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; -4.5*9.81];
u_data  = u_data - u_bar; % Adjust for unmeasured input

nu = size(u_data,1);


N = length(t); % Total number of data samples
ny = size(y_data,1);

% Create empty results table
VariableTypes = {'int16', 'double'}; % id, q, MAE_mean
VariableNames = {'k',     'MAE_mean'};
for i = 1:ny % Mae column for each measured state
    VariableNames = [VariableNames, strcat('MAE_', num2str(i))];
    VariableTypes = [VariableTypes, 'double'];
end
Size = [N, length(VariableTypes)];

results = table('Size',Size,'VariableTypes',VariableTypes,'VariableNames',VariableNames);
emptry_row = 1; % Keep track of next empty row to insert results 

for k = 1:length(t)
    % Data to test with, starting at sample k           
    y_run = y_data(:, k + (0:N_test-1));
    u_run = u_data(:, k + (0:N_test-1));
    t_run =      t(:, k + (0:N_test-1));

    % Initial condition
    y_hat_0 = y_data(:,k);

    % Initial delay coordinates
    y_delays = zeros((q-1)*ny, 1);
    j = k; % index of y_data
    for i = 1:ny:ny*(q-1) % index of y_delays
        j = j - 1; % previos index of y_data
        if j < 1 % Assume for time steps before k=1, y = initial condition
            j = 1;
        end        
        y_delays(i:(i+ny-1)) = y_data(:,j); % Insert into delay vector
    end

    % Run model predictions
    switch model_type
        case 'delay_A'
            y_hat = zeros(q*ny, N_test); % Empty estimated Y
            y_hat(:,1) = [y_hat_0; y_delays]; % Initial condition
            for j = 1:N_test-1
                y_hat(:,j+1) = A*y_hat(:,j) + B*u_run(:,j);
            end
            y_hat = y_hat(1:ny, :); % Only keep non-delay rows
            
        case 'delay_B'
            y_hat = zeros(ny,N_test); % Empty estimated Y
            y_hat(:,1) = y_hat_0; % Initial condition
            for j = 1:N_test-1
                upsilon = [y_delays; u_run(:,j)]; % Concat delays and control for use with B
                y_hat(:,j+1) = A*y_hat(:,j) + B*upsilon;
                if q ~= 1 % Special case if no delays
                    y_delays = [y_hat(:,j); y_delays(1:(end-ny),:)]; % Add y(k) to y_delay for next step [y(k); y(k-1); ...]
                end
            end
            
        otherwise
            error("Use either 'delay_A' or 'delay_B' ")
            
    end % switch

    % Vector of Mean Absolute Error on testing data
    MAE = sum(abs(y_hat - y_run), 2)./N_test; % For each measured state

    % Save results
    results(emptry_row,:) = [{k, mean(MAE)}, num2cell(MAE')]; % add to table of results
    emptry_row = emptry_row + 1; 

    % Plot and pause at every timestep
    if plot_and_pause
        if mod(k, 100) == 0
            figure(1)
            hold on;
            plot(t_run, y_run)
            threshold = 2*max(max(abs(y_data))); % for if y_hat is unstable
            y_hat(y_hat > threshold) = threshold;
            y_hat(y_hat < -threshold) = -threshold;
            plot(t_run, y_hat, '--', 'LineWidth', 1)
            hold off;
            disp('Pausing... Press enter to continue')
            pause
            clf('reset')
            disp('Continuing...')
        end
    end
    
end % k

if plot_results
    figure
    semilogy(t(results.k), results.MAE_mean)
    title("Model error over time")
    
    figure
    plot(t(results.k), results.MAE_mean*100)
    hold on
    title("Model error vs states")
    plot(t, y_data)
    legend('mean MAE scaled', 'x', 'z', 'theta')
    hold off
end


