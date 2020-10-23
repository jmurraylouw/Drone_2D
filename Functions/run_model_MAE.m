% Load simulation data
simulation_data_file = 'No_payload_data_6';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Desired sampling time
Ts = Ts_mpc;

% Resample time series to needed sample time
x_resamp = resample(out.x, 0:Ts:out.x.Time(end));  
u_resamp = resample(out.F_r, 0:Ts:out.x.Time(end));  

% Extract data
u_data  = u_resamp.Data';
x_data  = x_resamp.Data';
y_data  = x_data(1:3,:); % Measurement data (x, z, theta)
t       = x_resamp.Time'; % Time

% Adjust for constant disturbance / mean control values
u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; -4.5*9.81];
u_data  = u_data - u_bar; % Adjust for unmeasured input

nu = size(u_data,1);

N_test = 50;
pause_and_plot = 1;
plot_results = 1;

% NB: Ensure u_data and System matrices use same sample time
% results = model_MAE_accross_data(y_data, u_data, t, A_mpc, B_mpc(:, 1:nu), q, N_test, 'delay_A', pause_and_plot, plot_results);
results = model_MAE_accross_data(y_data, u_data, t, A_resamp, B_resamp, q, N_test, 'delay_B', pause_and_plot, plot_results);

% ???? Model and data is with different time scale


