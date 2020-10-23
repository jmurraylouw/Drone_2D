% Load simulation data
simulation_data_file = 'No_payload_data_5';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Extract data
u_data  = out.F_r.Data';
x_data  = out.x.Data';
y_data  = x_data(1:3,:); % Measurement data (x, z, theta)
t       = out.x.Time'; % Time

% Adjust for constant disturbance / mean control values
u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; -4.5*9.81];
u_data  = u_data - u_bar; % Adjust for unmeasured input


N_test = 50;
pause_and_plot = 1;
plot_results = 1;

results = model_MAE_accross_data(y_data, u_data, t, A, B, N_test, 'delays_A', pause_and_plot, plot_results);




