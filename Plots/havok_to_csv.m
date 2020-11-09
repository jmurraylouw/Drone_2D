%% Save data
description = 'HAVOK with different params. step x = 1 ';
m_sim = m; % m used in simulation
l_sim = l; % l used in simulation
cbeta_sim = cbeta;
save('Plots/havok_diff_params.mat', 'description', 'out', 'm_sim', 'l_sim', 'cbeta_sim')

% Write data to csv files
Ts_csv = 0.01;

%% HAVOK
% load('Plots/havok_noise.mat');

t = 0:Ts_csv:39;

X_pid = resample(out.x_pid, t);% Resample time series to desired sample time and period  
X_pid = X_pid.Data;

T_pid = resample(out.T_pid, t); % Thrust of motors
T_pid = T_pid.Data;

X_havok = resample(out.x_mpc, t);
X_havok = X_havok.Data;

T_havok = resample(out.T_mpc, t); % Thrust of motors
T_havok = T_havok.Data;

% Convert radians to degress
angle_columns = [3,4,7,8];
X_pid(:, angle_columns) = (180/pi).*X_pid(:, angle_columns);
X_havok(:, angle_columns) = (180/pi).*X_havok(:, angle_columns);


%% CSV
num_samples = length(t);
num_vars = 1 + (size(X_pid,2) + size(T_pid,2))*2;

table_data = array2table([t', X_pid, T_pid, X_havok, T_havok]);

table_data.Properties.VariableNames = {'t',  'x_pid', 'z_pid', 'th_pid', 'b_pid', 'dx_pid', 'dz_pid', 'dth_pid', 'db_pid', 'T1_pid', 'T2_pid','x_havok', 'z_havok', 'th_havok', 'b_havok', 'dx_havok', 'dz_havok', 'dth_havok', 'db_havok', 'T1_havok', 'T2_havok'};
%     , 'x_dmd', 'z_dmd', 'th_dmd', 'b_dmd', 'dx_dmd', 'dz_dmd', 'dth_dmd', 'db_dmd', 'T1_dmd', 'T2_dmd'};
writetable(table_data,'Plots/havok_diff_params.csv')

