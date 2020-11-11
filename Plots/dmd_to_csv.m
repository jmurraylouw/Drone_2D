%% Save data
% description = 'DMD with diff_params. step x = 1 ';
% m_sim = m; % m used in simulation
% l_sim = l; % l used in simulation
% cbeta_sim = cbeta;
% save('Plots/dmd_diff_params.mat', 'description', 'out', 'm_sim', 'l_sim', 'cbeta_sim')

% Write data to csv files
Ts_csv = 0.1;

%% DMD
load('Plots/dmd_diff_params.mat')

t = 0:Ts_csv:39;

X_pid = resample(out.x_pid, t);% Resample time series to desired sample time and period  
X_pid = X_pid.Data;

T_pid = resample(out.T_pid, t); % Thrust of motors
T_pid = T_pid.Data;

X_dmd = resample(out.x_mpc, t);
X_dmd = X_dmd.Data;

T_dmd = resample(out.T_mpc, t); % Thrust of motors
T_dmd = T_dmd.Data;

% Convert radians to degress
angle_columns = [3,4,7,8];
X_pid(:, angle_columns) = (180/pi).*X_pid(:, angle_columns);
X_dmd(:, angle_columns) = (180/pi).*X_dmd(:, angle_columns);


%% CSV
num_samples = length(t);
num_vars = 1 + (size(X_pid,2) + size(T_pid,2))*2;

table_data = array2table([t', X_pid, T_pid, X_dmd, T_dmd]);

table_data.Properties.VariableNames = {'t',  'x_pid', 'z_pid', 'th_pid', 'b_pid', 'dx_pid', 'dz_pid', 'dth_pid', 'db_pid', 'T1_pid', 'T2_pid','x_dmd', 'z_dmd', 'th_dmd', 'b_dmd', 'dx_dmd', 'dz_dmd', 'dth_dmd', 'db_dmd', 'T1_dmd', 'T2_dmd'};
%     , 'x_dmd', 'z_dmd', 'th_dmd', 'b_dmd', 'dx_dmd', 'dz_dmd', 'dth_dmd', 'db_dmd', 'T1_dmd', 'T2_dmd'};
writetable(table_data,'Plots/dmd_diff_params.csv')

