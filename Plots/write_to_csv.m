% Write data to csv files
Ts_csv = 0.1;

%% HAVOK
load('Plots/MPC_havok_vs_pid_step_1.mat');

t = 0:Ts_csv:50;

X_pid = resample(out.x_pid, t);% Resample time series to desired sample time and period  
X_pid = X_pid.Data;

T_pid = resample(out.T_pid, t); % Thrust of motors
T_pid = T_pid.Data;

X_havok = resample(out.x_mpc, t);
X_havok = X_havok.Data;

T_havok = resample(out.T_mpc, t); % Thrust of motors
T_havok = T_havok.Data;

%% DMD
load('Plots/MPC_dmd_vs_pid_step_1.mat');

X_dmd = resample(out.x_mpc, t);
X_dmd = X_dmd.Data;

T_dmd = resample(out.T_mpc, t); % Thrust of motors
T_dmd = T_dmd.Data;

%% CSV
num_samples = length(t);
num_vars = 1 + (size(X_pid,2) + size(T_pid,2))*3;

dmd_havok_pid = array2table([t', X_pid, T_pid, X_havok, T_havok, X_dmd, T_dmd]);

dmd_havok_pid.Properties.VariableNames = {'t',  'x_pid', 'z_pid', 'th_pid', 'b_pid', 'dx_pid', 'dz_pid', 'dth_pid', 'db_pid', 'T1_pid', 'T2_pid','x_havok', 'z_havok', 'th_havok', 'b_havok', 'dx_havok', 'dz_havok', 'dth_havok', 'db_havok', 'T1_havok', 'T2_havok', 'x_dmd', 'z_dmd', 'th_dmd', 'b_dmd', 'dx_dmd', 'dz_dmd', 'dth_dmd', 'db_dmd', 'T1_dmd', 'T2_dmd'};
writetable(dmd_havok_pid,'Plots/pos_step_pid_vs_mpc.csv')

