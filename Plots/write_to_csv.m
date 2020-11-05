% Write data to csv files
Ts_csv = 0.1;

t = 0:Ts_csv:50;
X_nl = resample(out.x_no_load, t);% Resample time series to desired sample time and period  
X_nl = X_nl.Data;

X_wl = resample(out.x_with_load, t);% Resample time series to desired sample time and period  
X_wl = X_wl.Data;

num_samples = length(t);
num_vars = 1 + size(X_nl,2) + size(X_wl,2);

pos_step_pid = array2table([t', X_nl, X_wl]);
% % t = time for both nl and wl
% % nl = no load
% % wl = with load
% % th = theta (pitch angle)
% % b = beta (payload angle)
pos_step_pid.Properties.VariableNames = {'t', 'x_nl', 'z_nl', 'th_nl', 'dx_nl', 'dz_nl', 'dth_nl', 'x_wl', 'z_wl', 'th_wl', 'b_wl', 'dx_wl', 'dz_wl', 'dth_wl', 'db_wl'};
writetable(pos_step_pid,'Plots/pos_step_pid.csv')

