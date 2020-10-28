% Simulated mpc predictions with for loop following simulation data
% Always run .slx simulation first
% Used to analyse prediction of mpc when designing

% Resample time series to mpc sample time
% x_resamp = resample(out.x_mpc, 0:Ts_mpc:(out.x_mpc.Time(end)) );  
% mo_resamp   = resample(out.mo,    0:Ts_mpc:(out.x_mpc.Time(end)) );  
% mv_resamp   = resample(out.mv,    0:Ts_mpc:(out.x_mpc.Time(end)) );  
% ref_resamp  = resample(out.ref,   0:Ts_mpc:(out.x_mpc.Time(end)) );  
% 
% % Extract data
% x_data = xmpc_resamp.Data';
% mo_data   = mo_resamp.Data';
% mv_data   = mv_resamp.Data';
% ref_data  = ref_resamp.Data';

% Extract data
mo_data   = out.mo.Data';
mv_data   = out.mv.Data';
ref_data  = out.ref.Data';

t = xmpc_resamp.Time';
N = size(mo_data,2); % Number of data samples

ph = mpc_drone_2d.PredictionHorizon; % Prediction Horizon
x_mpc = mpcstate(mpc_drone_2d); % Current state of mpc
v = []; % No measured distrubances

for k = 1:N % every timestep k
    ym = mo_data(:, k);
    r = ref_data(:, k);
    [mv, info] = mpcmove(mpc_drone_2d, x_mpc, ym, r, v);
    if mod(k,0.03/Ts_mpc) == 0 && (k*Ts_mpc > 26)
        subplot(2,1,1)
        plot(info.Topt + t(k), info.Yopt(:,y_rows));
        ylim([-15, 10])
        hold on;
        plot(info.Topt + t(k), ref_data(y_rows,(0:ph)+k)')
        plot(info.Topt + t(k), mo_data(y_rows,(0:ph)+k)', ':', 'LineWidth', 2)
        hold off;
        
        subplot(2,1,2)
        plot(info.Topt + t(k), info.Uopt)
        hold on;
        plot(info.Topt + t(k), mv_data(:,(0:ph)+k)', ':', 'LineWidth', 2)
        hold off;
        pause
    end
end
