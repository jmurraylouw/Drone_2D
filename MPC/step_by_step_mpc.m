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

load('Data/MPC_sim_2.mat')

% Extract data
mo_data   = out.mo.Data';
mv_data   = out.mv.Data';
ref_data  = out.ref.Data';

t = out.mo.Time';
N = size(mo_data,2); % Number of data samples

ph = mpc_drone_2d.PredictionHorizon; % Prediction Horizon
x_mpc = mpcstate(mpc_drone_2d); % Current state of mpc
v = []; % No measured distrubances

% x_mpc = x_mpc_19; % State of mpc at 19.5 s

for k = 1/Ts_mpc:N % every timestep k
    k*Ts_mpc
    ym = mo_data(:, k);
    r = ref_data(:, k);
    [mv, info] = mpcmove(mpc_drone_2d, x_mpc, ym, r, v);
    if mod(k, 0.2/Ts_mpc) == 0 && (k*Ts_mpc > 4.5)
        for state = y_rows
            figure(state)
            ylabel(state)
            ylim([-2.5, 2])
            plot(info.Topt + t(k), info.Yopt(:,state));
            ylim([-15, 10])
            hold on;
            plot(info.Topt + t(k), ref_data(state,(0:ph)+k)')
            plot(info.Topt + t(k), mo_data(state,(0:ph)+k)', ':', 'LineWidth', 2)
            legend('ref', 'prediction', 'actual')
            ylim([-2.5, 2.5])
            hold off;
        end
        
        figure(state + 1)
        plot(info.Topt + t(k), info.Uopt)
        hold on;
        plot(info.Topt + t(k), mv_data(:,(0:ph)+k)', ':', 'LineWidth', 2) % Actual input given
        hold off;
        pause
    end
end
