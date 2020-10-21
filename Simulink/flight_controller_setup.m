%% PID controllers
load('Data/Drone_2D_control_params.mat'); % Load controller gain values

% Dimensions
nx = 6; % Number of states
ny = 3; % Number of measurements
nu = 2; % Number of inputs

% Initial conditions
x0 = zeros(nx,1); % Initial state
y0 = zeros(ny,1); % Initial measurements
u0 = [0; 0]; % Initial input

% Model parameters
M     = 4.5; % Mass of drone body (at fulcrum)
I_yy  = 0.235; % Moment of inertia of drone body about body x axis
r     = 0.49*1/sqrt(2); % Distance from each rotor force to COM of drone
g     = -9.81; % Acceleration due to gravity (always negative)
C_Dx  = 0.2 ; % Damping coef. of drone through air in x direction (f = C_Dx*xdot)
C_Dz  = 0.2; % Damping coef. of drone in z direction (f = cy*zdot)
rho   = 1.225; % Air density (kg/m^3)
tau   = 0.07; % Motor time constant

% Mixing Matrix
MM = [-1, 1;
       1, 1]; % [T1; T2] = MM*[delta_E; delta_T]

% Forgetting factor
lambda = 0.985; % Exponential forgetting factor of moving average to smooth out spikes in input for SVD performance

% Way points
num_waypoints = 100; % Number of waypoints included in command
point_time_interval = 0; % Initial interval between commands

waypoints = table('Size', [(num_waypoints+1)*2, 3], 'VariableTypes', ["double", "double", "double"]);
waypoints.Properties.VariableNames = {'point_time', 'x_coord', 'z_coord'};

x_coord = 0;
z_coord = -5;
waypoints(1,:) = table(0,                   x_coord, z_coord); % Initial point
waypoints(2,:) = table(point_time_interval, x_coord, z_coord); % Initial point for 6 seconds

x_min        = -10;     x_max        = 10; % (m) minimum and maximum coordinates for waypoints
z_min        = -25;     z_max        = -5;
interval_min = 2;       interval_max = 8;  % (s) minimum and maximum time interval between commands

point_time = point_time_interval;
rng(0); % Initialise random number generator for repeatability
for i = 1:num_waypoints
    point_time_interval = (interval_max - interval_min).*rand() + interval_min; % (s) random time interval between commands
    point_time = point_time + point_time_interval;

    waypoints(2*i,  :) = table(point_time, x_coord, z_coord); % Previous point    
    x_coord    = (x_max - x_min).*rand() + x_min; % x coordinate of next waypoint
    z_coord    = (z_max - z_min).*rand() + z_min; % z coordinate of next waypoint   
    waypoints(2*i+1,:) = table(point_time, x_coord, z_coord); % Next point
end
i = i+1;
waypoints(2*i,  :) = table(point_time+interval_max, x_coord, z_coord); % Add time to reach final point

waypoints_ts = timeseries([waypoints.x_coord, waypoints.z_coord], waypoints.point_time); % timeseries object for From Workspace block

% Run with A and x
run_and_plot = 0;
if run_and_plot
    
    k_start = 500-q;
    window = 10*2; % Number of data points to run and test for

    % Initial condition
    y_hat_0 = [];
    for row = 0:q-1 % First column of spaced Hankel matrix
        k = k_start + row;
        y_hat_0 = [y_hat_0; y_data(:,k)];
    end
    
    % Start at end of initial condition k
    Y_data = y_data(:,k:k+window);
    U      = u_data(:,k:k+window);
    T      = t(:,k:k+window);
    N      = length(Y_data);
    
    % Run model
    Y_hat = zeros(length(y_hat_0),N); % Empty estimated Y
    Y_hat(:,1) = y_hat_0; % Initial condition
    for k = 1:N-1
        Y_hat(:,k+1) = A_hat*Y_hat(:,k) + B_hat*U(:,k);
    end

    y_hat = Y_hat(end-ny+1:end, :); % Extract only non-delay time series (last m rows)

%     % Vector of Mean Absolute Error on testing data
%     MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state

    % Plot data vs model
    figure;
    plot(T, Y_data);
    hold on;
%     plot(t_test, y_test);

    plot(T, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
    title('Training and Testing data vs Model');
    hold off;
end % run_and_plot

% DMD parameters
% Ts_dmd, ny, nu, x0, u0, N_train, q, model_intervals
load('Data/MPC_initial_plant.mat'); % load A, B, Ts_dmd
N_train = 30/Ts_dmd; % Num of data samples for training 
q = 6; 
model_intervals = 10; 

% Initial LTI system 
A_dmd = [A,       B(:, 1:end-nu);
         eye((q-1)*ny),   zeros((q-1)*ny,ny)];
B_dmd = [B(:, end-nu+1:end); zeros((q-1)*ny, nu)];
C_dmd = eye(q*ny);
D_dmd = zeros(q*ny, nu);
dmd_sys = ss(A_dmd,B_dmd,C_dmd,D_dmd,Ts_dmd); % LTI system

% Resample model to MPC sample time
Ts_mpc = 0.1;
mpc_sys = d2d(dmd_sys, Ts_mpc);
[A_mpc,B_mpc,C_mpc,D_mpc,~] = ssdata(mpc_sys);

% Initital conditions for extended measurment vector for MPC
% All delay states are also at y0
y_ext_0 = zeros(q*ny, 1); % Allocate space
for row = 0:q-1 % First column of spaced Hankel matrix
        y_ext_0(row*ny+1:(row+1)*ny, 1) = y0;
end

delays_0 = []; % Initial delay vector
for i = 1:q-1
    delays_0 = [delays_0; y0];
end

% MPC object
old_status = mpcverbosity('off'); % No display messages
mpc_sys.InputGroup.MV = 1:nu; % Munipulated Variable indices
mpc_sys.OutputGroup.MO = 1:q*ny; % Measured Output

tuning_weight = 1; % Smaller = robust, Larger = aggressive
mpc_drone_2d = mpc(mpc_sys,Ts_mpc);

t_s = 4; % Settling time (s)
mpc_drone_2d.PredictionHorizon  = t_s/Ts_mpc; % Prediction horizon (samples), initial guess according to MATLAB: Choose Sample Time and Horizons
mpc_drone_2d.ControlHorizon     = 1; % Control horizon (samples)

mpc_drone_2d.Weights.OutputVariables        = [1, 1, 0, zeros(1, (q-1)*ny)]*tuning_weight;
mpc_drone_2d.Weights.ManipulatedVariables   = 1e-4*[1, 1]*tuning_weight; % Weights of delay coordinates to 0
mpc_drone_2d.W.ManipulatedVariablesRate     = 1e-4*[1, 1]/tuning_weight;

% Output bounds
theta_min = -30*(pi/180);
theta_max = abs(theta_min); % Anton's pitch command constraint

mpc_drone_2d.OV(3).Min = theta_min;
mpc_drone_2d.OV(3).Max = theta_max;

% Input bounds

% Normalised
F_r_z_min = -150;
F_r_z_max = 150; % ??? Assume designed for 2:1 power to weight ratio

F_r_x_max = F_r_z_max;
F_r_x_min = -F_r_x_max;

mpc_drone_2d.MV(1).Min = F_r_x_min;
mpc_drone_2d.MV(1).Max = F_r_x_max;

mpc_drone_2d.MV(2).Min = F_r_z_min;
mpc_drone_2d.MV(2).Max = F_r_z_max;

% Display
mpc_drone_2d 

% Nominal operating conditions for AMPC block
u_bar = [0; M*g];
U_nom = zeros(nu,1);
X_nom = zeros(q*ny,1);
Y_nom = zeros(q*ny,1);
DX_nom = zeros(q*ny,1);











