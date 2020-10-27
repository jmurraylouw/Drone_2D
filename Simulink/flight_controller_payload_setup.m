%% Setup workspace for simulation: flight_controller_payload.slx
% Called from InitFunc callback

% PID controllers
load('Data/Drone_2D_control_params.mat'); % Load controller gain values

dx_sp_max = 12; % (m/s) Max x velocity command (x and z)
dx_sp_min = -12; % (m/s) Min x velocity command (x and z)

dz_sp_max = 12; % (m/s) Max z velocity command (x and z)
dz_sp_min = -12; % (m/s) Min z velocity command (x and z)

theta_sp_max = 30*pi/180; % (rad) Max. pitch rate command
theta_sp_min = -30*pi/180; % Min. pitch rate command

dtheta_sp_max = 120*pi/180; % Max. pitch rate command
dtheta_sp_min = -120*pi/180; % Min. pitch rate command
% 
% sigma = 0.001; % Std deviation of measurement noise

% Dimentions
y_rows = 1:4;
nx = 8; % Number of states
ny = length(y_rows); % Number of measurements
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

m     = 2; % Mass of swinging payload (kg)
l     = 1; % Length of pendulum (m)
cbeta = 0.01; % Rotational damping coef of payload at connection

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
interval_min = 4;       interval_max = 8;  % (s) minimum and maximum TIME interval between commands

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

% DMD parameters
% Ts_dmd, ny, nu, x0, u0, N_train, q, model_intervals
% load('Data/MPC_initial_plant.mat'); % load A_dmd, B_dmd, q, Ts_dmd from a previous DMD run
% model_intervals = 10; 

% Sample time of MPC:
t_s = 5; % For guidance, minimum desired settling time (s)
Ts_mpc = 0.05; % Guide: between 10% to 25% of desired response time

simulation_data_file = 'With_payload_data_1';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data for dmd model
start_time = 30;
end_time = 60;
q = 2;
sigma = 0;
plot_prediction = 1;
u_bar = [0; (M+m)*g];
 
[A_dmd, B_dmd] = model_DMD(out, start_time, end_time, Ts_mpc, q, y_rows, sigma, plot_prediction, u_bar);

% Disturbance model to account for model uncertainty (eliminate steady-state error)
CO = 2; % number of Controlled Outputs (x and z). theta is not controlled to a reference
dist_influence = 0; % 1e-5; % Disturbances include uncertainty in model
B_ud = dist_influence*[eye(CO); zeros(q*ny - CO, CO)]; % B of unmeasured disturbance, for distrubance force in x and z

% Change model structure so delays are included in A, not B 
A_mpc = [A_dmd,       B_dmd(:, 1:end-nu);
         eye((q-1)*ny),   zeros((q-1)*ny,ny)];

B_mpc = [[B_dmd(:, end-nu+1:end); zeros((q-1)*ny, nu)], B_ud];
C_mpc = eye(q*ny);
D_mpc = zeros(q*ny, nu + CO);
mpc_sys = ss(A_mpc,B_mpc,C_mpc,D_mpc,Ts_mpc); % LTI system

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
mpc_sys.OutputGroup.MO = 1:q*ny; % Measured Output

mpc_sys.InputGroup.MV = 1:nu; % Munipulated Variable indices
mpc_sys.InputGroup.UD = nu + (1:CO); % unmeasured disturbance indices, one for each 

tuning_weight = 1; % Smaller = robust, Larger = aggressive
mpc_drone_2d = mpc(mpc_sys,Ts_mpc);

% Guide: PH so PH*Ts == desired responce time
mpc_drone_2d.PredictionHorizon  = 25; %t_s/Ts_mpc; % Prediction horizon (samples), initial guess according to MATLAB: Choose Sample Time and Horizons
mpc_drone_2d.ControlHorizon     = 3; % Control horizon (samples)

mpc_drone_2d.Weights.OutputVariables        = [1, 1, 0, 0, zeros(1, (q-1)*ny)]*tuning_weight;
mpc_drone_2d.Weights.ManipulatedVariables   = 1e-3*[1, 1]*tuning_weight; % Weights of delay coordinates to 0
mpc_drone_2d.Weights.ManipulatedVariablesRate     = 1e-1*[5, 1]/tuning_weight;

% Output bounds
theta_min = -30*(pi/180);
theta_max = abs(theta_min); % Anton's pitch command constraint

mpc_drone_2d.OV(3).Min = theta_min;
mpc_drone_2d.OV(3).Max = theta_max;

% Input bounds

% Normalised
F_r_z_min = -100;
F_r_z_max = 100; % ??? Assume designed for 2:1 power to weight ratio

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












