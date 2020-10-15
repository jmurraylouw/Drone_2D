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

% Initital conditions for extended measurment vector for MPC
% All delay states are also at x0
y_ext_0 = zeros(q*ny, 1); % Allocate space
for row = 0:q-1 % First column of spaced Hankel matrix
        y_ext_0(row*ny+1:(row+1)*ny, 1) = y0;
end

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
point_time_interval = 6; % Initial interval between commands

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

% figure;
% title('Waypoint route')
% plot(waypoints.x_coord, waypoints.z_coord)
% xlabel('x');
% ylabel('z');

%% LTI system
Ts_mpc = 0.1; % MPC sampling time

dmd_sys = ss(A_hat,B_hat,eye(q*ny),zeros(q*ny,nu),Ts); % LTI system
dmd_sys = d2d(dmd_sys,Ts_mpc,'zoh'); % Resample to match MPC

[A_hat2,B_hat2,C_hat2,D_hat2,Ts] = ssdata(dmd_sys); % Extract resampled matrixes
dmd_sys = ss(A_hat,B_hat,C_hat2,D_hat2,Ts_mpc); % LTI system with new Ts 

%% MPC object
Ts_mpc = 0.1; % Sample time of MPC (s)
old_status = mpcverbosity('off'); % No display messages
dmd_sys.InputGroup.MV = 2; % Munipulated Variable
dmd_sys.OutputGroup.MO = 6; % Measured Variable

mpc_drone_2d = mpc(dmd_sys,Ts_mpc);
mpc_drone_2d.Weights.OutputVariables = [zeros(1,(q-1)*ny), 1, 1, 0]; % Set weights of delay coordinates to 0, so they do not follow reference

%% Nominal operating conditions for AMPC block
U_nom = u_bar;
X_nom = zeros(q*ny,1);
Y_nom = zeros(q*ny,1);
DX_nom = zeros(q*ny,1);








