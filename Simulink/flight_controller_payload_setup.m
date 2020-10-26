%% Setup workspace for simulation: flight_controller_payload.slx
% Called from InitFunc callback

%% PID controllers
load('Data/Drone_2D_control_params.mat'); % Load controller gain values

dx_sp_max = 12; % (m/s) Max x velocity command (x and z)
dx_sp_min = -12; % (m/s) Min x velocity command (x and z)

dz_sp_max = 12; % (m/s) Max z velocity command (x and z)
dz_sp_min = -12; % (m/s) Min z velocity command (x and z)

theta_sp_max = 30*pi/180; % (rad) Max. pitch rate command
theta_sp_min = -30*pi/180; % Min. pitch rate command

dtheta_sp_max = 120*pi/180; % Max. pitch rate command
dtheta_sp_min = -120*pi/180; % Min. pitch rate command

sigma = 0.001; % Std deviation of measurement noise

%% Initial conditions
n = 8; % Number of states
x0 = zeros(n,1); 
u0 = [0; 0];

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

% figure;
% title('Waypoint route')
% plot(waypoints.x_coord, waypoints.z_coord)
% xlabel('x');
% ylabel('z');

