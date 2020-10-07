%% Setup workspace for simulation: flight_controller_payload.slx
% Called from InitFunc callback

%% PID controllers
load('Data/Drone_2D_control_params.mat'); % Load controller gain values

%% Initial conditions
n = 6; % Number of states
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
cbeta = 0.1; % Rotational damping coef of payload at connection

% Mixing Matrix
MM = [-1, 1;
       1, 1]; % [T1; T2] = MM*[delta_E; delta_T]

% Forgetting factor
lambda = 0.985; % Exponential forgetting factor of moving average to smooth out spikes in input for SVD performance

% Way points
num_waypoints = 100; % Number of waypoints included in command
point_time_interval = 4; % (s) time interval between points

waypoints = table('Size', [(num_waypoints+1)*2, 3], 'VariableTypes', ["double", "double", "double"]);
waypoints.Properties.VariableNames = {'point_time', 'x_coord', 'z_coord'};

x_coord = 0;
z_coord = -5;
waypoints(1,:) = table(0,                   x_coord, z_coord); % Initial point
waypoints(2,:) = table(point_time_interval, x_coord, z_coord); % Initial point

x_min = -10; x_max = 10; % minimum and maximum coordinates for waypoints
z_min = -25; z_max = -5;

rng(0); % Initialise random number generator for repeatability
for i = 1:num_waypoints
    point_time = point_time_interval*i;
    waypoints(2*i,  :) = table(point_time_interval*i, x_coord, z_coord); % Previous point
    
    x_coord    = (x_max - x_min).*rand() + x_min; % x coordinate of next waypoint
    z_coord    = (z_max - z_min).*rand() + z_min; % z coordinate of next waypoint   
    waypoints(2*i+1,:) = table(point_time_interval*i, x_coord, z_coord); % Next point
end
i = i+1;
waypoints(2*i,  :) = table(point_time_interval*i, x_coord, z_coord); % Add time to reach final point

waypoints_ts = timeseries([waypoints.x_coord, waypoints.z_coord], waypoints.point_time); % timeseries object for From Workspace block

% figure;
% title('Waypoint route')
% plot(waypoints.x_coord, waypoints.z_coord)
% xlabel('x');
% ylabel('z');

