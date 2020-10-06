%% PID controllers
load('Data/Drone_2D_control_params.mat'); % Load controller gain values

%% Initial conditions
x0 = zeros(6,1); 
u0 = [0; 0];

% Model parameters
M = 4.5; % Mass of drone body (at fulcrum)
I_yy = 0.235; % Moment of inertia of drone body about body x axis
r = 0.49*1/sqrt(2); % Distance from each rotor force to COM of drone
g = -9.81; % Acceleration due to gravity (always negative)
C_Dx = 0.2 ;% Damping coef. of drone through air in x direction (f = C_Dx*xdot)
C_Dz = 0.2; % Damping coef. of drone in z direction (f = cy*zdot)
rho = 1.225; % Air density (kg/m^3)
tau = 0.07; % Motor time constant

% Mixing Matrix
MM = [-1, 1;
       1, 1]; % [T1; T2] = MM*[delta_E; delta_T]
   
% Way points
num_waypoints = 10; % Number of waypoints included in command
waypoints = table('Size', [num_waypoints, 3], 'VariableTypes', ["double", "double", "double"]);
waypoints.Properties.VariableNames = {'point_time', 'x_coord', 'z_coord'};
waypoints(1,:) = table(0, 0, -5); % Initial point
point_time_interval = 10; % (s) time interval between points
x_min = -10; x_max = 10; % minimum and maximum coordinates for waypoints
z_min = -25; z_max = -5;

rng(0); % Initialise random number generator for repeatability
for i = 1:2:num_waypoints*2+1
    x_coord    = (x_max - x_min).*rand() + x_min; % x coordinate of next waypoint
    z_coord    = (z_max - z_min).*rand() + z_min; % z coordinate of next waypoint
    waypoints(i+1,:) = table(point_time_interval*(i),   x_coord, z_coord);
    waypoints(i+2,:) = table(point_time_interval*(i+1), x_coord, z_coord);
end
waypoints
waypoints_ts = timeseries([waypoints.x_coord, waypoints.z_coord], waypoints.point_time); % timeseries object for import into simulink
waypoints_ts = setinterpmethod(waypoints_ts,'zoh'); % Set ZOH as interpolation method

% figure;
% title('Waypoint route')
% plot(waypoints.x_coord, waypoints.z_coord)
% xlabel('x');
% ylabel('z');

