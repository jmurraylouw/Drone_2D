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

% Anton controller gain values
kp_dtheta = 0.089;
ki_dtheta = 0.02;
kd_dtheta = 0.003;
N_dtheta = 500; % Ignore

kp_theta = 3;

kp_dx = 0.048;
ki_dx = 0.008;
kd_dx = 0.002;
N_dx = 5; % ignore

kp_x = 0.35;

kp_dz = 0.1;
ki_dz = 0.01;
kd_dz = 0;
N_dz = 5;

kp_z = 0.9;


% Mixing Matrix.
MM = [-1, 1;
       1, 1]; % [T1; T2] = MM*[delta_E; delta_T]