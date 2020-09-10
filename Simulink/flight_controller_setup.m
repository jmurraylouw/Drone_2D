%% PID controllers
load('drone_2D_PID_controllers.mat');

% Pitch Angular rate controller
[kp_dtheta, ki_dtheta, kd_dtheta, Tf_dtheta] = piddata(PID_dtheta);
N_dtheta = 1/Tf_dtheta; % Derivative filter constant

% Pitch Angle controller
kp_theta = PID_theta.Kp;

% Z Velocity controller
[kp_dz, ki_dz, kd_dz, Tf_dz] = piddata(PID_dz);
N_dz = 1/Tf_dz; % Derivative filter constant

% X Velocity controller
kp_dx = 10;
ki_dx = 0.2;
kd_dx = 1.3;
N_dx = 100;

%% Save exported controllers from Simulink Tuner
% For when a new controller is designed and exported
save('drone_2D_PID_controllers.mat', 'PID_dtheta', 'PID_theta', 'PID_dz')

%% Initial conditions
x0 = zeros(6,1); 
% x0(3) = 0.1;
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

% Mixing Matrix.
MM = [-1, 1;
       1, 1]; % [F1; F2] = MM*[delta_E; delta_T]