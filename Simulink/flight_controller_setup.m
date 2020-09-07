% Initial conditions
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
MM = [-1, -1;
       1, -1]; % [F1; F2] = MM*[delta_E; delta_T]