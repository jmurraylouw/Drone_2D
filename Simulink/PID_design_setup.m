%% Design for angular rate controller

% s = tf('s'); % transfer function variable
% G_Omega_y = (2*r/(tau*I_yy))/(s + (s + 1/tau)); % Omega_y/delta_E_r(s). Transfer function from elevator reference to angular rate

syms s % Laplace variable
syms delta_E
syms I_yy % Moment of inertia of drone body about body x axis (kg.m^3)
syms tau % Time constant for motors (s)
syms r % Distance from COM to rotor thrust (m)
syms theta_sp % Pitch set-point (rad)
syms dx M F_x C_x_lin
syms g
syms F_x_r

% Motor mixing algorithm
T1_r = -delta_E; % Thrust reference 1
T2_r = delta_E; % Thrust reference 2

% Motor model
G_motor = (1/tau)/(s + 1/tau);
T1 = T1_r*G_motor; % Actual thrust
T2 = T2_r*G_motor;

% Body moment
M_By = -T1*r + T2*r; % Momement around y-axis

% Newtons law
% M_By = I_yy*ddtheta
ddtheta = M_By/I_yy;

% Integrator
dtheta = (1/s)*ddtheta; % Angular velocity is integral of Angular accelration

% Transfer function from delta_E to dtheta
G_dtheta = dtheta/delta_E;

% Design requirements:

% Reject disturbances
% Zero steady-state error
% overshoot = 4.2%
% bandwidth = 12.52 rad/s
% 2% settling time = 0.6 s

PO = 4.2; % Percentage Overshoot (%)
wb = 13; % Bandwidth (rad/s)
ts = 0.6; % 2% settling time (s)

% Percentage overshoot:
zeta = sqrt( (log(PO/100))^2 / (pi^2 + (log(PO/100))^2) );  % Damping ratio
theta = atan(sqrt(1 - zeta^2) / zeta); % Max angle from real axis to dominant pole

% Settling time:
sigma = -log(0.02)/ts % Real part limit of dominant pole

% First implement PI controller
% D_PI = Kp*(s + z_c) / s
% Use I controller to reject disturbances (closed loop becomes type 2)
% Use P controller to place in performance envelope

% Pole of D_PI is at origin, so let zero be close to origin: z_c = 0.1
z_c = 0.1;

% Draw root locus
figure;
hold on;
rlocus();

% Plot requirement limits
plot([1, 1]*sigma, ylim); % Settling time requirement limit


% Transfer function of dtheta PID controller 
D_dtheta = kp_dtheta + ki_dtheta*(1/s) + kd_dtheta*(s/(1/N_dtheta*s + 1));
G_dtheta_cl = D_dtheta*G_dtheta/(1 + D_dtheta*G_dtheta); % Closed loop tf with PID control for dtheta

% TF from dtehta_sp to theta (seen by angular rate controller)
theta = (1/s)*dtheta;
% dtheta = G_dtheta_cl*theta_sp;
G_theta = G_dtheta_cl*(1/s);

% Transfer function with theta P controller 
D_theta = kp_theta;
G_theta_cl = D_theta*G_theta/(1 + D_theta*G_theta); % Closed loop tf with PID control for theta

%% F_x_r to theta_sp
% F_x_r = -delta_T*sin(theta_sp)
% Small angle approx: sin(theta_sp) == theta_sp
% Therefore: F_x_r = delta_T*(theta_sp)

% F_z = -delta_T*cos(theta)
% Small angle approx: cos(theta) == 1
% Therefore: F_z = -delta_T
% Linearise condition: F_z = Mg
delta_T = -M*g;
eqn = (F_x_r == -delta_T*(theta_sp));
theta_sp = solve(eqn, theta_sp);
G_theta_sp = theta_sp/F_x_r

%% theta to dx
% F_x = -delta_T*sin(theta)
% Small angle approx: sin(theta_sp) == theta_sp
F_x = -delta_T*(theta);
eqn = (F_x - C_x_lin*dx == s*M*dx); % Equation of Newton 2nd law in x direction
dx = solve(eqn, dx); % Solve for dx according to F_x from Newton 2nd law
G_th_dx = dx/theta; % TF from theta to dx

%% F_x_r to dx
G_dx = G_theta_sp*G_theta_cl*G_th_dx;

% Model parameters
M = 4.5; % Mass of drone body (at fulcrum)
I_yy = 0.235; % Moment of inertia of drone body about body x axis
r = 0.49*1/sqrt(2); % Distance from each rotor force to COM of drone
g = -9.81; % Acceleration due to gravity (always negative)
C_Dx = 0.2 ;% Damping coef. of drone through air in x direction (f = C_Dx*xdot)
C_Dz = 0.2; % Damping coef. of drone in z direction (f = cy*zdot)
rho = 1.225; % Air density (kg/m^3)
tau = 0.07; % Motor time constant
dx_bar = 5; % Average x velocity to linearise (m/s)
C_x_lin = C_Dx*dx_bar; % Drag coef of linearised drag with average velocity

% Substitute paramater values
G_dtheta = subs(G_dtheta);
G_dx = subs(G_dx);

% Convert to tf object
G_dtheta = sym2tf(G_dtheta)
G_dx = sym2tf(subs(G_dx));
G_th_dx = sym2tf(subs(G_th_dx))

% %% Display values from angular rate controller design
% D_Omega_y.Kp
% D_Omega_y.Ki
% D_Omega_y.Kd
% 
% save('Data/PID_tuned_values.mat', 'D_Omega_y', 'D_theta') % Save exported P structure from Simulink design
% 
% %% Display values from angular rate controller design
% D_theta.Kp
% 
% save('Data/PID_tuned_values.mat', 'D_Omega_y', 'D_theta') % Save exported P structure from Simulink design

function TF = sym2tf(sym_TF)
% Converts symbolic representation of transfer function
% to a transfer function object
% from: https://www.mathworks.com/matlabcentral/answers/310042-how-to-convert-symbolic-expressions-to-transfer-functions
    ExpFun = matlabFunction(sym_TF);
    ExpFun = str2func(regexprep(func2str(ExpFun), '\.([/^\\*])', '$1'));
    TF = ExpFun(tf('s'));
end








