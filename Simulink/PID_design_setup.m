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

%% From F_x_r to theta_sp
% F_x_r = delta_T*sin(theta_sp)
% Small angle approx: sin(theta_sp) == theta_sp
% Therefore: F_x_r = delta_T*(theta_sp)

% F_z = -delta_T*cos(theta)
% Small angle approx: cos(theta) == theta
% Therefore: F_z = -delta_T
% Linearise condition: F_z = Mg
delta_T = -M*g;
eqn = (F_x_r == delta_T*(theta_sp));
theta_sp = solve(eqn, theta_sp);
G_theta_sp = theta_sp/F_x_r;

%% From F_x to dx
eqn = (F_x - C_x_lin*dx == s*M*dx); % Equation of Newton 2nd law in x direction
dx = solve(eqn, dx); % Solve for dx according to F_x from Newton 2nd law
G_dx = dx/F_x; % TF from F_x to dx

%%
% Substitute paramater values
I_yy = 0.235; % Moment of inertia of drone body about body x axis (kg.m^3)
tau = 0.07; % Time constant for motors (s)
r = 0.49*1/sqrt(2); % Distance from COM to rotor thrust (m)
G_dtheta = subs(G_dtheta);
C_Dx = 0.2; % Body Aerodynamic Coefficient in x (m^2)
C_x_lin = C_Dx*dx_bar; % Drag coef of linearised drag with average velocity

% Convert to tf object
G_dtheta = sym2tf(G_dtheta)

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








