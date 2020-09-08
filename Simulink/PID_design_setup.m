%% Design for angular rate controller

% s = tf('s'); % transfer function variable
% G_Omega_y = (2*r/(tau*I_yy))/(s + (s + 1/tau)); % Omega_y/delta_E_r(s). Transfer function from elevator reference to angular rate

syms s % Laplace variable
syms delta_E
syms I_yy % Moment of inertia of drone body about body x axis (kg.m^3)
syms tau % Time constant for motors (s)
syms r % Distance from COM to rotor thrust (m)


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
Omega_y = (1/s)*ddtheta; % Angular velocity is integral of Angular accelration

% Transfer function from delta_E to Omega_y
G_Omega_y = Omega_y/delta_E;

% Substitute paramater values
I_yy = 0.235; % Moment of inertia of drone body about body x axis (kg.m^3)
tau = 0.07; % Time constant for motors (s)
r = 0.49*1/sqrt(2); % Distance from COM to rotor thrust (m)
G_Omega_y = subs(G_Omega_y);

% Convert to tf object
G_Omega_y = sym2tf(G_Omega_y)

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








