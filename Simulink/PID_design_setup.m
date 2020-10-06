%% Design of cascaded loop controllers for 2D drone

close all

syms s % Laplace variable
syms delta_E % Virtual elevator
syms I_yy % Moment of inertia of drone body about body x axis (kg.m^3)
syms tau % Time constant for motors (s)
syms r % Distance from COM to rotor thrust (m)
syms theta_sp % Pitch set-point (rad)
syms dx
syms F_x 
syms C_x_lin
syms g
syms F_x_r

% Model parameters
M = 4.5; % Mass of drone body (at fulcrum)
I_yy = 0.235; % Moment of inertia of drone body about body x axis
r = 0.49*1/sqrt(2); % Distance from each rotor force to COM of drone
g = -9.81; % Acceleration due to gravity (always negative)
C_Dx = 0.2;% Damping coef. of drone through air in x direction (f = C_Dx*xdot)
C_Dz = 0.2; % Damping coef. of drone in z direction (f = cy*zdot)
rho = 1.225; % Air density (kg/m^3)
tau = 0.07; % Motor time constant
dx_bar = 1; % Average x velocity to linearise at (m/s)
dz_bar = 1; % Average z velocity to linearise at (m/s)

% % Anton controller gain values
% kp_dtheta = 0.089;
% ki_dtheta = 0.02;
% kd_dtheta = 0.003;
% N_dtheta = 500; % Ignore
% 
% kp_theta = 3;
% 
% kp_dx = 0.048;
% ki_dx = 0.008;
% kd_dx = 0.002;
% N_dx = 5; % ignore
% 
% kp_x = 0.35;
% 
% kp_dz = 0.1;
% ki_dz = 0.01;
% kd_dz = 0;
% N_dz = 5;
% 
% kp_z = 0.9;

% Motor mixing algorithm
T1_r = -delta_E; % Thrust reference 1
T2_r = delta_E; % Thrust reference 2

% Motor model
G_motor = (1/tau)/(s + 1/tau);
T1 = T1_r*G_motor; % Actual thrust
T2 = T2_r*G_motor;
motor_wb = bandwidth(sym2tf(G_motor));

% Body moment
M_By = -T1*r + T2*r; % Momement around y-axis

% Newtons law
% M_By = I_yy*ddtheta
ddtheta = M_By/I_yy;

% Integrator
dtheta = (1/s)*ddtheta; % Angular velocity is integral of Angular accelration

%%-------------------------------------------------------------
%% Angular rate controller
%%-------------------------------------------------------------

% Plant Transfer Function
% Transfer function from delta_E to dtheta
G_dtheta = dtheta/delta_E;
G_dtheta = subs(G_dtheta); % Substitute paramater values
G_dtheta_tf = sym2tf(G_dtheta);

%% Design requirements:

% Reject disturbances
% Zero steady-state error
PO = 4.2; % Percentage Overshoot (%)
wb_buffer = 1.9; % Add buffer because wb drops when kd and LPF is added
wb_inner = motor_wb; % Inner control loop bandwidth is from motor controller (rad/s)
wb = 12.52 + wb_buffer; % Desired dandwidth (rad/s)
ts = 0.6; % 2% settling time (s)

% Percentage overshoot:
zeta = sqrt( (log(PO/100))^2 / (pi^2 + (log(PO/100))^2) );  % Damping ratio
zeta = 1/sqrt(2); % Optimally damped
theta_pole = atan(sqrt(1 - zeta^2) / zeta); % Max angle from real axis to dominant pole

% Settling time:
sigma_p = log(0.02)/ts; % Real part limit of dominant pole, p for pole to avoid confusion with noise

% Plot requirement limits
figure;
hold on;
y_limits = [-8, 8];
plot([1, 1]*sigma_p, y_limits, '--') % Settling time requirement limit
x_theta = max(y_limits)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(y_limits), '--')

% Root locus with P controller
rlocus(sym2tf(G_dtheta));
title('G(dtheta) with P controller root locus varied by kp')

%% PI controller for dtheta:

% D_PI = Kp*(s + z_c) / s
% Use I controller to reject disturbances (closed loop becomes type 2)
% Use P controller to place in performance envelope

% Desired poles to meet zeta and sigma_p requirement
omega_p = abs(sigma_p*tan(theta_pole)); % Imag part of desired pole
desired_pole = sigma_p + 1i*omega_p;

% TF of PI controller
% D_pi = kp*(s + ki/kp) / s; % transfer function of Pi controller
% D_pi = kp*(s + z_c) / s;

% Place zero of D_pi to reach desired pole
syms z_c_angle % Angle from controller zero to desired pole
pole_angles = sum(angle(desired_pole - [pole(G_dtheta_tf); 0])); % sum of angles from open-loop poles to desired pole
% Note: controller adds open-loop pole at origin
% Note: angle from point A to B = angle(B-A)

zero_angles = sum(angle(desired_pole - zero(G_dtheta_tf))); % sum of angles from open-loop zeros to desired pole
% Note: controller adds open-loop zero at z_c, with angle z_c_angle to
% desired pole

% Angle from zero to desired pole needed for Angle Criteria
z_c_angle = (pole_angles - zero_angles - pi);

% Zero needed to provide angle to complete Angle Criteria
syms z_c % Zero of controller. z_c = ki/kp
z_c_angle_eqn = ( tan(z_c_angle) == (imag(desired_pole))/(real(desired_pole) - (z_c)) );
z_c = double(solve(z_c_angle_eqn, z_c));

% Use Magnitude criteria to determine kp 
syms kp_dtheta
D_pi = (s - z_c) / s; % Controller without kp
magnitude = abs(double(subs(D_pi*G_dtheta, s, desired_pole))); % substitute desired pole into s
% mag_criteria = (magnitude == 1);
kp_dtheta = 1/magnitude;
% D_pi = kp_dtheta*D_pi; % Add kp to PI controller

% Calculate ki from z_c and kp
ki_dtheta = -kp_dtheta*z_c;

% Calculate Kp needed for desired bandwidth
% % kp_min = 0.001;
% % kp_max = 10;
% % wb_tol = 0.001; % tolerance on bandwidth frequency
% % % % kp_dtheta = kp_for_bandwidth(D_pi*G_dtheta,wb,wb_tol,kp_min,kp_max);

% Root locus with PI controller
figure;
rlocus(sym2tf(D_pi*G_dtheta));
title('G(dtheta) with PI controller root locus varied by kp')
hold on;

% Plot current poles for kp needed for bandwidth
current_pole = rlocus(sym2tf(D_pi*G_dtheta), kp_dtheta);
plot(real(current_pole), imag(current_pole), 'rs', 'Markersize', 7); % Plot current pole locatiosn
ylim([min(imag(current_pole)), max(imag(current_pole))]*1.2); % adjust ylim to include plotted poles

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit
x_theta = max(ylim)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(ylim), '--');

% Update Pi controller with final kp value
D_pi = kp_dtheta*D_pi;

%% PID controller for dtheta:

N_dtheta = 2*wb; % Frequency of Low Pas Filter (Manually tune)

% Manual root locus plot of closed loop system with PID controller
% % syms kd_dtheta
D_pid = @(kd_dtheta) kp_dtheta + ki_dtheta*(1/s) + kd_dtheta*s;
G_dtheta_cl = @(kd_dtheta) D_pid(kd_dtheta)*G_dtheta/(1 + D_pid(kd_dtheta)*G_dtheta); % Closed loop tf with P control for dtheta

figure;
hold on;
title('G(dtheta) with PID controller root locus varied by kd');

% for kd_dtheta = 0:0.01:0.5
%     poles = pole(sym2tf(G_dtheta_cl(kd_dtheta)));
%     plot(3real(poles), imag(poles), 'k.'); % Plot pole of current k
% end

% Starting poles
poles = pole(sym2tf(G_dtheta_cl(0)));
plot(real(poles), imag(poles), 'bx'); % Plot pole of current k

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit
x_theta = max(ylim)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(ylim), '--');
grid on;

% Choose kd to place poles within desired region
kd_dtheta = 0.095;
kd_dtheta = 0;
D_pid = D_pid(kd_dtheta);
G_dtheta_cl = G_dtheta_cl(kd_dtheta);

% Current poles
poles = pole(sym2tf(G_dtheta_cl));
plot(real(poles), imag(poles), 'rs', 'MarkerSize', 7); % Plot pole of current k

% Step responce
t_dist = 10;
controller_step_responce(G_dtheta, [D_pid, D_pi, kp_dtheta], {'PID', 'PI', 'P'}, t_dist)
title('Step responce of dtheta controllers')

% Performance parameters
sys = sym2tf(G_dtheta_cl);
dtheta_performance = stepinfo(sys);
wb = bandwidth(sys);
dtheta_performance.Bandwidth = wb;
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
dtheta_performance.slow_factor = slow_factor;
dtheta_performance

%%-------------------------------------------------------------
%% Angle controller
%%-------------------------------------------------------------

%% Design requirements:
% Zero steady-state error
% Overdamped
% Timescale seperation from inner loop

PO = 0; % Percentage Overshoot (%)
wb_inner = dtheta_performance.Bandwidth;
wb = 4.41; % Desired dandwidth (rad/s).
ts = 1.95; % 2% settling time (s)

%% Plant Transfer Function
% TF from dtehta_sp to theta (seen by angular rate controller)
theta = (1/s)*dtheta;
dtheta = G_dtheta_cl*theta_sp;
G_theta = simplifyFraction(G_dtheta_cl*(1/s));
G_theta_tf = sym2tf(G_theta);

% Calculate Kp needed for desired bandwidth (Binary search)
wb_tol = 0.001;
kp_min = 0.001;
kp_max = 10;
kp_theta = kp_for_bandwidth(G_theta,wb,wb_tol,kp_min,kp_max);

% Transfer function inclucding P controller
D_p = kp_theta;
G_theta_cl = D_p*G_theta/(1 + D_p*G_theta); % Closed loop tf with PID control for theta
% G_theta_cl = theta/theta_sp

% % Bode of closed loop plant with Kp
% figure;
% bode(sym2tf(G_theta_cl));
% title('G(theta) closed-loop with Kp for desired bandwidth');
% grid on;

% Root locus of plant with P controller
figure;
rlocus(sym2tf(G_theta));
title('G(theta) with P controller root locus varied by kp')
hold on;

% Plot current poles for kp needed for bandwidth
current_pole = rlocus(sym2tf(kp_theta*G_theta), kp_theta);
plot(real(current_pole), imag(current_pole), 'rs', 'Markersize', 7); % Plot current pole locatiosn

% Settling time:
sigma_p = log(0.02)/ts; % Real part limit of dominant pole, p for pole to avoid confusion with noise

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit

% Step responce
t_dist = 10;
controller_step_responce(G_theta, D_p, {'P'}, t_dist)
title('Step responce of theta controllers')

% Performance parameters
sys = sym2tf(D_p*G_theta/(1 + D_p*G_theta));
theta_performance = stepinfo(sys);
wb = bandwidth(sys);
theta_performance.Bandwidth = wb;
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
theta_performance.slow_factor = slow_factor;
theta_performance

%%-------------------------------------------------------------
%% X Velocity controller
%%-------------------------------------------------------------

%% Design requirements:

% Reject disturbances
% Zero steady-state error
PO = 12; % Percentage Overshoot (%)
wb_inner = theta_performance.Bandwidth;
wb = 2.166; % Desired dandwidth (rad/s). Slower than previous wb by a factor
ts = 11.6; % 2% settling time (s)

% Plant Transfer Function:

% F_x_r to theta_sp:
% ------------------
% F_x_r = -delta_T*sin(theta_sp)
% Small angle approx: sin(theta_sp) == theta_sp
% Therefore: F_x_r = delta_T*(theta_sp)

% F_z = -delta_T*cos(theta)
% Small angle approx: cos(theta) == 1
% Therefore: F_z = -delta_T
% Linearise condition: F_z = Mg
delta_T = -M*g;
syms theta_sp F_x_r
eqn = (F_x_r == -delta_T*(theta_sp));
theta_sp = solve(eqn, theta_sp);
G_F_x_r = theta_sp/F_x_r;

% theta to dx:
% ------------
% F_x = -delta_T*sin(theta)
% Small angle approx: sin(theta_sp) == theta_sp
syms theta dx
F_x = -2*delta_T*(theta);
C_x_lin = C_Dx*dx_bar; % Linearised drag coef. at average velocity
eqn = (F_x - C_x_lin*rho*dx == s*M*dx); % Equation of Newton 2nd law in x direction
dx = solve(eqn, dx); % Solve for dx according to F_x from Newton 2nd law
G_th_dx = dx/theta; % TF from theta to dx

% F_x_r to dx:
G_dx = G_F_x_r*G_theta_cl*G_th_dx; % dx/F_xr
G_dx_tf = sym2tf(G_dx);

% PI controller for dx:

% Percentage overshoot:
zeta = sqrt( (log(PO/100))^2 / (pi^2 + (log(PO/100))^2) );  % Damping ratio
theta_pole = atan(sqrt(1 - zeta^2) / zeta); % Max angle from real axis to dominant pole

% Settling time:
sigma_p = log(0.02)/ts; % Real part limit of dominant pole, p for pole to avoid confusion with noise

% Now implement PI controller
% D_PI = Kp*(s + z_c) / s
% Use I controller to reject disturbances (closed loop becomes type 2)
% Use P controller to place in performance envelope

% Pole of D_PI is at origin, so let zero be close to origin: z_c = 0.1
z_c = 0.1; % z_c = ki/kp
D_pi = (s + z_c) / s; % transfer function of Pi controller without kp

% Place kp for bandwidth
kp_min = 0.001;
kp_max = 30;
wb_tol = 0.001; % tolerance on bandwidth frequency
kp_dx = kp_for_bandwidth(D_pi*G_dx,wb,wb_tol,kp_min,kp_max);

% Calculate ki from z_c and kp
ki_dx = kp_dx*z_c;
D_pi = kp_dx + ki_dx*(1/s); % PI controller TF

% Root locus of plant with PI controller
figure;
rlocus(sym2tf(D_pi*G_dx));
title('G(dx) with PI controller root locus varied by kp')
hold on;

% Plot current poles for kp needed for bandwidth
current_pole = rlocus(sym2tf(D_pi*G_dx), kp_dx);
plot(real(current_pole), imag(current_pole), 'rs', 'Markersize', 7); % Plot current pole locatiosn
ylim([min(imag(current_pole)), max(imag(current_pole))]*1.2); % adjust ylim to include plotted poles

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit
x_theta = max(ylim)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(ylim), '--');

% PID controller for dx:

N_dx = wb*2; % Frequency of Low Pas Filter (Manually tune)

% Manual root locus plot of closed loop system with PID controller
% % syms kd_dx
D_pid = @(kd_dx) kp_dx + ki_dx*(1/s) + kd_dx*s; % Without LPF
% D_pid = @(kd_dx) kp_dx + ki_dx*(1/s) + kd_dx*(N_dx/(1 + N_dx/s)); % Controller TF including low pass filter
G_dx_cl = @(kd_dx) D_pid(kd_dx)*G_dx/(1 + D_pid(kd_dx)*G_dx); % Closed loop tf with P control for dtheta

figure;
hold on;
for kd_dx = 0:0.1:2
    poles = pole(sym2tf(G_dx_cl(kd_dx)));
    plot3(real(poles), imag(poles), kd_dx*(real(poles)./real(poles)), 'k.'); % Plot pole of current k
end

% Starting poles
poles = pole(sym2tf(G_dx_cl(0)));
plot(real(poles), imag(poles), 'bx'); % Plot pole of current k

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit
x_theta = max(ylim)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(ylim), '--');

% Choose kd to place poles within desired region
kd_dx = 0; % Manually adjust

% Insert final parameter into TFs
D_pid = D_pid(kd_dx);
G_dx_cl = G_dx_cl(kd_dx);

% Current poles:
poles = pole(sym2tf(G_dx_cl));
plot3(real(poles), imag(poles), kd_dx*(real(poles)./real(poles)), 'rs', 'MarkerSize', 7); % Plot pole of current k

% Step responce
t_dist = 10;
controller_step_responce(G_dx, [D_pid, D_pi, kp_dx], {'PID', 'PI', 'P'}, t_dist)
title('Step responce of dx controllers')

% Performance parameters
sys = sym2tf(G_dx_cl);
dx_performance = stepinfo(sys);
wb = bandwidth(sys);
dx_performance.Bandwidth = wb;
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
dx_performance.slow_factor = slow_factor;
dx_performance


%%-------------------------------------------------------------
%% X Position controller
%%-------------------------------------------------------------

%% Design requirements:
% Zero steady-state error
% Overdamped
% Timescale seperation from inner loop

PO = 0; % Percentage Overshoot (%)
wb_inner = dx_performance.Bandwidth;
wb = 0.59; % Desired dandwidth (rad/s).
ts = 11.51; % 2% settling time (s)

%% Plant Transfer Function
% TF from x_sp to x (seen by position controller)
% x = (1/s)*dx;
% dx = G_dx_cl*x_sp;
G_x = simplifyFraction(G_dx_cl*(1/s));
G_x_tf = sym2tf(G_x);

% Calculate Kp needed for desired bandwidth (Binary search)
wb_tol = 0.001;
kp_min = 0.001;
kp_max = 10;
kp_x = kp_for_bandwidth(G_x,wb,wb_tol,kp_min,kp_max);

% Transfer function inclucding P controller
D_p = kp_x;
G_x_cl = D_p*G_x/(1 + D_p*G_x); % Closed loop tf with PID control for x
% G_x_cl = x/x_sp

% % Bode of closed loop plant with Kp
% figure;
% bode(sym2tf(G_theta_cl));
% title('G(theta) closed-loop with Kp for desired bandwidth');
% grid on;

% Root locus of plant with P controller
figure;
rlocus(sym2tf(G_x));
title('G(x) with P controller root locus varied by kp')
hold on;

% Plot current poles for kp needed for bandwidth
current_pole = rlocus(sym2tf(kp_x*G_x), kp_x);
plot(real(current_pole), imag(current_pole), 'rs', 'Markersize', 7); % Plot current pole locatiosn

% Settling time:
sigma_p = log(0.02)/ts; % Real part limit of dominant pole, p for pole to avoid confusion with noise

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit

% Step responce
t_dist = 10;
controller_step_responce(G_x, D_p, {'P'}, t_dist)
title('Step responce of x controller')

% Performance parameters
sys = sym2tf(G_x_cl);
x_performance = stepinfo(sys);
wb = bandwidth(sys);
x_performance.Bandwidth = wb;
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
x_performance.slow_factor = slow_factor;
x_performance


%%-------------------------------------------------------------
%% Z Velocity controller
%%-------------------------------------------------------------

%% Design requirements:
close all
% Reject disturbances
% Zero steady-state error
PO = 12; % Percentage Overshoot (%)
% ??? correct other controllers to use performance struct for wb_inner
wb_inner = motor_wb;
wb = 2.9; % Desired dandwidth (rad/s). Slower than previous wb by a factor
ts = 11.6; % 2% settling time (s)

% Plant Transfer Function

% Linearise at: theta_sp = 0 -> theta = delta_E = 0

% F_z_r to delta_T:
% ------------------
% delta_T = -0.5*F_z_r/cos(theta);
% Small angle approx: cos(theta) == 1
% Therefore: F_z = -delta_T

syms F_z_r
delta_T = -0.5*F_z_r;

% delta_T to dz:
% ------------
% [T1; T2] = MM*[delta_E; delta_T]
% MM = [-1, 1;
%        1, 1];
% delta_E = 0
% Therefore: T1 = delta_T, T2 = delta_T
T1_r = delta_T; % Thrust reference 1
T2_r = delta_T; % Thrust reference 2

% Motor model
G_motor = (1/tau)/(s + 1/tau);
T1 = T1_r*G_motor; % Actual thrust
T2 = T2_r*G_motor;

% theta = 0, therefore F_z = sum of motor thrusts
F_z = - T1 - T2; % Applied force in z direction. Z is down, but T1,T2 positive up

syms dz
C_z_lin = C_Dz*dz_bar; % Linearised drag coef. at average velocity dz_bar

% Linearise at hover, therefore ignore gravity
eqn = (F_z - C_z_lin*rho*dz == s*M*dz); % Newton 2nd law in z direction
dz = solve(eqn, dz); % Solve for dz in terms of F_z_r
G_dz = dz/F_z_r; % TF from F_z_r to dz
G_dz_tf = sym2tf(G_dz);

% PI controller for dz:

% Percentage overshoot:
zeta = sqrt( (log(PO/100))^2 / (pi^2 + (log(PO/100))^2) );  % Damping ratio
theta_pole = atan(sqrt(1 - zeta^2) / zeta); % Max angle from real axis to dominant pole

% Settling time:
sigma_p = log(0.02)/ts; % Real part limit of dominant pole, p for pole to avoid confusion with noise

% Now implement PI controller
% D_PI = Kp*(s + z_c) / s
% Use I controller to reject disturbances (closed loop becomes type 2)
% Use P controller to place in performance envelope

% Pole of D_PI is at origin, so let zero be close to origin: z_c = 0.1
z_c = 1.4; % z_c = ki/kp
D_pi = (s + z_c) / s; % transfer function of Pi controller without kp

% Place kp for bandwidth
kp_min = 0.001;
kp_max = 30;
wb_tol = 0.001; % tolerance on bandwidth frequency
kp_dz = kp_for_bandwidth(D_pi*G_dz,wb,wb_tol,kp_min,kp_max);

% Calculate ki from z_c and kp
ki_dz = kp_dz*z_c;
D_pi = kp_dz + ki_dz*(1/s); % PI controller TF

% Root locus of plant with PI controller
figure;
rlocus(sym2tf(D_pi*G_dz));
title('G(dz) with PI controller root locus varied by kp')
hold on;

% Plot current poles for kp needed for bandwidth
current_pole = rlocus(sym2tf(D_pi*G_dz), kp_dz);
plot(real(current_pole), imag(current_pole), 'rs', 'Markersize', 7); % Plot current pole locatiosn

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit
x_theta = max(ylim)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(ylim), '--');

% PID controller for dz:

N_dz = wb*2; % Frequency of Low Pas Filter (Manually tune)

% Manual root locus plot of closed loop system with PID controller
syms kd_dz
% D_pid = @(kd_dz) kp_dz + ki_dz*(1/s) + kd_dz*s; % Without LPF
D_pid = @(kd_dz) kp_dz + ki_dz*(1/s) + kd_dz*(N_dz/(1 + N_dz/s)); % Controller TF including low pass filter
G_dz_cl = @(kd_dz) D_pid(kd_dz)*G_dz/(1 + D_pid(kd_dz)*G_dz); % Closed loop tf with P control for dtheta

figure;
title('G(dz) with PID controller root locus varied by kp')
hold on;
% for kd_dz = 0:0.05:1.1
%     poles = pole(sym2tf(G_dz_cl(kd_dz)));
%     plot3(real(poles), imag(poles), kd_dz*(real(poles)./real(poles)), 'k.'); % Plot pole of current k
% end

% Starting poles
poles = pole(sym2tf(G_dz_cl(0)));
plot(real(poles), imag(poles), 'bx'); % Plot pole of current k

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit
x_theta = max(ylim)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(ylim), '--');

% Choose kd to place poles within desired region
kd_dz = 0; % Manually adjust

% Insert final parameter into TFs
D_pid = D_pid(kd_dz);
G_dz_cl = simplifyFraction(G_dz_cl(kd_dz));

% Current poles:
poles = pole(sym2tf(G_dz_cl));
plot3(real(poles), imag(poles), kd_dz*(real(poles)./real(poles)), 'rs', 'MarkerSize', 7); % Plot pole of current k

% Step responce
t_dist = 10;
controller_step_responce(G_dz, [D_pid, D_pi, kp_dz], {'PID', 'PI', 'P'}, t_dist)
title('Step responce of dz controllers')

% Performance parameters
sys = sym2tf(G_dz_cl);
dz_performance = stepinfo(sys);
wb = bandwidth(sys);
dz_performance.Bandwidth = wb;
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
dz_performance.slow_factor = slow_factor;
dz_performance


%%-------------------------------------------------------------
%% Z Position controller
%%-------------------------------------------------------------

%% Design requirements:
% Zero steady-state error
% Overdamped
% Timescale seperation from inner loop

PO = 0; % Percentage Overshoot (%)
wb_inner = dz_performance.Bandwidth;
wb = 1; % Desired dandwidth (rad/s).
ts = 11.51; % 2% settling time (s)

% Plant Transfer Function
% TF from z_sp to z (seen by position controller)
% z = (1/s)*dz;
% dz = G_dz_cl*z_sp;
G_z = simplifyFraction(G_dz_cl*(1/s));
G_z_tf = sym2tf(G_z);

% Calculate Kp needed for desired bandwidth (Binary search)
wb_tol = 0.001;
kp_min = 0.001;
kp_max = 10;
kp_z = kp_for_bandwidth(G_z,wb,wb_tol,kp_min,kp_max);

% Transfer function inclucding P controller
D_p = kp_z;
G_z_cl = D_p*G_z/(1 + D_p*G_z); % Closed loop tf with PID control for z
% G_z_cl = z/z_sp

% Bode of closed loop plant with Kp
% figure;
% bode(sym2tf(G_theta_cl));
% title('G(theta) closed-loop with Kp for desired bandwidth');
% grid on;

% Root locus of plant with P controller
figure;
rlocus(sym2tf(G_z));
title('G(z) with P controller root locus varied by kp')
hold on;

% Plot current poles for kp needed for bandwidth
current_pole = rlocus(sym2tf(kp_z*G_z), kp_z);
plot(real(current_pole), imag(current_pole), 'rs', 'Markersize', 7); % Plot current pole locatiosn

% Settling time:
sigma_p = log(0.02)/ts; % Real part limit of dominant pole, p for pole to avoid confusion with noise

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit

% Step responce
t_dist = 10;
controller_step_responce(G_z, D_p, {'P'}, t_dist)
title('Step responce of z controller')

% Performance parameters
sys = sym2tf(G_z_cl);
z_performance = stepinfo(sys);
wb = bandwidth(sys);
z_performance.Bandwidth = wb;
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
z_performance.slow_factor = slow_factor;
z_performance

%% Save gain values
description = 'PID gain values and Linearised TF models from PID_design_setup.m for drone 2D';
save('Data/Drone_2D_control_params.mat', 'description', ...
'kp_dtheta', 'ki_dtheta', 'kd_dtheta', 'N_dtheta', ...
'kp_theta', ...
'kp_dx', 'ki_dx', 'kd_dx', 'N_dx', ...
'kp_x', ...
'kp_dz', 'ki_dz', 'kd_dz', 'N_dz', ...
'kp_z', ...
'G_dtheta_tf', 'G_theta_tf', 'G_dx_tf', 'G_x_tf', 'G_dz_tf', 'G_z_tf') % Linearised TF used to compare linear models with plant

function TF = sym2tf(sym_TF)
% Converts symbolic representation of transfer function
% to a transfer function object
% from: https://www.mathworks.com/matlabcentral/answers/310042-how-to-convert-symbolic-expressions-to-transfer-functions

    ExpFun = matlabFunction(sym_TF);
    ExpFun = str2func(regexprep(func2str(ExpFun), '\.([/^\\*])', '$1'));
    TF = ExpFun(tf('s'));
end

function kp = kp_for_bandwidth(G,wb,wb_tol,kp_min,kp_max)
% Find proportional feedback needed for specific bandwidth of transfer function
% Using Binary search
% kp = proportional feedback gain
% G = open loop transfer function in symbolic format, not tf object
% wb = desired bandwidth
% wb_tol = allowable tolerance on bandwidth
% kp_min = minimum kp
% kp_min = minimum kp

kp_min_param = kp_min; % Doesnt chnage during function
kp_max_param = kp_max;

    max_iterations = 200;
    kp_found = 0;
    for i = 1:max_iterations
        kp_mid = mean([kp_max,kp_min]);
        sys = sym2tf(kp_mid*G/(1 + kp_mid*G));
        wb_actual = bandwidth(sys); % Final bandwidth
        if abs(wb - wb_actual) < wb_tol
            kp_found = 1;
            break;
        elseif wb_actual > wb
            kp_max = kp_mid;
        else
            kp_min = kp_mid;
        end
    end
    kp = kp_mid;
    
    % If still not found, increment until wb_actual > wb
    if ~kp_found
        for kp = kp_min_param:0.01:kp_max_param
            sys = sym2tf(kp*G/(1 + kp*G));
            wb_actual = bandwidth(sys); % Final bandwidth
            if wb_actual>wb
                kp_found = 1;
                break;
            end
        end
        if ~kp_found
            error('kp_dtheta not found. Try changing search space')
        end
    end

    
end

function controller_step_responce(G, D_array, Legend, t_dist)
%% Step responce of P, Pi, and PID controllers with disturbance at t_dist
% G = symbolic tf of plant
% D_array = Array of controllers to simulate (symbolic tf of controllers)
% Legend = cell array of entries to use as legend
% t_dist = Time of step disturbance (s)
    

    Ts = 0.01;
    t = 0:Ts:10;
    u = ones(1,length(t)); % Input time series    
    u(:, t_dist/Ts:end) = 2; % Add step disturbance at t = t_dist

    % Plot
    figure;
    hold on
    grid on
    
    % Simulate controllers
    for i = 1:length(D_array)
        sys = sym2tf(D_array(i)*G/(1 + D_array(i)*G));
        lsim(sys,u,t)    
    end
    
    legend(Legend)

end





