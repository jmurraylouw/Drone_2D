%% Design of cascaded loop controllers for 2D drone

close all

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

%%-------------------------------------------------------------
%% Angular rate controller
%%-------------------------------------------------------------

%% Plant Transfer Function
% Transfer function from delta_E to dtheta
G_dtheta = dtheta/delta_E;
G_dtheta = subs(G_dtheta); % Substitute paramater values

%% Design requirements:

% Reject disturbances
% Zero steady-state error
PO = 4.2; % Percentage Overshoot (%)
wb_buffer = 1.9; % Add buffer because wb drops when kd and LPF is added
wb = 12.52 + wb_buffer; % Desired dandwidth (rad/s)
ts = 0.6; % 2% settling time (s)

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

% Calculate Kp needed for desired bandwidth
kp_min = 0.001;
kp_max = 10;
wb_tol = 0.001; % tolerance on bandwidth frequency
kp_dtheta = kp_for_bandwidth(D_pi*G_dtheta,wb,wb_tol,kp_min,kp_max);
D_pi = kp_dtheta*D_pi;

% Calculate ki from z_c and kp
ki_dtheta = kp_dtheta*z_c;

% Draw root locus
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

%% PID controller for dtheta:

N_dx = 2*wb; % Frequency of Low Pas Filter (Manually tune)

% Manual root locus plot of closed loop system with PID controller
syms kd_dx
% D_pid = @(kd_dx) kp_dx + ki_dx*(1/s) + kd_dx*s; % Without LPF
D_pid = @(kd_dx) kp_dx + ki_dx*(1/s) + kd_dx*s*(N_dx/(1 + N_dx)); % Controller TF including low pass filter
G_dx_cl = @(kd_dx) D_pid(kd_dx)*G_dx/(1 + D_pid(kd_dx)*G_dx); % Closed loop tf with P control for dtheta

% Manual root locus plot of closed loop system with PID controller
syms kd_dtheta
D_pid = @(kd_dtheta) kp_dtheta + ki_dtheta*(1/s) + kd_dtheta*s*(N_dtheta/(1 + N_dtheta));
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
dtheta_performance

%%-------------------------------------------------------------
%% Angle controller
%%-------------------------------------------------------------

%% Design requirements:
% Zero steady-state error
% Overdamped
% Timescale seperation from inner loop

PO = 0; % Percentage Overshoot (%)
wb_inner = wb;
wb = 4.41; % Desired dandwidth (rad/s).
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
ts = 1.95; % 2% settling time (s)

%% Plant Transfer Function
% TF from dtehta_sp to theta (seen by angular rate controller)
% theta = (1/s)*dtheta;
% dtheta = G_dtheta_cl*theta_sp;
G_theta = simplifyFraction(G_dtheta_cl*(1/s));

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
theta_performance.slow_factor = slow_factor;
theta_performance

%%-------------------------------------------------------------
%% X Velocity controller
%%-------------------------------------------------------------

%% Design requirements:

% Reject disturbances
% Zero steady-state error
PO = 12; % Percentage Overshoot (%)
wb_inner = wb;
wb = 2.166; % Desired dandwidth (rad/s). Slower than previous wb by a factor
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
ts = 11.6; % 2% settling time (s)

%% Plant Transfer Function

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

%% PI controller for dx:

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

%% PID controller for dx:

% ???? Low pass filter (need to decide how to pick filter value
N_dx = wb*2; % Frequency of Low Pas Filter (Manually tune)

% Manual root locus plot of closed loop system with PID controller
syms kd_dx
% D_pid = @(kd_dx) kp_dx + ki_dx*(1/s) + kd_dx*s; % Without LPF
D_pid = @(kd_dx) kp_dx + ki_dx*(1/s) + kd_dx*(N_dx/(1 + N_dx/s)); % Controller TF including low pass filter
G_dx_cl = @(kd_dx) D_pid(kd_dx)*G_dx/(1 + D_pid(kd_dx)*G_dx); % Closed loop tf with P control for dtheta

figure;
hold on;
% for kd_dx = 0:0.01:2
%     poles = pole(sym2tf(G_dx_cl(kd_dx)));
%     plot3(real(poles), imag(poles), kd_dx*(real(poles)./real(poles)), 'k.'); % Plot pole of current k
% end

% Starting poles
poles = pole(sym2tf(G_dx_cl(0)));
plot(real(poles), imag(poles), 'bx'); % Plot pole of current k

% Plot requirement limits
plot([1, 1]*sigma_p, ylim, '--'); % Settling time requirement limit
x_theta = max(ylim)/tan(theta_pole); % x to plot theta line
plot([-1, 0, -1]*x_theta, [1, 0, -1]*max(ylim), '--');

% Choose kd to place poles within desired region
kd_dx = 0.01; % Manually adjust

% Insert final parameter into TFs
D_pid = D_pid(kd_dx);
G_dx_cl = G_dx_cl(kd_dx);

% Current poles:
poles = pole(sym2tf(G_dx_cl));
plot(real(poles), imag(poles), 'rs', 'MarkerSize', 7); % Plot pole of current k

% Step responce
t_dist = 10;
controller_step_responce(G_dx, [D_pid, D_pi, kp_dx], {'PID', 'PI', 'P'}, t_dist)
title('Step responce of dx controllers')

% Performance parameters
sys = sym2tf(G_dx_cl);
dx_performance = stepinfo(sys);
wb = bandwidth(sys);
dx_performance.Bandwidth = wb;
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
slow_factor = wb_inner/wb; % Factor that outer controller is slower than inner controller
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
title('Step responce of theta controllers')

% Performance parameters
sys = sym2tf(G_x_cl);
x_performance = stepinfo(sys);
wb = bandwidth(sys);
x_performance.Bandwidth = wb;
x_performance.slow_factor = slow_factor;
x_performance


%% Convert sym to tf objects
% G_dtheta = sym2tf(G_dtheta);
% G_th_dx = sym2tf(subs(G_th_dx));
% G_dx = sym2tf(subs(G_dx));

%% Save gain values
description = 'PID gain values from PID_design_setup.m for drone 2D';
save('Data/Drone_2D_control_params.mat', 'description', ...
'kp_dtheta', 'ki_dtheta', 'kd_dtheta', ...
'kp_theta', ...
'kp_dx', 'ki_dx', 'kd_dx')


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

    max_iterations = 1000;
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
    if ~kp_found
        error('kp_dtheta not found. Try changing search space')
    end
    kp = kp_mid;
    
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





