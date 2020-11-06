%% Derive system equations for 2D drone with suspended payload
% Two vertical forces at distance, r, from COM represent the rotor forces
% North East Down axis system. Therefore z is down
clear all

% Symbolic variables
syms M % Mass of drone body (at fulcrum)
syms I_yy % Moment of inertia of drone body about body x axis
syms r % Distance from each rotor force to COM of drone
syms g % Acceleration due to gravity (always negative)
syms C_Dx % Damping coef. of drone through air in x direction (f = C_Dx*xdot)
syms C_Dz % Damping coef. of drone in z direction (f = cy*zdot)
syms rho % Air density
syms t % Time

syms T1 % Rotor thrust on drone on left of COM
syms T2 % Rotor thrust on drone on right of COM

syms x(t) % x position of drone
syms z(t) % z position of drone
syms theta(t) % Pitch angle of drone (horisontal = 0 rad)

syms dx % dx/dt of drone 
syms dz % dz/dt of drone 
syms dtheta % dtheta/dt of drone

n = 6; % number of states
X = sym('X', [n, 1]); % State vector [x; z; theta; dx; dz; dtheta]
states = [x; z; theta]; % non-rate states

% Rates
dx        = diff(x, t);
dz        = diff(z, t);
dtheta    = diff(theta, t);

% Drone body equations
KE_M = 0.5*M*(dx^2 + dz^2) + 0.5*I_yy*dtheta^2; % Kinetic energy of drone body (linear + rotational)
PE_M = M*g*z; % Potential energy of drone body

% Lagrangian
L = KE_M - PE_M;
L = simplify(L);

% Non-conservative Forces 
% Time-constant model for motor thrusts are modelled in simulink
% Qx = -(T1 + T2)*sin(theta) - 0.5*rho*dx*abs(dx)*C_Dx; % Aerodynamic drag from A. Erasmus thesis, (Chp 3.45)
% Qz = -(T1 + T2)*cos(theta) - 0.5*rho*dz*abs(dz)*C_Dz; % NB: z is down

% Linearise at small angle:
% sin(theta) = theta
% cos(theta) = 1

Qx = -(T1 + T2)*(theta); % Aerodynamic drag from A. Erasmus thesis, (Chp 3.45)
Qz = -(T1 + T2)*1; % NB: z is down

% Non-conservative Torques
% ?? no aerodynamic drag on rotation?
Qtheta = T2*r - T1*r; % Torques caused be rotor forces

% Lagrangian equations
eq_x     = euler_lag(L, x, Qx, t); 
eq_z     = euler_lag(L, z, Qz, t);
eq_theta = euler_lag(L, theta, Qtheta, t);

% Clear symbol connections
syms dx  dz  dtheta
syms ddx ddz ddtheta
dstates  = [dx;  dz;  dtheta];
ddstates = [ddx;  ddz;  ddtheta];

eqns = [eq_x; eq_z; eq_theta]; % VEquations to solve with

% Substitute symbols into derivatives
old = [diff(states,t); diff(states,t,t)];
new = [dstates;        ddstates];
eqns = subs(eqns, old, new);

% Solve
solution = solve(eqns, ddstates);
ddstates = struct2cell(solution); % Convert to cell from struct
ddstates = [ddstates{:}]; % Convert to normal syms array from cell

% Simplify
ddstates = simplifyFraction(ddstates);

% Substitute state variables with y
old = [states; dstates];
new = X;
ddstates = subs(ddstates, old, new);

%% Display to copy into script
for i = 1:n/2
    fprintf('dx(%d,1) = %s;\n',(i + n/2),ddstates(i))
end

%% Display continuous state space

% Model parameters
M     = 4.5; % Mass of drone body (at fulcrum)
I_yy  = 0.235; % Moment of inertia of drone body about body x axis
r     = 0.49*1/sqrt(2); % Distance from each rotor force to COM of drone
g     = -9.81; % Acceleration due to gravity (always negative)
C_Dx  = 0.2 ; % Damping coef. of drone through air in x direction
C_Dz  = 0.2; % Damping coef. of drone in z direction
rho   = 1.225; % Air density (kg/m^3)

m     = 2; % Mass of swinging payload (kg)
l     = 1; % Length of pendulum (m)
cbeta = 0.05; % Rotational damping coef of payload at connection

C_px = 0.05; % Damping coef. of drone through air in x direction
C_pz = 0.05; % Damping coef. of drone through air in z direction

params = [M; I_yy; r; g; C_Dx; C_Dz; rho; m; l; cbeta; C_px; C_pz]; % Vector of all parameters


% Motor constant
tau   = 0.07; % Motor time constant

% Mixing Matrix
MM = [-1, 1;
       1, 1]; % [T1; T2] = MM*[delta_E; delta_T]

% System dimensions
nx = 6;
nu = 2;

% Linearisation point / Hover condition
x_bar = zeros(nx,1);
u_bar = -0.5*(M)*g*[1; 1];

% Jacobians
f_x = @(x) lin_drone_nl(x, u_bar, params); % Function of only x, u is constant
A_lin = jaccsd(f_x,x_bar) % System matrix from Jacobian at hover condition

f_u = @(u) lin_drone_nl(x_bar, u, params); % Function of only u, x is constant
B_lin = jaccsd(f_u, u_bar) % Input matrix from Jacobian at hover condition


% %% Display pretty equations
% disp('')
% disp('Pretty')
% disp('------')
% for i = 1:n/2
%     disp(i + n/2)
%     pretty(ddstates(i))
%     disp("-------------------------------")
% end

function dx = lin_drone_nl(x,u,params) % linearised drone no payload
    %% Differential equations for drone in 2D without a payload
    % dx: x_dot, time derivative of state vector
    % x: state vector [x, z, theta]
    % u: input vector
    % params: parameters for drone model
    
    % Parameters:
    M     = params(1); % Mass of drone body (at fulcrum)
    I_yy  = params(2); % Moment of inertia of drone body
    r     = params(3); % Distance from each rot`or force to COM of drone
    g     = params(4); % Acceleration due to gravity (always negative)
    C_Dx  = params(5); % Damping coef. of drone through air in x direction (f = cx*xdot)
    C_Dz  = params(6); % Damping coef. of drone in z direction (f = cy*zdot)
    rho   = params(7); % Air density
    
    % States (extract notation so we can copy over from symbolic toolbox equations)
    X1 = x(1); % x
    X2 = x(2); % z
    X3 = x(3); % theta
    X4 = x(4); % dx
    X5 = x(5); % dz
    X6 = x(6); % dtheta
    
    % Inputs
    T1 = u(1);
    T2 = u(2);
    
    % System equations
    dx = zeros(6,1);
    
    dx(1,1) = x(4);
    dx(2,1) = x(5);
    dx(3,1) = x(6);
    
    c_vel = 0.027; % Damping according to velocity, not velocity squared, for linear model
    
    dx(4,1) = -(T1*X3 + T2*X3 + c_vel*X4)/(M);
    dx(5,1) = -(T1 + T2 + M*g + c_vel*X5)/(M);
    dx(6,1) = -(T1*r - T2*r)/I_yy;
end

function dx = nlin_drone_nl(x,u,params)
    %% Differential equations for drone in 2D WITHOUT a payload
    % dx: x_dot, time derivative of state vector
    % x: state vector [x, z, theta]
    % u: input vector
    % params: parameters for drone model
    
    % Parameters:
    M     = params(1); % Mass of drone body (at fulcrum)
    I_yy  = params(2); % Moment of inertia of drone body
    r     = params(3); % Distance from each rot`or force to COM of drone
    g     = params(4); % Acceleration due to gravity (always negative)
    C_Dx  = params(5); % Damping coef. of drone through air in x direction (f = cx*xdot)
    C_Dz  = params(6); % Damping coef. of drone in z direction (f = cy*zdot)
    rho   = params(7); % Air density
    
    % States (extract notation so we can copy over from symbolic toolbox equations)
    X1 = x(1); % x
    X2 = x(2); % z
    X3 = x(3); % theta
    X4 = x(4); % dx
    X5 = x(5); % dz
    X6 = x(6); % dtheta
    
    % Inputs
    T1 = u(1);
    T2 = u(2);
    
    % System equations
    dx = zeros(6,1);
    
    dx(1,1) = x(4);
    dx(2,1) = x(5);
    dx(3,1) = x(6);
    
    dx(4,1) = -(2*T1*sin(X3) + 2*T2*sin(X3) + C_Dx*rho*(X4)/(2*M)); %% Removed abs(X4 becasue non-differentiable)
    dx(5,1) = -(2*M*g + 2*T1*cos(X3) + 2*T2*cos(X3) + C_Dz*X5*rho)/(2*M);
    dx(6,1) = -(T1*r - T2*r)/I_yy;
end
