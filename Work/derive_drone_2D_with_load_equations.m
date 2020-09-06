%% Derive system equations for 2D drone with suspended payload
% Two vertical forces at distance, r, from COM represent the rotor forces
% North East Down axis system. Therefore z is down
clear all

% Define symbolic variables
syms M % Mass of drone body (at fulcrum)
syms m % Mass of swinging payload
syms I % Moment of inertia of drone body
syms l % Length of pendulum
syms r % Distance from each rotor force to COM of drone
syms g % Acceleration due to gravity (always negative)
syms cx % Damping coef. of drone through air in x direction (f = cx*xdot)
syms cz % Damping coef. of drone in z direction (f = cy*zdot)
syms ctheta % Damping coef. of drone in theta direction (f = ct*thetadot)
syms cbeta % Damping coef. of drone in theta direction (f = ct*thetadot)
syms t % Time

syms F1 % Rotor force on drone on left of COM
syms F2 % Rotor force on drone on right of COM

syms x(t) % x position of drone
syms z(t) % z position of drone
syms theta(t) % Pitch angle of drone (horisontal = 0 rad)
syms beta(t) % Suspended angle of payload cable (vertical down = 0 rad)
syms dx % dx/dt of drone 
syms dz % dz/dt of drone 
syms dtheta % dtheta/dt of drone
syms dbeta % dbeta/dt of payload cable

n = 8; % number of states
X = sym('X', [n, 1]); % State vector [x; z; theta; beta; dx; dz; dtheta; dbeta]
states = [x; z; theta; beta]; % non-rate states

% Rates
dx        = diff(x, t);
dz        = diff(z, t);
dtheta    = diff(theta, t);
dbeta     = diff(beta, t);

% Drone body equations
KE_M = 0.5*M*(dx^2 + dz^2) + 0.5*I*dtheta^2; % Kinetic energy of drone body (linear + rotational)
PE_M = M*g*z; % Potential energy of drone body

% Payload equations
x_m = x + l*sin(beta); % x position of payload
z_m = z + l*cos(beta); % z position of payload

KE_m = 0.5*m*( diff(x_m,t)^2 + diff(z_m,t)^2 ); % Kinetic energy of payload
PE_m = m*g*z_m; % Potential energy of payload

% Lagrangian
L = (KE_M + KE_m) - (PE_M + PE_m);
L = simplify(L);

% Non-conservative Forces 
Qx = -(F1 + F2)*sin(theta) - cx*dx; % ?? change cx air damping according to pitch
Qz = -(F1 + F2)*cos(theta) - cz*dz; % NB: z is down

% Non-conservative Torques
Qtheta = F2*r - F1*r - ctheta*dtheta; % Torques caused be rotor forces and air damping
Qbeta  = -cbeta*dbeta; % Torques caused air damping on rotation of cable

% Lagrangian equations
eq_x     = euler_lag(L, x, Qx, t); 
eq_z     = euler_lag(L, z, Qz, t);
eq_theta = euler_lag(L, theta, Qtheta, t);
eq_beta  = euler_lag(L, beta,  Qbeta, t);

% Clear symbol connections
syms dx  dz  dtheta  dbeta
syms ddx ddz ddtheta ddbeta
dstates  = [dx;  dz;  dtheta;  dbeta];
ddstates = [ddx;  ddz;  ddtheta;  ddbeta];

eqns = [eq_x; eq_z; eq_theta; eq_beta]; % VEquations to solve with

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

%% Display to copy
for i = 1:n/2
    disp(i + n/2)
    (ddstates(i))
    disp("-------------------------------")
end

%% Display pretty equations
disp('Pretty')
disp('------')
for i = 1:n/2
    disp(i + n/2)
    pretty(ddstates(i))
    disp("-------------------------------")
end


