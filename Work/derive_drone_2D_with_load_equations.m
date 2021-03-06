%% Derive system equations for 2D drone with suspended payload
% Two vertical forces at distance, r, from COM represent the rotor forces
% North East Down axis system. Therefore z is down
clear all

% Symbolic variables
syms M % Mass of drone body (at fulcrum)
syms I_yy % Moment of inertia of drone body about body x axis
syms r % Distance from each rotor force to COM of drone
syms g % Acceleration due to gravity (always negative)
syms C_Dx % Damping coef. of drone through air in x direction (f = 0.5*rho*dx*abs(dx)*C_Dx)
syms C_Dz % Damping coef. of drone in z direction
syms rho % Air density
syms C_px % Damping coef. of paylaod through air in x direction (f = 0.5*rho*dx*abs(dx)*C_Dx)
syms C_pz % Damping coef. of payload in z direction
syms t % Time

syms m % Mass of swinging payload
syms l % Length of pendulum
syms cbeta % Rotational damping coef of payload at connection

syms T1 % Rotor thrust on drone on left of COM
syms T2 % Rotor thrust on drone on right of COM

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
KE_M = 0.5*M*(dx^2 + dz^2) + 0.5*I_yy*dtheta^2; % Kinetic energy of drone body (linear + rotational)
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
% ???? Add damping of payload
% Time-constant model for motor thrusts are modelled in simulink
Qx = -(T1 + T2)*sin(theta) - 0.5*rho*dx*abs(dx)*C_Dx; % Aerodynamic drag from A. Erasmus thesis, (Chp 3.45)
Qz = -(T1 + T2)*cos(theta) - 0.5*rho*dz*abs(dz)*C_Dz; % NB: z is down

% Non-conservative Torques
% ?? no aerodynamic drag on rotation?
Qtheta = T2*r - T1*r; % Torques caused by rotor forces

F_wpx = -0.5*C_px*rho*diff(x_m,t)*abs(diff(x_m,t)); % Force due to drag of wind/air on payload in x direction
F_wpz = -0.5*C_pz*rho*diff(z_m,t)*abs(diff(z_m,t)); % Force due to drag of wind/air on payload in z direction
Qbeta  = -cbeta*dbeta + F_wpx*cos(beta)*l - F_wpz*sin(beta)*l; % Torques caused air damping on rotation of cable
% Qbeta = 0; % try no damping

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

eqns = [eq_x; eq_z; eq_theta; eq_beta]; % Equations to solve with

% Substitute symbols into derivatives
old = [diff(states,t); diff(states,t,t)];
new = [dstates;        ddstates];
eqns = subs(eqns, old, new);

% Solve
solution = solve(eqns, ddstates);
ddstates = struct2cell(solution); % Convert to cell from struct
ddstates = [ddstates{:}]; % Convert to normal syms array from cell

% Simplify
ddstates = simplify(ddstates);
ddstates = simplifyFraction(ddstates);

% Substitute state variables with y
old = [states; dstates];
new = X;
ddstates = subs(ddstates, old, new);

%% Display to copy into script
for i = 1:n/2
    fprintf('dx(%d,1) = %s;\n',(i + n/2),ddstates(i))
end

%% Display pretty equations
disp('Pretty')
disp('------')
for i = 1:n/2
    disp(i + n/2)
    pretty(ddstates(i))
    disp("-------------------------------")
end

