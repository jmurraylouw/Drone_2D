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
Qx = -(T1 + T2)*sin(theta) - 0.5*rho*dx*abs(dx)*C_Dx; % Aerodynamic drag from A. Erasmus thesis, (Chp 3.45)
Qz = -(T1 + T2)*cos(theta) - 0.5*rho*dz*abs(dz)*C_Dz; % NB: z is down

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

%% Display pretty equations
disp('')
disp('Pretty')
disp('------')
for i = 1:n/2
    disp(i + n/2)
    pretty(ddstates(i))
    disp("-------------------------------")
end


