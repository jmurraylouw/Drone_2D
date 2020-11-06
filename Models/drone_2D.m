function dx = drone_2D(x,u,params)
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
    
    dx(4,1) = -(2*T1*sin(X3) + 2*T2*sin(X3) + C_Dx*X4*rho*abs(X4))/(2*M);
    dx(5,1) = -(2*M*g + 2*T1*cos(X3) + 2*T2*cos(X3) + C_Dz*X5*rho*abs(X5))/(2*M);
    dx(6,1) = -(T1*r - T2*r)/I_yy;
end