function dx = drone_2D(x,u,params)
    %% Differential equations for drone in 2D without a payload
    % dx: x_dot, time derivative of state vector
    % x: state vector [x, z, theta]
    % u: input vector
    % params: parameters for drone model
    
    % Parameters:
    M = params(1); % Mass of drone body (at fulcrum)
    I = params(2); % Moment of inertia of drone body
    r = params(3); % Distance from each rotor force to COM of drone
    g = params(4); % Acceleration due to gravity (always negative)
    cx = params(5); % Damping coef. of drone through air in x direction (f = cx*xdot)
    cz = params(6); % Damping coef. of drone in z direction (f = cy*zdot)
    ctheta = params(7); % Damping coef. of drone in theta direction (f = ct*thetadot)
    
    % Once-off calculations to speed up computation 
    C3 = cos(x(3)); % cos(theta)
    S3 = sin(x(3)); % sin(theta)

    % System equations
    dx = zeros(6,1);
    
    dx(1,1) = x(4);
    dx(2,1) = x(5);
    dx(3,1) = x(6);
    
    dx(4,1) = -(x(4)*cx + u(1)*S3 + u(2)*S3)/M;
    dx(5,1) = -(M*g + x(5)*cz + u(1)*C3 + u(2)*C3)/M;
    dx(6,1) = -(u(1)*r - u(2)*r + x(6)*ctheta)/I;

end