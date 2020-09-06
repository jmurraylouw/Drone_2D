function dx = drone_2D(x,u,params)
    %% Differential equations for drone in 2D without a payload
    % dx: x_dot, time derivative of state vector
    % x: state vector
    % u: input vector
    % params: parameters for drone model
    
    % Parameters:
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
    
    % Once-off calculations to speed up computation
    C4 = cos(x(4)); % cos(beta)
    S4 = sin(x(4)); % sin(beta)
    C4_2 = C4^2;
    S4_2 = S4^2;
    
    C3 = cos(x(3)); % cos(theta)
    S3 = sin(x(3)); % sin(theta)

    
    dx = zeros(8,1);
    
    dx(1,1) = x(5);
    dx(2,1) = x(6);
    dx(3,1) = x(7);
    dx(4,1) = x(8);
        
    dx(5,1) = ;
    dx(6,1) = ;
    dx(7,1) = ;
    dx(8,1) = ;

end