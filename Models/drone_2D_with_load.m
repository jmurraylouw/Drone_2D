function dx = drone_2D_with_load(x,u,params)
    %% Differential equations for drone in 2D with suspended payload
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

    % Differential equations
    dx = zeros(8,1);
    
    dx(1,1) = x(5);
    dx(2,1) = x(6);
    dx(3,1) = x(7);
    dx(4,1) = x(8);
        
    dx(5,1) = -(x(5)*cx*l*m*C4_2 - x(8)*cbeta*m*C4 - M*x(8)*cbeta*C4 + F1*M*l*C4_2*S3 + F2*M*l*C4_2*S3 + F1*M*l*S3*S4_2 + F2*M*l*S3*S4_2 + F1*l*m*C4_2*S3 + F2*l*m*C4_2*S3 - M*x(8)^2*l^2*m*S4^3 + M*x(5)*cx*l*C4_2 + M*x(5)*cx*l*S4_2 - x(6)*cz*l*m*C4*S4 - M*x(8)^2*l^2*m*C4_2*S4 - F1*l*m*C3*C4*S4 - F2*l*m*C3*C4*S4)/((M + m)*(M*l*C4_2 + M*l*S4_2));
    dx(6,1) = -(M^2*g*l*S4_2 + M*x(8)*cbeta*S4 + x(8)*cbeta*m*S4 + M^2*g*l*C4_2 + M*g*l*m*C4_2 + M*g*l*m*S4_2 + x(6)*cz*l*m*S4_2 + F1*M*l*C3*C4_2 + F2*M*l*C3*C4_2 + F1*M*l*C3*S4_2 + F2*M*l*C3*S4_2 + F1*l*m*C3*S4_2 + F2*l*m*C3*S4_2 - M*x(8)^2*l^2*m*C4^3 + M*x(6)*cz*l*C4_2 + M*x(6)*cz*l*S4_2 - x(5)*cx*l*m*C4*S4 - M*x(8)^2*l^2*m*C4*S4_2 - F1*l*m*C4*S3*S4 - F2*l*m*C4*S3*S4)/((M + m)*(M*l*C4_2 + M*l*S4_2));
    dx(7,1) = -(F1*r - F2*r + x(7)*ctheta)/I;
    dx(8,1) = -(M*x(8)*cbeta + x(8)*cbeta*m + F1*l*m*C3*S4 - F1*l*m*C4*S3 + F2*l*m*C3*S4 - F2*l*m*C4*S3 - x(5)*cx*l*m*C4 + x(6)*cz*l*m*S4)/(M*m*l^2*C4_2 + M*m*l^2*S4_2);

end