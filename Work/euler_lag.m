function eq = euler_lag(L, q, Q, t)
% Returns eq from lagragian equation in form eq == 0
% q is variable being derived with respect to
% Q is non-conservative forces
% Equation L has to be in a form with all derivatives as 
% "diff()" not in "dot"-form:

syms q_ qdot_

% Substitute time-dependant vars for non-time-dependant vars
% (to avoid dir of function)
L_ = subs(L, diff(q,t), qdot_); 
L_ = subs(L_, q, q_);

% Lagrange differentials
dL_dqdot    = diff(L_, qdot_); % Derivative of L with respect to qdot
dL_dq       = diff(L_, q_); % Derivative of L with respect to q

% Reinstate variables as time dependant before applying time derivative
dL_dqdot = subs(dL_dqdot, qdot_, diff(q,t));
dL_dqdot = subs(dL_dqdot, q_, q);

dL_dq = subs(dL_dq, qdot_, diff(q,t));
dL_dq = subs(dL_dq, q_, q);

% Lagrangian equation
eq = ( diff(dL_dqdot, t) - dL_dq == Q ); 

