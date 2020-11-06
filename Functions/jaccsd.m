function J = jaccsd(f,x)
% JACCSD Jacobian through complex step differentiation
% By Yi Cao at Cranfield University, 02/01/2008
% f =   function to apply Jacobian to. 
%       Insert as fuction handle. i.e. @function_name
% J = f'(x,u) Jacobian matrix at x and u
% x = state vector at linearisation point

f_x = f(x);
n = numel(x);
m = numel(f_x);
J = zeros(m,n);
h = n*eps; % very small value
for k = 1:n
    x1 = x;
    x1(k) = x1(k)+ h*1i;
    J(:,k) = imag(f(x1))/h;
end

end