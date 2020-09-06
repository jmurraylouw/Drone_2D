x0 = zeros(6,1); % Initial conditions
x0(3) = 0.1;

MM = [-1, -1;
       1, -1]; % Mixing Matrix. [F1; F2] = MM*[delta_E; delta_T]