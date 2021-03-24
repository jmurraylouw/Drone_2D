
Ts_mpc = 0.03;
R = zeros(6,6);

for i = 1:5
    R(i+1,1) = Ts_mpc;
    T_y = 6;
    for j = 1:5
        R(1,j+1) = T_y;
        R(i+1,j+1) = T_y/Ts_mpc;
        T_y = T_y + 1;
    end
    Ts_mpc = Ts_mpc + 0.01;

end
R