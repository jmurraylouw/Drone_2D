load('MPC/toy_data.mat') % Simulation data of toy problem
load('MPC/toy_model.mat')

rng('default');
rng(1); % Repeatable random numbers
rand_sys = rss(3,3,2);
[A_plant, B_plant, ~, ~] = ssdata(rand_sys); % Continuous plant model

C_plant = eye(3);
D_plant = zeros(3,2);
x0 = zeros(3,1);

%% Implentation of HAVOK

% Extract data
u_data  = out.u.Data';
x_data  = out.y.Data';
y_data  = out.y.Data';
t       = out.tout'; % Time

% Adjust for constant disturbance / mean control values
% u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; -4.5*9.81];
% u_data  = u_data - u_bar; % Adjust for unmeasured input
u_bar = [0; 0];
u_data  = u_data - u_bar; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 2000; % Num of data samples for testing
x_test = x_data(:,end-N_test+1:end);
y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Data dimentions
nx = size(x_data,1); % number of states
ny = size(y_data,1); % number of measurements
nu = size(u_data,1); % number of inputs
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

% Add noise
rng('default');
rng(1); % Repeatable random numbers
sigma = 0.001; % Noise standard deviation
y_data_noise = y_data + sigma*randn(size(y_data));

% Training data - Last sample of training is first sample of testing
N_train = 5000; % Number of sampels in training data
y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
t_train = t(:,end-N_test-N_train+2:end-N_test+1);

% Parameters
q = 3; % Override
p = 1; % Override

w = N_train - q + 1; % num columns of Hankel matrix
D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)

% Create Hankel matrix with measurements
Y = zeros(q*ny,w); % Augmented state with delay coordinates [...; Y(k-2); Y(k-1); Y(k)]
for row = 0:q-1 % Add delay coordinates
    Y(row*ny+1:(row+1)*ny, :) = y_train(:, row + (0:w-1) + 1);
end

Upsilon = u_train(:, q:end); % Leave out last time step to match V_til_1
YU_bar = [Y; Upsilon];

% % SVD of the Hankel matrix
% [U1,S1,V1] = svd(YU_bar, 'econ');
% figure, semilogy(diag(S1), 'x'), hold on;
% title('Singular values of Omega, showing p truncation')
% plot(p,S1(p,p), 'ro'), hold off;
% 
% % Truncate SVD matrixes
% U_tilde = U1(:, 1:p); 
% S_tilde = S1(1:p, 1:p);
% V_tilde = V1(:, 1:p);
% 
% % Setup V2 one timestep into future from V1
% V_til_2 = V_tilde(2:end  , :)'; % Turnd on side (wide short matrix)
% V_til_1 = V_tilde(1:end-1, :)';
% 
% % DMD on V
% AB_tilde = V_til_2*pinv(V_til_1); % combined A and B matrix, side by side
% AB_tilde = stabilise(AB_tilde,3);
% 
% % convert to x coordinates
% AB_bar = (U_tilde*S_tilde)*AB_tilde*pinv(U_tilde*S_tilde);
% A_bar = AB_bar(1:q*m, 1:q*m);
% B_bar = AB_bar(1:q*m, q*m+1:end);
% % A_bar = stabilise(A_bar,10);

% DMD of Y
Y2 = Y(:, 2:end  );
Y1 = Y(:, 1:end-1);

YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y and U, above and below
AB = Y2*pinv(YU); % combined A and B matrix, side by side
A_hat  = AB(:,1:q*ny); % Extract A matrix
B_hat  = AB(:,(q*ny+1):end);

A_hat; % Estimated discrete A matrix of system
B_hat; % Estimated discrete B matrix of system

n_hat = size(A_hat, 1);
l_hat = nu;
m_hat = n_hat;

C_hat = eye(m_hat);
D_hat = zeros(m_hat, l_hat);

x0_hat = zeros(n_hat,1);

%% LTI system
Ts_mpc = 0.1; % MPC sampling time

dmd_sys = ss(A_hat,B_hat,eye(q*ny),zeros(q*ny,nu),Ts); % LTI system
dmd_sys = d2d(dmd_sys,Ts_mpc,'zoh'); % Resample to match MPC

[A_hat2,B_hat2,C_hat2,D_hat2,Ts] = ssdata(dmd_sys); % Extract resampled matrixes
dmd_sys = ss(A_hat,B_hat,C_hat2,D_hat2,Ts_mpc); % LTI system with new Ts 

%% MPC object
old_status = mpcverbosity('off');
dmd_sys.InputGroup.MV = 2;
dmd_sys.OutputGroup.MO = 6;
MPCobj = mpc(dmd_sys,Ts_mpc);

% Set weights of delay coordinates to 0
MPCobj.Weights.OutputVariables = [zeros(1,(q-1)*ny), 1, 1, 0];
%% Nominal operating conditions for AMPC block
U_nom = u_bar;
X_nom = zeros(n_hat,1);
Y_nom = zeros(m_hat,1);
DX_nom = zeros(n_hat,1);

%% Way points
num_waypoints = 100; % Number of waypoints included in command
point_time_interval = 6; % Initial interval between commands

waypoints = table('Size', [(num_waypoints+1)*2, 4], 'VariableTypes', ["double", "double", "double", "double"]);
waypoints.Properties.VariableNames = {'point_time', 'x_coord', 'z_coord' , 'th_coord'};

x_coord = 0;
z_coord = -5;
th_coord = 0;

waypoints(1,:) = table(0,                   x_coord, z_coord, th_coord); % Initial point
waypoints(2,:) = table(point_time_interval, x_coord, z_coord, th_coord); % Initial point for 6 seconds

x_min        = -10;     x_max        = 10; % (m) minimum and maximum coordinates for waypoints
z_min        = -25;     z_max        = -5;
interval_min = 2;       interval_max = 8;  % (s) minimum and maximum time interval between commands

point_time = point_time_interval;
rng(0); % Initialise random number generator for repeatability
for i = 1:num_waypoints
    point_time_interval = (interval_max - interval_min).*rand() + interval_min; % (s) random time interval between commands
    point_time = point_time + point_time_interval;

    waypoints(2*i,  :) = table(point_time, x_coord, z_coord, th_coord); % Previous point    
    x_coord    = (x_max - x_min).*rand() + x_min; % x coordinate of next waypoint
    z_coord    = (z_max - z_min).*rand() + z_min; % z coordinate of next waypoint   
    th_coord    = (x_max - x_min).*rand() + x_min; % theta coordinate of next waypoint
    
    waypoints(2*i+1,:) = table(point_time, x_coord, z_coord, th_coord); % Next point
end
i = i+1;
waypoints(2*i,  :) = table(point_time+interval_max, x_coord, z_coord, th_coord); % Add time to reach final point

% Add zero reference for delay coordinates
waypoints_ts = timeseries([zeros(size(waypoints.x_coord,1), (q-1)*ny), waypoints.x_coord, waypoints.z_coord, waypoints.th_coord], waypoints.point_time); % timeseries object for From Workspace block

%% Initital conditions for extended state vector
% All delay states are also at x0
x_ext_0 = zeros(q*ny, 1); % Allocate space
for row = 0:q-1 % First column of spaced Hankel matrix
        x_ext_0(row*ny+1:(row+1)*ny, 1) = x0;
end

%% Run with A and x
run_and_plot = 0;
if run_and_plot
    
    % Initial condition
    y_hat_0 = zeros(q*ny,1);
    for row = 0:q-1 % First column of spaced Hankel matrix
        y_hat_0(row*ny+1:(row+1)*ny, 1) = y_train(:, end - ((q-1)+1) + row + 1);
    end

    % Run model
    Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
    Y_hat(:,1) = y_hat_0; % Initial condition
    for k = 1:N_test-1
        Y_hat(:,k+1) = A_hat*Y_hat(:,k) + B_hat*u_test(:,k);
    end

    y_hat = Y_hat(end-ny+1:end, :); % Extract only non-delay time series (last m rows)

    % Vector of Mean Absolute Error on testing data
    MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state

    %% Plot data vs model
    figure;
    plot(t_train, y_train);
    hold on;
    plot(t_test, y_test);

    plot(t_test, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
    plot((D + t(N-N_test-N_train)).*[1, 1], ylim, 'r');
    plot(t(N-N_test-N_train).*[1,1], ylim, 'k');
    plot(t(N-N_test).*[1,1], ylim, 'k');
    title('Training and Testing data vs Model');
    % legend('','','','','','', 'x hat','z hat','theta hat', 'x hat bar','z hat bar','theta hat bar');
    hold off;
end % run_and_plot


%% Save data
save('MPC/toy_model.mat', 'u_input', 'waypoints_ts', ...
'A_hat',   'B_hat',   'C_hat',   'D_hat',   'x0_hat', ...
'A_plant', 'B_plant', 'C_plant', 'D_plant', 'x0', ...
'MPCobj')

function A = stabilise(A_unstable,max_iterations)
    % If some eigenvalues are unstable due to machine tolerance,
    % Scale them to be stable
    A = A_unstable;
    count = 0;
    while (sum(abs(eig(A)) > 1) ~= 0)       
        [Ve,De] = eig(A);
        unstable = abs(De)>1; % indexes of unstable eigenvalues
        De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-14 + count*2); % Normalize all unstable eigenvalues (set abs(eig) = 1)
        A = Ve*De/(Ve); % New A with margininally stable eigenvalues
        A = real(A);
        count = count+1;
        if(count > max_iterations)
            break
        end
    end

end






