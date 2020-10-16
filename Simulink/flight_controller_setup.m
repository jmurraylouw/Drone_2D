%% PID controllers
load('Data/Drone_2D_control_params.mat'); % Load controller gain values

% Parameters
% q = 90; % Override
% p = 33; % Override

% Dimensions
nx = 6; % Number of states
ny = 3; % Number of measurements
nu = 2; % Number of inputs

% Initial conditions
x0 = zeros(nx,1); % Initial state
y0 = zeros(ny,1); % Initial measurements
u0 = [0; 0]; % Initial input

% Initital conditions for extended measurment vector for MPC
% All delay states are also at x0
y_ext_0 = zeros(q*ny, 1); % Allocate space
for row = 0:q-1 % First column of spaced Hankel matrix
        y_ext_0(row*ny+1:(row+1)*ny, 1) = y0;
end

% Model parameters
M     = 4.5; % Mass of drone body (at fulcrum)
I_yy  = 0.235; % Moment of inertia of drone body about body x axis
r     = 0.49*1/sqrt(2); % Distance from each rotor force to COM of drone
g     = -9.81; % Acceleration due to gravity (always negative)
C_Dx  = 0.2 ; % Damping coef. of drone through air in x direction (f = C_Dx*xdot)
C_Dz  = 0.2; % Damping coef. of drone in z direction (f = cy*zdot)
rho   = 1.225; % Air density (kg/m^3)
tau   = 0.07; % Motor time constant

% Mixing Matrix
MM = [-1, 1;
       1, 1]; % [T1; T2] = MM*[delta_E; delta_T]

% Forgetting factor
lambda = 0.985; % Exponential forgetting factor of moving average to smooth out spikes in input for SVD performance

% Way points
num_waypoints = 100; % Number of waypoints included in command
point_time_interval = 6; % Initial interval between commands

waypoints = table('Size', [(num_waypoints+1)*2, 3], 'VariableTypes', ["double", "double", "double"]);
waypoints.Properties.VariableNames = {'point_time', 'x_coord', 'z_coord'};

x_coord = 0;
z_coord = -5;
waypoints(1,:) = table(0,                   x_coord, z_coord); % Initial point
waypoints(2,:) = table(point_time_interval, x_coord, z_coord); % Initial point for 6 seconds

x_min        = -10;     x_max        = 10; % (m) minimum and maximum coordinates for waypoints
z_min        = -25;     z_max        = -5;
interval_min = 2;       interval_max = 8;  % (s) minimum and maximum time interval between commands

point_time = point_time_interval;
rng(0); % Initialise random number generator for repeatability
for i = 1:num_waypoints
    point_time_interval = (interval_max - interval_min).*rand() + interval_min; % (s) random time interval between commands
    point_time = point_time + point_time_interval;

    waypoints(2*i,  :) = table(point_time, x_coord, z_coord); % Previous point    
    x_coord    = (x_max - x_min).*rand() + x_min; % x coordinate of next waypoint
    z_coord    = (z_max - z_min).*rand() + z_min; % z coordinate of next waypoint   
    waypoints(2*i+1,:) = table(point_time, x_coord, z_coord); % Next point
end
i = i+1;
waypoints(2*i,  :) = table(point_time+interval_max, x_coord, z_coord); % Add time to reach final point

waypoints_ts = timeseries([waypoints.x_coord, waypoints.z_coord], waypoints.point_time); % timeseries object for From Workspace block

% figure;
% title('Waypoint route')
% plot(waypoints.x_coord, waypoints.z_coord)
% xlabel('x');
% ylabel('z');

%% Implentation of HAVOK
% 
% Extract data
% u_data  = out.F_r.Data';
% x_data  = out.x.Data';
% y_data  = x_data(1:3,:);
% t       = out.tout'; % Time
% 
% Adjust for constant disturbance / mean control values
% % u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; M*g];
% u_data  = u_data - u_bar; % Adjust for unmeasured input
% 
% Testing data - Last 50 s is for testing and one sample overlaps training 
% N_test = 2000; % Num of data samples for testing
% x_test = x_data(:,end-N_test+1:end);
% y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
% u_test = u_data(:,end-N_test+1:end);
% t_test = t(:,end-N_test+1:end);
% 
% Data dimentions
% nx = size(x_data,1); % number of states
% ny = size(y_data,1); % number of measurements
% nu = size(u_data,1); % number of inputs
% Ts = t(2)-t(1);     % Sample time of data
% N  = length(t);     % Number of data samples
% 
% Add noise
% rng('default');
% rng(1); % Repeatable random numbers
% sigma = 0.001; % Noise standard deviation
% y_data_noise = y_data + sigma*randn(size(y_data));
% 
% Training data - Last sample of training is first sample of testing
% N_train = 5000; % Number of sampels in training data
% y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
% u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
% t_train = t(:,end-N_test-N_train+2:end-N_test+1);
% 
% 
% w = N_train - q + 1; % num columns of Hankel matrix
% D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)
% 
% Create Hankel matrix with measurements
% Y = zeros(q*ny,w); % Augmented state with delay coordinates [...; Y(k-2); Y(k-1); Y(k)]
% for row = 0:q-1 % Add delay coordinates
%     Y(row*ny+1:(row+1)*ny, :) = y_train(:, row + (0:w-1) + 1);
% end
% 
% Upsilon = u_train(:, q:end); % Leave out last time step to match V_til_1
% YU_bar = [Y; Upsilon];
% 
% DMD of Y
% Y2 = Y(:, 2:end  );
% Y1 = Y(:, 1:end-1);
% 
% YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y and U, above and below
% AB = Y2*pinv(YU); % combined A and B matrix, side by side
% A_hat  = AB(:,1:q*ny); % Extract A matrix
% B_hat  = AB(:,(q*ny+1):end);
% 
% A_hat; % Estimated discrete A matrix of system
% B_hat; % Estimated discrete B matrix of system
% 
% n_hat = size(A_hat, 1);
% l_hat = nu;
% m_hat = n_hat;
% 
% C_hat = eye(m_hat);
% D_hat = zeros(m_hat, l_hat);
% 
% x0_hat = zeros(n_hat,1);


% Run with A and x
run_and_plot = 1;
if run_and_plot
    
    k_start = 500-q;
    window = 10*2; % Number of data points to run and test for

    % Initial condition
    y_hat_0 = [];
    for row = 0:q-1 % First column of spaced Hankel matrix
        k = k_start + row;
        y_hat_0 = [y_hat_0; y_data(:,k)];
    end
    
    % Start at end of initial condition k
    Y_data = y_data(:,k:k+window);
    U      = u_data(:,k:k+window);
    T      = t(:,k:k+window);
    N      = length(Y_data);
    
    % Run model
    Y_hat = zeros(length(y_hat_0),N); % Empty estimated Y
    Y_hat(:,1) = y_hat_0; % Initial condition
    for k = 1:N-1
        Y_hat(:,k+1) = A_hat*Y_hat(:,k) + B_hat*U(:,k);
    end

    y_hat = Y_hat(end-ny+1:end, :); % Extract only non-delay time series (last m rows)

%     % Vector of Mean Absolute Error on testing data
%     MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state

    % Plot data vs model
    figure;
    plot(T, Y_data);
    hold on;
%     plot(t_test, y_test);

    plot(T, y_hat, '--', 'LineWidth', 1); % Plot only non-delay coordinate
    title('Training and Testing data vs Model');
    hold off;
end % run_and_plot

% 
%

% % View only non-delay coordinates
% C_hat = [zeros(ny,(q-1)*ny), eye(ny)];


% LTI system
Ts_mpc = 0.1; % MPC sampling time

% Resample to correct sample time
dmd_sys = ss(A_hat,B_hat,eye(q*ny),zeros(q*ny,nu),Ts); % LTI system
mpc_sys = d2d(dmd_sys,Ts_mpc,'zoh'); % Resample to match MPC

[A_mpc,B_mpc,C_mpc,D_mpc,~] = ssdata(mpc_sys); % Extract resampled matrixes
mpc_sys = ss(A_mpc,B_mpc,C_mpc,D_mpc,Ts_mpc); % LTI system with new Ts 

% MPC object
Ts_mpc = 0.1; % Sample time of MPC (s)
old_status = mpcverbosity('off'); % No display messages
mpc_sys.InputGroup.MV = 2; % Munipulated Variable
mpc_sys.OutputGroup.MO = 6; % Measured Variable

mpc_drone_2d = mpc(mpc_sys,Ts_mpc);
mpc_drone_2d.Weights.OutputVariables = [zeros(1,(q-1)*ny), 1, 1, 0]; % Set weights of delay coordinates to 0, so they do not follow reference

% Nominal operating conditions for AMPC block
U_nom = u_bar;
X_nom = zeros(q*ny,1);
Y_nom = zeros(q*ny,1);
DX_nom = zeros(q*ny,1);










