%% Setup workspace for simulation: flight_controller_payload.slx
% Called from InitFunc callback

% PID controllers
% _2 is most like Anton.
% _4 is most damped
load('Data/Drone_2D_control_params_2.mat'); % Load controller gain values

% Low Pass Filter cut-off frequency
N_dtheta = Inf;
N_dx = 5*2*pi; % Rad/s
N_dz = 5*2*pi; % Rad/s

dx_sp_max = 12; % (m/s) Max x velocity command (x and z)
dx_sp_min = -12; % (m/s) Min x velocity command (x and z)

dz_sp_max = 12; % (m/s) Max z velocity command (x and z)
dz_sp_min = -12; % (m/s) Min z velocity command (x and z)

theta_sp_max = 30*pi/180; % (rad) Max. pitch rate command
theta_sp_min = -30*pi/180; % Min. pitch rate command

dtheta_sp_max = 120*pi/180; % Max. pitch rate command
dtheta_sp_min = -120*pi/180; % Min. pitch rate command
% 
% sigma = 0.001; % Std deviation of measurement noise

% Model parameters
M     = 4.5; % Mass of drone body (at fulcrum)
I_yy  = 0.235; % Moment of inertia of drone body about body x axis
r     = 0.49*1/sqrt(2); % Distance from each rotor force to COM of drone
g     = -9.81; % Acceleration due to gravity (always negative)
C_Dx  = 0.2 ; % Damping coef. of drone through air in x direction
C_Dz  = 0.2; % Damping coef. of drone in z direction
rho   = 1.225; % Air density (kg/m^3)
tau   = 0.07; % Motor time constant

m     = 2; % Mass of swinging payload (kg)
l     = 1; % Length of pendulum (m)
cbeta = 0.03; % Rotational damping coef of payload at connection

C_px = 0.01; % Damping coef. of drone through air in x direction
C_pz = 0.01; % Damping coef. of drone through air in z direction

% Noise parameters
Ts_noise = 0.01; % Sampling time of noise
omega_b_noise = 6e-8; % Anglular velocity. Noise power of bandwidth limited white noise
quat_noise = 6e-8; % Angles
vel_e_noise = 4e-8; % Linear Velocity
pos_e_noise = 4e-7; % Linear position

% Mixing Matrix
MM = [-1, 1;
       1, 1]; % [T1; T2] = MM*[delta_E; delta_T]

% Forgetting factor
lambda = 0.985; % Exponential forgetting factor of moving average to smooth out spikes in input for SVD performance

% DMD parameters
% Ts_dmd, ny, nu, x0, u0, N_train, q, model_intervals
% load('Data/MPC_initial_plant.mat'); % load A_dmd, B_dmd, q, Ts_dmd from a previous DMD run
% model_intervals = 10; 

% Dimentions
y_rows = 1:4;
nx = 8; % Number of states
ny = length(y_rows); % Number of measurements
nu = 2; % Number of inputs

% Initial conditions
x0 = zeros(nx,1); % Initial state
y0 = zeros(ny,1); % Initial measurements
u0 = -0.5*(M+m)*g*[1; 1]; % Initial input

x0_no_load = zeros(6,1); % Initial state
y0_no_load = zeros(3,1); % Initial measurements
Ts_csv = 0.1; % Sample time of To Worspace blocks for csv data

% Excitement signal
Ts_excite = 0; % Sample time
var_excite = 0; % Variance of excitement signal (pulse train)

simulation_data_file = 'With_payload_data_2';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

% Disturbance model to account for model uncertainty (eliminate steady-state error)
CO = 2; % number of Controlled Outputs (x and z). theta is not controlled to a reference
dist_influence = 0; % Disturbances include uncertainty in model

% Internal plant model
model = 'havok'; % Choose which model to use for MPC
switch model
    case 'dmd'
%         start_time = 20;
%         end_time = 100;
%         y_rows = 1:4;
%         q = 6;
%         sigma = 0;
%         plot_prediction = 0;
% 
%         [A_dmd, B_dmd] = model_DMD(out, start_time, end_time, Ts_mpc, q, y_rows, sigma, plot_prediction);
%         
        % Change model structure so delays are included in A, not B 
        A_mpc = [A_dmd,       B_dmd(:, 1:end-nu);
                 eye((q-1)*ny),   zeros((q-1)*ny,ny)];

        % B_mpc = [B_dmd(:, end-nu+1:end); zeros((q-1)*ny, nu)];
        B_ud = dist_influence*[eye(CO); zeros(q*ny - CO, CO)]; % B of unmeasured disturbance, for distrubance force in x and z
        B_mpc = [[B_dmd(:, end-nu+1:end); zeros((q-1)*ny, nu)], B_ud];
    
    case 'havok'
%         load('Data/havoc_model_5.mat')
        A_mpc = A_havok;
        B_ud = dist_influence*[eye(CO); zeros(q*ny - CO, CO)]; % B of unmeasured disturbance, for distrubance force in x and z
        B_mpc = [B_havok, B_ud];
        
    otherwise
        error("Choose only 'dmd' or 'havok' ")
end

% Sample time of MPC:
Ts_mpc = Ts_havok;

C_mpc = eye(q*ny);
D_mpc = zeros(q*ny, nu + CO);

mpc_sys = ss(A_mpc,B_mpc,C_mpc,D_mpc,Ts_mpc); % LTI system

% Initital conditions for extended measurment vector for MPC
% All delay states are also at y0
y_ext_0 = zeros(q*ny, 1); % Allocate space
for row = 0:q-1 % First column of spaced Hankel matrix
        y_ext_0(row*ny+1:(row+1)*ny, 1) = y0;
end

delays_0 = []; % Initial delay vector
for i = 1:q-1
    delays_0 = [delays_0; y0];
end

% MPC object
old_status = mpcverbosity('off'); % No display messages
mpc_sys.OutputGroup.MO = 1:q*ny; % Measured Output

mpc_sys.InputGroup.MV = 1:nu; % Munipulated Variable indices
mpc_sys.InputGroup.UD = nu + (1:CO); % unmeasured disturbance indices, one for each 

tuning_weight = 1; % Smaller = robust, Larger = aggressive
mpc_drone_2d = mpc(mpc_sys,Ts_mpc);

% Manually set covariance
x_mpc = mpcstate(mpc_drone_2d); % Initial state
covariance = zeros(size(x_mpc.Covariance));
covariance(y_rows, y_rows) = diag([2e-3, 1e-3, 1e-5, 1e-4]);
x_mpc = mpcstate(mpc_drone_2d, [], [], [], [], covariance);

Ty = 6; % Prediction period, For guidance, minimum desired settling time (s)
Tu = 3; % Control period, desired control settling time
mpc_drone_2d.PredictionHorizon  = floor(Ty/Ts_mpc); %t_s/Ts_mpc; % Prediction horizon (samples), initial guess according to MATLAB: Choose Sample Time and Horizons
mpc_drone_2d.ControlHorizon     = floor(Tu/Ts_mpc); % Control horizon (samples)
mpc_drone_2d.Weights.OutputVariables        = [1, 1, 0, 10, zeros(1, (q-1)*ny)]*tuning_weight;
mpc_drone_2d.Weights.ManipulatedVariables   = 5e-1*[1, 1]*tuning_weight; % Weights of delay coordinates to 0
mpc_drone_2d.Weights.ManipulatedVariablesRate     = 1e-2*[1, 1]/tuning_weight;

% Output bounds
% theta_min = -30*(pi/180);
% theta_max = abs(theta_min); % Anton's pitch command constraint
% 
% mpc_drone_2d.OV(3).Min = theta_min;
% mpc_drone_2d.OV(3).Max = theta_max;
% 
% % Input bounds
% 
% % Normalised
% F_r_z_min = -150;
% F_r_z_max = abs((M + m)*g);
% 
% F_r_x_max = F_r_z_max;
% F_r_x_min = -F_r_x_max;
% 
% mpc_drone_2d.MV(1).Min = F_r_x_min;
% mpc_drone_2d.MV(1).Max = F_r_x_max;
% 
% mpc_drone_2d.MV(2).Min = F_r_z_min;
% mpc_drone_2d.MV(2).Max = F_r_z_max;

% Display
mpc_drone_2d 

% Nominal operating conditions for AMPC block
u_bar = [0; M*g];
U_nom = zeros(nu,1);
X_nom = zeros(q*ny,1);
Y_nom = zeros(q*ny,1);
DX_nom = zeros(q*ny,1);

% Way points
num_waypoints = 100; % Number of waypoints included in command

clear table
waypoints = table('Size', [(num_waypoints+1)*2, 3], 'VariableTypes', ["double", "double", "double"]);
waypoints.Properties.VariableNames = {'point_time', 'x_coord', 'z_coord'};

waypoint_opt = 'random xz'; % waypoint option
switch waypoint_opt
    case 'random xz'
        x_coord = 0;
        z_coord = 0;
        waypoints(1,:) = table(0,                   x_coord, z_coord); % Initial point

        x_min   = 0.5;     x_max   = 1.5; % (m) minimum and maximum step size for waypoints
        z_min   = 0.5;     z_max   = 1.5;
        interval_min = 4;       interval_max = 10;  % (s) minimum and maximum TIME interval between commands

        rng(0); % Initialise random number generator for repeatability
        point_time = 0; % Currently at time zero
        next_point = 1; % Index of next point
        for i = 1:num_waypoints
            % Step x only
            time_interval = (interval_max - interval_min).*rand() + interval_min; % (s) random time interval between commands
            point_time = point_time + time_interval;

            waypoints(next_point,  :) = table(point_time, x_coord, z_coord); % Previous point    
            next_point = next_point + 1;
            x_step = ((x_max - x_min).*rand() + x_min)*sign(randn()); % x step of next waypoint (size)*(direction)
            x_coord = x_coord + x_step;
            waypoints(next_point,:) = table(point_time, x_coord, z_coord); % Next point
            next_point = next_point + 1; 
            
            % Step z only
            time_interval = (interval_max - interval_min).*rand() + interval_min; % (s) random time interval between commands
            point_time = point_time + time_interval;

            waypoints(next_point,  :) = table(point_time, x_coord, z_coord); % Previous point    
            next_point = next_point + 1;
            z_step = ((z_max - z_min).*rand() + z_min)*sign(randn()); % z step of next waypoint (size)*(direction)
            z_coord = z_coord + z_step;
            waypoints(next_point,:) = table(point_time, x_coord, z_coord); % Next point
            next_point = next_point + 1;            
        end
        i = i+1;
        waypoints(2*i,  :) = table(point_time+interval_max, x_coord, z_coord); % Add time to reach final point

    case 'random x'
        x_coord = 0;
        z_coord = 0; % constant z
        waypoints(1,:) = table(0, x_coord, z_coord); % Initial point

        x_min        = -5;     x_max         = 5; % (m) minimum and maximum coordinates for waypoints
        interval_min = 2;     interval_max = 5;  % (s) minimum and maximum TIME interval between commands

        point_time = 0;
        rng(0); % Initialise random number generator for repeatability
        for i = 1:num_waypoints
            time_interval = (interval_max - interval_min).*rand() + interval_min; % (s) random time interval between commands
            point_time = point_time + time_interval;
            if i == 1 % set first interval
                point_time = 5;
            end
            waypoints(2*i,  :) = table(point_time, x_coord, z_coord); % Previous point    
            x_coord    = (x_max - x_min).*rand() + x_min; % x coordinate of next waypoint
            waypoints(2*i+1,:) = table(point_time, x_coord, z_coord); % Next point
        end
        i = i+1;
        waypoints(2*i,  :) = table(point_time+interval_max, x_coord, z_coord); % Add time to reach final point
        
    case 'regular x'
        time_interval = 4; % (s) interval between commands
        step_size = 2;
        x_coord = step_size;
        z_coord = 0;
        point_time = 0;
        waypoints(1,:) = table(0, x_coord, z_coord); % Initial point
        for i = 1:num_waypoints
        
            point_time = point_time + time_interval;
        
            waypoints(2*i,  :) = table(point_time, x_coord, z_coord); % Previous point    
            x_coord    = x_coord + step_size; % x coordinate of next waypoint
            z_coord    = 0; % z coordinate of next waypoint   
            waypoints(2*i+1,:) = table(point_time, x_coord, z_coord); % Next point
        end
        waypoints(end,:) = table(200, x_coord, z_coord); % Final point
        
     otherwise
        error("Unknown waypoint option")
end

waypoints_ts = timeseries([waypoints.x_coord, waypoints.z_coord], waypoints.point_time); % timeseries object for From Workspace block
% plot(waypoints_ts.Time, waypoints_ts.Data)
