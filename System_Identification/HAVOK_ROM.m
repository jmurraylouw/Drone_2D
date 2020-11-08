%% Implentation of Hankel Alternative View Of Koopman for 2D Drone
close all;

% Extract data
simulation_data_file = 'With_payload_data_12';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

Ts = 0.03;     % Desired sample time
Ts_havok = Ts;
y_rows = 1:4;

% Adjust for constant disturbance / mean control values
% u_bar = mean(out.u.Data,1); % Input needed to keep at a fixed point
u_bar = [0, (m + M)*g];
out.u.Data  = out.u.Data - u_bar; % Adjust for unmeasured input

% Training data
train_time = 0:Ts:200;
x_train = resample(out.x, train_time );% Resample time series to desired sample time and training period  
u_train = resample(out.u, train_time );  
t_train = x_train.Time';
N_train = length(t_train);

x_train = x_train.Data';
y_train = x_train(y_rows,:);
u_train = u_train.Data';

% Testing data
test_time = 200:Ts:300;
x_test = resample(out.x, test_time );  
u_test = resample(out.u, test_time );  
t_test = x_test.Time';
N_test = length(t_test); % Num of data samples for testing

x_test = x_test.Data';
y_test = x_test(y_rows,:); % One sample of testing data overlaps for initial condition
u_test = u_test.Data';

% Data dimentions
nx = size(x_train,1); % number of states
ny = size(y_train,1); % number of measurements
nu = size(u_train,1); % number of inputs

% % Add noise
% rng('default');
% rng(1); % Repeatable random numbers
sigma = 0; % Noise standard deviation
% y_data_noise = y_data + sigma*randn(size(y_data));

comment = 'no_angle_measure'; % Extra comment to differentiate this run

% Read previous results
sigma = 0;
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
results_file = ['Data/havok_results_', comment, '_', simulation_data_file, '_sig=', sig_str, '.mat'];

try
    load(results_file);
    results(~results.q,:) = []; % remove empty rows
    
    % Parameters
    best_row = find(results.MAE_mean == min(results.MAE_mean));
    best_results = results(best_row,:);
    q = double(best_results.q);
    p = double(best_results.p);
    
    only_q_Ts = 0; % Try best result for specific q
    if only_q_Ts
        q = 3;
        q_results = results((results.q == q & results.Ts == Ts),:);
        best_row = find(q_results.MAE_mean == min(q_results.MAE_mean));
        best_results = q_results(best_row,:)
        p = double(best_results.p);
    end
    
catch
    disp('No saved results file')  
end

% Override parameters:
q = 8;
p = 23;

q
p

w = N_train - q + 1; % num columns of Hankel matrix
D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)

% Create Hankel matrix with measurements
Y = zeros((q)*ny,w); % Augmented state Y[k] at top
for row = 0:q-1 % Add delay coordinates
    Y((end - ny*(row+1) + 1):(end - ny*row), :) = y_train(:, row + (1:w));
end

Upsilon = u_train(:, q:end); % Leave out last time step to match V_til_1
YU_bar = [Y; Upsilon];

% SVD of the Hankel matrix
[U1,S1,V1] = svd(YU_bar, 'econ');
figure, semilogy(diag(S1), 'x'), hold on;
title('Singular values of Omega, showing p truncation')
plot(p, S1(p,p), 'ro'), hold off;

% Truncate SVD matrixes
U_tilde = U1(:, 1:p); 
S_tilde = S1(1:p, 1:p);
V_tilde = V1(:, 1:p);

% Superposition of Y and U influence in V space:
% (U_tilde*S_tilde)*V_tilde     \approx =                       [Y, U]
%                               V_tilde = pinv(U_tilde*S_tilde)*[Y; U] 
% T = pinv(U_tilde*S_tilde) = [TY; TU]
% V_tilde = [TY, TU]*[Y; U] = TY*Y + TU*U
%                           = VY   + VU    

% Transformation matrices
T = pinv(U_tilde*S_tilde); % Transformation matrix from YU_bar to V
TY = T(:, 1:q*ny);
TU = T(:, q*ny+1:end);

% Map X and U to V-space
Yv = TY*Y; % Y matrix in V-space (Y contribution to V_tilde)
Uv = TU*Upsilon; % Upslion matrix in V-space (Upsilon contribution to V_tilde)

% Setup for DMD regression
Yv2 = Yv(:, 2:end  ); % One timestep ahead of VY1
Yv1 = Yv(:, 1:end-1);
Uv1 = Uv(:, 1:end-1);
% Yv2 = Av*Yv1 + Bv*Yv1    <- reduced order model format
% Yv2 = [Av, Bv]*[Yv1;Uv1]

% Regression
ABv = Yv2*pinv([Yv1; Uv1]); % [A,B] for V-space
Av = ABv(:, 1:p);
Bv = ABv(:, p+1:end);

% Setup V2 one timestep into future from V1
V_til_2 = V_tilde(2:end  , :)'; % Turnd on side (wide short matrix)
V_til_1 = V_tilde(1:end-1, :)';

% DMD on V
AB_tilde = V_til_2*pinv(V_til_1); % combined A and B matrix, side by side
AB_tilde = stabilise(AB_tilde,3);

% Convert to x coordinates
AB_havok = (U_tilde*S_tilde)*AB_tilde*pinv(U_tilde*S_tilde);

% System matrixes from HAVOK
A_havok = AB_havok(1:q*ny, 1:q*ny);
B_havok = AB_havok(1:q*ny, q*ny+1:end);
% A_havok = stabilise(A_havok,10);

% Make matrix sparse
A_havok(ny+1:end, :) = [eye((q-1)*ny), zeros((q-1)*ny, ny)]; % Add Identity matrix to carry delays over to x(k+1)
B_havok(ny+1:end, :) = zeros((q-1)*ny, nu); % Input has no effect on delays

% DMD of Y
Y2 = Y(:, 2:end  );
Y1 = Y(:, 1:end-1);

YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y and U, above and below
AB = Y2*pinv(YU); % combined A and B matrix, side by side

% System matrixes from DMD
A  = AB(:,1:q*ny); % Extract A matrix
B  = AB(:,(q*ny+1):end);

% A = stabilise(A,10);

% Compare to testing data

%% Run with HAVOK with Reduced Order Modelling and Control

% Initial condition (last entries of training data)
y_hat_0 = zeros(q*ny,1); % Y[k] at top
for row = 0:q-1 % First column of spaced Hankel matrix
    y_hat_0(row*ny+1:(row+1)*ny, 1) = y_test(:,q-row);
end

% Convert to V-space:
yv_0 = TY*y_hat_0; % Initial condintion in V-space
Uv_test = TU*u_test; % Control inputs in V-space

% Run model
Yv_hat = zeros(length(yv_0),N_test); % Empty estimated Y
Yv_hat(:,q) = yv_0; % Initial condition
for k = q:N_test-1
    Yv_hat(:,k+1) = Av*Yv_hat(:,k) + Bv*Uv_test(:,k);
end

% Convert to Y-space (original)
Y_hat_v = pinv(TY)*Yv_hat;
y_hat_v = Y_hat_v(1:ny, :); % Extract only non-delay time series

%% Plot data vs model
figure;
plot(t_train, y_train);
hold on;
plot(t_test, y_test);

% plot(t_test, y_hat, 'k--', 'LineWidth', 1); % Plot only non-delay coordinate
plot(t_test, y_hat_v, 'r--', 'LineWidth', 1); % Plot only non-delay coordinate
title('Training and Testing data vs ROM HAVOK');
% legend('','','','','','', 'x hat','z hat','theta hat', 'x hat bar','z hat bar','theta hat bar');
hold off;

MAE_v = sum(abs(y_hat_v - y_test), 2)./N_test % For each measured state

%% Run with HAVOK (A_havok, B_havok and x)
% figure;
% plot(U1(:,1:5))
% title('First 5 modes of SVD')

% Compare to testing data
% Initial condition (last entries of training data)
y_hat_0 = zeros(q*ny,1); % Y[k] at top
for row = 0:q-1 % First column of spaced Hankel matrix
    y_hat_0(row*ny+1:(row+1)*ny, 1) = y_test(:,q-row);
end

% Run model
Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
Y_hat(:,q) = y_hat_0; % Initial condition
for k = q:N_test-1
    Y_hat(:,k+1) = A_havok*Y_hat(:,k) + B_havok*u_test(:,k);
end

y_hat_bar = Y_hat(1:ny, :); % Extract only non-delay time series

% Vector of Mean Absolute Error on testing data
MAE = sum(abs(y_hat_bar - y_test), 2)./N_test % For each measured state


% %% Run with DMD (A and x)
% 
% % Initial condition
% y_hat_0 = zeros(q*ny,1);
% for row = 0:q-1 % First column of spaced Hankel matrix
%     y_hat_0(row*ny+1:(row+1)*ny, 1) = y_train(:, end - ((q-1)+1) + row + 1);
% end
%             
% % Run model
% Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
% Y_hat(:,1) = y_hat_0; % Initial condition
% for k = 1:N_test-1
%     Y_hat(:,k+1) = A*Y_hat(:,k) + B*u_test(:,k);
% end
% 
% y_hat = Y_hat(end-ny+1:end, :); % Extract only non-delay time series (last m rows)
% 
% % Vector of Mean Absolute Error on testing data
% MAE = sum(abs(y_hat - y_test), 2)./N_test % For each measured state
% 
% % Compare MAE and MAE_til
% MAE_error_percent = (MAE - MAE_bar)./MAE_bar*100

%% Plot data vs model
% figure;
% plot(t_train, y_train);
% hold on;
% plot(t_test, y_test);
% 
% % plot(t_test, y_hat, 'k--', 'LineWidth', 1); % Plot only non-delay coordinate
% plot(t_test, y_hat_bar, 'r--', 'LineWidth', 1); % Plot only non-delay coordinate
% title('Training and Testing data vs Model (red = HAVOK, black = DMD)');
% % legend('','','','','','', 'x hat','z hat','theta hat', 'x hat bar','z hat bar','theta hat bar');
% hold off;

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