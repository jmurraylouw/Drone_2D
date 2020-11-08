%% Implentation of Hankel Alternative View Of Koopman for 2D Drone
% Grid search of parameters
% Saves all the results for different parameter combinations

close all;

total_timer = tic; % Start timer for this script

% Search space
q_min = 2; % Min value of q in grid search
q_max = 10; % Max value of q in grid search
q_increment = 1; % Increment value of q in grid search

p_min = 2; % Min value of p in grid search
p_max = 300; % Max value of p in grid search
p_increment = 1; % Increment value of p in grid search

q_search = q_min:q_increment:q_max; % List of q parameters to search in
% p_search defined before p for loop

comment = 'no_angle_measure'; % Extra comment to differentiate this run

% Extract data
simulation_data_file = 'With_payload_data_12';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

Ts = 0.03;     % Desired sample time
y_rows = 1:4;
MAE_weight = [1; 1; 1]; % Weighting of error of each state when calculating mean

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
x_test = resample(out.x, train_time );  
u_test = resample(out.u, train_time );  
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
% % sigma = 0.001; % Noise standard deviation
% y_data_noise = y_data + sigma*randn(size(y_data));

% Create empty results table
VariableTypes = {'double', 'int16',   'int16', 'int16', 'double'}; % id, q, p, MAE
VariableNames = {'Ts',     'N_train', 'q',     'p',     'MAE_mean'};
for i = 1:ny % Mae column for each measured state
    VariableNames = [VariableNames, strcat('MAE_', num2str(i))];
    VariableTypes = [VariableTypes, 'double'];
end
Size = [length(q_search)*length(p_min:p_increment:p_max), length(VariableTypes)];

% Read previous results
sigma = 0;
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
results_file = ['Data/havok_results_', comment, '_', simulation_data_file, '_sig=', sig_str, '.mat'];

try
    load(results_file);
    results(~results.q,:) = []; % remove empty rows
    results = [results; table('Size',Size,'VariableTypes',VariableTypes,'VariableNames',VariableNames)];
    
catch
    disp('No saved results file')  
    
    results = table('Size',Size,'VariableTypes',VariableTypes,'VariableNames',VariableNames);
    emptry_row = 1; % Keep track of next empty row to insert results 
end

% Grid search
format compact % display more compact
for q = q_search
        q_is_new = 1; % 1 = first time using this q this session
        q
        tic;
        
        p_max_new = min([p_max, q*ny]); % Max p to avoid out of bounds 
        p_search = p_min:p_increment:p_max_new; % List of p to search, for every q
        for p = p_search
            p_is_new = 1; % 1 = first time using this p this session
            
            if ~isempty(find(results.q == q & results.p == p & results.Ts == Ts & results.N_train == N_train, 1)) 
                continue % continue to next p if this combo has been searched before
            end
            
            if q_is_new % Do this only when q is seen first time
                q_is_new = 0; % q is no longer new
            
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

            end
            
            % Truncate SVD matrixes
            U_tilde = U1(:, 1:p); 
            S_tilde = S1(1:p, 1:p);
            V_tilde = V1(:, 1:p);

            % Setup V2 one timestep into future from V1
            V_til_2 = V_tilde(2:end  , :)'; % Turnd on side (wide short matrix)
            V_til_1 = V_tilde(1:end-1, :)';

            % DMD on V
            AB_tilde = V_til_2*pinv(V_til_1); % combined A and B matrix, side by side
            AB_tilde = stabilise(AB_tilde,3);
            
            % convert to x coordinates
            AB = (U_tilde*S_tilde)*AB_tilde*pinv(U_tilde*S_tilde);
            A = AB(1:q*ny, 1:q*ny);
            B = AB(1:q*ny, q*ny+1:end);            

            % Make matrix sparse
            A(ny+1:end, :) = [eye((q-1)*ny), zeros((q-1)*ny, ny)]; % Add Identity matrix to carry delays over to x(k+1)
            B(ny+1:end, :) = zeros((q-1)*ny, nu); % Input has no effect on delays
            
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
                Y_hat(:,k+1) = A*Y_hat(:,k) + B*u_test(:,k);
            end

            y_hat = Y_hat(1:ny, :); % Extract only non-delay time series

            % Vector of Mean Absolute Error on testing data
            MAE = sum(abs(y_hat - y_test), 2)./N_test; % For each measured state
            
            % Save results
            results(emptry_row,:) = [{Ts, N_train, q, p, mean(MAE.*MAE_weight)}, num2cell(MAE')]; % add to table of results
            emptry_row = emptry_row + 1; 
            
        end % p
        
        save(results_file, 'results', 'emptry_row')
        toc;
end % q
format short % back to default/short display

% Save results
results(~results.q,:) = []; % remove empty rows
save(results_file, 'results', 'emptry_row')

best_results_overall = results((results.MAE_mean == min(results.MAE_mean)),:)

%% Only for this Ts:
results_Ts = results((results.Ts == Ts),:);
best_results_Ts = results_Ts((results_Ts.MAE_mean == min(results_Ts.MAE_mean)),:)

total_time = toc(total_timer); % Display total time taken

%% For one q:
results_q = results((results.q == 5),:);
figure
semilogy(results_q.p, results_q.MAE_1, 'r.')
hold on
semilogy(results_q.p, results_q.MAE_mean, 'k.')
hold off

%% Plot spread of results for this Ts
plot_results = 1;
if plot_results
    semilogy(results_Ts.q, results_Ts.MAE_mean, '.')
    ylim([1e-3, 1e0])
end

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
