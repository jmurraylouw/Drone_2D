%% Implentation of Hankel Alternative View Of Koopman for 2D Drone
% Grid search of parameters
% Saves all the results for different parameter combinations

close all;
clear all;

total_timer = tic; % Start timer for this script

% Search space
q_min = 30; % Min value of q in grid search
q_max = 150; % Max value of q in grid search
q_increment = 1; % Increment value of q in grid search

p_min = 30; % Min value of p in grid search
p_max = 300; % Max value of p in grid search
p_increment = 1; % Increment value of p in grid search

q_search = q_min:q_increment:q_max; % List of q parameters to search in
% p_search defined before p for loop

% Extract data
simulation_data_file = 'No_payload_data_2';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

u_data  = out.F_r.Data';
x_data  = out.x.Data';
measured_states = [1,2,3];
y_data  = x_data(measured_states,:); % Measurement data (x, z, theta)
t       = out.tout'; % Time

% Adjust for constant disturbance / mean control values
u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; M*g];
u_data  = u_data - u_bar; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 2000; % Num of data samples for testing
x_test = x_data(:,end-N_test+1:end);
y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Data dimentions
n = size(x_data,1); % number of states
m = size(y_data,1); % number of measurements
l = size(u_data,1); % number of inputs
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
    
% Create empty results table
VariableTypes = {'int16','int16','double'}; % id, q, p, MAE
VariableNames = {'q', 'p', 'MAE_mean'};
for i = 1:m % Mae column for each measured state
    VariableNames = [VariableNames, strcat('MAE_', num2str(i))];
    VariableTypes = [VariableTypes, 'double'];
end
Size = [length(q_search)*length(p_min:p_increment:p_max), length(VariableTypes)];

% Read previous results
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
results_file = ['Data/havok_results_', simulation_data_file, '_sig=', sig_str, '.mat'];

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
        
        p_max_new = min([p_max, q*m]); % Max p to avoid out of bounds 
        p_search = p_min:p_increment:p_max_new; % List of p to search, for every q
        for p = p_search
            p_is_new = 1; % 1 = first time using this p this session
            
            if ~isempty(find(results.q == q & results.p == p, 1)) 
                continue % continue to next p if this combo has been searched before
            end
            
            if q_is_new % Do this only when q is seen first time
                q_is_new = 0; % q is no longer new
            
                w = N_train - q + 1; % num columns of Hankel matrix
                D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)

                % Create Hankel matrix with measurements
                Y = zeros(q*m,w); % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
                for row = 0:q-1 % Add delay coordinates
                    Y(row*m+1:(row+1)*m, :) = y_train(:, row + (0:w-1) + 1);
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
            AB_bar = (U_tilde*S_tilde)*AB_tilde*pinv(U_tilde*S_tilde);
            A_bar = AB_bar(1:q*m, 1:q*m);
            B_bar = AB_bar(1:q*m, q*m+1:end);            

            % Compare to testing data
            % Initial condition (last entries of training data)
            y_hat_0 = zeros(q*m,1);
            for row = 0:q-1 % First column of spaced Hankel matrix
                y_hat_0(row*m+1:(row+1)*m, 1) = y_train(:, end - ((q-1)+1) + row + 1);
            end

            % Run model
            Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
            Y_hat(:,1) = y_hat_0; % Initial condition
            for k = 1:N_test-1
                Y_hat(:,k+1) = A_bar*Y_hat(:,k) + B_bar*u_test(:,k);
            end

            y_hat = Y_hat(end-m+1:end, :); % Extract only non-delay time series (last m rows)

            % Vector of Mean Absolute Error on testing data
            MAE = sum(abs(y_hat - y_test), 2)./N_test; % For each measured state
        
            % Save results
            results(emptry_row,:) = [{q, p, mean(MAE)}, num2cell(MAE')]; % add to table of results
            emptry_row = emptry_row + 1; 
            
        end % p
        
        save(results_file, 'results', 'emptry_row')
        toc;
end % q
format short % back to default/short display

% Save results
results(~results.q,:) = []; % remove empty rows
save(results_file, 'results', 'emptry_row')

best_row = find(results.MAE_mean == min(results.MAE_mean));
best_results = results(best_row,:)

total_time = toc(total_timer); % Display total time taken


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
