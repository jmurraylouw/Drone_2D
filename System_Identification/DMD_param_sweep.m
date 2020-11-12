%% Implentation of DMD
% Grid search of parameters
% Saves all the results for different parameter combinations

% close all;

total_timer = tic; % Start timer for this script

% % Search space
q_min = 10; % Min value of q in grid search
q_max = 40; % Max value of q in grid search
q_increment = 1; % Increment value of q in grid search

p_min = 2; % Min value of p in grid search
p_max = q_max*4; % Max value of p in grid search
p_increment = 1; % Increment value of p in grid search

q_search = q_min:q_increment:q_max; % List of q parameters to search in
% % p_search defined before p for loop
% 
% comment = ''; % Extra comment to differentiate this run
% 
% % Extract data
% % simulation_data_file = 'With_payload_and_noise_data_3';
% simulation_data_file = 'With_payload_and_noise_scale_0';
% load(['Data/', simulation_data_file, '.mat']) % Load simulation data
% 
% Ts = 0.03;     % Desired sample time
% Ts_dmd = Ts;
% y_rows = 1:4;
% MAE_weight = [1; 1; 1; 1]; % Weighting of error of each state when calculating mean
% output_scale = [1; 1; 1; 1]; % Weighting of error of each state when calculating mean
% 
% % Adjust for constant disturbance / mean control values
% % u_bar = mean(out.u.Data,1); % Input needed to keep at a fixed point
% u_bar = [0, (m + M)*g];
% out.u.Data  = out.u.Data - u_bar; % Adjust for unmeasured input
% 
% % Training data
% train_time = 0:Ts:200;
% x_train = resample(out.x, train_time );% Resample time series to desired sample time and training period  
% u_train = resample(out.u, train_time );  
% t_train = x_train.Time';
% N_train = length(t_train);
% 
% x_train = x_train.Data';
% y_train = output_scale.*x_train(y_rows,:);
% u_train = u_train.Data';
% 
% % Testing data
% test_time = 300:Ts:400;
% x_test = resample(out.x, test_time );  
% u_test = resample(out.u, test_time );  
% t_test = x_test.Time';
% N_test = length(t_test); % Num of data samples for testing
% 
% x_test = x_test.Data';
% y_test = output_scale.*x_test(y_rows,:); % One sample of testing data overlaps for initial condition
% u_test = u_test.Data';
% 
% % Data dimentions
% nx = size(x_train,1); % number of states
% ny = size(y_train,1); % number of measurements
% nu = size(u_train,1); % number of inputs  

% % Add noise
% rng('default');
% rng(1); % Repeatable random numbers
% sigma = 0.01; % Noise standard deviation
% y_train = y_train + sigma*randn(size(y_train));

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
results_file = ['Data/dmd_results_', comment, simulation_data_file, '_sig=', sig_str, '.mat'];

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
                
                % Hankel matrix with delay measurements
                if q == 1 % Special case if no delay coordinates
                    Upsilon = u_train(:, q:end);
                else
                    Upsilon = zeros((q-1)*ny,w); % Augmented state Y[k] at top
                    for row = 0:q-2 % Add delay coordinates
                        Upsilon((end - ny*(row+1) + 1):(end - ny*row), :) = y_train(:, row + (1:w));
                    end

                    Upsilon = [Upsilon; u_train(:, q:end)]; % Leave out last time step to match V_til_1
                end

                % Matrix with time series of states
                Y = y_train(:, q-1 + (1:w));

                % DMD of Y
                Y2 = Y(:, 2:end  );
                Y1 = Y(:, 1:end-1);

                YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y above and U below

                % SVD of the Hankel matrix
                [U1,S1,V1] = svd(YU, 'econ');

            end

            % Truncate SVD matrixes
            U_tilde = U1(:, 1:p); 
            S_tilde = S1(1:p, 1:p);
            V_tilde = V1(:, 1:p);

            % YU = \approx U_tilde*S_tilde*V_tilde'
            AB = Y2*pinv(U_tilde*S_tilde*V_tilde'); % combined A and B matrix, side by side
            % AB = Y2*pinv(YU); % combined A and B matrix, side by side

            % System matrixes from DMD
            A_dmd  = AB(:,1:ny); % Extract A matrix
            B_dmd  = AB(:,(ny+1):end);
            % A = stabilise(A,1);
            
            % Initial condition
            y_hat_0 = y_test(:,q);

            % Initial delay coordinates
            y_delays = zeros((q-1)*ny,1);
            k = q; % index of y_data
            for i = 1:ny:ny*(q-1) % index of y_delays
                k = k - 1; % previosu index of y_data
                y_delays(i:(i+ny-1)) = y_test(:,k);
            end

            % Run model
            y_hat = zeros(ny,N_test); % Empty estimated Y
            y_hat(:,1) = y_hat_0; % Initial condition
            for k = 1:N_test-1
                upsilon = [y_delays; u_test(:,k)]; % Concat delays and control for use with B
                y_hat(:,k+1) = A_dmd*y_hat(:,k) + B_dmd*upsilon;
                if q ~= 1
                    y_delays = [y_hat(:,k); y_delays(1:(end-ny),:)]; % Add y(k) to y_delay for next step [y(k); y(k-1); ...]
                end
            end

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

best_results = results((results.MAE_mean == min(results.MAE_mean)),:)

%% Only for this Ts:
% results_Ts = results((results.Ts == Ts),:);
% best_results_Ts = results_Ts((results_Ts.MAE_mean == min(results_Ts.MAE_mean)),:)
% 
% total_time = toc(total_timer); % Display total time taken
% 
%% For one q:
results_q = results((results.q == best_results.q),:);
figure
% semilogy(results_q.p, results_q.MAE_1, 'r.')
% hold on
semilogy(results_q.p, results_q.MAE_mean, 'k.')
% hold off

%% Plot results
plot_results = 1;
if plot_results
    figure
    semilogy(results.q, results.MAE_mean, '.')
    y_limits = [5e-1, 1e0];
%     ylim(y_limits)
    xlim([min(results.q), max(results.q)])
    title('DMD')
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
