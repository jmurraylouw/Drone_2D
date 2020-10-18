function sfunc_DMD(block)
% Applies Dynamic Mode Decomposition with Control on a 
% Moving Window of past input data to obtain the state space system matrixes

setup(block);
% Get dimentions of state and input vectors

%endfunction

function setup(block)

% Dialog parameters
block.NumDialogPrms     = 8;

Ts      = block.DialogPrm(1).Data;
ny      = block.DialogPrm(2).Data;
nu      = block.DialogPrm(3).Data;
y0      = block.DialogPrm(4).Data;
u0      = block.DialogPrm(5).Data;
N_train = block.DialogPrm(6).Data;
q       = block.DialogPrm(7).Data;
model_intervals = block.DialogPrm(8).Data; 

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 2;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% INPUT port properties
% y (measurement vector)
block.InputPort(1).Dimensions        = ny;
block.InputPort(1).DatatypeID        = 0;  % double
block.InputPort(1).Complexity        = 'Real';

% u (input vector)
block.InputPort(2).Dimensions        = nu;
block.InputPort(2).DatatypeID        = 0;  % double
block.InputPort(2).Complexity        = 'Real';


% OUTPUT port properties
% A (Flattened system matrix)
block.OutputPort(1).DatatypeID       = 0; % double
block.OutputPort(1).Complexity       = 'Real';
block.OutputPort(1).SamplingMode     = 'Sample';
block.OutputPort(1).Dimensions       = [ny, ny];

% B (Input matrix)
block.OutputPort(2).DatatypeID       = 0; % double
block.OutputPort(2).Complexity       = 'Real';
block.OutputPort(2).SamplingMode     = 'Sample';
block.OutputPort(2).Dimensions       = [ny, nu];

% Sample time
block.SampleTimes = [Ts, 0]; % Set sample time

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

% Internal registry for block methods
block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required
block.RegBlockMethod('SetInputPortSamplingMode', @SetInputPortSamplingMode);
  
%end setup

function DoPostPropSetup(block)   

    block.NumDworks = 2;
    
    % Dialog parameters
    Ts      = block.DialogPrm(1).Data;
    ny      = block.DialogPrm(2).Data;
    nu      = block.DialogPrm(3).Data;
    y0      = block.DialogPrm(4).Data;
    u0      = block.DialogPrm(5).Data;
    N_train = block.DialogPrm(6).Data;
    q       = block.DialogPrm(7).Data;
    model_intervals = block.DialogPrm(8).Data; 
  
    w = N_train - q + 1; % num columns of Hankel matrix
    
    % Matrix of window of state vectors
    % Needs to reshape because Dwork cannot store matrix, only vector
    block.Dwork(1).Name            = 'Y';
    block.Dwork(1).Dimensions      = ny*w;
    block.Dwork(1).DatatypeID      = 0;      % double
    block.Dwork(1).Complexity      = 'Real'; % real
    block.Dwork(1).UsedAsDiscState = true;

    % Matrix of window of input vectors
    % Needs to reshape because Dwork cannot store matrix
    block.Dwork(2).Name            = 'Upsilon';
    block.Dwork(2).Dimensions      = ((q-1)*ny + nu)*w; % height(Delays + u)*width
    block.Dwork(2).DatatypeID      = 0;      % double
    block.Dwork(2).Complexity      = 'Real'; % real
    block.Dwork(2).UsedAsDiscState = true;
        
function InitializeConditions(block)

%end InitializeConditions

function Start(block)
    % Dialog parameters
    Ts      = block.DialogPrm(1).Data;
    ny      = block.DialogPrm(2).Data;
    nu      = block.DialogPrm(3).Data;
    y0      = block.DialogPrm(4).Data;
    u0      = block.DialogPrm(5).Data;
    N_train = block.DialogPrm(6).Data;
    q       = block.DialogPrm(7).Data;
    model_intervals = block.DialogPrm(8).Data; 
    
    w = N_train - q + 1; % num columns of Hankel matrix

    % Training data with only initial conditions
    y_train = y0 + zeros(ny,N_train); % Assume all data before simulation is at init conditions
    u_train = u0 + zeros(nu,N_train);
    
    % Inititialise Hankel matrix with delay measurements
    Delays = zeros((q-1)*ny,w); % Augmented state with delay coordinates [...; Y(k-2); Y(k-1); Y(k)]
    for row = 0:q-2 % Add delay coordinates
        Delays((end - ny*(row+1) + 1):(end - ny*row), :) = y_train(:, row + (1:w));
    end
    Upsilon = [Delays; u_train(:, q:end)]; % Leave out last time step to match V_til_1

    % Matrix with time series of states
    Y = y_train(:, q:end);
    
    % Reshape to vector and Assign to memory/Dwork
    block.Dwork(1).Data = reshape(Y,       1, ny*w );
    block.Dwork(2).Data = reshape(Upsilon, 1, ((q-1)*ny + nu)*w );
    
%end Start

function Outputs(block)

    % Dialog parameters
    Ts      = block.DialogPrm(1).Data;
    ny      = block.DialogPrm(2).Data;
    nu      = block.DialogPrm(3).Data;
    y0      = block.DialogPrm(4).Data;
    u0      = block.DialogPrm(5).Data;
    N_train = block.DialogPrm(6).Data;
    q       = block.DialogPrm(7).Data;
    model_intervals = block.DialogPrm(8).Data; 
    
    w = N_train - q + 1; % num columns of Hankel matrix
    
    % Dwork Memory
    Y_dwork = block.Dwork(1).Data;
    U_dwork = block.Dwork(2).Data;

    Y       = reshape(Y_dwork, ny,              w); % Reshape vector into matrix
    Upsilon = reshape(U_dwork, ((q-1)*ny + nu), w); % Reshape vector into matrix
    
    % Current timestep datadata
    u       = block.InputPort(1).Data;
    x       = block.InputPort(2).Data;
    
    % Update Upsilon
    new_row = [Y(:,end); Upsilon(1:(end - ny - nu), end); u];
    Upsilon = [Upsilon(:,2:end), new_row]; % Add new row of updated values

    % Update Y
    Y = [Y(:,2:end), y]; % Forget oldest entry, add y_k
        
    % DMD of Y
    Y2 = Y(:, 2:end  ); 
    Y1 = Y(:, 1:end-1); % One time-step behin of Y2

    YU = [Y1; Upsilon(:,1:end-1)]; % Combined matrix of Y and U, above and below
    AB = Y2*pinv(YU); % combined A and B matrix, side by side

    % System matrixes from DMD
    A  = AB(:,1:ny); % Extract A matrix
    B  = AB(:,(ny+1):end);

    % Adjust slightly unstable eigenvalues
    A = stabilise(A,3);
    
    % Output
    block.OutputPort(1).Data = A;
    block.OutputPort(2).Data = B;
    
    % Reshape to vector and update memory/Dwork
    block.Dwork(1).Data = reshape(Y,       1, ny*w );
    block.Dwork(2).Data = reshape(Upsilon, 1, ((q-1)*ny + nu)*w );
   
%end Outputs

function SetInputPortSamplingMode(block, port, mode)
    block.InputPort(port).SamplingMode = mode;
    
%end SetInputPortSamplingMode

function Terminate(block)

%end Terminate

