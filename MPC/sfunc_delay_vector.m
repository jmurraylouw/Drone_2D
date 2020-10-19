function sfunc_delay_vector(block)
% Outputs a vector of delay coordinates of measurements
% e.g. with q = 3 delays:
% delays = [x(k-2); x(k-1)]
% Used with an MPC for a state-space model created with DMD/HAVOK

	setup(block);
  
%endfunction

function setup(block)
  
    block.NumDialogPrms  = 4;

    %% Register number of input and output ports
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 1;

    %% Setup functional port properties to dynamically
    %% inherited.
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    %% Extract Dialog params
    delays_0  = block.DialogPrm(1).Data; % For reference sake
    ny        = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;
    Ts_mpc    = block.DialogPrm(4).Data;

    %% Port dimentions
    block.InputPort(1).Dimensions        = ny;
    block.InputPort(1).DirectFeedthrough = false;

    block.OutputPort(1).Dimensions       = (q-1)*ny;

    %% Set block sample time to same as MPC
    block.SampleTimes = [Ts_mpc 0];

    %% Set the block simStateCompliance to default (i.e., same as a built-in block)
    block.SimStateCompliance = 'DefaultSimState';

    %% Register methods
    block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
    block.RegBlockMethod('InitializeConditions',    @InitConditions);  
    block.RegBlockMethod('Outputs',                 @Output);  
    block.RegBlockMethod('Update',                  @Update);  
  
%endfunction

function DoPostPropSetup(block)

    %% Extract Dialog params
    delays_0  = block.DialogPrm(1).Data; % For reference sake
    ny        = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;
    Ts_mpc    = block.DialogPrm(4).Data;

    %% Setup Dwork
    block.NumDworks = 1;
    block.Dwork(1).Name = 'delays_next'; % Previous x_ext
    block.Dwork(1).Dimensions      = ny*(q-1);
    block.Dwork(1).DatatypeID      = 0;
    block.Dwork(1).Complexity      = 'Real';
    block.Dwork(1).UsedAsDiscState = true;

%endfunction

function InitConditions(block)

    %% Extract Dialog params
    delays_0  = block.DialogPrm(1).Data; % For reference sake
    
    %% Initialize Dwork
    block.Dwork(1).Data = delays_0; % Initial delays
  
%endfunction

function Output(block)
    
    %% Extract Dialog params
    delays_0  = block.DialogPrm(1).Data; % For reference sake
    ny        = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;
    Ts_mpc    = block.DialogPrm(4).Data;

    delays = block.Dwork(1).Data; % Current delays vector, defined at previous timestep
    block.OutputPort(1).Data = delays;
  
%endfunction

function Update(block)
    %% Extract Dialog params
    delays_0  = block.DialogPrm(1).Data; % For reference sake
    ny        = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;
    Ts_mpc    = block.DialogPrm(4).Data;
    
    y = block.InputPort(1).Data; % Measurement vector
    delays = block.Dwork(1).Data; % Current delays vector, defined at previous timestep
    
    % Save delays vector for next iteration
    delays_next = [y; delays(1:length(delays)-ny)];
    block.Dwork(1).Data = delays_next; 
  
%endfunction

