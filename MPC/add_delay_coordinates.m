function add_delay_coordinates(block)
% Creates an extended state vector with current state and delay coordinates
% e.g. with q = 3 delays:
% x_ext = [x(k-2); x(k-1); x(k)]
% Used with an MPC for a state-space model created with HAVOK

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
    %   x_ext_0   = block.DialogPrm(1).Data; % For reference sake
    ny         = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;
    Ts_mpc    = block.DialogPrm(4).Data;

    %% Port dimentions
    block.InputPort(1).Dimensions        = ny;
    block.InputPort(1).DirectFeedthrough = true;

    block.OutputPort(1).Dimensions       = q*ny;

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
    ny        = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;

    %% Setup Dwork
    block.NumDworks = 1;
    block.Dwork(1).Name = 'x_ext_prev'; % Previous x_ext
    block.Dwork(1).Dimensions      = q*ny;
    block.Dwork(1).DatatypeID      = 0;
    block.Dwork(1).Complexity      = 'Real';
    block.Dwork(1).UsedAsDiscState = true;

%endfunction

function InitConditions(block)

    %% Initialize Dwork
    block.Dwork(1).Data = block.DialogPrm(1).Data; % x_ext_0
  
%endfunction

function Output(block)
    %% Extract Dialog params
    ny        = block.DialogPrm(2).Data;

    x = block.InputPort(1).Data; % Measurement vector
    x_ext_prev = block.Dwork(1).Data; % Previous xtended measurement vector
    x_ext = [x_ext_prev((ny+1):length(x_ext_prev)); x] % Extended measurement vector for output, discard oldest delay, add new measurement
    
    block.OutputPort(1).Data = x_ext;
    
    block.Dwork(1).Data = x_ext; % Save current x_ext for next iteration
  
%endfunction

function Update(block)

  
%endfunction

