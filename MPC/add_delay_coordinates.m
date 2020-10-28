function add_delay_coordinates(block)
% Creates an extended state vector with current state and delay coordinates
% e.g. with q = 3 delays:
% x_ext = [x(k); x(k-1); x(k-2)]
% Used with an MPC for a state-space model created with HAVOK
% Parameters: y_ext_0, ny, q, Ts_mpc

	setup(block);
  
%endfunction

function setup(block)
  
    block.NumDialogPrms  = 5;

    %% Register number of input and output ports
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 1;

    %% Setup functional port properties to dynamically
    %% inherited.
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    %% Extract Dialog params
    %   y_ext_0   = block.DialogPrm(1).Data; % For reference sake
    y_rows    = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;
    Ts_mpc    = block.DialogPrm(4).Data;
    nx        = block.DialogPrm(5).Data;
    
    ny        = length(y_rows);
    
    %% Port dimentions
    block.InputPort(1).Dimensions        = nx;
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
    y_rows    = block.DialogPrm(2).Data;
    q         = block.DialogPrm(3).Data;
    ny        = length(y_rows);

    %% Setup Dwork
    block.NumDworks = 1;
    block.Dwork(1).Name = 'y_ext_prev'; % Previous x_ext
    block.Dwork(1).Dimensions      = q*ny;
    block.Dwork(1).DatatypeID      = 0;
    block.Dwork(1).Complexity      = 'Real';
    block.Dwork(1).UsedAsDiscState = true;

%endfunction

function InitConditions(block)

    %% Initialize Dwork
    block.Dwork(1).Data = block.DialogPrm(1).Data; % y_ext_0
  
%endfunction

function Output(block)
    %% Extract Dialog params
    y_rows    = block.DialogPrm(2).Data;
    ny        = length(y_rows);

    x = block.InputPort(1).Data; % Measurement vector
    y = x(y_rows);
    
    y_ext_prev = block.Dwork(1).Data; % Previous xtended measurement vector
    y_ext = [y; y_ext_prev(1:length(y_ext_prev) - ny)]; % Extended measurement vector for output, discard oldest delay, add new measurement
    
    block.OutputPort(1).Data = y_ext;
    
    block.Dwork(1).Data = y_ext; % Save current x_ext for next iteration
  
%endfunction

function Update(block)

  
%endfunction

