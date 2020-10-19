function add_filler_reference(block)
% Adds zeros to create correct dimension reference signal for MPC block
% Connect input to y_reference and output to ref port of MPC
% e.g. with q = 3 delays:
% ref_ext = [ref(k); 0; 0]
% Used with an MPC for a state-space model created with DMD
% Parameters: ny, q, Ts_mpc

	setup(block);
  
%endfunction

function setup(block)
  
    block.NumDialogPrms  = 3;

    %% Register number of input and output ports
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 1;

    %% Setup functional port properties to dynamically
    %% inherited.
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    %% Extract Dialog params
    ny        = block.DialogPrm(1).Data;
    q         = block.DialogPrm(2).Data;
    Ts_mpc    = block.DialogPrm(3).Data;

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
%endfunction

function InitConditions(block)
%endfunction

function Output(block)
    %% Extract Dialog params
    ny        = block.DialogPrm(1).Data;
    q         = block.DialogPrm(2).Data;
    
    ref = block.InputPort(1).Data; % Reference vector
    ref_ext = [ref; zeros((q-1)*ny,1)]; % Extended reference for output
    
    block.OutputPort(1).Data = ref_ext;
  
%endfunction

function Update(block) 
%endfunction

