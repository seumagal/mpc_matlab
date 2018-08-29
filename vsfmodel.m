function  vsfodel(block)
	setup(block);

function setup(block)
   block.NumDialogPrms     = 10; 
   block.SetPreCompInpPortInfoToDynamic;
   block.SetPreCompOutPortInfoToDynamic;

	%Inputs
	block.NumInputPorts = 2;

   % In1: Front wheels steer angle
   block.InputPort(1).Complexity         = 'Real'; 
   block.InputPort(1).DataTypeId         = 0; %double
   block.InputPort(1).SamplingMode       = 'Sample';
   block.InputPort(1).Dimensions         = 1;
   block.InputPort(1).DirectFeedthrough  = true;
	% In2: Ext Yaw Moment
   block.InputPort(2).Complexity         = 'Real'; 
   block.InputPort(2).DataTypeId         = 0;
   block.InputPort(2).SamplingMode       = 'Sample';
   block.InputPort(2).Dimensions         = 1;
   block.InputPort(2).DirectFeedthrough  = true;

   	% In2: Ext Yaw Moment
   block.InputPort(2).Complexity         = 'Real'; 
   block.InputPort(2).DataTypeId         = 0;
   block.InputPort(2).SamplingMode       = 'Sample';
   block.InputPort(2).Dimensions         = 1;
   block.InputPort(2).DirectFeedthrough  = true;

   	% In2: Wheel torques
   block.InputPort(2).Complexity         = 'Real'; 
   block.InputPort(2).DataTypeId         = 0;
   block.InputPort(2).SamplingMode       = 'Sample';
   block.InputPort(2).Dimensions         = 1;
   block.InputPort(2).DirectFeedthrough  = true;

	
	%Outputs
	block.NumOutputPorts = 6;
   % Out1: Vehicle states
   block.OutputPort(1).complexity        = 'real'; 
   block.OutputPort(1).datatypeid        = 0;
   block.OutputPort(1).samplingmode      = 'sample';
   %block.OutputPort(1).name              = 'states';                          
   block.OutputPort(1).dimensions        = 5;
	% Out2: Wheels load
   block.OutputPort(2).complexity        = 'real'; 
   block.OutputPort(2).datatypeid        = 0;
   block.OutputPort(2).samplingmode      = 'sample';
   %block.OutputPort(2).name              = 'wheels load' ;
   block.OutputPort(2).dimensions        = 4;
	% Out3: Wheels lat force
   block.OutputPort(3).complexity        = 'real'; 
   block.OutputPort(3).datatypeid        = 0;
   block.OutputPort(3).samplingmode      = 'sample';
   %block.OutputPort(3).name              = 'wheels lat force'; 
   block.OutputPort(3).dimensions        = 4;
	% Out4: Wheels long force
   block.OutputPort(4).complexity        = 'real'; 
   block.OutputPort(4).datatypeid        = 0;
   block.OutputPort(4).samplingmode      = 'sample';
   %block.OutputPort(4).name              = 'wheels long force';
   block.OutputPort(4).dimensions        = 4;
	% Out5: Wheels slip angle
   block.OutputPort(5).complexity        = 'real'; 
   block.OutputPort(5).datatypeid        = 0;
   block.OutputPort(5).samplingmode      = 'sample';
   %block.OutputPort(5).name              = 'wheels slip angle';
   block.OutputPort(5).dimensions        = 4;

   block.OutputPort(6).complexity        = 'real'; 
   block.OutputPort(6).datatypeid        = 0;
   block.OutputPort(6).samplingmode      = 'sample';
   %block.OutputPort(6).name              = 'wheels slip angle';
   block.OutputPort(6).dimensions        = 5;
	
   block.NumContStates = 5;

   %% Set block sample time to variable sample time
   block.SampleTimes = [-1 0];

   %% Set the block simStateCompliance to default (i.e., same as a built-in block)
   block.SimStateCompliance = 'DefaultSimState';

   %% Register methods
	  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
	  block.RegBlockMethod('Outputs',                 @Output);  
	  block.RegBlockMethod('Derivatives',             @Derivative);  
	  block.RegBlockMethod('Update',                  @Update);
	  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup); 
	%end	

function DoPostPropSetup(block)
	block.NumDworks=4;
	% wheel load
	  block.Dwork(1).Name            = 'wheel_load'; 
	  block.Dwork(1).Dimensions      = 4;
	  block.Dwork(1).DatatypeID      = 0;
	  block.Dwork(1).Complexity      = 'Real';
	  block.Dwork(1).UsedAsDiscState = true;
	% wheel lat
	  block.Dwork(2).Name            = 'wheel_lat'; 
	  block.Dwork(2).Dimensions      = 4;
	  block.Dwork(2).DatatypeID      = 0;
	  block.Dwork(2).Complexity      = 'Real';
	  block.Dwork(2).UsedAsDiscState = true;
	% wheel long
	  block.Dwork(3).Name            = 'wheel_long'; 
	  block.Dwork(3).Dimensions      = 4;
	  block.Dwork(3).DatatypeID      = 0;
	  block.Dwork(3).Complexity      = 'Real';
	  block.Dwork(3).UsedAsDiscState = true;
	% wheel slip
	  block.Dwork(4).Name            = 'wheel_slip'; 
	  block.Dwork(4).Dimensions      = 4;
	  block.Dwork(4).DatatypeID      = 0;
	  block.Dwork(4).Complexity      = 'Real';
	  block.Dwork(4).UsedAsDiscState = true;
      
function InitConditions(block)
	block.ContStates.Data = zeros(5,1);
	block.ContStates.Data(2) = block.DialogPrm(11).Data; %long speed

    
function Update(block)

	MFA = block.DialogPrm(1).Data;
	MFB = block.DialogPrm(2).Data;
	z   = block.DialogPrm(3).Data;
	zf0 = block.DialogPrm(4).Data;
	zr0 = block.DialogPrm(5).Data;
	a   = block.DialogPrm(6).Data;
	b   = block.DialogPrm(7).Data;
	tf2  = block.DialogPrm(8).Data;
    tr2 = block.DialogPrm(9).Data;
    long_slip = 0.001*ones(1,4);
	
	w_slip_angle = [block.InputPort(1).Data, block.InputPort(1).Data, 0, 0];
	num 	     = block.ContStates.Data(1) + block.ContStates.Data(3)*[a,a,-b,-b];
	den 	     = block.ContStates.Data(2) + block.ContStates.Data(3)*[-tf2,tf2,-tr2,tr2];
	w_slip_angle = w_slip_angle - atan2(num,den);
    
	if isempty(block.Derivatives.Data)
		lat_accel = 0;
		long_accel = 0;
	else
		lat_accel = block.Derivatives.Data(1);
		long_accel  = block.Derivatives.Data(2);
	end

	load_f  = [zf0,zf0,zr0,zr0] +...
              [-z(2),z(2),-z(7),z(7)]*lat_accel +...
              [-1,-1,1,1]*z(1)*long_accel +...
              [-z(3),z(3),-z(5),z(5)]*block.ContStates.Data(5) +...
              [-z(4),z(4),-z(6),z(6)]*block.ContStates.Data(4);
          
    block.Dwork(1).Data = load_f;
	load_f = load_f/1000;
    % Forca longitudinal das rodas
	aux1 				  = MFB(2)*load_f + MFB(3);
	aux2                   = (((MFB(4)*load_f + MFB(5)).*exp(-MFB(6)*abs(load_f)))./(MFB(1)*aux1)).*long_slip;
    aux3 					  = load_f.*(MFB(7)*load_f + MFB(8)) + MFB(9);
	long_f                = load_f.*aux1.*sin( MFB(1)*atan( (1-aux3).*aux2 + aux3.*atan( aux2 ) ) );

	% Forca lateral das rodas
	aux1 			      = MFA(2)*load_f + MFA(3);
	aux2 				  = MFA(7)*load_f+MFA(8);
	aux3 				  = (MFA(4)*sin(2*atan2(load_f,MFA(5)))./(MFA(1)*load_f.*aux1)).*w_slip_angle;
	lat_f                 = aux1.*load_f.* sin( MFA(1)*atan( (1-aux2).*aux3 + aux2.*atan( aux3 ) ) );
    
	block.Dwork(2).Data = lat_f;
	block.Dwork(3).Data = long_f;
	block.Dwork(4).Data = w_slip_angle;

function Derivative(block)
	a           = block.DialogPrm(6).Data;
	b           = block.DialogPrm(7).Data;
	tf2         = block.DialogPrm(8).Data;
    tr2         = block.DialogPrm(9).Data;
	d           = block.DialogPrm(10).Data;
	
	long_f 		= block.Dwork(3).Data;
	lat_f  		= block.Dwork(2).Data;	
	    
	rotate 		= [cos(block.InputPort(1).Data), -sin(block.InputPort(1).Data); sin(block.InputPort(1).Data), cos(block.InputPort(1).Data)];
	auxf 		= rotate*[long_f(1);lat_f(1)];
	long_f(1) 	= auxf(1,1);
	lat_f(1) 	= auxf(2,1);
	auxf 		= rotate*[long_f(2);lat_f(2)];
	long_f(2) 	= auxf(1,1);
	lat_f(2) 	= auxf(2,1);
    
	Fy=sum(lat_f);
	Fx=sum(long_f);
	Mz = [a,a,-b,-b]*lat_f + [-tf2,tf2,-tr2,tr2]*long_f + block.InputPort(2).Data;
	
	dx(5) = block.ContStates.Data(4);
	dx(4) = d(1)*Fy - d(2)*Mz + d(3)*sin(block.ContStates.Data(5)) - d(4)*block.ContStates.Data(5) - d(5)*block.ContStates.Data(4);
	dx(3) = d(6)*Mz - d(7)*dx(4);
	dx(2) = d(8)*Fx + d(9)*block.ContStates.Data(3)*block.ContStates.Data(1);
	dx(1) = d(10)*Fy + d(11)*dx(4) - d(12)*block.ContStates.Data(3)*block.ContStates.Data(2);

	block.Derivatives.Data = dx;

function Output(block)
	%States
	block.OutputPort(1).Data = block.ContStates.Data(1:5);
	%Load
	block.OutputPort(2).Data = block.Dwork(1).Data;
	%Lat
	block.OutputPort(3).Data = block.Dwork(2).Data;
	%long
	block.OutputPort(4).Data = block.Dwork(3).Data;
	%wheel slip	
	block.OutputPort(5).Data = block.Dwork(4).Data;	

    
	% Força normal nas rodas
	if isempty(block.Derivatives.Data)
		block.OutputPort(6).Data = zeros(5,1);
	else
		block.OutputPort(6).Data = block.Derivatives.Data(1:5);
    end
    