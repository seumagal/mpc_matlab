	clear 
    %% Par�metros do modelo
	a                   =  1.035;  %dist�ncia entre cg e eixo frontal
	b                   =  1.655;  %dist�ncia entre cg e eixo traseiro
    tf                  =  1.52;   %comprimento do eixo frontal
    tr                  =  1.54;   %comprimento do eixo traseiro
	hg                  =  0.542;  %altura do cg
	m                   =  1704.7; %massa
	hs                  =  0.451;  %altura do centro de rolagem
	ms                  =  1526.9; %massaa sobre centro de rolagem
	kf					=  47298;  %coeficiente de rigidez frontal a rolagem
	kr 					=  47311;  %coeficiente de rigidez traseira a rolagem
	cf					=  2823;   %coeficiente de amortecimento frontal a rolagem
	cr 					=  2652;   %coeficiente de amortecimental traseiro a rolagem
	roll_damping        = cf+cr;   %coeficiente de rigidez a rolagem
	roll_stiffness      = kf+kr;   %coeficiente de amortecimento a rolagem
	Izz                 = 3048.1;  %momento de in�rcia do eixo de guinada
	Ixx                 = 744;     %momento de in�rcia do eixo de rolagem
	Ixz                 = 50;      %produto de in�rcia dos eixos de rolagem e guinada
	g                   = 9.81;    %acelera��o gravitacional
    l = a+b;                       %dist�ncia entre os eixos frontal e traseiro       
	u = 80/3.6;                    %velocidade de lineariza��o
    C_alpha = 90624*[1,1,1,1]/2;   %Conering stiffness
    
    %Coeficiente utilizado para calcular a taxa de guinada desejada
	Ku = m*((a/C_alpha(3))-(b/C_alpha(1)))/(l*l);
    Ku = abs(Ku)
    
    %% Modelo linear Mx'= A1x + B1u + E1steer
    M1 = [      m*u,    0,  -ms*hs, 0;
	              0,  Izz,    -Ixz, 0; 
		   -ms*hs*u, -Ixz,    Ixx,  roll_damping;
		          0,    0,      0,  1 ];

	A1 = [-sum(C_alpha)         , -C_alpha*[a;a;-b;-b]/u - m*u  , 0             ,0;
	      -C_alpha*[a;a;-b;-b]  , -C_alpha*[a*a;a*a;b*b;b*b], 0             ,0;
			0                   , ms*hs*u                     , 0 ,(ms*hs*g-roll_stiffness);
			0                   , 0                           , 1             ,0];
	
    B1 = [0; 1; 0; 0];  
    
    E1 = [ C_alpha*[1;1;0;0]; C_alpha*[a;a;0;0];0;0 ];
    
    %% Modelo espa�o de estado: x` = Ax + Bu + Esteer y'= Cx + Du
	A = M1\A1;
	B = M1\B1;
	C = eye(size(A));
    E = M1\E1;
    D = zeros(size([B,E]));

    %% Modelo discreto
    sys = ss( A,[B,E]  ,C,D);
    d_sys = c2d(sys,8e-4);
     
    %% Pondera��es dos 
	Q = diag([4,1,0.001,0.001]);
	R = 1e-12;
	
    %% Restri��es
    lb = [ -pi/36; -pi; -pi; -pi/36 ]; %Valor m�nimo dos estados
    ub = -lb;                          %Valor m�ximo dos estados
    slew_lb = -1e4;                    %Varia��o m�nima do comando
    slew_ub = 1e4;                     %Varia��o m�xima do comando
    cmd_lb = -1e5;                     %Valpr m�nimo do comando
    cmd_ub = 1e5;                      %Valor m�ximo do comando
    
    %% Criando um descrito do controle preditivo
    mpcobj = ClassLinearMPC.settings(d_sys.a, d_sys.b(:,1), d_sys.c, ....
             d_sys.d(:,1), R,Q, 10, lb, ub, slew_lb, slew_ub, cmd_lb, cmd_ub,...
             d_sys.b(:,2), -pi, pi );
    
    %% Carrega os sinais para o ester�amento das rodas dianteiras
    load simin;
    % simin   fishhook
    % simin2  fishhook sem retornar ester�amento para zero graus
    % simin3  doublelane
    
    steer = simin(1:8:300001,:);
    steer(:,2) = steer(:,2)*pi/16;
    
    %% Estado inicial do sistema
    mpcobj.u = 0;
    mpcobj.u_d = 0;
    mpcobj.x = zeros(5,1);