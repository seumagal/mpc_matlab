%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Este script executa uma simulação para o controle de estabilidade lateral 
  veicular mediante aplicação do controle preditivo para sistemas LTI,
  com perturbação mensurável (entrada não controlada) e restrições
  para os estados e sinal de controle.
%}

clear 
%% Configuração

    %Horizonte de predição
    simlinear = [];
    openloop = [];
    Horizon = 101;
    Ts = 8e-5; %Período de amostragem
    
    %param_coeff =   [1 2 5 10 20 ]; %Parametrização geral Nij
    param_coeff = [ 1000; 4; 100 ]; %Parametrização exponencial alpha; ne 
    % Ponderações da função custo
    Q_Y = diag(2); % Ponderação do erro em y
    Q_U = 1e-12;   % Ponderação do custo do controle

    % Restrições 
    C_C =  [ 1 0 0 0 ]; %Restringe apenas o primeiro estado (deslizamento)

    lb = -pi/18;        %Val   or mínimo do estado regulado (slip angle)
    ub = -lb;           %Valor máximo do estado regulado
    rate_ub = 10*pi*Ts/180; %Restricao superior do slip rate
    rate_lb = -rate_ub; %Restricao superior do slip rate
    slew_lb = -1200;     %Diferença mínima entre duas amostras consecutivas do comando
    slew_ub = 1200;      %Diferença máxima entre duas amostras consecutivas do comando
    cmd_lb = -2000;      %Valor máximo do comandos
    cmd_ub = 2000;       %Valor mínimo do comando
    rate_factor = 14;   %Fator utilizado na atualizacao dinamica da restricao do slip rate
    
    rate_th =  pi*Ts/180; %Limiar para nao restringir a variacao do slip_rate em zero.
     
    csvwrite('cmd_ub.csv',cmd_ub);
    csvwrite('cmd_lb.csv',cmd_lb);
    % Manobra de+ teste
    maneuver = 'doublelane'; %'doublelane'
    max_steer = 10*pi/180;

%% Parâmetros do modelo
    a                   =  1.1;  %distância entre cg e eixo frontal
    b                   =  1.3;  %distância entre cg e eixo traseiropl
    tf                  =  1.4;   %comprimento do eixo frontal
    tr                  =  1.41;   %comprimento do eixo traseiro
    hg                  =  0.6;  %altura do cg
    m                   =  1070; %massa
    hs                  =  0.55;  %altura do centro de rolagem
    ms                  =  900; %massaa sobre centro de rolagem
    kf					=  65590;  %coeficiente de rigidez frontal a rolagem
    kr 					=  65590;  %coeficiente de rigidez traseira a rolagem
    cf					=  2100;   %coeficiente de amortecimento frontal a rolagem
    cr 					=  2100;   %coeficiente de amortecimental traseiro a rolagem
    roll_damping        =  cf+cr;   %coeficiente de rigidez a rolagem
    roll_stiffness      =  kf+kr;   %coeficiente de amortecimento a rolagem
    Izz                 =  2100;  %momento de inércia do eixo de guinada
    Ixx                 =  500;     %momento de inércia do eixo de rolagem
    Ixz                 =  47.5;      %produto de inércia dos eixos de rolagem e guinada
    g                   =  9.80665 ;    %aceleração gravitacional
    l = a+b;                       %distância entre os eixos frontal e traseiro       
    uspeed = 80/3.6;               %velociade longitudinal de linearização
    MFA                 = [ 1.6  ,-34 , 1250 ,  2320  , 12.8 ,     0 , -0.0053 , 0.1925         ];
 	MFB                 = [ 1.55 ,0   , 1000 ,   60  , 300  , 0.17 , 0      , 0      ,0.2 ];
	Reff                 = 0.307;

    C_alpha = MFA(4)*sind( 2*atand( ([b,b,a,a]*m*g/(2000*l))/MFA(5) ) )*180/pi;    
    Ku = m*(b*C_alpha(3) - a*C_alpha(2))/(l*C_alpha(1)*C_alpha(3));
    %Ku = 0.06;  %Coeficiente utilizado para calcular a taxa de guinada desejada

    l = a+b;
	pneumatic_trail     = 0;
	u = 80/3.6;
	friction_coefficient = 0.9;
    J = 1;
    Ts = 8e-5;
    
	eta = m*Ixx*Izz - m*Ixz*Ixz - ms*ms*hs*hs*Izz;
    

	d(1) = ms*hs*Izz;
	d(2) = m*Ixz;
	d(3) = m*Izz*ms*hs*g;
	d(4) = m*Izz*roll_stiffness;
	d(5) = m*Izz*roll_damping;
	d(1:5) = d(1:5)/eta;
	d(6) = 1/Izz;
	d(7) = Ixz/Izz;
	d(8) = 1/m;
    d(9) = 1;
    d(10) = ms*hs/m;
	d(11) = 1/m;
	d(12) = ms*hs/m;
	d(13) = 1;

	zf0 = m*g*b/(2*l);
	zr0 = m*g*a/(2*l);
	z(1) = m*hg*b/(2*l);
	z(2) = m*hg/(l*tf);
	z(3) = kf/tf;
	z(4) = cf/tf;
	z(5) = kr/tr;
	z(6) = cr/tr;
	z(7) = m*hg*a/(2*l);
	tf2 = tf/2;
    tr2 = tr/2;
    C_1 = 4;
    M_1 = 0.5;
    M_2 = 0.5;
    GAMMA = 8;    
    
    %states_max = [14.55 23.80 0.35 0.01 0.02]
    %states_min = [-12.31 6.18 -0.26 -0.01 -0.02];
    dac_max = [3  ,17 , 0.5 , 0.2 , 0.2 ,1.5]
    dac_min = [-7 ,5  ,-1.5 ,-0.4 ,-0.4 ,-3];

	params = [MFA,MFB,z,zf0,zr0,a,b,tf2,tr2,d,Reff,J];
    
%% Modelo simlinear Mx'= A1x + B1u + E1steer
    M1 = [      m*uspeed,    0,  -ms*hs,  0;
                       0,  Izz,    -Ixz,  0;
                -ms*hs*uspeed, -Ixz,     Ixx,  roll_damping;
                       0,    0,       0,  1 ];

    A1 = [-sum(C_alpha)         , -(C_alpha*[a;a;-b;-b]/uspeed) - m*uspeed  , 0  , 0;
          -C_alpha*[a;a;-b;-b]  , -C_alpha*[a*a;a*a;b*b;b*b]/uspeed              , 0  , 0;
            0                   ,  ms*hs*uspeed                           , 0  , (ms*hs*g-roll_stiffness);
            0                   , 0                                       , 1  , 0 ];

    B1 = [0; 1; 0; 0];  

    E1 = [ C_alpha*[1;1;0;0]; C_alpha*[a;a;0;0];0;0 ];

%% Modelo espaço de estado: x` = Ax + Bu + Esteer y'= Cx + Du
    A = M1\A1;
    B = M1\B1;
    C = eye(size(A));
    E = M1\E1;
    D = zeros(size([B,E]));

%% Modelo discreto
    sys = ss( A,[B,E]  ,C,D);
    d_sys = c2d(sys,8e-4);
    
    csvwrite('dA.csv',d_sys.a);
    csvwrite('dB.csv',d_sys.b);
    
    abs(eig(d_sys.a))
    
    %codegen -report model_update.m -args {zeros(size(d_sys.a)), zeros(size(d_sys.b)), zeros(4,1),0,0}

%% Criando um descritor do controle preditivo
    dA = d_sys.a;
    dB = d_sys.b;
    dC = d_sys.c;
    dD = d_sys.d;
    
    slew_sub = 0;
    slew_slb = 0;
    mpcobj = ClassLinearMPC.settings(d_sys.a, d_sys.b(:,1), d_sys.c(2,:), d_sys.d(:,1),...
             Q_U,Q_Y, C_C, Horizon, ...
             lb, ub, ... %slew_sub, slew_slb, ...
             slew_lb, slew_ub, cmd_lb, cmd_ub,...
             d_sys.b(:,2), ...
             'exponencial', param_coeff, 8e-4 );
    G1 = mpcobj.G1;
    G2 = mpcobj.G2;
    G3 = mpcobj.G3;
    G4 = mpcobj.G4;
    F1 = mpcobj.F1;
    F2 = mpcobj.F2;
    F3 = mpcobj.F3;
    np = mpcobj.np;
    nc = mpcobj.nc;
    nu = mpcobj.nu;
    N = Horizon;
    
    AINEQ = mpcobj.AINEQ;
    H = mpcobj.H;
    PI = mpcobj.PI_E;
    if isempty( simlinear )
        clear mpcobj
        clear sys
        clear d_sys
    end
%% Carrega os sinais para o esterçamento das rodas dianteiras
    load simin;
    load dlcinput;
    switch maneuver 
        case 'fishhook'
        disp('fishhook')
        steer = simin(1:8:300001,:);

        case 'doublelane'
        disp('doublelane')
        steer = dlcinput;
        steer(:,1) = steer(:,1);

        otherwise
            disp('entrada nao implementada')

    end

    steer(:,2) = steer(:,2)*max_steer;

    csvwrite('steer.csv',steer);
    
%% Estado inicial do sistema
    mpcobj.u = 0;
    mpcobj.u_d = 0;
    mpcobj.x = zeros(4,1);

    x = zeros(4,size(steer,1));
    desired_yaw_rate = zeros(size(steer(:,1)));
    u = zeros(size(steer(:,1)));

    strlen = 0; 
    
    save modelParam
    
    if isempty(simlinear)
        return;
    end
%% Laço que executa as iterações da simulação
for i=1:size(steer,1)
   

    desired_yaw_rate(i) = ...
        steer(i,2)*uspeed/(l + Ku*l*uspeed*uspeed);

    TRACK = desired_yaw_rate(i)*ones(Horizon,1);
    

    sub = atan((( ub - mpcobj.x(1) )/(ub - lb))*14)*rate_ub;
    lub = atan((( mpcobj.x(1) - lb )/(ub - lb))*14)*rate_lb;
    sub( sub < pi*Ts/180 ) = pi*Ts/180;
    lub( lub > -pi*Ts/180 ) = -pi*Ts/180;
    ssub(i)=sub;
    slub(i)=lub;
    if isempty(openloop)
        [u(i), QP] = mpcobj.write_next_command(TRACK, steer(i,2), lub, sub );
    else
        u(i) = 0;
    end
    x(:,i) =  model_update_mex(d_sys.a,d_sys.b,mpcobj.x,u(i),steer(i,2));
    mpcobj.x = x(:,i);
    
end
qpOASES_sequence('c',QP);
figure()
plotresult
save exe_teste
