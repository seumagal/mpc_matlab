%{
---------------------------------------------------------------------------
  Criador: Zo� Magalh�es (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matr�cula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: Andr� Murilo
  
  Este script executa uma simula��o para o controle de estabilidade lateral 
  veicular mediante aplica��o do controle preditivo para sistemas LTI,
  com perturba��o mensur�vel (entrada n�o controlada) e restri��es
  para os estados e sinal de controle.
%}

clear 
%% Configura��o

    %Horizonte de predi��o 
    Horizon = 101;

    %param_coeff =   [1 2 5 10 20 ]; %Parametriza��o geral Nij
    param_coeff = [ 1000; 3 ]; %Parametriza��o exponencial alpha; ne 
    % Pondera��es da fun��o custo
    Q_Y = diag(2); % Pondera��o do erro em y
    Q_U = 1e-12;   % Pondera��o do custo do controle

    % Restri��es 
    C_C = [ 1 0 0 0 ]; %Restringe apenas o primeiro estado (deslizamento)

    lb = -pi/36;        %Valor m�nimo do estado regulado
    ub = -lb;           %Valor m�ximo do estado regulado
    slew_lb = -2e4;     %Diferen�a m�nima entre duas amostras consecutivas do comando
    slew_ub = 2e4;      %Diferen�a m�xima entre duas amostras consecutivas do comando
    cmd_lb = -1e4;      %Valor m�ximo do comando
    cmd_ub = 1e4;       %Valor m�nimo do comando

    % Manobra de teste
    maneuver = 'doublelane'; %'doublelane'
    max_steer = pi/18;

%% Par�metros do modelo
    a                   =  1.035;  %dist�ncia entre cg e eixo frontal
    b                   =  1.655;  %dist�ncia entre cg e eixo traseiro
    tf                  =  1.52;   %comprimento do eixo frontal
    tr                  =  1.54;   %comprimento do eixo traseiro
    hg                  =  0.542;  %altura do cg
    m                   =  1704.7; %massa
    hs                  =  0.451;  %altura do centro de rolagem
    ms                  =  1526.9; %massaa sobre centro de rolagem
    kf					=  47298/4;  %coeficiente de rigidez frontal a rolagem
    kr 					=  47311/4;  %coeficiente de rigidez traseira a rolagem
    cf					=  2823;   %coeficiente de amortecimento frontal a rolagem
    cr 					=  2652;   %coeficiente de amortecimental traseiro a rolagem
    roll_damping        =  cf+cr;   %coeficiente de rigidez a rolagem
    roll_stiffness      =  kf+kr;   %coeficiente de amortecimento a rolagem
    Izz                 =  3048.1;  %momento de in�rcia do eixo de guinada
    Ixx                 =  744;     %momento de in�rcia do eixo de rolagem
    Ixz                 =  0;      %produto de in�rcia dos eixos de rolagem e guinada
    g                   =  9.80665 ;    %acelera��o gravitacional
    l = a+b;                       %dist�ncia entre os eixos frontal e traseiro       
    uspeed = 80/3.6;               %velociade longitudinal de lineariza��o
    MFA                 = [ 1.6  ,-34 , 1250 ,  2320  , 12.8 ,     0 , -0.0053 , 0.1925         ];

    C_alpha = MFA(4)*sin(2*atan(([b,b,a,a]*m*g/(2000*l))/MFA(5)));    
    Ku = 0.06;  %Coeficiente utilizado para calcular a taxa de guinada desejada

%% Modelo linear Mx'= A1x + B1u + E1steer
    M1 = [      m*uspeed,    0,  -ms*hs,  0;
                       0,  Izz,    -Ixz,  0;
           -ms*hs*uspeed, -Ixz,     Ixx,  roll_damping;
                       0,    0,       0,  1 ];

    A1 = [-sum(C_alpha)         , -(C_alpha*[a;a;-b;-b]/uspeed) - m*uspeed  , 0  , 0;
          -C_alpha*[a;a;-b;-b]  , -C_alpha*[a*a;a*a;b*b;b*b]              , 0  , 0;
            0                   , ms*hs*uspeed                            , 0  , (ms*hs*g-roll_stiffness);
            0                   , 0                                       , 1  , 0 ];

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
    abs(eig(d_sys.a))
    
    %codegen -report model_update.m -args {zeros(size(d_sys.a)), zeros(size(d_sys.b)), zeros(4,1),0,0}

%% Criando um descritor do controle preditivo

    mpcobj = ClassLinearMPC.settings(d_sys.a, d_sys.b(:,1), d_sys.c(2,:), ....
             d_sys.d(:,1), Q_U,Q_Y, C_C, Horizon, lb, ub, slew_lb, slew_ub, cmd_lb, cmd_ub,...
             d_sys.b(:,2), 'exponencial', param_coeff, 8e-4 );



%% Carrega os sinais para o ester�amento das rodas dianteiras
    load simin;

    switch maneuver 
        case 'fishhook'
        disp('fishhook')
        steer = simin(1:8:300001,:);

        case 'doublelane'
        disp('doublelane')
        steer = simin3(1:8:300001,:);

        otherwise
            disp('entrada nao implementada')

    end

    steer(:,2) = steer(:,2)*max_steer;

%% Estado inicial do sistema
    mpcobj.u = 0;
    mpcobj.u_d = 0;
    mpcobj.x = zeros(4,1);

    x = zeros(4,size(steer,1));
    desired_yaw_rate = zeros(size(steer(:,1)));
    u = zeros(size(steer(:,1)));

    strlen = 0; 
    

 
%% La�o que executa as itera��es da simula��o
for i=1:size(steer,1)
   
    %tic
    %disp(steer(i,1))
    
    %{
    fprintf(repmat('\b',1,strlen));
    str = sprintf('\n%2.4f', steer(i,1) );
    fprintf(str);
    strlen = numel(str);
    %}
    %disp('3')
    %tic
    % Calcula a taxa de guinada desejada 
    desired_yaw_rate(i) = ...
        steer(i,2)*uspeed/(l + Ku*l*uspeed*uspeed);
    %toc
    % Calcula a trajet�ria desejada 
    TRACK = desired_yaw_rate(i)*ones(Horizon,1);
    %disp('2')
    %tic
    % Atualiza o sinal de comando
    [u(i), QP] = mpcobj.write_next_command(TRACK, steer(i,2) );
   %toc;
    % Atualiza os estados
   % disp('1')
    %tic
    x(:,i) =  model_update_mex(d_sys.a,d_sys.b,mpcobj.x,u(i),steer(i,2));
    %toc
    mpcobj.x = x(:,i);
    %toc
end
qpOASES_sequence('c',QP);
figure()
plotresult
save exe_teste
