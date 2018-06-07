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

  Histórico de modificações:
---------------------------------------------------------------------------
- 24/05/2018 - Zoé Magalhães
- Início do controle de versão.
---------------------------------------------------------------------------
%}

clear 
%% Configuração


%Horizonte de predição 
Horizon = 10
% Ponderações da função custo
Q_Y = diag([4,1,0.001,0.001]); % Ponderação do erro em y
Q_U = 1e-12;                   % Ponderação do custo do controle

% Restrições - defina como [] para testar sem restrições
C_C = [ 1 0 0 0 ];

lb = -pi/18; %Valor mínimo dos estados
ub = -lb;                          %Valor máximo dos estados
slew_lb = -1000;                    %Variação mínima do comando
slew_ub = 1000;                     %Variação máxima do comando
cmd_lb = -0.8e5;                     %Valor mínimo do comando
cmd_ub = 0.8e5;                      %Valor máximo do comando

% Manobra de teste
maneuver = 'doublelane'; %'doublelane'
max_steer = pi/8;

%% Parâmetros do modelo
a                   =  1.035;  %distância entre cg e eixo frontal
b                   =  1.655;  %distância entre cg e eixo traseiro
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
Izz                 = 3048.1;  %momento de inércia do eixo de guinada
Ixx                 = 744;     %momento de inércia do eixo de rolagem
Ixz                 = 50;      %produto de inércia dos eixos de rolagem e guinada
g                   = 9.81;    %aceleração gravitacional
l = a+b;                       %distância entre os eixos frontal e traseiro       
uspeed = 80/3.6;               %velocidade longitudinal de linearização
C_alpha = 90624*[1,1,1,1]/2;   %Conering stiffness

%Coeficiente utilizado para calcular a taxa de guinada desejada
Ku = m*((a/C_alpha(3))-(b/C_alpha(1)))/(l*l);
Ku = abs(Ku);

%% Modelo linear Mx'= A1x + B1u + E1steer
M1 = [      m*uspeed,    0,  -ms*hs,  0;
                   0,  Izz,    -Ixz,  0; 
       -ms*hs*uspeed, -Ixz,     Ixx,  roll_damping;
                   0,    0,       0,  1 ];

A1 = [-sum(C_alpha)         , -C_alpha*[a;a;-b;-b]/uspeed - m*uspeed  , 0  , 0;
      -C_alpha*[a;a;-b;-b]  , -C_alpha*[a*a;a*a;b*b;b*b]              , 0  , 0;
        0                   , ms*hs*uspeed                            , 0  , (ms*hs*g-roll_stiffness);
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


%% Criando um descrito do controle preditiv

mpcobj = ClassLinearMPC.settings(d_sys.a, d_sys.b(:,1), d_sys.c, ....
         d_sys.d(:,1), Q_U,Q_Y, C_C, Horizon, lb, ub, slew_lb, slew_ub, cmd_lb, cmd_ub,...
         d_sys.b(:,2) );



%% Carrega os sinais para o esterçamento das rodas dianteiras
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
%% Laço que executa as iterações da simulação
for i=1:size(steer,1)
    
    fprintf(repmat('\b',1,strlen));
    str = sprintf('%2.4f', steer(i,1) );
    fprintf(str);
    strlen = numel(str);
    
    % Calcula a taxa de guinada desejada 
    desired_yaw_rate(i) = ...
        steer(i,2)*uspeed/(l + Ku*l*uspeed*uspeed);
    if abs(desired_yaw_rate(i)) > 0.4 
        desired_yaw_rate(i) = 0.4*sign(steer(i,2));
    end
    
    % Calcula a trajetória desejada 
    TRACK = zeros(4,mpcobj.n);
    TRACK(2,:) = desired_yaw_rate(i);
    
    % Atualiza o sinal de comando
    u(i) = mpcobj.write_next_command(TRACK, steer(i,2) );
    
    % Atualiza os estados
    mpcobj.x = d_sys.a*mpcobj.x + d_sys.b*[u(i);steer(i,2)];
    
    %registra os estados para plotar depois
    x(:,i) = mpcobj.x;
end

save exe_teste