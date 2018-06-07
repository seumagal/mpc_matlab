%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Este script implementa a função de configuração do descritor de controle
  preditivo.

  Histórico de modificações:
---------------------------------------------------------------------------
- 24/05/2018 - Zoé Magalhães
- Início do controle de versão.
- Contempla sistemas com e sem restrição
- Contempla sistemas com e sem perturbações mensuráveis. 
---------------------------------------------------------------------------
%}
%%

%{
@brief Função para configurar o controlador
@arg arg_A é a matriz A do espaço de estado discreto
@arg arg_B é a matriz B do espaço de estado discreto
@arg arg_C é a matriz C do espaço de estado discreto
@arg arg_D é a matriz D do espaço de estado discreto
@arg arg_Q_Y é a matriz de ponderação dos erros das saídas
             controladas na função custo
@arg arg_Q_U é a mamtriz de ponderação do vetor de controle na
             função custo
@arg arg_n é o horizonte de predição em número de amostras
@arg arg_lb é o vetor de restrição inferior dos estados
@arg arg_ub é o vetor de restrição superior dos estados
@arg arg_slew_lb é o vetor de restrição inferior da variação entre
                 duas amostras do vetor de comando.
@arg arg_slew_ub é o vetor de restrição superior da variação entre
                 duas amostras do vetor de comando
@arg arg_cmd_lb é o vetor de restrição inferior do vetor de
                comando
@arg arg_cmd_ub é o vetor de restrição superior do vetor de
                comando
@arg_E é a matriz E de ganho das perturbações na atualização dos
       estados. x[n] = A x[n-1] + Bu + Ev, em que v é o vetor
       de perturbações
%}
function obj = settings( arg_A, arg_B, arg_C, arg_D, arg_Q_U, arg_Q_Y, ...
                         arg_C_C, arg_n, arg_lb, arg_ub, ...
                         arg_slew_lb, arg_slew_ub, ...
                         arg_cmd_lb, arg_cmd_ub,...
                         arg_E )
% Verifica consitencia dos argumentos                  
    if size( arg_A, 1 ) ~= size ( arg_B, 1 )
        disp('erro: size( arg_A, 1 ) ~= size ( arg_B, 1 )');
        return
    end

    if size( arg_C, 2 ) ~= size( arg_Q_Y, 1 )
        disp('erro: size( arg_C, 2 ) ~= size( arg_Q_Y, 1 )');
        return
    end

    if size( arg_C, 2 ) ~= size( arg_A, 2 )
        disp('erro: size( arg_C, 2 ) ~= size( arg_A, 2 )');
        return
    end

    if size( arg_B, 2 ) ~= size( arg_Q_U, 1 )
        disp('erro: size( arg_B, 2 ) ~= size( arg_Q_U, 1 )');
        return
    end
    
    if size(arg_lb,1) ~= size(arg_ub,1)
        disp('erro: size(arg_lb,1) ~= size(arg_ub,1)');
        return
    end
    
    if size(arg_lb,1) ~= size(arg_C_C,1) && size(arg_lb,1) ~= 0
        disp('erro: size(arg_lb,1) ~= size(arg_C_C,1) && size(arg_lb,1) ~= 0');
        return
    end
    
    if size(arg_slew_lb,1) ~= size(arg_slew_ub,1)
        disp('erro: size(arg_slew_lb,1) ~= size(arg_slew_ub,1)');
        return
    end
    
    if size(arg_slew_lb,1) ~= size(arg_B,2) && size(arg_slew_lb,1) ~= 0
        disp('erro: size(arg_slew_lb,1) ~= size(arg_B,2) && size(arg_slew_lb,1) ~= 0');
        return
    end
    
    if size(arg_cmd_lb,1) ~= size(arg_cmd_ub,1)
        disp('erro: size(arg_cmd_lb,1) ~= size(arg_cmd_ub,1)');
        return
    end
    
    if size(arg_cmd_lb,1) ~= size(arg_B,2)  && size(arg_cmd_lb,1) ~= 0
        disp('erro: size(arg_cmd_lb,1) ~= size(arg_B,2)');
        return
    end
    
    if arg_n <= 0 
        disp('erro: arg_n <= 0');
        return
    end

 % Configura o descritor   
    obj = ClassLinearMPC;
    obj.A = arg_A;
    obj.B = arg_B;
    obj.Q_U = arg_Q_U;
    obj.Q_Y = arg_Q_Y;
    obj.C = arg_C;
    obj.D = arg_D;
    obj.n = arg_n;
    obj.c_PHI = cell(0);
    obj.c_PSI = cell(0);
    obj.C_C = arg_C_C;
    
    obj.nu = size(obj.B,2);
    obj.nx = size(obj.A,2);
    obj.ny = size(obj.C,2);
    obj.nv = size(arg_E,2);
    obj.nc = size(obj.C_C,1);
    
    
    % Se existe perturbação
    if size(arg_E,1) == size(arg_B,1)
        obj.B = [ obj.B, arg_E ];
        obj.Q_U = [ obj.Q_U, zeros(obj.nu,obj.nv);
                    zeros(obj.nv, obj.nu + obj.nv) ];
        obj.nu = obj.nu + obj.nv;
        obj.D = [ obj.D, zeros(obj.ny, obj.nv); zeros(obj.nv, obj.nu ) ];
    end
    
    obj.H  = zeros(obj.nu*obj.n );
    obj.F1 = zeros(obj.nu*obj.n, obj.nx);
    obj.F2 = zeros(obj.nu*obj.n, obj.ny*obj.n);
    obj.F3 = zeros(obj.nu*obj.n, obj.nu);
    
    obj.AINEQ = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu , obj.nu*obj.n );
    obj.G1 = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu, obj.nx );
    obj.G2 = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu, obj.nu - obj.nv )
    obj.G3 = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu, 1 );

    for i = 1:obj.n
        obj.c_PHI{i} = obj.A^i;
        obj.c_PSI{i} = zeros(obj.nx,obj.nu*obj.n);
        for j = 1:i
            sel = obj.sel_matrix(j, obj.nu);
            obj.c_PSI{i} = obj.c_PSI{i} + (obj.A^(i-j)*obj.B*sel);
        end

        sel = obj.sel_matrix(i,obj.nu);
        obj.H = obj.H + obj.c_PSI{i}'*obj.C'*obj.Q_Y*obj.C*obj.c_PSI{i} + sel'*obj.Q_U*sel;
        obj.F1 = obj.F1 + obj.c_PSI{i}'*obj.C'*obj.Q_Y*obj.C*obj.c_PHI{i};
        obj.F3 = obj.F3 + sel'*obj.Q_U;
        sel = obj.sel_matrix(i,obj.ny);
        obj.F2 = obj.F2 + obj.c_PSI{i}'*obj.C'*obj.Q_Y*sel;
        
        obj.AINEQ((i-1)*obj.nc + 1: i*obj.nc, : ) =  obj.C_C*obj.c_PSI{i};
        obj.AINEQ(2*obj.n*obj.nc+(i-1)*obj.nu+1:2*obj.n*obj.nc+i*obj.nu - obj.nv,(i-1)*obj.nu+1:i*obj.nu - obj.nv )= eye(obj.nu - obj.nv);
        if i > 1
            obj.AINEQ(2*obj.n*obj.nc+(i-1)*obj.nu+1:2*obj.n*obj.nc+i*obj.nu - obj.nv,(i-2)*obj.nu+1:(i-1)*obj.nu - obj.nv )= ...
                -eye(obj.nu - obj.nv);
        end
        obj.G1((i-1)*obj.nc + 1: i*obj.nc, : ) = - obj.C_C*obj.c_PHI{i};
    end

    obj.AINEQ( obj.n*obj.nc + 1: 2*obj.n*obj.nc , : ) = - obj.AINEQ( 1: obj.n*obj.nc , : );
    obj.AINEQ( 2*obj.n*obj.nc + obj.n*obj.nu + 1 : 2*obj.n*obj.nc + 2*obj.n*obj.nu, : ) = ...
        - obj.AINEQ( 2*obj.n*obj.nc + 1 : 2*obj.n*obj.nc + obj.n*obj.nu, : );
    obj.G1( obj.n*obj.nc + 1: 2*obj.n*obj.nc , : ) = - obj.G1( 1: obj.n*obj.nc , : );
    
    nu =  obj.nu - obj.nv; 
    obj.G2( 2*obj.n*obj.nc + 1 : 2*obj.n*obj.nc + nu, 1 : nu  ) = eye(nu);
    obj.G2(2*obj.n*obj.nc+obj.n*obj.nu+1:2*obj.n*obj.nc+obj.n*obj.nu+nu,1:nu)=-eye(nu);

    obj.H = 2*obj.H;
    obj.F1 = 2*obj.F1;
    obj.F2 = -2*obj.F2;
    obj.F3 = 2*obj.F3; 
    PI =  obj.sel_matrix( 1, obj.nu );
    AUX = PI/obj.H;
    obj.K =  AUX*obj.F1;
    obj.G = -AUX*obj.F2;
    obj.L = -AUX*obj.F3;
    obj.constrained = 0;
    obj.u = zeros(obj.nu,1);

    if ( size(arg_lb,1) == obj.nc ) && ...
       ( size ( arg_ub,1 ) == obj.nc ) && ...
       ( size(arg_slew_lb,1) == obj.nu - obj.nv ) && ...
       ( size(arg_slew_ub,1) == obj.nu - obj.nv )
            obj.cmd_lb = arg_cmd_lb;
            obj.cmd_ub = arg_cmd_ub;
            state_ub = arg_ub;
            state_lb = arg_lb;
            obj.G3(1:obj.n*obj.nc,1) = repmat(state_ub,obj.n,1);
            obj.G3(obj.n*obj.nc+1:2*obj.n*obj.nc,1) = - repmat(state_lb,obj.n,1);
            obj.G3(2*obj.n*obj.nc + 1: 2*obj.n*obj.nc + obj.n*obj.nu ) = repmat([arg_slew_ub; zeros(obj.nv,1) ],obj.n,1);
            obj.G3(2*obj.n*obj.nc+obj.n*obj.nu+1:2*obj.n*obj.nc+2*obj.n*obj.nu)=-repmat([arg_slew_lb; zeros(obj.nv,1)],obj.n,1);

            obj.constrained = 1;
    end
end   