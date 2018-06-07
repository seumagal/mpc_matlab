%{
---------------------------------------------------------------------------
  Criador: Zo� Magalh�es (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matr�cula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: Andr� Murilo
  
  Este script implementa a fun��o de configura��o do descritor de controle
  preditivo.

  Hist�rico de modifica��es:
---------------------------------------------------------------------------
- 24/05/2018 - Zo� Magalh�es
- In�cio do controle de vers�o.
- Contempla sistemas com e sem restri��o
- Contempla sistemas com e sem perturba��es mensur�veis. 
---------------------------------------------------------------------------
%}
%%

%{
@brief Fun��o para configurar o controlador
@arg arg_A � a matriz A do espa�o de estado discreto
@arg arg_B � a matriz B do espa�o de estado discreto
@arg arg_C � a matriz C do espa�o de estado discreto
@arg arg_D � a matriz D do espa�o de estado discreto
@arg arg_Q_Y � a matriz de pondera��o dos erros das sa�das
             controladas na fun��o custo
@arg arg_Q_U � a mamtriz de pondera��o do vetor de controle na
             fun��o custo
@arg arg_n � o horizonte de predi��o em n�mero de amostras
@arg arg_lb � o vetor de restri��o inferior dos estados
@arg arg_ub � o vetor de restri��o superior dos estados
@arg arg_slew_lb � o vetor de restri��o inferior da varia��o entre
                 duas amostras do vetor de comando.
@arg arg_slew_ub � o vetor de restri��o superior da varia��o entre
                 duas amostras do vetor de comando
@arg arg_cmd_lb � o vetor de restri��o inferior do vetor de
                comando
@arg arg_cmd_ub � o vetor de restri��o superior do vetor de
                comando
@arg_E � a matriz E de ganho das perturba��es na atualiza��o dos
       estados. x[n] = A x[n-1] + Bu + Ev, em que v � o vetor
       de perturba��es
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
    
    
    % Se existe perturba��o
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