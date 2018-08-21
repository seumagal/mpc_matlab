%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Este script implementa a função de configuração do descritor de controle
  preditivo.
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
@arg arg_E é a matriz E de ganho das perturbações na atualização dos
       estados. x[n] = A x[n-1] + Bu + Ev, em que v é o vetor
       de perturbações
@arg arg_parameterization é uma string informando se a parâmetrização é
- 'none' : trivial
- 'general': atualização do comando nos instantes definidos pelos coeficientes
- 'exponencial' : comando formado pela combinação de  exponenciais.
@arg arg_param_coefficient são os coeficientes da parâmetrização.
@arg arag_sampling_time é o período de amostragem
%}
function obj = settings( arg_A, arg_B, arg_C, arg_D, arg_Q_U, arg_Q_Y, ...
                         arg_C_C, arg_n, arg_lb, arg_ub, ...
                         arg_slew_lb, arg_slew_ub, ...
                         arg_cmd_lb, arg_cmd_ub,...
                         arg_E,  arg_parameterization, arg_param_coefficient, ...
                         arg_sampling_time )   
% Verifica consitencia dos argumentos                  
    if size( arg_A, 1 ) ~= size ( arg_B, 1 )
        disp('erro: size( arg_A, 1 ) ~= size ( arg_B, 1 )');
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
    obj.nx = size(obj.A,1);
    obj.ny = size(obj.C,1);
    obj.nv = size(arg_E,2);
    obj.nc = size(obj.C_C,1);
    
    alpha = 1000;
    
    % Se existe perturbação
    if size(arg_E,1) == size(arg_B,1)
        obj.B = [ obj.B, arg_E ];
        obj.Q_U = [ obj.Q_U, zeros(obj.nu,obj.nv);
                    zeros(obj.nv, obj.nu + obj.nv) ];
        obj.nu = obj.nu + obj.nv;
        obj.D = [ obj.D, zeros(obj.nx, obj.nv); zeros(obj.nv, obj.nu ) ];
    end
    
    
    switch arg_parameterization
       
        case 'general'
            coeff = arg_param_coefficient( arg_param_coefficient < obj.n & ...
                                           arg_param_coefficient > 0 );
            coeff = [ 0, coeff ]; 
            nr = length(coeff);
            np = nr*obj.nu;
            Pi_r = zeros(obj.nu*obj.n, np);
            for i=0:obj.n-1     
                if i < coeff(end)
                    [ Nji, ji ] = max(coeff(coeff<=i)); 
                    delta = (i - Nji)/(coeff(ji+1) - Nji);
                    selThis = obj.sel_matrix(ji,obj.nu, nr);
                    selNext = obj.sel_matrix(ji+1,obj.nu, nr);
                
                    Pi_r( i*obj.nu+1:(i+1)*obj.nu, : ) = ...
                        (eye( obj.nu )*(1 - delta))*selThis + delta*selNext;
                else
                    Pi_r( i*obj.nu+1:(i+1)*obj.nu, : ) = obj.sel_matrix(nr,obj.nu,nr);
                end
             
            end
            
        case 'exponencial'
            lambda = arg_param_coefficient( 1, : );
            ne = arg_param_coefficient( 2, : );
            ne = [ ne, ones(1,obj.nv)];
            lambda = [lambda, zeros(1,obj.nv)];
            Pi_e =[];
            for i=0:obj.n-1
                Mi=[];
                for j=1:obj.nu
                    ell = 1:ne(j);
                    Mji = exp( -lambda(j)*i*arg_sampling_time./((ell-1)*alpha+1));
                    Mi=blkdiag(Mi,Mji);
                end
                Pi_e = [Pi_e;Mi];    
            end
            
        case 'none'
            disp('parametrizacao trivial');
        otherwise
            disp('parametrizacao trivial');
    end 
    
    obj.H  = zeros(obj.nu*obj.n );
    obj.F1 = zeros(obj.nu*obj.n, obj.nx);
    obj.F2 = zeros(obj.nu*obj.n, obj.ny*obj.n);
    obj.F3 = zeros(obj.nu*obj.n, obj.nu);
    
    obj.AINEQ = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu , obj.nu*obj.n );
    obj.G1 = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu, obj.nx );
    obj.G2 = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu, obj.nu - obj.nv );
    obj.G3 = zeros( 2*obj.n*obj.nc + 2*obj.n*obj.nu, 1 );

    for i = 1:obj.n
        obj.c_PHI{i} = obj.A^i;
        obj.c_PSI{i} = zeros(obj.nx,obj.nu*obj.n);
        for j = 1:i
            sel = obj.sel_matrix(j, obj.nu, obj.n);
            obj.c_PSI{i} = obj.c_PSI{i} + (obj.A^(i-j)*obj.B*sel);
        end

        sel = obj.sel_matrix(i,obj.nu, obj.n);
        obj.H = obj.H + obj.c_PSI{i}'*obj.C'*obj.Q_Y*obj.C*obj.c_PSI{i} + sel'*obj.Q_U*sel;
        
        obj.F1 = obj.F1 + obj.c_PSI{i}'*obj.C'*obj.Q_Y*obj.C*obj.c_PHI{i};
        obj.F3 = obj.F3 + sel'*obj.Q_U;
        sel = obj.sel_matrix(i,obj.ny,obj.n);
        obj.F2 = obj.F2 + obj.c_PSI{i}'*obj.C'*obj.Q_Y*sel;
        
        obj.AINEQ((i-1)*obj.nc + 1: i*obj.nc, : ) =  obj.C_C*obj.c_PSI{i};
        obj.AINEQ(2*obj.n*obj.nc+(i-1)*obj.nu+1:2*obj.n*obj.nc+i*obj.nu - obj.nv,(i-1)*obj.nu+1:i*obj.nu - obj.nv)= eye(obj.nu - obj.nv);
        if i > 1
            obj.AINEQ(2*obj.n*obj.nc+(i-1)*obj.nu+1:2*obj.n*obj.nc+i*obj.nu - obj.nv ,(i-2)*obj.nu+1:(i-1)*obj.nu - obj.nv )= ...
                -eye(obj.nu - obj.nv);
        end
        obj.G1((i-1)*obj.nc + 1: i*obj.nc, : ) = -obj.C_C*obj.c_PHI{i};
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
    %PI =  obj.sel_matrix( 1, obj.nu, obj.n );
    %AUX = PI/obj.H;
    %obj.K =  AUX*obj.F1;
    %obj.G = -AUX*obj.F2;
    %obj.L = -AUX*obj.F3;
    obj.constrained = 0;
    obj.u = zeros(obj.nu,1);

    if ( size(arg_lb,1) == obj.nc ) &&...
       ( size ( arg_ub,1 ) == obj.nc ) && ...
       ( size(arg_slew_lb,1) == obj.nu - obj.nv ) && ...
       ( size(arg_slew_ub,1) == obj.nu - obj.nv )
            obj.cmd_lb = arg_cmd_lb;
            obj.cmd_ub = arg_cmd_ub;
                obj.G3(1:obj.n*obj.nc,1) = repmat(arg_ub,obj.n,1);
                obj.G3(obj.n*obj.nc+1:2*obj.n*obj.nc,1) = - repmat(arg_lb,obj.n,1);
                obj.G3(2*obj.n*obj.nc + 1: 2*obj.n*obj.nc + obj.n*obj.nu ) = repmat([arg_slew_ub; zeros(obj.nv,1) ],obj.n,1);
                obj.G3(2*obj.n*obj.nc+obj.n*obj.nu+1:2*obj.n*obj.nc+2*obj.n*obj.nu)=-repmat([arg_slew_lb; zeros(obj.nv,1)],obj.n,1);              
            obj.constrained = 1;
    end
    
    switch arg_parameterization
   
        case 'general'
            obj.H = Pi_r'*obj.H*Pi_r;
            obj.PI_R = Pi_r;
            obj.AINEQ = [ obj.AINEQ*Pi_r; -Pi_r; Pi_r];

            obj.np = nr*obj.nu;
            
        case 'exponencial'
            obj.H = Pi_e'*obj.H*Pi_e;
            obj.PI_E = Pi_e;
            obj.AINEQ = [ obj.AINEQ*Pi_e; -Pi_e; Pi_e];                
            obj.np = sum(ne);
        case 'none'
            obj.PI_R = [];
            obj.PI_E = [];
        otherwise
            obj.PI_R = [];
            obj.PI_E = [];
            obj.np = obj.nu * obj.n;


    end

    obj.H = (obj.H + obj.H')/2;
    
    csvwrite('H.csv',obj.H);
    csvwrite('PI_E.csv',obj.PI_E);
    csvwrite('A_CONSTRAINT.csv',obj.AINEQ);
    csvwrite('G1.csv',obj.G1);
    csvwrite('G2.csv',obj.G2);
    csvwrite('G3.csv',obj.G3);
    csvwrite('F1.csv',obj.F1);
    csvwrite('F2.csv',obj.F2);
    csvwrite('F3.csv',obj.F3);
    
    %Configura o solver da QP
    %{
    hmax0 = norm(obj.H,2);
    hmaxg = norm(obj.AINEQ'*obj.AINEQ,2);
    rho0 = 1*hmax0/hmaxg;
    rho_max = 10;
    n_rho = 10*obj.np;
    beta_plus = 1.1;
    beta_minus = 0.4;
    p0 = zeros( obj.np,1);
    gam_min = 2;
    N_iter = 50;
    obj.QP =  ClassPGE.settings( rho0, rho_max, n_rho, N_iter , p0, gam_min, ...
                                 beta_plus, beta_minus, size(obj.H), [obj.np,1],...
                                 size(obj.AINEQ),[size(obj.AINEQ,1),1],obj.np);
    %}
end   
