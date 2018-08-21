%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Função de configuração do algoritmo PGE  

---------------------------------------------------------------------------
%}
%%

%{
    @biref Configura o algoritmo PGE para a progração quadrática 
    J = 0.5p' H p +F'p com restrições A p<=B, p_min <= p <= p_max
    @param[in] arg_rho0 valor inicial da variavel rho do algoritmo
    @param[in] arg_rho_max limiar superior de rho
    @param[in] arg_n_rho comprimento do vetor rho
    @param[in] arg_N_iter máximo de interações para encontrar a solução
    @param[in] arg_p0 solução inicial
    @param[in] arg_gam_min limite inferior de gama
    @param[in] arg_beta_plus valor do coeficiente  beta_plus do algoritmo
    @param[in] arg_size_H tamanho da matriz H da funcao custo J = 0.5p'Hp + Fp 
    @param[in] arg_size_F tamanho da matriz F da funcao custo
    @param[in] arg_size_A tamanho da matriz A da restrição Ap<B
    @param[in] arg_size_B tamanho da matriz B da restrição Ap<B
    @param[in] arg_size_np quantidade de variáveis de decisão da QP.
%}
function obj = settings( arg_rho0, arg_rho_max, arg_n_rho, arg_N_iter, ...
                         arg_p0, arg_gam_min, arg_beta_plus, arg_beta_minus, ...
                         arg_size_H, arg_size_F, arg_size_A, arg_size_B, ...
                         arg_np )
    obj = ClassPGE;
    
    obj.config.rho0 = arg_rho0;
    obj.config.rho_max = arg_rho_max;
    obj.config.N_iter = arg_N_iter;
    obj.config.p0 = arg_p0;
    obj.config.gam_min = arg_gam_min;
    obj.config.beta_plus = arg_beta_plus;
    obj.config.beta_minus = arg_beta_minus;
    obj.config.n_rho = arg_n_rho;
    obj.config.rho = 0;

    tPGE = coder.typeof(obj.config);
    
    codegen -report ./@ClassPGE/solver_run.m -args {tPGE,zeros(arg_size_H),zeros(arg_size_F),zeros(arg_size_A),zeros(arg_size_B),zeros(arg_np,1),zeros(arg_np,1)}
    obj.solver = @solver_run_mex
    
end   
