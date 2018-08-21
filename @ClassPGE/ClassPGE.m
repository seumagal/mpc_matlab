%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Classe que representa uma unidade de execução do algoritmo PGE para
  programação quadrática
---------------------------------------------------------------------------
%}
%%
classdef ClassPGE 
    properties
        config
        J
        sol
        solver
    end
    
    methods(Static)        


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
            @param[in] arg_beta_minus valor do coeficiente beta_minus do algoritmo
        %}
        obj = settings( arg_rho0, arg_rho_max, arg_n_rho, arg_N_inter, ...
                        arg_p0, arg_gam_min, arg_beta_plus, arg_beta_minus,...
                        arg_size_H, arg_size_F, arg_size_A, arg_size_B, arg_np )
    end 
    
    methods
        %{
            @brief Função que executa a programação quadrática pelo algoritmo
            PGE para minimizar a a função custo J = p'Hp + arg_F'p com as restrições
            Ap <= B, arg_p_min <= p <= arg_p_max;
            @param[in] arg_H é a matriz H da função custo
            @param[in] arg_F é a matriz arg_F da função custo
            @param[in] arg_A é matriz A da inequação de restrição
            @param[in] arg_B é a matriz B da inequação de restrição
            @param[in] arg_p_min é o limite inferior do espaço de busca
            @param[in] arg_p_max é o limite superior do espaço de bsuca

        
            solver(obj, arg_H,arg_F,arg_A,arg_B,arg_p_mi,arg_p_max)
        %}
    end
    
end
