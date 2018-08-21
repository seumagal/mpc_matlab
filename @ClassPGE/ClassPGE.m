%{
---------------------------------------------------------------------------
  Criador: Zo� Magalh�es (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matr�cula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: Andr� Murilo
  
  Classe que representa uma unidade de execu��o do algoritmo PGE para
  programa��o quadr�tica
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
            @biref Configura o algoritmo PGE para a progra��o quadr�tica 
            J = 0.5p' H p +F'p com restri��es A p<=B, p_min <= p <= p_max
            @param[in] arg_rho0 valor inicial da variavel rho do algoritmo
            @param[in] arg_rho_max limiar superior de rho
            @param[in] arg_n_rho comprimento do vetor rho
            @param[in] arg_N_iter m�ximo de intera��es para encontrar a solu��o
            @param[in] arg_p0 solu��o inicial
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
            @brief Fun��o que executa a programa��o quadr�tica pelo algoritmo
            PGE para minimizar a a fun��o custo J = p'Hp + arg_F'p com as restri��es
            Ap <= B, arg_p_min <= p <= arg_p_max;
            @param[in] arg_H � a matriz H da fun��o custo
            @param[in] arg_F � a matriz arg_F da fun��o custo
            @param[in] arg_A � matriz A da inequa��o de restri��o
            @param[in] arg_B � a matriz B da inequa��o de restri��o
            @param[in] arg_p_min � o limite inferior do espa�o de busca
            @param[in] arg_p_max � o limite superior do espa�o de bsuca

        
            solver(obj, arg_H,arg_F,arg_A,arg_B,arg_p_mi,arg_p_max)
        %}
    end
    
end
