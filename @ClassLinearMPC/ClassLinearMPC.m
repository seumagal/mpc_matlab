%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Este script define a classe que representa um descritor do alogritmo
  de controle preditivo.

---------------------------------------------------------------------------
%}
%%
classdef ClassLinearMPC < handle
    properties
        % Matriz da representação em espaço de estados discreto
        % x[n] = Ax[n-1] + Bu[n]
        % y[n] = Cx[n] + Du[n]
        A       
        B
        C       
        D
        % Ponderações da função custo
        Q_U   % Ponderação do comando    
        Q_Y   % Ponderação do erro nas saídas controladas
        % Matrizes da função custo
        % J = (1/2)tilde_u' H tilde_u + (F1x[n] + F_2 tilde_y_D + F_3 u_d)'tilde_u
        H
        F1
        F2
        F3
        
        % Células de matrizes auxiliares usadas no cálculo das matrizes
        % da função custo e da inequação que representa as restrições
        c_PHI
        c_PSI
        
        % Matrizes da lei de controle obtida da minimização da função
        % custo sem restrições
        % u = -K x[n] + G tilde_y_d + L u_d
        K
        G
        L

        x     % Último estado medido 
        u     % Último vetor de comando calculado
        u_d   % Vetor de comado desejado para ss
        n     % Horizonte de predição
        nu    % Comprimento do vetor de controle
        nx    % Comprimento do vetor de estados
        ny    % Comprimento do vetor de saída
        nc    % Comrpimento do vetor de saídas com restrição
        nv    % Comprimento do vetor de perturbações mensuráveis (entradas não controladas)
        constrained %flag que indica se o modelo tem restrições(=1)
        np    % Comprimento do comando parâmetrizado.
        
        %Matrizes da inequação que representa as restrições
        AINEQ 
        G1
        G2
        G3
        BINEQ
        C_C
        
        % Restrições no comando
        cmd_lb
        cmd_ub
        
        % Parametrização
        PI_R
        PI_E
        
        %
        QP
        
    end
    
    methods(Static)        
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
         @arg arg_C_C é a matriz Cc da equação da saída com restrição y_c = Cc x + Dc u
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

        obj = settings( arg_A, arg_B, arg_C, arg_D, arg_Q_U, ...
                        arg_Q_Y, arg_C_C, arg_n, arg_lb, arg_ub, ...
                        arg_slew_lb, arg_slew_ub, ...
                        arg_cmd_lb, arg_cmd_ub, ...
                        arg_E, arg_parameterization, arg_param_coefficient,...
                        arg_sampling_time, solver )         
    end 
    
    methods
        %{
         @brief Retorna a matriz de seleção PI, tal que se tilde_s for
         o vetor coluna resultante da concatenação das amostras do vetor 
         coluna s no horizonte de predição N: 
         tilde_s = [s[n+1]; s[n+2]; ... ; s[n + N] ] , 
         então: PI * tilde_s = s[n+arg_sel].
         @arg arg_obj é o objeto que invoca o método
         @arg arg_sel indica a amostra selecionada
         @arg arg_n é o comprimento das amostras
         @arg arg_N é o numero de amostras concatenadas
         @return out_sel_matrix é a matriz de seleção PI.
        %}
        out_sel_matrix = sel_matrix( obj, arg_sel, arg_n, arg_N )
        
        %{
         @brief 
         @arg arg_obj é o objeto que invoca o método
         @arg arg_TRACK é a trajetória a ser rastreada composta
         pela concatenação das amostras vetor coluna de saída desejada y_d no
         horizonte de predição.
         arg_TRACK = [y_d[n+1],y_d[n+2], ..., y_d[n+N]].
         @arg disturbance é o vetor de perturbação atual, considerado
         constante durante todo o horizonte de predição.
         @return out_next_command é o vetor de comando atualizado.
         @return QP contexto da programação quadrática
         @note além de retorna esta função também registra em obj.u
               o vetor de controle calculado.
        %}
        [out_next_command, QP] = write_next_command( obj, arg_TRACK, arg_disturbance )
    end
    
end
