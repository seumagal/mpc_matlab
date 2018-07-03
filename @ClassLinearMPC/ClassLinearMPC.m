%{
---------------------------------------------------------------------------
  Criador: Zo� Magalh�es (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matr�cula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: Andr� Murilo
  
  Este script define a classe que representa um descritor do alogritmo
  de controle preditivo.

---------------------------------------------------------------------------
%}
%%
classdef ClassLinearMPC < handle
    properties
        % Matriz da representa��o em espa�o de estados discreto
        % x[n] = Ax[n-1] + Bu[n]
        % y[n] = Cx[n] + Du[n]
        A       
        B
        C       
        D
        % Pondera��es da fun��o custo
        Q_U   % Pondera��o do comando    
        Q_Y   % Pondera��o do erro nas sa�das controladas
        % Matrizes da fun��o custo
        % J = (1/2)tilde_u' H tilde_u + (F1x[n] + F_2 tilde_y_D + F_3 u_d)'tilde_u
        H
        F1
        F2
        F3
        
        % C�lulas de matrizes auxiliares usadas no c�lculo das matrizes
        % da fun��o custo e da inequa��o que representa as restri��es
        c_PHI
        c_PSI
        
        % Matrizes da lei de controle obtida da minimiza��o da fun��o
        % custo sem restri��es
        % u = -K x[n] + G tilde_y_d + L u_d
        K
        G
        L

        x     % �ltimo estado medido 
        u     % �ltimo vetor de comando calculado
        u_d   % Vetor de comado desejado para ss
        n     % Horizonte de predi��o
        nu    % Comprimento do vetor de controle
        nx    % Comprimento do vetor de estados
        ny    % Comprimento do vetor de sa�da
        nc    % Comrpimento do vetor de sa�das com restri��o
        nv    % Comprimento do vetor de perturba��es mensur�veis (entradas n�o controladas)
        constrained %flag que indica se o modelo tem restri��es(=1)
        np    % Comprimento do comando par�metrizado.
        
        %Matrizes da inequa��o que representa as restri��es
        AINEQ 
        G1
        G2
        G3
        BINEQ
        C_C
        
        % Restri��es no comando
        cmd_lb
        cmd_ub
        
        % Parametriza��o
        PI_R
        PI_E
        
        %
        QP
        
    end
    
    methods(Static)        
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
         @arg arg_C_C � a matriz Cc da equa��o da sa�da com restri��o y_c = Cc x + Dc u
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

        obj = settings( arg_A, arg_B, arg_C, arg_D, arg_Q_U, ...
                        arg_Q_Y, arg_C_C, arg_n, arg_lb, arg_ub, ...
                        arg_slew_lb, arg_slew_ub, ...
                        arg_cmd_lb, arg_cmd_ub, ...
                        arg_E, arg_parameterization, arg_param_coefficient,...
                        arg_sampling_time, solver )         
    end 
    
    methods
        %{
         @brief Retorna a matriz de sele��o PI, tal que se tilde_s for
         o vetor coluna resultante da concatena��o das amostras do vetor 
         coluna s no horizonte de predi��o N: 
         tilde_s = [s[n+1]; s[n+2]; ... ; s[n + N] ] , 
         ent�o: PI * tilde_s = s[n+arg_sel].
         @arg arg_obj � o objeto que invoca o m�todo
         @arg arg_sel indica a amostra selecionada
         @arg arg_n � o comprimento das amostras
         @arg arg_N � o numero de amostras concatenadas
         @return out_sel_matrix � a matriz de sele��o PI.
        %}
        out_sel_matrix = sel_matrix( obj, arg_sel, arg_n, arg_N )
        
        %{
         @brief 
         @arg arg_obj � o objeto que invoca o m�todo
         @arg arg_TRACK � a trajet�ria a ser rastreada composta
         pela concatena��o das amostras vetor coluna de sa�da desejada y_d no
         horizonte de predi��o.
         arg_TRACK = [y_d[n+1],y_d[n+2], ..., y_d[n+N]].
         @arg disturbance � o vetor de perturba��o atual, considerado
         constante durante todo o horizonte de predi��o.
         @return out_next_command � o vetor de comando atualizado.
         @return QP contexto da programa��o quadr�tica
         @note al�m de retorna esta fun��o tamb�m registra em obj.u
               o vetor de controle calculado.
        %}
        [out_next_command, QP] = write_next_command( obj, arg_TRACK, arg_disturbance )
    end
    
end
