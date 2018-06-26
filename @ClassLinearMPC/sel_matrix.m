%{
---------------------------------------------------------------------------
  Criador: Zo� Magalh�es (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matr�cula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: Andr� Murilo
  
  Este script implementa o m�todo do descritor de controle preditivo
  que fornece a matriz de sele��o de uma amostra no horizonte
  de predi��o.

  Hist�rico de modifica��es:
---------------------------------------------------------------------------
- 24/05/2018 - Zo� Magalh�es
- In�cio do controle de vers�o.
---------------------------------------------------------------------------
%}

%{
 @brief Retorna a matriz de sele��o PI, tal que se tilde_s for
 o vetor coluna resultante da concatena��o das amostras do vetor 
 coluna s no horizonte de predi��o N: 
 tilde_s = [s[n+1]; s[n+2]; ... ; s[n + N] ] , 
 ent�o: PI * tilde_s = s[n+arg_sel].
 @arg arg_obj � o objeto que invoca o m�todo
 @arg arg_sel indica a amostra selecionada
 @arg arg_n � o comprimento das amostras
 @arg arg_N � o n�mero de amostras concatenadas
 @return out_sel_matrix � a matriz de sele��o PI.
%}
function out_sel_matrix = sel_matrix( obj, arg_sel, arg_n, arg_N )
    out_sel_matrix = zeros( arg_n ,arg_n*arg_N);
    out_sel_matrix( :, (arg_sel - 1)*arg_n + 1:arg_sel*arg_n) = eye(arg_n);
end