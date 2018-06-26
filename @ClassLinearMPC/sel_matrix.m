%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Este script implementa o método do descritor de controle preditivo
  que fornece a matriz de seleção de uma amostra no horizonte
  de predição.

  Histórico de modificações:
---------------------------------------------------------------------------
- 24/05/2018 - Zoé Magalhães
- Início do controle de versão.
---------------------------------------------------------------------------
%}

%{
 @brief Retorna a matriz de seleção PI, tal que se tilde_s for
 o vetor coluna resultante da concatenação das amostras do vetor 
 coluna s no horizonte de predição N: 
 tilde_s = [s[n+1]; s[n+2]; ... ; s[n + N] ] , 
 então: PI * tilde_s = s[n+arg_sel].
 @arg arg_obj é o objeto que invoca o método
 @arg arg_sel indica a amostra selecionada
 @arg arg_n é o comprimento das amostras
 @arg arg_N é o número de amostras concatenadas
 @return out_sel_matrix é a matriz de seleção PI.
%}
function out_sel_matrix = sel_matrix( obj, arg_sel, arg_n, arg_N )
    out_sel_matrix = zeros( arg_n ,arg_n*arg_N);
    out_sel_matrix( :, (arg_sel - 1)*arg_n + 1:arg_sel*arg_n) = eye(arg_n);
end