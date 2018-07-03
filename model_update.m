%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  

  Função que atualiza o estado com base no estado atual, sinal
  de comando e perturbação

%}

%{
    @brief out_x = atg_A*arg_x + arg_B*[arg_u;arg_e]
%}
function [out_x] = model_update( arg_A, arg_B, arg_x, arg_u, arg_e )
    out_x = arg_A*arg_x + arg_B*[arg_u;arg_e];
end
