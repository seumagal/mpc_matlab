%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Este script implementa o método do descritor de controle preditivo
  que atualiza o vetor de comando

  Histórico de modificações:
---------------------------------------------------------------------------
- 07/06/2018 - Zoé Magalhães
- Início do controle de versão.
- Contempla restrições nos estados, no comando e na variação do comando 
  perturbações mensuráveis.
---------------------------------------------------------------------------
%}
%%
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

 @note além de retorna esta função também registra em obj.u
       o vetor de controle calculado.
%}
function out_next_command = write_next_command( obj, arg_TRACK, arg_disturbance )
    

    %Converte a trajetória para um vetor da concatenacao das saidas na
    %trajetoria
    TRACK = reshape(arg_TRACK,[],1);
   
    
    if obj.constrained == 0 % Sem restrição
        out_next_command = -obj.K*obj.x + obj.G*TRACK + obj.L*obj.u_d;
    else %Com restrição
        
        B_INEQ = obj.G1*obj.x + obj.G2*obj.u + obj.G3;
        FINEQ = obj.F1*obj.x + obj.F2*TRACK + obj.F3*[obj.u_d; arg_disturbance];
        track_lb = repmat([obj.cmd_lb;arg_disturbance],obj.n,1);
        track_ub = repmat([obj.cmd_ub;arg_disturbance],obj.n,1);
        options = optimset('display','off');

        
        [utrack, fval, exitflag] = quadprog(obj.H, FINEQ, obj.AINEQ, B_INEQ,...
                                             [], [], track_lb, track_ub, [],options );
        if isempty(utrack) %Condição inserida para facilitar debug, o projeto deve estar bem feito para que a otimização tenha solução.
            out_next_command = obj.u(1:end-obj.nv);
        else
            out_next_command =  utrack(1:obj.nu-obj.nv,1);
        end
    end


    obj.u = out_next_command;
end        
