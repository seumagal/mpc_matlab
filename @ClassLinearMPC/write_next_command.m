%{
---------------------------------------------------------------------------
  Criador: Zo� Magalh�es (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matr�cula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: Andr� Murilo
  
  Este script implementa o m�todo do descritor de controle preditivo
  que atualiza o vetor de comando

  Hist�rico de modifica��es:
---------------------------------------------------------------------------
- 07/06/2018 - Zo� Magalh�es
- In�cio do controle de vers�o.
- Contempla restri��es nos estados, no comando e na varia��o do comando 
  perturba��es mensur�veis.
---------------------------------------------------------------------------
%}
%%
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

 @note al�m de retorna esta fun��o tamb�m registra em obj.u
       o vetor de controle calculado.
%}
function out_next_command = write_next_command( obj, arg_TRACK, arg_disturbance )
    

    %Converte a trajet�ria para um vetor da concatenacao das saidas na
    %trajetoria
    TRACK = reshape(arg_TRACK,[],1);
   
    
    if obj.constrained == 0 % Sem restri��o
        out_next_command = -obj.K*obj.x + obj.G*TRACK + obj.L*obj.u_d;
    else %Com restri��o
        
        B_INEQ = obj.G1*obj.x + obj.G2*obj.u + obj.G3;
        FINEQ = obj.F1*obj.x + obj.F2*TRACK + obj.F3*[obj.u_d; arg_disturbance];
        track_lb = repmat([obj.cmd_lb;arg_disturbance],obj.n,1);
        track_ub = repmat([obj.cmd_ub;arg_disturbance],obj.n,1);
        options = optimset('display','off');

        
        [utrack, fval, exitflag] = quadprog(obj.H, FINEQ, obj.AINEQ, B_INEQ,...
                                             [], [], track_lb, track_ub, [],options );
        if isempty(utrack) %Condi��o inserida para facilitar debug, o projeto deve estar bem feito para que a otimiza��o tenha solu��o.
            out_next_command = obj.u(1:end-obj.nv);
        else
            out_next_command =  utrack(1:obj.nu-obj.nv,1);
        end
    end


    obj.u = out_next_command;
end        
