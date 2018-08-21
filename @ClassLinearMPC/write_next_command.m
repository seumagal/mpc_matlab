%{
---------------------------------------------------------------------------
  Criador: Zo� Magalh�es (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matr�cula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: Andr� Murilo
  
  Este script implementa o m�todo do descritor de controle preditivo
  que atualiza o vetor de comando

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
function [out_next_command, out_QP] = write_next_command( obj, arg_TRACK, arg_disturbance )
    persistent QP;
    %Converte a trajet�ria para um vetor da concatenacao das saidas na
    %trajetoria
    TRACK = reshape(arg_TRACK,[],1);
    
    if obj.constrained == 0 % Sem restri��o
        out_next_command = -obj.K*obj.x + obj.G*TRACK + obj.L*obj.u_d;
    else %Com restri��o
            
        B_INEQ = obj.G1*obj.x + obj.G2*obj.u + obj.G3;
        FINEQ = obj.F1*obj.x + obj.F2*TRACK + obj.F3*[obj.u_d; arg_disturbance];
        cmd_lb = repmat([obj.cmd_lb;arg_disturbance],obj.n,1);    
        cmd_ub = repmat([obj.cmd_ub;arg_disturbance],obj.n,1);
    
        options = qpOASES_options('reliable','initialStatusBounds',0,'maxIter',20);

%        options = optimset('display','off','MaxIter',1000,'TolCon',eps);

        if ~isempty(obj.PI_R)
            FINEQ = obj.PI_R'*FINEQ;
            B_INEQ = [B_INEQ;-cmd_lb;cmd_ub];
            
            
            [p, fval, exitflag] = quadprog((obj.H + obj.H')/2, FINEQ, obj.AINEQ, B_INEQ,...
                                                     [], [], cmd_lb(1:obj.np), cmd_ub(1:obj.np),[],options );
        

                if isempty(p) 
                    out_next_command = obj.u;
                else
                    out_next_command = obj.PI_R(1:obj.nu-obj.nv,:)*p;
                end
            

%            p = obj.QP.solver(obj.QP.config, obj.H,FINEQ,obj.AINEQ, B_INEQ, cmd_lb(1:obj.np), cmd_ub(1:obj.np));
%            out_next_command = obj.PI_R(1:obj.nu-obj.nv,:)*p;
        else
            if ~isempty(obj.PI_E)
                FINEQ = obj.PI_E'*FINEQ;
                B_INEQ = [B_INEQ;-cmd_lb;cmd_ub];
                

                %[p, fval, exitflag] = quadprog( obj.H, FINEQ, obj.AINEQ, B_INEQ,...
                %                                [], [], [ ones(4,1)*cmd_lb(1); cmd_lb(2)],[ ones(4,1)*cmd_ub(1); cmd_ub(2)], [],options );
                 
                if isempty(QP)
                    [QP,p,fval,exitflag,iter] = qpOASES_sequence( 'i', obj.H,FINEQ,obj.AINEQ,[ ones(obj.np-1,1)*cmd_lb(1); cmd_lb(2)],[ ones(obj.np-1,1)*cmd_ub(1); cmd_ub(2)],[], B_INEQ, options );
                else
                    [p,fval,exitflag,iter] = qpOASES_sequence( 'h', QP,FINEQ,[ ones(obj.np-1,1)*cmd_lb(1); cmd_lb(2)],[ ones(obj.np-1,1)*cmd_ub(1); cmd_ub(2)],[], B_INEQ, options );
                end
                out_QP = QP;
                if isempty(p) %Condi��o inserida para facilitar debug, o projeto deve estar bem feito para que a otimiza��o tenha solu��o.
                    out_next_command = obj.u;
                else
                    out_next_command = obj.PI_E(1:obj.nu-obj.nv,:)*p;
                end                                                
                
                                 
                %tic
 %                p = obj.QP.solver(obj.QP.config, obj.H,FINEQ,obj.AINEQ, B_INEQ, [ ones(4,1)*cmd_lb(1); cmd_lb(2)],[ ones(4,1)*cmd_ub(1); cmd_ub(2)]);
 %               out_next_command = obj.PI_E(1:obj.nu-obj.nv,:)*p;
                %toc
            else  
                
                [p, fval, exitflag] = quadprog( obj.H, FINEQ, obj.AINEQ, B_INEQ,...
                                                [], [], cmd_lb, cmd_ub, [],options );
        
                if isempty(p)
                    out_next_command = obj.u;
                else
                    out_next_command =  p(1:obj.nu-obj.nv,1);
                end                                                 
                                               
            
             %   p = obj.QP.solver(obj.QP.config, obj.H,FINEQ,obj.AINEQ, B_INEQ, cmd_lb, cmd_ub);
             %   out_next_command =  p(1:obj.nu-obj.nv,1);
            end
        end
    end
    
    obj.u = out_next_command;
end        
