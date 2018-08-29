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
 @arg arg_disturbance � o vetor de perturba��o atual, considerado
 constante durante todo o horizonte de predi��o.
 @arg arg_state_slew_lb � a taxa m�nima de varia��o dos estados restritos
 @arg arg_state_slew_ub � a taxa m�zima de varia��o dos estados restritos
 @return out_next_command � o vetor de comando atualizado.

 @note al�m de retorna esta fun��o tamb�m registra em obj.u
       o vetor de controle calculado.
%}
function [out_next_command, out_QP] = write_next_command( obj, arg_TRACK, arg_disturbance, arg_state_slew_lb, arg_state_slew_ub )
    persistent QP;
    %Converte a trajet�ria para um vetor da concatenacao das saidas na
    %trajetoria
    TRACK = reshape(arg_TRACK,[],1);
    
    if obj.constrained == 0 % Sem restri��o
        x = [obj.x; arg_disturbance ];
        out_next_command = -obj.K*x + obj.G*TRACK + obj.L*obj.u_d;
        out_QP = 0;
    else %Com restri��o
        G5 = zeros(size(obj.G4,1),1);
        rows_slew = 2*obj.n*obj.nc + 2*obj.n*obj.nu + 1 : 3*obj.n*obj.nc + 2*obj.n*obj.nu;
        
        G5( rows_slew, : ) = repmat(arg_state_slew_ub,obj.n,1);
        G5( rows_slew + obj.n*obj.nc, : ) = -repmat(arg_state_slew_lb,obj.n,1);
        
        B_INEQ = ( obj.G1 + obj.G4) *obj.x + obj.G2*obj.u + obj.G3 + G5;
        F = obj.F1*obj.x + obj.F2*TRACK + obj.F3*[obj.u_d; arg_disturbance];
        cmd_lb = repmat([obj.cmd_lb;arg_disturbance],obj.n,1);    
        cmd_ub = repmat([obj.cmd_ub;arg_disturbance],obj.n,1);
    
        options = qpOASES_options('maxIter',50);

%        options = optimset('display','off','MaxIter',1000,'TolCon',eps);

        if ~isempty(obj.PI_R)
            F = obj.PI_R'*F;
            B_INEQ = [B_INEQ;-cmd_lb;cmd_ub];
            
            
            [p, fval, exitflag] = quadprog((obj.H + obj.H')/2, F, obj.AINEQ, B_INEQ,...
                                                     [], [], cmd_lb(1:obj.np), cmd_ub(1:obj.np),[],options );
        

                if isempty(p) 
                    out_next_command = obj.u;
                else
                    out_next_command = obj.PI_R(1:obj.nu-obj.nv,:)*p;
                end
            

%            p = obj.QP.solver(obj.QP.config, obj.H,F,obj.AINEQ, B_INEQ, cmd_lb(1:obj.np), cmd_ub(1:obj.np));
%            out_next_command = obj.PI_R(1:obj.nu-obj.nv,:)*p;
        else
            if ~isempty(obj.PI_E)
                F = obj.PI_E'*F;
                B_INEQ = [B_INEQ(1:606,:);-cmd_lb;cmd_ub];
                

                %[p, fval, exitflag] = quadprog( obj.H, F, obj.AINEQ, B_INEQ,...
                %                                [], [], [ ones(4,1)*cmd_lb(1); cmd_lb(2)],[ ones(4,1)*cmd_ub(1); cmd_ub(2)], [],options );
                 
                if isempty(QP)
                    options = qpOASES_options('maxIter',1000);
                    [QP,p,fval,exitflag,iter] = qpOASES_sequence( 'i', obj.H,F,obj.AINEQ([1:606,809:end],:),[ ones(obj.np-1,1)*cmd_lb(1); cmd_lb(2)],[ ones(obj.np-1,1)*cmd_ub(1); cmd_ub(2)],[], B_INEQ, options );
                else
                    [p,fval,exitflag,iter] = qpOASES_sequence( 'h', QP,F,[ ones(obj.np-1,1)*cmd_lb(1); cmd_lb(2)],[ ones(obj.np-1,1)*cmd_ub(1); cmd_ub(2)],[], B_INEQ, options );
                end

                out_QP = QP;
                if isempty(p) %Condi��o inserida para facilitar debug, o projeto deve estar bem feito para que a otimiza��o tenha solu��o.
                    out_next_command = obj.u;
                else
                    out_next_command = obj.PI_E(1:obj.nu-obj.nv,:)*p;
                end                                                
                
                                 
                %tic
 %                p = obj.QP.solver(obj.QP.config, obj.H,F,obj.AINEQ, B_INEQ, [ ones(4,1)*cmd_lb(1); cmd_lb(2)],[ ones(4,1)*cmd_ub(1); cmd_ub(2)]);
 %               out_next_command = obj.PI_E(1:obj.nu-obj.nv,:)*p;
                %toc
            else  
                
                [p, fval, exitflag] = quadprog( obj.H, F, obj.AINEQ, B_INEQ,...
                                                [], [], cmd_lb, cmd_ub, [],options );
        
                if isempty(p)
                    out_next_command = obj.u;
                else
                    out_next_command =  p(1:obj.nu-obj.nv,1);
                end                                                 
                                                           
             %   p = obj.QP.solver(obj.QP.config, obj.H,F,obj.AINEQ, B_INEQ, cmd_lb, cmd_ub);
             %   out_next_command =  p(1:obj.nu-obj.nv,1);
            end
        end
    end
    
    obj.u = out_next_command;
end        
