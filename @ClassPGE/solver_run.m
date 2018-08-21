%{
---------------------------------------------------------------------------
  Criador: Zoé Magalhães (zr.magal@gmail.com)
  Mestrando do PPMEC-Unb, matrícula 170172767
  Disciplina: Controle Preditivo 01/2018
  Professor: André Murilo
  
  Função que executa a programação quadrática pelo algoritmo PGE.
%}


%{
    @brief Função que executa a programação quadrática pelo algoritmo
    PGE para minimizar a a função custo J = p'Hp + arg_F'p com as restrições
    Ap <= B, arg_p_min <= p <= arg_p_max;
    @param[in] arg_H é a matriz H da função custo
    @param[in] arg_F é a matriz arg_F da função custo
    @param[in] arg_A é matriz A da inequação de restrição
    @param[in] arg_B é a matriz B da inequação de restrição
    @param[in] arg_p_min é o limite inferior do espaço de busca
    @param[in] arg_p_max é o limite superior do espaço de bsuca
    
    @return out_sol é o vetor com a solução da QP
%}

function [out_sol] = solver_run(obj, arg_H,arg_F,arg_A, arg_B,arg_p_min,arg_p_max)
    hmax0 = norm(arg_H,2);
    hmaxg = norm(arg_A'*arg_A,2);
    np = size(arg_H,1);
    indrho = 1;
    lesp = zeros(np,obj.N_iter);
    lesp(:,1) = obj.p0;
    beta_plus = obj.beta_plus;
    beta_minus = obj.beta_minus;
    gam_min = obj.gam_min;
    gam = gam_min;
    obj.rho = obj.rho0;
    
    for i=1:obj.N_iter -1,
        hmax=hmax0+2*obj.rho*hmaxg;
        if(indrho==obj.n_rho),
            obj.rho=min(obj.rho_max,obj.beta_plus*obj.rho);
            indrho=1;
        end
        inter = arg_A*lesp(:,i)-arg_B;
        G = arg_H*lesp(:,i)+arg_F+2*obj.rho*arg_A'*max(0,inter);
        p1=lesp(:,i)-1/(hmax)*G;
        p2=lesp(:,i)-gam/(hmax)*G;
        inter1=arg_A*p1-arg_B;
        inter2=arg_A*p2-arg_B;
        J1=0.5*p1'*arg_H*p1+arg_F'*p1+obj.rho*(norm(max(0,inter1)))^2;
        J2=0.5*p2'*arg_H*p2+arg_F'*p2+obj.rho*(norm(max(0,inter2)))^2;
        if ( J1<J2 )
            lesp(:,i+1)=p1;
            gam=max(gam_min,beta_minus*gam);
        else
            lesp(:,i+1)=p2;
            gam=beta_plus*gam;
        end
        
        for j=1:np
            lesp(j,i+1)=...
            min(arg_p_max(j),max(arg_p_min(j),lesp(j,i+1)));
        end
        indrho=indrho+1;
        
        %disp(obj.J);
    end
    obj.rho0=obj.rho;
    obj.p0 = lesp(:,end);
    out_sol = lesp(:,end);
end
