function [GLL,GRR,GB]=get_surfgreen(hamiltonian,hopping_r,nbands,nrpts,kpoint,Np,a,omega,eta)
    accuracy=1e-16;
    [h00,h01]=MTB.ham.get_slab_h00(hamiltonian,hopping_r,nbands,nrpts,kpoint,Np,a);
    epsiloni=h00;
    epsilons=h00;
    epsilons_t=h00;
    alphai=h01;
    betai=h01';
    omegac=omega+1j*eta;
    G=zeros(nbands*Np,nbands*Np,3);
    for iter=1:1000
        g0 = inv(omegac*eye(nbands*Np,nbands*Np)-epsiloni);

        % a_{i-1}*(w-e_{i-1})^{-1}
        mat1=alphai*g0;

        % b_{i-1}*(w-e_{i-1})^{-1}
        mat2=betai*g0;

        g0=mat1*betai;
        % a_{i-1}*(w-e_{i-1})^{-1}*b_{i-1}
        epsiloni = epsiloni + g0;
        
        %es_i=es_{i-1}+a_{i-1}*(w-e_{i-1})^{-1}*b_{i-1}
        epsilons = epsilons + g0;
        
        g0=mat2*alphai;
        %b_{i-1}*(w-e_{i-1})^{-1}*a_{i-1}
        epsiloni=epsiloni + g0;

        %es_i=es_{i-1}+b_{i-1}*(w-e_{i-1})^{-1}*a_{i-1}
        epsilons_t=epsilons_t + g0;

        %a_i=a_{i-1}(w-e_{i-1})*a_{i-1}
        alphai=mat1*alphai;

        %b_i=b_{i-1}(w-e_{i-1})*b_{i-1}
        betai=mat2*betai;

        real_temp=sum(abs(alphai),'all');

        if (real_temp<accuracy)
            break;
        end
    end

    GLL=inv(omegac*eye(nbands*Np,nbands*Np) - epsilons);
    GRR=inv(omegac*eye(nbands*Np,nbands*Np) - epsilons_t);
    GB=inv(omegac*eye(nbands*Np,nbands*Np) - epsiloni);    
end

