function [h00,h01]=get_slab_h00(hamiltonian,hopping_r,nbands,nrpts,kpoint,Np,a)
    h00=zeros(nbands*Np,nbands*Np);
    h01=zeros(nbands*Np,nbands*Np);
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij=ham_qlayer2qlayer(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
    for i=1:Np
        for j=1:Np
            if (abs(i-j) <= ijmax)
                h00(nbands*(i-1)+1:nbands*i,nbands*(j-1)+1:nbands*j) = Hij(:,:,(j-i)+abs(ijmin)+1);
            end
        end
    end

    for i=1:Np
        for j=Np+1:Np*2
            if ((j-i) <= ijmax)
                h01(nbands*(i-1)+1:nbands*i,nbands*(j-1-Np)+1:nbands*(j-Np))=Hij(:,:,(j-i)+abs(ijmin)+1);
            end
        end
    end

    for i=1:Np*nbands
        for j=1:Np*nbands
            if(abs(h00(i,j)-conj(h00(j,i)))>1e-4)
                fprintf('there is something wrong with h00')
            end
        end
    end
end


function Hij=ham_qlayer2qlayer(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij = zeros(nbands,nbands,ijmax-ijmin+1);
    for j=1:nrpts
        if  ((ijmin<= hopping_r(j,3)) && (hopping_r(j,3) <= ijmax))
            kdotr = dot(kpoint,hopping_r(j,1:2)*a);
            ratio = exp(1j*kdotr);
            Hij(:,:,hopping_r(j,3)+abs(ijmin)+1) = Hij(:,:,hopping_r(j,3)+abs(ijmin)+1)+hamiltonian(:,:,j)*ratio;
        end
    end
end

