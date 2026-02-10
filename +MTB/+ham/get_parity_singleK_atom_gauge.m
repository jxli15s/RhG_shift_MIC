function [Energy,Psik,hk]=get_parity_singleK_atom_gauge(hamiltonian,hopping_r,wpos,nbands,nrpts,kpoint,a,b)

        tau=zeros(nbands,nbands,3);
        for i = 1:nbands
            for j = 1:nbands
                tau(i,j,:)=wpos(i,:)-wpos(j,:);
            end
        end

        kpoint=kpoint*b;
        hk=get_hk(hamiltonian,hopping_r,tau,nbands,nrpts,kpoint,a);
        % ishermitian(hk)
        [V,D]=eig(hk);
        [Energy,ind]=sort(diag(D));
        Psik=V(:,ind);
end
        
% 
% function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
%     hk=zeros(nbands,nbands);
%     for j=1:nrpts
%         temp=hamiltonian(:,:,j)*...
%         exp(1j*dot(kpoint,(hopping_r(j,:)*a))); 
%         hk=hk+temp;
%     end
% end

function hk=get_hk(hamiltonian,hopping_r,tau,nbands,nrpts,kpoint,a)
    tphase=zeros(nbands,nbands);
    for i = 1:nbands
        for j = 1:nbands
            r=reshape(tau(i,j,:),1,3);
            tphase(i,j)=exp(1j*r*kpoint');
        end
    end
    hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    hk=hamiltonian*exp(1j*((hopping_r*a)*kpoint'));
    hk=reshape(hk,nbands,nbands);
    hk=hk.*tphase;
    hk=(hk+hk')/2.0;
end

            
            
            
            