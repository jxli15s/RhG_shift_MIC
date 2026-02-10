function [Energy,Psik,hk]=get_parity_singleK(hamiltonian,hopping_r,nbands,nrpts,kpoint,a,b)

        kpoint=kpoint*b;
        hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
        
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

function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    hk=hamiltonian*exp(1j*((hopping_r*a)*kpoint'));
    hk=reshape(hk,nbands,nbands);
    hk=(hk+hk')/2.0;
end

            
            
            
            