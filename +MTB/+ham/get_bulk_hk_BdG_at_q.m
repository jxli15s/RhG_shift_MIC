function [Energy,Psik,hk]=get_bulk_hk_BdG_at_q(hamiltonian,hopping_r,nbands,nrpts,kpoint,a,b,mu,delta)
    kpoint=kpoint*b;
    %%solve the Hamiltonian(k)
    % hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    h_onsite=eye(nbands,nbands)*mu;
    hk_e=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
    hk_e=hk_e-h_onsite;
    hk_h=get_hk(hamiltonian,hopping_r,nbands,nrpts,-kpoint,a);
    hk_h=-hk_h.';
    hk_h=hk_h+h_onsite;
    hk_delta=get_hk_delta(nbands,delta);
    hk=[hk_e,zeros(size(hk_e));zeros(size(hk_h)),hk_h];
    hk=hk+hk_delta;
    % ishermitian(hk)
    % hk=(hk+hk')/2.0;
    [V,D]=eig(hk);
    [Energy,ind]=sort(diag(D));
    Psik=V(:,ind);
end

function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    hk=hamiltonian*exp(1j*((hopping_r*a)*kpoint'));
    hk=reshape(hk,nbands,nbands);
    hk=(hk+hk')/2.0;
end

function hk_delta=get_hk_delta(nbands,delta)
    hk_delta=zeros(2*nbands,2*nbands);
    sigma_x=[0, 1; 1, 0];
    sigma_y=[0, -1j; 1j, 0];
    sigma_z=[1, 0; 0, -1];
    h_delta=1j*sigma_y*delta;
    h_sigma_z=eye(nbands/2);
    hk_delta_up=kron(h_sigma_z,h_delta);
    hk_delta=kron(1j*sigma_y,hk_delta_up);
end

% function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
%     phase=exp(1j*((hopping_r*a)*kpoint'));
%     phase=kron(phase,ones(nbands^2,1));
%     phase=reshape(phase,nbands,nbands,nrpts);
%     hk=hamiltonian.*phase;
%     hk=sum(hk,3);
% end

function hk=get_hk_new(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    phase=reshape(kron(exp(1j*((hopping_r*a)*kpoint')),ones(nbands^2,1)),nbands,nbands,nrpts);
    hk=sum(hamiltonian.*phase,3);
end
        
            
            
            
            