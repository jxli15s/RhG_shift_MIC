function Energy=get_bulk_bands_BdG_sparse_at_q(hamiltonian,hopping_r,nbands,nrpts,kpoint,a,b,mu,h_delta,numEigs)
    %%solve the Hamiltonian(k)
    % hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    h_onsite=speye(nbands)*mu;
    kpoint=kpoint*b;
    hk_e=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
    hk_e=hk_e-h_onsite;
    hk_h=get_hk(hamiltonian,hopping_r,nbands,nrpts,-kpoint,a);
    hk_h=-hk_h.';
    hk_h=hk_h+h_onsite;
    hk=[hk_e,h_delta;h_delta', hk_h];
    % ishermitian(hk)
    hk=(hk+hk')/2;
    vals=sort(eigs(hk,numEigs,'smallestabs'));
    Energy=vals;
end

function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
hk = sparse(nbands, nbands);         % 初始化为 norb x norb 的全零 sparse 矩阵
for j=1:nrpts
    temp=hamiltonian{j}*...
        exp(1j*dot(kpoint,(hopping_r(j,:)*a)));
    hk=hk+temp;
end
hk=(hk+hk')/2.0;
end
        
            
            
            
            