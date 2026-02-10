function Energy=get_slab_mu_E_sparse_BdG_del(hamiltonian,hopping_r,del_index,nslab,nbands,numEigs,nrpts,kpoint,a,b,mus,delta)
    %%solve the Hamiltonian(k)
    Energy=zeros(numEigs,length(mus));

    %hbar = parfor_progressbar_v1(length(kpoints),'Computing...');
    parfor i=1:length(mus)
        fprintf('Current running k_i = %2d \n',i)
        mu=mus(i)
        Energy(:,i)=MTB.ham.get_slab_q_sparse_BdG_del(hamiltonian,hopping_r,del_index,nslab,nbands,numEigs,nrpts,kpoint,a,b,mu,delta)
    end
   % close(hbar);
end