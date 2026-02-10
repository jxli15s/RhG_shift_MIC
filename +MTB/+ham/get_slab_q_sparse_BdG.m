function Energy=get_slab_q_sparse_BdG(hamiltonian,hopping_r,nslab,nbands,numEigs,nrpts,kpoint,a,b,mu,delta)
        % kpoint=kpoint*b;
        hk=get_slab_hk_BdG(hamiltonian,hopping_r,nslab,nbands,nrpts,kpoint,a,mu,delta);
        vals=sort(real(eigs(hk,numEigs,'smallestabs')));
        Energy=vals;
end

function hk=get_slab_hk_BdG(hamiltonian,hopping_r,nslab,nbands,nrpts,kpoint,a,mu,delta)
    hk_e=sparse(nbands*nslab,nbands*nslab);
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij_e=ham_qlayer2qlayer(hamiltonian,hopping_r,nrpts,kpoint,a);
    for i1=1:nslab
        for i2=1:nslab
            if (ijmin<=(i2-i1)) && ((i2-i1) <= ijmax)
                hk_e((i2-1)*nbands+1:(i2-1)*nbands+nbands, ...
                    (i1-1)*nbands+1:(i1-1)*nbands+nbands) ...
                = Hij_e{i1-i2+abs(ijmin)+1};
            end
        end
    end
    h_onsite=sparse(eye(nbands*nslab,nbands*nslab)*mu);
    hk_e=hk_e-h_onsite;

    hk_h=sparse(nbands*nslab,nbands*nslab);
    kpoint=-kpoint;
    Hij_h=ham_qlayer2qlayer(hamiltonian,hopping_r,nrpts,kpoint,a);
    for i1=1:nslab
        for i2=1:nslab
            if (ijmin<=(i2-i1)) && ((i2-i1) <= ijmax)
                hk_h((i2-1)*nbands+1:(i2-1)*nbands+nbands, ...
                    (i1-1)*nbands+1:(i1-1)*nbands+nbands) ...
                = Hij_h{i1-i2+abs(ijmin)+1};
            end
        end
    end
    hk_h=-hk_h.';
    hk_h=hk_h+h_onsite;
    hk_delta=get_hk_delta(nbands,nslab,delta);
    hk=[hk_e,sparse(zeros(size(hk_e)));sparse(zeros(size(hk_h))),hk_h];
    hk=hk+hk_delta;
    hk=(hk+hk')/2.0;
end


function hk_delta=get_hk_delta(nbands,nslab,delta)
    hk_delta=zeros(2*nbands*nslab,2*nbands*nslab);
    sigma_x=[0, 1; 1, 0];
    sigma_y=[0, -1j; 1j, 0];
    sigma_z=[1, 0; 0, -1];
    h_delta=1j*sigma_y*delta;
    h_sigma_z=kron(sigma_z,eye(nbands*nslab/4));
    hk_delta_up=kron(h_sigma_z,h_delta);
    hk_delta=kron(1j*sigma_y,hk_delta_up);
    hk_delta=sparse(hk_delta);
end

function Hij=ham_qlayer2qlayer(hamiltonian,hopping_r,nrpts,kpoint,a)
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij = repmat({sparse(size(hamiltonian,1),size(hamiltonian,2))},1,ijmax-ijmin+1);
    for j=1:nrpts
        if  ((ijmin<= hopping_r(j,3)) && (hopping_r(j,3)<= ijmax))
            kdotr = dot(kpoint,hopping_r(j,1:2)*a);
            ratio = exp(1j*kdotr);
            Hij{hopping_r(j,3)+abs(ijmin)+1} = Hij{hopping_r(j,3)+abs(ijmin)+1} + sparse(hamiltonian(:,:,j)*ratio);
        end
    end
end

