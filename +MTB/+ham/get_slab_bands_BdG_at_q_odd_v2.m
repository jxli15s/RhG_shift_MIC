function [Energy,Psik]=get_slab_bands_BdG_at_q_odd_v2(hamiltonian,hopping_r,del_index,nslab,nbands,nrpts,kpoint,a,b,mu,delta)
        % This is different form get_slab_bands_BdG_at_q and get_slab_bands_BdG_at_q_odd.
        % This is a general version to deleate orbitals in BdG Hamiltonian
        % kpoint=kpoint*b;
        % hk_delta=get_hk_delta(nbands,nslab,del_index,delta);
        hk=get_slab_hk_BdG(hamiltonian,hopping_r,del_index,nslab,nbands,nrpts,kpoint,a,mu,delta);
        [V,D]=eig(hk);
        [Energy,ind]=sort(diag(D));
        Psik=V(:,ind);
end

function hk=get_slab_hk_BdG(hamiltonian,hopping_r,del_index,nslab,nbands,nrpts,kpoint,a,mu,delta)
    hk_e=zeros(nbands*nslab,nbands*nslab);
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij_e=ham_qlayer2qlayer(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
    for i1=1:nslab
        for i2=1:nslab
            if (ijmin<=(i2-i1)) && ((i2-i1) <= ijmax)
                hk_e((i2-1)*nbands+1:(i2-1)*nbands+nbands, ...
                    (i1-1)*nbands+1:(i1-1)*nbands+nbands) ...
                = Hij_e(:,:,i1-i2+abs(ijmin)+1);
            end
        end
    end
    h_onsite=eye(nbands*nslab-size(del_index,1),nbands*nslab-size(del_index,1))*mu;
    hk_e(del_index,:)=[];
    hk_e(:,del_index)=[];
    hk_e=hk_e-h_onsite;
    % % %% add zeeman
    % % sigma_x=[0, 1; 1, 0];
    % % zeeman=kron(eye(size(hk_e)/2),sigma_x*delta/2);
    % % hk_e=hk_e+zeeman;
    


    hk_h=zeros(nbands*nslab,nbands*nslab);
    kpoint=-kpoint;
    Hij_h=ham_qlayer2qlayer(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
    for i1=1:nslab
        for i2=1:nslab
            if (ijmin<=(i2-i1)) && ((i2-i1) <= ijmax)
                hk_h((i2-1)*nbands+1:(i2-1)*nbands+nbands, ...
                    (i1-1)*nbands+1:(i1-1)*nbands+nbands) ...
                = Hij_h(:,:,i1-i2+abs(ijmin)+1);
            end
        end
    end


    hk_h=-hk_h.';
    hk_h(del_index,:)=[];
    hk_h(:,del_index)=[];
    hk_h=hk_h+h_onsite;

    % % %add zeeman
    % % hk_h=hk_h-zeeman;
   
    hk_delta=get_hk_delta(nbands,nslab,del_index,delta);
    hk=[hk_e,zeros(size(hk_e));zeros(size(hk_h)), hk_h];
    % size(hk)
    hk=hk+hk_delta;
    % hk=(hk+hk')/2.0;
end


function hk_delta=get_hk_delta(nbands,nslab,del_index,delta)
    hk_delta=zeros(2*nbands*nslab-2*size(del_index,1),2*nbands*nslab-2*size(del_index,1));
    sigma_x=[0, 1; 1, 0];
    sigma_y=[0, -1j; 1j, 0];
    sigma_z=[1, 0; 0, -1];
    h_delta=1j*sigma_y*delta;
    % hk_delta_up=eye(size(hk_delta)/4);
    % hk_delta_up=kron(hk_delta_up,h_delta);
    hk_delta_up=zeros(size(hk_delta)/2);
    % size(hk_delta_up)
    hk_delta_0 = kron(eye(nbands*(nslab-1)/4),h_delta);
    % size(hk_delta_0)
    hk_delta_pi = -kron(eye(nbands*(nslab-1)/4),h_delta);
    hk_delta_up(1:size(hk_delta_0,1),1:size(hk_delta_0,1))=hk_delta_0;
    hk_delta_up(end-size(hk_delta_pi,1)+1:end,end-size(hk_delta_pi,1)+1:end)=hk_delta_pi;
    hk_delta=kron(1j*sigma_y,hk_delta_up);
end

function Hij=ham_qlayer2qlayer(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij = zeros(nbands,nbands,ijmax-ijmin+1);
    for j=1:nrpts
        if  (ijmin<= hopping_r(j,3) <= ijmax)
            kdotr = dot(kpoint,hopping_r(j,1:2)*a);
            ratio = exp(1j*kdotr);
            Hij(:,:,hopping_r(j,3)+abs(ijmin)+1) = Hij(:,:,hopping_r(j,3)+abs(ijmin)+1)+hamiltonian(:,:,j)*ratio;
        end
    end
end

