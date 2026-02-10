function Energy=get_slab_q_sparse_BdG_del(hamiltonian,hopping_r,del_index,nslab,nbands,numEigs,nrpts,kpoint,a,b,mu,delta)
        % kpoint=kpoint*b;
        hk=get_slab_hk_BdG(hamiltonian,hopping_r,del_index,nslab,nbands,nrpts,kpoint,a,mu,delta);
        vals=sort(real(eigs(hk,numEigs,'smallestabs')));
        Energy=vals;
end

function hk=get_slab_hk_BdG(hamiltonian,hopping_r,del_index,nslab,nbands,nrpts,kpoint,a,mu,delta)
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
    hk_e(del_index,:)=[];
    hk_e(:,del_index)=[];


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
    hk_h(del_index,:)=[];
    hk_h(:,del_index)=[];
    
    hk_delta=get_hk_delta_sparse(nbands,nslab,del_index,delta);
    n1=size(hk_e,1);
    n2=size(hk_h,2);
    hk=[hk_e,sparse(n1,n2);sparse(n2,n1),hk_h];
    hk=hk+hk_delta;
    hk=(hk+hk')/2.0;
end


function hk_delta = get_hk_delta_sparse(nbands, nslab, del_index, delta)
% 构造稀疏 BdG pairing 矩阵 hk_delta
% 输入:
%   nbands     - 每层的轨道数
%   nslab      - 层数
%   del_index  - 被删掉的态索引（影响最终大小）
%   delta      - pairing 振幅
% 输出:
%   hk_delta   - 2*(nbands*nslab - del_num) × 2*(nbands*nslab - del_num) 稀疏 pairing 矩阵

    del_num = size(del_index, 1);
    N_total = nbands * nslab;
    N_eff = N_total - del_num;

    % Pauli matrices
    sigma_y = sparse([0, -1i; 1i, 0]);

    % Pairing on sublattice
    h_delta = 1i * sigma_y * delta;  % 2×2 matrix

    % 构造上下三角块（相对粒子-空穴子空间）
    up_size = N_eff;
    hk_delta_up = spalloc(up_size, up_size, up_size);  % 预留空间

    % block_size = size(h_delta, 1);
    n_block = (nbands * (nslab - 1)) / 4;

    % % 安全检查是否整除
    % if mod(n_block,1) ~= 0
    %     error("nbands*(nslab-1)/4 must be an integer.");
    % end

    % 构造两个 pairing 块
    hk_delta_0 = kron(speye(n_block), h_delta);
    hk_delta_pi = -kron(speye(n_block), h_delta);

    % 插入到 hk_delta_up
    hk_delta_up(1:size(hk_delta_0,1), 1:size(hk_delta_0,2)) = hk_delta_0;
    hk_delta_up(end-size(hk_delta_pi,1)+1:end, end-size(hk_delta_pi,2)+1:end) = hk_delta_pi;

    % 构造最终 BdG pairing 稀疏矩阵（粒子-空穴空间）
    hk_delta = kron(1i * sigma_y, hk_delta_up);  % 2N_eff × 2N_eff 稀疏矩阵
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

