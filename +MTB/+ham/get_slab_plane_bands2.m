function [Unk,Enk]=get_slab_plane_bands2(obj,del_index,Kx,Ky,nslab)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Kx,Ky: the k-mesh we will solve
    % dim_H: num of the bands/orbits
    % Unk: Norbits*Nands*knum*knum
    % Enk: knum*knum*Nbands
    knum=size(Kx,1);
    % Unk=zeros(dim_H,dim_H,knum,knum);
    % Enk=zeros(knum,knum,dim_H);
    itotal=knum^2;
    dim_H=size(obj.ham,1);
    nrpts=size(obj.hopr2,1);
    unktem=zeros(dim_H*nslab-size(del_index,1),dim_H*nslab-size(del_index,1),itotal);
    size(unktem)
    enktem=zeros(itotal,dim_H*nslab-size(del_index,1));
    size(enktem)
    parfor ll=1:itotal
          [i,j]=ind2sub([knum,knum],ll); 
          k=[Kx(i,j),Ky(i,j)];
          fprintf("Processing on %d %d of %d KPOINTS\n", i, j, itotal);
          [Etem, psi] = slovek(obj,del_index,k,nslab,dim_H,nrpts);
          enktem(ll,:)=Etem;
          unktem(:,:,ll)=psi;
    end
    Unk=reshape(unktem,[dim_H*nslab-size(del_index,1),dim_H*nslab-size(del_index,1),knum,knum]);
    Enk=reshape(enktem,[knum,knum,dim_H*nslab-size(del_index,1)]);
    % parfor i=1:knum
    %     for j=1:knum
    %         k=[Kx(i,j),Ky(i,j)];
    %         % [Etem, psi] = obj.solve_kpoint(k);
    %         [Etem, psi] = slovek(obj,k);
    %         Unk(:,:,i,j)=psi;
    %         Enk(i,j,:)=Etem;
    %     end
    % end 
end

function [E, Psik]=slovek(obj,del_index,k,nslab,dim_H,nrpts)
    hk=MTB.ham.get_slab_hk_v2(obj.ham,obj.hopr2,del_index,nslab,dim_H,nrpts,k,obj.a2);
    hk=(hk+hk')./2.0;
    [V,D]=eig(hk);
    [E,ind]=sort(diag(D));
    Psik=V(:,ind);
end

% function hk=get_hk(obj,kpoint)
% %    kpoint=kpoint*obj.b
%     hk=zeros(size(obj.ham,1),size(obj.ham,1));
%     nrpts=size(obj.ham,3);
%     for j=1:nrpts
%         temp=obj.ham(:,:,j)*...
%         exp(1j*dot(kpoint,(obj.hopr(j,:)*obj.a))); 
%         hk=hk+temp;
%     end
% end
