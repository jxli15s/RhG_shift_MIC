function [Unk,Enk]=get_bulk_plane_bands(obj,Kx,Ky,Kz)
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
    unktem=zeros(dim_H,dim_H,itotal);
    enktem=zeros(itotal,dim_H);
    parfor ll=1:itotal
          [i,j]=ind2sub([knum,knum],ll); 
          k=[Kx(i,j),Ky(i,j),Kz(i,j)];
          [Etem, psi] = slovek(obj,k);
          enktem(ll,:)=Etem;
          unktem(:,:,ll)=psi;
    end
    Unk=reshape(unktem,[dim_H,dim_H,knum,knum]);
    Enk=reshape(enktem,[knum,knum,dim_H]);
end

function [E, Psik]=slovek(obj,k)
    % hk=get_hk(obj,k);
    hk=obj.get_hk(k);
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

% function hk=get_hk(obj,kpoint)
% hamiltonian=reshape(obj.ham,[],size(obj.hopr,1));
% hk=hamiltonian*exp(1j*((obj.hopr*obj.a)*kpoint'));
% hk=reshape(hk,size(obj.ham,1),[]);
% hk=(hk+hk')/2.0;
% end
