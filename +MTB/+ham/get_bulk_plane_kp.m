function [Unk,Enk]=get_bulk_plane_kp(pars,get_hk,nbands,Kx,Ky,Kz)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Kx,Ky: the k-mesh we will solve
    % dim_H: num of the bands/orbits
    % Unk: Norbits*Nands*knum*knum
    % Enk: knum*knum*Nbands
    
    % Unk=zeros(dim_H,dim_H,knum,knum);
    % Enk=zeros(knum,knum,dim_H);
    knum=size(Kx,1);
    itotal=knum^2;
    dim_H=nbands;
    unktem=zeros(dim_H,dim_H,itotal);
    enktem=zeros(itotal,dim_H);
    parfor ll=1:itotal
          [i,j]=ind2sub([knum,knum],ll); 
          k=[Kx(i,j),Ky(i,j),Kz(i,j)];
          [Etem, psi] = slovek(k,pars,get_hk);
          enktem(ll,:)=Etem;
          unktem(:,:,ll)=psi;
    end
    Unk=reshape(unktem,[dim_H,dim_H,knum,knum]);
    Enk=reshape(enktem,[knum,knum,dim_H]);
end

function [E, Psik]=slovek(k,pars,get_hk)
    [hk,~,~,~]=get_hk(k,pars);
    [V,D]=eig(hk);
    [E,ind]=sort(diag(D));
    Psik=V(:,ind);
end


