function [Unk,Enk,vxk,vyk,vzk]=get_bulk_plane_kp_velocity(pars,get_hk,nbands,Kx,Ky,Kz)
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
    vxtem=zeros(dim_H,dim_H,itotal);
    vytem=zeros(dim_H,dim_H,itotal);
    vztem=zeros(dim_H,dim_H,itotal);
    parfor ll=1:itotal
          [i,j]=ind2sub([knum,knum],ll); 
          k=[Kx(i,j),Ky(i,j),Kz(i,j)];
          [Etem, psi, vx, vy, vz] = slovek(k,pars,get_hk);  
          enktem(ll,:)=Etem;
          unktem(:,:,ll)=psi;
          vxtem(:,:,ll)=vx;
          vytem(:,:,ll)=vy;
          vztem(:,:,ll)=vz;
    end
    Unk=reshape(unktem,[dim_H,dim_H,knum,knum]);
    Enk=reshape(enktem,[knum,knum,dim_H]);
    vxk=reshape(vxtem,[dim_H,dim_H,knum,knum]);
    vyk=reshape(vytem,[dim_H,dim_H,knum,knum]);
    vzk=reshape(vztem,[dim_H,dim_H,knum,knum]);
end

function [E, Psik,vx,vy,vz]=slovek(k,pars,get_hk)
    [hk,hkx,hky,hkz]=get_hk(k,pars);
    [V,D]=eig(hk);
    [E,ind]=sort(diag(D));
    Psik=V(:,ind);
    Psikp = Psik';
    vx = Psikp  * hkx * Psik;
    vy = Psikp  * hky * Psik;
    vz = Psikp  * hkz * Psik;
end


