function [Unk,Enk]=get_bulk_plane_bands_add_electric(obj,Electric_field_in_evpA,Kx,Ky,Kz)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Kx,Ky,Kz: the k-mesh we will solve
    % dim_H: num of the bands/orbits
    % Unk: Norbits*Nands*knum*knum
    % Enk: knum*knum*Nbands
    knum=size(Kx,1);
    % Unk=zeros(dim_H,dim_H,knum,knum);
    % Enk=zeros(knum,knum,dim_H);
    itotal=knum^2;
    dim_H=size(obj.ham,1);
    nrpts=size(obj.ham,3);
    unktem=zeros(dim_H,dim_H,itotal);
    enktem=zeros(itotal,dim_H);
    hamiltonian=obj.ham;
    hopping_r=obj.hopr;
    wpos=obj.wpos;
    a=obj.a;
    parfor ll=1:itotal
        [i,j]=ind2sub([knum,knum],ll);
        kpoint=[Kx(i,j),Ky(i,j),Kz(i,j)];
        % fprintf("Processing on %d %d of %d KPOINTS\n", i, j, itotal);
        hk=get_hk(hamiltonian,hopping_r,wpos,Electric_field_in_evpA,dim_H,nrpts,kpoint,a);
        [V,D]=eig(hk);
        [Etem,ind]=sort(diag(D));
        psi=V(:,ind);
        enktem(ll,:)=Etem;
        unktem(:,:,ll)=psi;
    end
    Unk=reshape(unktem,[dim_H,dim_H,knum,knum]);
    Enk=reshape(enktem,[knum,knum,dim_H]);
end

function hk=get_hk(hamiltonian,hopping_r,wpos,Electric_field_in_evpA,nbands,nrpts,kpoint,a)
    hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    hk=hamiltonian*exp(1j*((hopping_r*a)*kpoint'));
    hk=reshape(hk,nbands,nbands);
    hke=hk_add_elec(wpos,nbands,Electric_field_in_evpA);
    hk=hk+hke;
    hk=(hk+hk')./2;
end


function hke=hk_add_elec(wpos,dim_H,Electric_field_in_evpA)
    %wpos=round(wpos);
    hke=zeros(dim_H,dim_H);
    minrz=min(wpos(:,3));
    maxrz=max(wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    wpos(:,3)=wpos(:,3)-rz;
    for i = 1:dim_H
        hke(i,i)=wpos(i,3)*Electric_field_in_evpA;
    end 
end

% function hke=hk_add_elec(wpos,Electric_field_in_evpA)
%     wpos=(wpos(:,3)-(min(wpos(:,3))+max(wpos(:,3)))/2.0)*Electric_field_in_evpA;
%     hke=diag(wpos);
% end