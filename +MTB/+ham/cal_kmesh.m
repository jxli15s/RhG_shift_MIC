function [Unk,Enk]=cal_kmesh(obj,Kx,Ky,Kz,dim_H)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Kx,Ky: the k-mesh we will solve
    % dim_H: num of the bands/orbits
    % Unk: Norbits*Nands*knum*knum
    % Enk: knum*knum*Nbands
    knum=size(Kx,1);
    % Unk=zeros(dim_H,dim_H,knum,knum);
    % Enk=zeros(knum,knum,dim_H);
    slovek=@(obj,k) obj.solve_kpoint(k);
    itotal=knum^2;
    unktem=zeros(dim_H,dim_H,itotal);
    enktem=zeros(itotal,dim_H);
    parfor ll=1:itotal
          [i,j]=ind2sub([knum,knum],ll); 
          k=[Kx(i,j),Ky(i,j)];
          [Etem, psi] = slovek(obj,k);
          enktem(ll,:)=Etem;
          unktem(:,:,ll)=psi;
    end
    Unk=reshape(unktem,[dim_H,dim_H,knum,knum]);
    Enk=reshape(enktem,[knum,knum,dim_H]);
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

function [Etem,psi]=solve_kpoint(k)
    
end