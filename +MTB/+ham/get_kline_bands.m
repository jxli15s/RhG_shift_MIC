function Enk=get_kline_bands(obj,kline)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % kline: the k-line knum*3
    % dim_H: num of the bands/orbits
    % Unk: Norbits*Nands*knum*knum
    % Enk: knum*knum*Nbands
    knum=size(kline,1);
    dim_H=size(obj.ham,1);
    % unktem=zeros(dim_H,dim_H,itotal);
    Enk=zeros(knum,dim_H);
    parfor ik=1:knum
          k=kline(ik,:);
          [Etem, ~] = slovek(obj,k);
          Enk(ik,:)=Etem;
    end
end

function [E, Psik]=slovek(obj,k)
    % hk=get_hk(obj,k);
    hk=obj.get_hk(k);
    [V,D]=eig(hk);
    [E,ind]=sort(diag(D));
    Psik=V(:,ind);
end