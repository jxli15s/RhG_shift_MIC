function [Eaxis,bcd,Omega_k_df]=get_bcd(obj,Omega_k,Enk,Enk_d,Enum,Emin,Emax,tem)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Omega_k: Thek Berry curvature knum*knum*Nbands
    % Enk: knum*knum*Nbands
    G=25812; %h/e^2
    e=1.6*10^-19;
    tau=1*10^-12;
    Ey=1000;
    Dxz=1*10^-10;
    hbar=1.054572*10^-34;
    %coef=tau*Dxz*e*Ey*pi/hbar%tau*D_xz*e*Ey*pi/hbar*(1/G)
    coef=1;
    %coef=1/G*10^8; %(Ohm cm)^-1 for 3D A to cm
    knum=size(Omega_k,1);
    dkx=obj.b(1,:)./knum;
    dkx=norm(dkx);
    Eaxis=linspace(Emin,Emax,Enum);
    bcd=zeros(1,Enum);
    kb=8.61733*10^-5; % eV/K Boltzmann constant
    for i=1:Enum
        ef=Eaxis(i);
        % sigma(j)=sum(Omega_k(Enk<E));
        fi=get_chern_at_ef(Enk,ef,kb,tem);
        fi_d=get_chern_at_ef(Enk_d,ef,kb,tem);
        bcd(i)=sum((fi_d-fi)./dkx.*Omega_k,"all");
        Omega_k_df=(fi_d-fi)./dkx.*Omega_k;
        % bcd(i)=sum((fi_d-fi)./dkx.*Omega_k,"all");
        % Omega_k_df=(fi_d-fi)./dkx.*Omega_k;        
    end
    bcd=bcd.*coef;
end

function fi=get_chern_at_ef(Enk,ef,kb,tem)
    fi=1.0/(exp((Enk-ef)/kb/tem)+1.0);%Fermi-Dirac statistics
end
