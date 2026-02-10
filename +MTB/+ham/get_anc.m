function [Eaxis,sigma]=get_anc(Omega_k,Enk,Enum,Emin,Emax,tem)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Omega_k: Thek Berry curvature knum*knum*Nbands
    % Enk: knum*knum*Nbands
    G=25812; %h/e^2
    e=1.6*10^-19; 
    % coef=1/G*10^10;%/e/(6.2415*10^18); %(A/m/K) eV*S/m/K for 3D
    coef=1/G;%/e/(6.2415*10^18); %(A/m/K) eV*S/m/K for 2D there is no length
    Eaxis=linspace(Emin,Emax,Enum);
    sigma=zeros(1,Enum);
    kb=8.61733*10^-5; % eV/K Boltzmann constant
    parfor i=1:Enum
        ef=Eaxis(i);
        fprintf('The Current Ef is %12.6f \n',i)
        % sigma(j)=sum(Omega_k(Enk<E));
        sigma(i)=get_chern_at_ef(Omega_k,Enk,ef,kb,tem);
    end
    sigma=sigma.*coef/tem;
end

function sigma=get_chern_at_ef(Omega_k,Enk,ef,kb,tem)
    E=Enk-ef;
    f=1.0/(exp(E/kb/tem)+1.0);
    phase=(E.*f+kb*tem*log(1+exp(-(E/kb/tem))));
    tmp = find((-(E/kb/tem))>30);
    phase(tmp) = E(tmp).*f(tmp)-E(tmp);
    sigma=sum(Omega_k.*phase,'all');
end