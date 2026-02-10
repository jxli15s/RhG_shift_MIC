function [Eaxis,sigma]=get_ahc(Omega_k,Enk,Enum,Emin,Emax,tem)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Omega_k: Thek Berry curvature knum*knum*Nbands
    % Enk: knum*knum*Nbands
    G=25812; %h/e^2
    e=1.6*10^-19;
    % coef=1/G*10^8; %(Ohm cm)^-1 for 3D A to cm
    %coef=1/G; %(Ohm cm)^-1 for 2D no length
    coef=1;
    Eaxis=linspace(Emin,Emax,Enum);
    sigma=zeros(1,Enum);
    kb=8.61733*10^-5; % eV/K Boltzmann constant
    for i=1:Enum
        ef=Eaxis(i);
        % sigma(j)=sum(Omega_k(Enk<E));
        sigma(i)=get_chern_at_ef(Omega_k,Enk,ef,kb,tem);
    end
    sigma=sigma.*coef;
end

function sigma=get_chern_at_ef(Omega_k,Enk,ef,kb,tem)
    fi=1.0/(exp((Enk-ef)/kb/tem)+1.0);%Fermi-Dirac statistics
    sigma=sum(Omega_k.*fi,'all');
end