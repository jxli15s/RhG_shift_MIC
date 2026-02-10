function [Eaxis,Dos,TDos]=get_dos_FermiDirac(En,Tem,Enum,Emin,Emax,nk,plottap)
% 2023-10-19
            %{
                 c: the broadening of Gaussian broadening
                 En:knum*knum*Nband the energy spectrum of the bands
                 Enum: number of energy we will calculate

            %}
            En=En(:);
            Eaxis=linspace(Emin,Emax,Enum);
            Dos=zeros(1,Enum);
            TDos=zeros(1,Enum);
            Occupied=numel(En(En<0))/nk/nk;
            kb=8.617333262*10^-5;
            parfor i=1:Enum
                E=Eaxis(i);
                Dos(i)=sum(1 ./ (4 * kb * Tem) .* sech((En - E) ./ (2 * kb * Tem)).^2)/nk/nk
                % Dos(i)=sum(1./(exp((En-E)/kb/Tem)+1))/nk/nk %Fermi-Dirac distribution
                % Dos(i)=1/(eps*sqrt(2*pi))*sum(exp((-(En-E).^2)./(2*eps*eps)))/nk/nk;%Gaussian distribution
                %Dos(i)=1/pi*sum(eps/(eps^2 + (En-E).^2))
                TDos(i)=abs(Occupied-numel(En(En<E))/nk/nk);
            end
            %dE=(Emax-Emin)/(Enum-1);
            %C=sum(Dos);
            %C=sum(Dos)*dE;
            Nband=TDos(1)+TDos(end);

            %Dos=Dos./C.*Nband;

            if plottap==1
            figure('Color','white')
%            plot(Eaxis,Dos,'k-')
            plot(Eaxis,Dos,'Linestyle','-','Color','#4DA1D7','LineWidth',2)
            hold on;
            % plot(Eaxis,TDos,'Linestyle','-','LineWidth',2)
            xlabel('E(eV)')
            ylabel('Dos')
            set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)
            end
            
        end
