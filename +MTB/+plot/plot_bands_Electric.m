function plot_bands_Electric(Energy,nbands,efermi,kpath,labels,kk,figname,E_gate)
linesize=2;
xrange=[kpath(1),kpath(end)];
figure();
for i=1:length(Energy(:,1))
    plot(kpath,Energy(i,:)-efermi,'Color','black','LineWidth',linesize);
    hold on
end
plot(kpath,zeros(1,length(kpath)),'--black','LineWidth',2)
for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[min(min(Energy))-efermi-1 max(max(Energy))-efermi-1],'--k','LineWidth',linesize)
end
grid off
box on
xlim(xrange)
xticks(kk)
xticklabels(labels)
ylim([-1 1])
% yticks(-1::1)
ylabel('Energy (eV)','FontSize',24)
%title(num2str(E_gate)+'meV','Interpreter','latex','FontSize',24,'Color',[0.5 0 0.5])
titleE=E_gate+" mV/nm";
title(titleE,'FontSize',24,'Color',[0.5 0 0.5]);
ax=gca;
ax.LineWidth=1;
ax.XColor=[0.5 0 0.5];
ax.YColor=[0.5 0 0.5];
ax.YAxis.FontSize=18;
ax.XAxis.FontSize=18;
% ax.FontWeight='bold'
ax.Layer='bottom';
ax.PlotBoxAspectRatio=[1.5 1 1];

%set(gca,'LineWidth',1,'Color',[1 1 1],'XColor',[0.5 0 0.5],'YColor',[0.5 0 0.5],'Layer','top','PlotBoxAspectRatio',[1.5 1 1])
%set(gca,'LineWidth',1,'Color',[1 1 1],'XColor',[0.5 0 0.5],'YColor',[0.5 0 0.5],'Layer','top')
print(figname+'band_wannier90','-dpng','-r600')
end