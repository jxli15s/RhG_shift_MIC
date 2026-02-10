function wx=get_wilsonloop_slab_BdG_v4_del(obj,k_mesh_dim,nslab,band1,band2,del_index,mu,delta)
[nbands,~,nrpts]=size(obj.ham);
[KX,KY]=K_mesh(k_mesh_dim);
numsub=(nslab*nbands-size(del_index,1))*2;
knumx=k_mesh_dim(1);
knumy=k_mesh_dim(2);

nband=band2-band1+1;
wx=zeros(knumy,nband);
parfor i=1:knumy
    unk=zeros(numsub,numsub,knumx);
    tic;
    for j=1:knumx
        fprintf("Processing on %d %d of %d KPOINTS\n", i, j, knumx*knumy);
        kpoint=[KX(i,j),KY(i,j)];
        % [~,unk(:,:,j)]=MTB.ham.get_slab_bands_BdG_at_q(obj.ham,obj.hopr2,nslab,nbands,nrpts,kpoint,obj.a2,obj.b2,mu,delta)
        [~,unk(:,:,j)]=MTB.ham.get_slab_bands_BdG_at_q_odd(obj.ham,obj.hopr2,del_index,nslab,nbands,nrpts,kpoint,obj.a2,obj.b2,mu,delta);
    end
    unk(:,:,knumx)=unk(:,:,1);%.*(exp(-1i.*G1*tau'))';
    tem=eye(nband);
    for j=1:knumx-1
        tem=tem*unk(:,band1:band2,j)'*unk(:,band1:band2,j+1);
    end
    wx(i,:)=sort(angle(eig(tem))) / pi;
    toc;
end


%Plot
%kx=linspace(0,1,knumx);
%figure('Color','white')
%plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
%xticks([0,1/2,1])
%xticklabels({'0','\pi','2\pi'})
%ylim([-1,1])
%ylabel('Wilson loop bands')
%set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

end


function [Kx,Ky]=K_mesh(k_mesh_dim)
kx=linspace(0,1,k_mesh_dim(1));
ky=linspace(0,0.12,k_mesh_dim(2));
[Kx,Ky]=meshgrid(kx,ky);
end
