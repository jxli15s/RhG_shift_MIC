function wx=get_wilsonloop_slab_BdG_v2_del(obj,knum,nslab,band1,band2,del_index,mu,delta)
[nbands,~,nrpts]=size(obj.ham);
[~,KX,KY]=K_mesh(obj,knum);
numsub=(nslab*nbands-size(del_index,1))*2;

% tau=obj.wpos(:,1:2);
G1=obj.b2(1,1:2);
G2=obj.b2(2,1:2);
nband=band2-band1+1;
wx=zeros(knum,nband);
parfor i=1:knum
    unk=zeros(numsub,numsub,knum);
    for j=1:knum
        fprintf("Processing on %d %d of %d KPOINTS\n", i, j, knum*knum);
        kpoint=[KX(i,j),KY(i,j)];
        % [~,unk(:,:,j)]=MTB.ham.get_slab_bands_BdG_at_q(obj.ham,obj.hopr2,nslab,nbands,nrpts,kpoint,obj.a2,obj.b2,mu,delta)
        [~,unk(:,:,j)]=MTB.ham.get_slab_bands_BdG_at_q_odd(obj.ham,obj.hopr2,del_index,nslab,nbands,nrpts,kpoint,obj.a2,obj.b2,mu,delta);
    end
    unk(:,:,knum)=unk(:,:,1);%.*(exp(-1i.*G1*tau'))';
    tem=eye(nband);
    for j=1:knum-1
        tem=tem*unk(:,band1:band2,j)'*unk(:,band1:band2,j+1);
        % tem=tem*unk(:,band1:band2,i,j)'*unk(:,band1:band2,i,j+1);
    end
    wx(i,:)=sort(angle(eig(tem))) / pi;
end


%Plot
kx=linspace(0,1,knum);
figure('Color','white')
plot(kx,wx,'.','Color','#007EC9','MarkerSize',15)
xticks([0,1/2,1])
xticklabels({'0','\pi','2\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

end


function [Kpts,KX,KY]=K_mesh(obj,knum)
Gm1=obj.b2(1,1:2);
Gm2=obj.b2(2,1:2);
kx=linspace(0,1,knum);
ky=linspace(0,1,knum);
kx=kx(1:knum);
ky=ky(1:knum);
[Kx,Ky]=meshgrid(kx,ky);
num_kpts=knum^2;
Kpts=zeros(num_kpts,2);
tem=0;
for j=1:knum
    for i=1:knum
        tem=tem+1;
        Kpts(tem,:)=Kx(i,j).*Gm1+Ky(i,j).*Gm2;
    end
end
KX=Kx;
KY=Ky;
end
