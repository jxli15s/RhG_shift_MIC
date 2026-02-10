function [wx,unk]=get_wilsonloop_xy(obj,knum,band1,band2)
[~,KX,KY]=K_mesh(obj,knum);
numsub=size(obj.ham,1);
unk=zeros(numsub,numsub,knum,knum);
tau=obj.wpos(:,1:2);
parfor i=1:knum
    i
    for j=1:knum
        % k=[KX(i,j),KX(i,j),KY(i,j)];
        % k=[0,KX(i,j),KY(i,j)]; % for Mirror x
        k=[KX(i,j),KX(i,j),KY(i,j)]; % for (110) plane
        [~,unk(:,:,i,j)]=obj.solve_kpoint(k);
        % phase=(exp(-1i.*k*tau'))';
        % unk(:,:,i)=phase.*unk(:,:,i);
    end
end

G1=obj.b(1,1:2);
G2=obj.b(2,1:2);
unk(:,:,knum,:)=unk(:,:,1,:);%.*(exp(-1i.*G1*tau'))';
% unk(:,:,:,knum)=unk(:,:,:,1);%.*(exp(-1i.*G2*tau'))';
nband=band2-band1+1;
wx=zeros(knum,nband);
for i=1:knum
    tem=eye(nband);
    for j=1:knum-1
        tem=tem*unk(:,band1:band2,j,i)'*unk(:,band1:band2,j+1,i);
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
Gm1=obj.b(1,1:2);
Gm2=obj.b(2,1:2);
%kx=linspace(0,1,knum+1);
kx=linspace(-0.5,0.5,knum+1);
ky=linspace(0,1,knum+1);
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
KX=Kx.*Gm1(1)+Ky.*Gm2(1);
KY=Kx.*Gm1(2)+Ky.*Gm2(2);

end