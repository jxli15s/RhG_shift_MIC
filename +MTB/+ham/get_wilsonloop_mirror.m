function [wx,unk1]=get_wilsonloop_mirror(obj,M,knum,band1,band2)
[~,KX,KY]=K_mesh(obj,knum);
numsub=size(obj.ham,1)/2;
unk1=zeros(numsub,numsub,knum,knum);
unk2=zeros(numsub,numsub,knum,knum);
tau=obj.wpos(:,1:2);

parfor i=1:knum
    i
    for j=1:knum
        % k=[KX(i,j),KX(i,j),KY(i,j)];
        k=[0,KX(i,j),KY(i,j)]; % for Mirror_x
        % k=[KX(i,j),KX(i,j),KY(i,j)]; % for Mirror_xy
        [~,~,unk1(:,:,i,j),~]=solve_kpoint(obj,M,k);
        [~,~,~,unk2(:,:,i,j)]=solve_kpoint(obj,M,k);
        % phase=(exp(-1i.*k*tau'))';
        % unk(:,:,i)=phase.*unk(:,:,i);
    end
end

unk1(:,:,knum,:)=unk1(:,:,1,:);%.*(exp(-1i.*G1*tau'))';
% unk(:,:,:,knum)=unk(:,:,:,1).*(exp(-1i.*G2*tau'))';
nband=band2-band1+1;
wx1=zeros(knum,nband);
for i=1:knum
    tem=eye(nband);
    for j=1:knum-1
        tem=tem*unk1(:,band1:band2,j,i)'*unk1(:,band1:band2,j+1,i); % along ky give us kx evolution
        % tem=tem*unk(:,band1:band2,i,j)'*unk(:,band1:band2,i,j+1); 
    end
    wx1(i,:)=sort(angle(eig(tem))) / pi;
end

unk2(:,:,knum,:)=unk2(:,:,1,:);%.*(exp(-1i.*G1*tau'))';
% unk(:,:,:,knum)=unk(:,:,:,1).*(exp(-1i.*G2*tau'))';
nband=band2-band1+1;
wx2=zeros(knum,nband);
for i=1:knum
    tem=eye(nband);
    for j=1:knum-1
        tem=tem*unk2(:,band1:band2,j,i)'*unk2(:,band1:band2,j+1,i);
        % tem=tem*unk(:,band1:band2,i,j)'*unk(:,band1:band2,i,j+1);
    end
    wx2(i,:)=sort(angle(eig(tem))) / pi;
end

%Plot
kx=linspace(0,1,knum);
figure('Color','white')
plot(kx,wx1,'.','Color','#007EC9','MarkerSize',15)
hold on;
plot(kx,wx2,'.','Color','#c92f00','MarkerSize',15)
xticks([0,1/2,1])
xticklabels({'-\pi','0','\pi'})
ylim([-1,1])
ylabel('Wilson loop bands')
set(gca,'Fontsize',20,'FontName','Times New Roman','linewidth',0.8)

wx=[wx1,wx2];

end


function [Kpts,KX,KY]=K_mesh(obj,knum)
% Gm1=obj.b(1,1:2);
% Gm2=obj.b(2,1:2);
Gm1=obj.b(2,2:3);
Gm2=obj.b(3,2:3);
%kx=linspace(0,1,knum+1);
kx=linspace(-0.5,0.5,knum+1);
% kx=linspace(-0.17,0.17,knum+1);
ky=linspace(0,1,knum+1);
kx=kx(1:knum);
ky=ky(1:knum);
[Kx,Ky]=meshgrid(kx,ky);
num_kpts=knum^2;
Kpts=zeros(num_kpts,2);
tem=0;
fprintf("here")
for j=1:knum
    for i=1:knum
        tem=tem+1;
        Kpts(tem,:)=Kx(i,j).*Gm1+Ky(i,j).*Gm2;
    end
end
% KX=Kx.*Gm1(1)+Ky.*Gm2(1);
% KY=Kx.*Gm1(2)+Ky.*Gm2(2);
KX=Kx;
KY=Ky;
end

function [e1,e2,psik1,psik2]=solve_kpoint(obj,M,k)
    %slove_kpoint: this function give the energy level and eigenstate at the given
    %k-point.
    k=k*obj.b;
    % for i =1:size(obj.wpos,1)
    %     M(i,:)=exp(-2j*(obj.wpos(i,1))*k(1)).*M(i,:);
    %     % m(i,:)=exp(-2j*(g.wpos(i,:))*kpoint').*m(i,:);
    % end
    HK=obj.get_hk(k);
    HK=inv(M)*HK*M;
    hk1=HK(1:size(HK,1)/2,1:size(HK,2)/2);
    [V1,D1]=eig(hk1);
    [e1,ind1]=sort(diag(real(D1)));
    psik1=V1(:,ind1);
    hk2=HK(size(HK,1)/2+1:end,size(HK,2)/2+1:end);
    [V2,D2]=eig(hk2);
    [e2,ind2]=sort(diag(real(D2)));
    psik2=V2(:,ind2);
    % [E,ind]=sort(diag(D));
    E=[e1;e2];
    % E=real(E);
    Psik=blkdiag(psik1,psik2);
end
