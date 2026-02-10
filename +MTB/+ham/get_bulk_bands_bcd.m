function [Energy,Omega_k,kpath,kk]=get_bulk_bands_bcd(obj,hkpoints,nk,bandindex)
    a=obj.a;
    b=obj.b;
    hamiltonian=obj.ham;
    hopping_r=obj.hopr;
    nbands=size(hamiltonian,1);
    nrpts=size(hopping_r,1);
    kTrans=@(x)x(1)*b(1,:)+x(2)*b(2,:)+x(3)*b(3,:); %This can be done b'*hk'
    kPoints=@(startPoint,endPoint,numPoints)...
        [linspace(startPoint(1),endPoint(1),numPoints)',...
        linspace(startPoint(2),endPoint(2),numPoints)',...
        linspace(startPoint(3),endPoint(3),numPoints)'...
        ];
    hkpoints2=cell(1,length(hkpoints));
    for i=1:length(hkpoints)
        hkpoints2(i)={kTrans(hkpoints{i})};
    end

    kpoints=[];
    for i=1:length(hkpoints2)-1
        kpoints=[kpoints;kPoints(hkpoints2{i},hkpoints2{i+1},nk)];
        kpoints(end,:)=[];
    end
    kpoints=[kpoints;kTrans(hkpoints{end})];
    
    kpath=zeros(1,length(kpoints));
    kpath(1)=0;
    for i=2:length(kpoints)
        temp=kpath(i-1);
        kpath(i)=norm(kpoints(i,:)-kpoints(i-1,:))+temp;
    end
    kk=zeros(1,length(hkpoints2));
    for i=1:length(hkpoints2)
        temp=(i-1)*(nk-1)+1;
        kk(i)=kpath(temp);
    end

    %%solve the Hamiltonian(k)
    % hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    enktem=zeros(length(kpoints),nbands);
    unktem=zeros(nbands,nbands,length(kpoints));
    enktem_dk=zeros(length(kpoints),nbands);
    Omega=zeros(length(kpoints),size(bandindex,2));
    for i=1:length(kpoints)
        fprintf("Processing on %d of %d KPOINTS\n", i, length(kpoints));
        knum=301;
        dkx=obj.b(1,:)/knum;
        dky=obj.b(2,:)/knum;
        dS=norm(dkx)*norm(dky);

        kpoint=kpoints(i,:)+dkx;
        hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
        % ishermitian(hk)
        knum=301;
        [V,D]=eig(hk);
        [Etem,ind]=sort(diag(D));
        psi=V(:,ind);
        enktem_dk(i,:)=Etem;

        kpoint=kpoints(i,:);
        hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
        % ishermitian(hk)

        [V,D]=eig(hk);
        [Etem,ind]=sort(diag(D));
        psi=V(:,ind);
        Energy(:,i)=Etem;
        enktem(i,:)=Etem;
        unktem(:,:,i)=psi;



        nsband=size(bandindex,2);
        itotal=knum;
        qm=@MTB.ham.get_quan_metric;
        delta=0.001;
        for j=1:nsband
            band=bandindex(j);
            omega(i,j)=-2.*imag(qm(obj,kpoint,dkx,dky,unktem(:,:,i),Energy(:,i),band,delta));
        end
    end
    Omega_k=omega%.*dS;
end
% 
% function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
%     hk=zeros(nbands,nbands);
%     for j=1:nrpts
%         temp=hamiltonian(:,:,j)*...
%         exp(1j*dot(kpoint,(hopping_r(j,:)*a))); 
%         hk=hk+temp;
%     end
% end

function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    hamiltonian=reshape(hamiltonian,nbands*nbands,nrpts);
    hk=hamiltonian*exp(1j*((hopping_r*a)*kpoint'));
    hk=reshape(hk,nbands,nbands);
    hk=(hk+hk')/2.0;
end
