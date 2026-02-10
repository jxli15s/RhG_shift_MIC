function [Enk,Unk,kpath,kk]=get_bulk_bands_psi_add_electric(hamiltonian,hopping_r,wpos,Electric_field_in_evpA,nbands,nrpts,hkpoints,nk,a,b)
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
    Enk=zeros(nbands,length(kpoints));
    Unk=zeros(nbands,nbands,length(kpoints));
    parfor i=1:length(kpoints)
        hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoints(i,:),a,wpos,Electric_field_in_evpA);
        % vals=sort(real(eig(hk)));
        [V,D]=eig(hk);
        [Etem,ind]=sort(diag(D));
        psi=V(:,ind);
        Enk(:,i)=Etem;
        Unk(:,:,i)=psi;
    end
end

function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a,wpos,Electric_field_in_evpA)
    hk=zeros(nbands,nbands);
    for j=1:nrpts
        temp=hamiltonian(:,:,j)*...
        exp(1j*dot(kpoint,(hopping_r(j,:)*a))); 
        hk=hk+temp;
    end
    hke=hk_add_elec(wpos,nbands,Electric_field_in_evpA);
    hk=hk+hke;
    hk=(hk+hk')/2.0;
end

function hke=hk_add_elec(wpos,nbands,Electric_field_in_evpA)
%    wpos=round(wpos)
    hke=zeros(nbands,nbands);
    minrz=min(wpos(:,3));
    maxrz=max(wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    wpos(:,3)=wpos(:,3)-rz;
    for i = 1:nbands
        hke(i,i)=wpos(i,3)*Electric_field_in_evpA;
    end 
end