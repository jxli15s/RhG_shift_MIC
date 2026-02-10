function [Enk,Unk,kpath,kk]=get_bulk_bands_full(hamiltonian,hopping_r,nbands,nrpts,hkpoints,nk,a,b)
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
    num_k=length(kpoints);
    Enk=zeros(nbands,num_k);
    Unk=zeros(nbands,nbands,num_k);
    parfor i=1:length(kpoints)
        hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoints(i,:),a);
        % hk=get_hk_new(hamiltonian,hopping_r,nbands,nrpts,kpoints(i,:),a);
        % ishermitian(hk)
        [V,D]=eig(hk);
        [E,ind]=sort(diag(real(D)));
        E=real(E);
        psik=V(:,ind);
        Enk(:,i)=E;
        Unk(:,:,i)=psik;
    end
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

% function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
%     phase=exp(1j*((hopping_r*a)*kpoint'));
%     phase=kron(phase,ones(nbands^2,1));
%     phase=reshape(phase,nbands,nbands,nrpts);
%     hk=hamiltonian.*phase;
%     hk=sum(hk,3);
% end

function hk=get_hk_new(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    phase=reshape(kron(exp(1j*((hopping_r*a)*kpoint')),ones(nbands^2,1)),nbands,nbands,nrpts);
    hk=sum(hamiltonian.*phase,3);
end
        
            
            
            
            