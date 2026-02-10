function [Wk,Psik,u_k_super]=get_bulk_unfolding_bands(hamiltonian,hopping_r,nbands,nrpts,hkpoints,nk,g,gs)
    a=gs.a;
    b=gs.b;
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
    parfor i=1:length(kpoints)
        hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoints(i,:),a);
        % hk=get_hk_new(hamiltonian,hopping_r,nbands,nrpts,kpoints(i,:),a);
        % ishermitian(hk)
        % vals=sort(eig(hk));
        % Energy(:,i)=vals;
        [V,D]=eig(hk);
        [E,ind]=sort(diag(D));
        Energy(:,i)=E;
        Psik(:,:,i)=V(:,ind);
    end
    u_k_super=construct_supercell_wavefunction(kpoints,Psik, g,gs);
    % Wk=compute_projection(u_k_super,Psik);  
    Wk=compute_projection_v2(kpoints,Psik,g,gs);
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

function u_k_super=construct_supercell_wavefunction(kpoints,Psik, g,gs)
    R=floor((gs.wpos)*inv(g.a))*g.a;
    [nbands,~,numk]=size(Psik);
    [nbands_u,~,~]=size(g.ham);
    % u_k_super=zeros(nbands,numk);
    for ik=1:numk
        Us2u=exp(1j*R*kpoints(ik,:)');
        % u_k_super(:,:,ik)=kron(ones(1,nbands_u),Us2u);
        u_k_super(:,ik)=Us2u;

    end
end

function Wk=compute_projection(u_k_super,Psik)
    [nbands,~,numk]=size(Psik);
    Wk=zeros(nbands,numk);
    for ik=1:numk
        % Wk(:,ik) = sum(abs(Psik(:,:,ik)' * u_k_super(:,ik)).^2, 2);
        Wk(:,ik) = abs(Psik(:,:,ik)' * u_k_super(:,ik)).^2;
        % Wk(:,ik)=abs(diag(Psik(:,:,ik)' * (diag(u_k_super(:,ik))*Psik(:,:,ik))));
        % Wk(:,ik)=sum(abs(Psik(:,:,ik)' * (diag(u_k_super(:,ik))*Psik(:,:,ik))).^2,2);
    end
end

function Wk=compute_projection_v2(kpoints,Psik,g,gs)
    [nbands,~,numk]=size(Psik);
    Wk=zeros(nbands,numk);
    R=floor((gs.wpos)*inv(g.a))*g.a;
    num_u=size(g.wpos,1);
    parfor ik=1:numk
        % Wk(:,ik) = sum(abs(Psik(:,:,ik)' * u_k_super(:,ik)).^2, 2);
        for ib=1:nbands
            for i=1:nbands
                for j=1:nbands
                    if mod(i,num_u)==mod(j,num_u)
                        Wk(ib,ik)=Wk(ib,ik)+...
                            exp(-1j*(R(i,:)-R(j,:))*kpoints(ik,:)')*Psik(i,ib,ik)*Psik(j,ib,ik)';
                    end
                end
            end
        end
    end
end
            
            
            