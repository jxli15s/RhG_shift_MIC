function [Energy,kpath,kk]=get_slab_bands_BdG_del(hamiltonian,hopping_r,del_index,nslab,nbands,nrpts,hkpoints,nk,a,b,mu,delta)

    kPoints=@(startPoint,endPoint,numPoints)...
        [linspace(startPoint(1),endPoint(1),numPoints)',...
        linspace(startPoint(2),endPoint(2),numPoints)'...
        ];

%    hkpoints2=zeros(length(hkpoints),2);
    hkpoints=reshape(cell2mat(hkpoints),2,[])';
    hkpoints2=hkpoints*b;

    kpoints=[];
    for i=1:length(hkpoints2)-1
        kpoints=[kpoints;kPoints(hkpoints2(i,:),hkpoints2(i+1,:),nk)];
        kpoints(end,:)=[];
    end

    kpoints=[kpoints;hkpoints2(end,:)];

    kpath=zeros(1,length(kpoints));
    kpath(1)=0;
    for i=2:length(kpoints)
        temp = norm(kpoints(i,:)-kpoints(i-1,:));
        kpath(i)=kpath(i-1)+temp;
    end

    kk=zeros(1,length(hkpoints2));
    for i=1:length(hkpoints2)
        temp=(i-1)*(nk-1)+1;
        kk(i)=kpath(temp);
    end

    %%solve the Hamiltonian(k)
    % Ek=zeros(nbands*nslab*2,length(kpoints));
    % T=zeros(1,length(kpoints));
    % hbar = parfor_progressbar_v1(length(kpoints),'Computing...');
    parfor i=1:length(kpoints)
        fprintf('Current running k_i = %2d \n',i)
        [vals,~]=MTB.ham.get_slab_bands_BdG_at_q_odd_v2(hamiltonian,hopping_r,del_index,nslab,nbands,nrpts,kpoints(i,:),a,b,mu,delta);
        Energy(:,i)=vals;
    end
    % close(hbar);
end


