function [Energy,kpath,kk]=get_slab_bands(hamiltonian,hopping_r,nslab,nbands,nrpts,hkpoints,nk,a,b)

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
    Ek=zeros(nbands,length(kpoints));
%    hbar = parfor_progressbar_v1(length(kpoints),'Computing...');
    parfor i=1:length(kpoints)
        fprintf('Current running k_i = %2d \n',i)
        hk=MTB.ham.get_slab_hk(hamiltonian,hopping_r,nslab,nbands,nrpts,kpoints(i,:),a);
%        hbar.iterate(1);
%        fprintf('Job on k_i = %2d time: %8.4f s\n',i,T(i))
        vals=sort(eig(hk));
        Energy(:,i)=vals;
    end
%    close(hbar);
end

function hk=get_slab_hk(hamiltonian,hopping_r,nslab,nbands,nrpts,kpoint,a)
    hk=zeros(nbands*nslab,nbands*nslab);
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij=ham_qlayer2qlayer(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
    for i1=1:nslab
        for i2=1:nslab
            if (ijmin<=(i2-i1)) && ((i2-i1) <= ijmax)
                hk((i2-1)*nbands+1:(i2-1)*nbands+nbands, ...
                    (i1-1)*nbands+1:(i1-1)*nbands+nbands) ...
                = Hij(:,:,i1-i2+abs(ijmin)+1);
            end
        end
    end
    hk=(hk+hk')/2.0;
end

function Hij=ham_qlayer2qlayer(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
    ijmax=max(hopping_r(:,3));
    ijmin=min(hopping_r(:,3));
    Hij = zeros(nbands,nbands,ijmax-ijmin+1);
    for j=1:nrpts
        if  (ijmin<= hopping_r(j,3) <= ijmax)
            kdotr = dot(kpoint,hopping_r(j,1:2)*a);
            ratio = exp(1j*kdotr);
            Hij(:,:,hopping_r(j,3)+abs(ijmin)+1) = Hij(:,:,hopping_r(j,3)+abs(ijmin)+1)+hamiltonian(:,:,j)*ratio;
        end
    end
end

