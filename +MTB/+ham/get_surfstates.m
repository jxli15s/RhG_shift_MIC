function [dos_l,dos_r,dos_bulk,kk,kpath]=get_surfstates(hamiltonian,hopping_r,nbands,nrpts,hkpoints,nk,Np,a,b,omegamax,omegamin,omeganum)

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
    
    omegas=linspace(omegamin,omegamax,omeganum);
    eta=(omegamax-omegamin)/omeganum*3.0;

    knum=size(kpoints,1);
    dimG=[omeganum,knum];
    itotal=prod(dimG);
    dos_l=zeros(1,itotal);
    dos_r=zeros(1,itotal);
    dos_bulk=zeros(1,itotal);

    parfor l=1:itotal
        [omega_j,k_i]=ind2sub(dimG,l);
        fprintf('Processing on K_i %3d Omega_j %3d\n',k_i,omega_j)
        [GLL,GRR,GB]=MTB.ham.get_surfgreen(hamiltonian,hopping_r,nbands,nrpts,kpoints(k_i,:),Np,a,omegas(omega_j),eta);
        % dos_l(l)=sum(-imag(diag(GLL)));
        for i=1:size(GLL,1)
            dos_l(l)=dos_l(l)-imag(GLL(i,i))
        end
        dos_r(l)=sum(-imag(diag(GRR)));
        dos_bulk(l)=sum(-imag(diag(GB)));
    end

    dos_l=reshape(dos_l,dimG);
    dos_r=reshape(dos_r,dimG);
    dos_bulk=reshape(dos_bulk,dimG);
end



 
