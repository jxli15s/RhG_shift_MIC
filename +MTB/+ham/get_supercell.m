function gs=get_supercell(g,n1,n2)
% Slove the Hamiltonian on the K-mesh 
    % obj: The object of the continuum model system
    % Kx,Ky: the k-mesh we will solve
    % dim_H: num of the bands/orbits
    % Unk: Norbits*Nands*knum*knum
    % Enk: knum*knum*Nbands
    gs = MTB.geometry("supercell");
    gs.a=[n1,0,0;0,n2,0;0,0,1]*g.a;
    gs.b=inv(gs.a')*2*pi;
    pvolm=det(g.a); % primitive cell volum
    svolm=det(gs.a); % supercell volum
    % c=fix(svolm/pvolm)+1;
    norbitalc = size(g.wpos,1); % number of orbitals in a cell
    norbitals = size(g.wpos,1)*n1*n2; %total number of orbitals
    natomc = size(g.atoms,1); % number of atoms in a cell
    natoms = natomc*n1*n2; %total number of atoms
    gs.atoms=zeros(natoms,3);
    gs.wpos=zeros(norbitals,3);
    nx=max(g.hopr(:,1));
    ny=max(g.hopr(:,2));
    hop_dim=[ceil(nx/n1),ceil(ny/n2)];
    nhoprs=(2*hop_dim(1)+1)*(2*hop_dim(2)+1);
    cellindex=zeros(n1*n2,3);
    kk=1;
    ll=1;
    for j=1:n2
        for i=1:n1
            cellindex(n1*(j-1)+i,:)=[i-1,j-1,0]*g.a/gs.a;
            for k=1:natomc
                gs.atoms(kk,:) = (g.atoms(k,:)+[i-1,j-1,0])*g.a/gs.a;
                kk=kk+1;
            end
            for l=1:norbitalc
                gs.wpos(ll,:) = g.wpos(l,:)+[i-1,j-1,0]*g.a;
                ll=ll+1;
            end
            
        end
    end

    gs.ham=zeros(norbitals,norbitals,nhoprs);
    mm=1;
    for i=-hop_dim(1):hop_dim(1)
        for j=-hop_dim(2):hop_dim(2)
            gs.hopr(mm,:)=[i,j,0];
            gs.ham(:,:,mm)=get_t(gs,g,i,j,cellindex,norbitalc,norbitals);
            mm=mm+1;
        end
    end
   
end
    
function ham=get_t(gs,g,i,j,cellindex,norbitalc,norbitals)
    ham=zeros(norbitals,norbitals);
    for l=1:size(cellindex,1)
        for m=1:size(cellindex,1)
            r1=cellindex(l,:);
            r2=cellindex(m,:)+[i,j,0];
            d_frac_r=round((r2-r1)*gs.a/g.a);
            temp=find(ismember(g.hopr,d_frac_r,'rows'));
           if temp
               ham((l-1)*norbitalc+1:l*norbitalc,(m-1)*norbitalc+1:m*norbitalc)=g.ham(:,:,temp);
           end
        end
    end
end

% function hk=get_hk(obj,kpoint)
% %    kpoint=kpoint*obj.b
%     hk=zeros(size(obj.ham,1),size(obj.ham,1));
%     nrpts=size(obj.ham,3);
%     for j=1:nrpts
%         temp=obj.ham(:,:,j)*...
%         exp(1j*dot(kpoint,(obj.hopr(j,:)*obj.a))); 
%         hk=hk+temp;
%     end
% end
