function [hamiltonian,hopping_r] = read_hr(filename1,filename2)
    %read wannier90_hr.dat
    %Return:
    %hamiltonian: real space hopping parameters
    %hopping_r: real space hopping coordinates
    data1=importdata(filename1);
    numData=data1.data;
    nbands=numData(1);
    nrpts=numData(2);
    degrpts=numData(3:2+nrpts);
    data2=readmatrix(filename2);
    hopping_r=data2(1:nbands*nbands:end,1:3);
    hamiltonian=data2(:,6)+1j*data2(:,7);
    degrpts=kron(degrpts,ones(nbands*nbands,1));
    hamiltonian=hamiltonian./degrpts;
    hamiltonian=reshape(hamiltonian,[nbands,nbands,nrpts]);
    % for i=1:nrpts
    %     hamiltonian(:,:,i)=hamiltonian(:,:,i)/degrpts(i);
    % end
end