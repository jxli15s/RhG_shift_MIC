function Energy= get_votex_sparse_at_q(hamiltonian,hopping_r,nbands,nrpts,kpoint,a,b,mu,h_delta,numEigs)
% Build s-wave vortex BdG Hamiltonian from xy open boundary tight-binding
%
% Inputs:
%   gs : structure from get_xy_open_wannier
%   mu : chemical potential
%   delta0 : pairing amplitude (can be real)
%   center : [x0, y0] vortex center (optional, default = system center)
%
% Output:
%   Hbdg : sparse BdG Hamiltonian

% Number of orbitals
h_onsite=speye(nbands)*mu;
kpoint=kpoint*b;
hk_e=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a);
hk_e=hk_e-h_onsite;
hk_h=get_hk(hamiltonian,hopping_r,nbands,nrpts,-kpoint,a);
hk_h=-hk_h.';
hk_h=hk_h+h_onsite;
hk=[hk_e,h_delta;h_delta', hk_h];
% ishermitian(hk)
hk=(hk+hk')/2;
% ishermitian(hk)
vals=sort(eigs(hk,numEigs,'smallestabs'));
Energy=vals;
end       
            
function hk=get_hk(hamiltonian,hopping_r,nbands,nrpts,kpoint,a)
hk = sparse(nbands, nbands);         % 初始化为 norb x norb 的全零 sparse 矩阵
for j=1:nrpts
    temp=hamiltonian{j}*...
        exp(1j*dot(kpoint,(hopping_r(j,:)*a)));
    hk=hk+temp;
end
hk=(hk+hk')/2.0;
end