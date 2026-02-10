function ph = precompute_dtau_phases(mesh, tau)
%PRECOMPUTE_DTAU_PHASES  Precompute exp(±i p·(τ_a-τ_b)) for all orbital pairs.
%
% In atom gauge, the density vertex / orbital phase enters the Fock term as:
%   exp[-i (k-k')·(τ_a-τ_b)] = exp[-i q·dtau_ab]
%
% For FFT convolution we use local p-grid and q-grid on the patch.
% We precompute:
%   phase_plus  = exp(+i p·dtau_ab)
%   phase_minus = exp(-i p·dtau_ab)
% in centered order (p=0 at center).

Norb = size(tau,1);
Nky  = mesh.Nky; Nkx = mesh.Nkx;
P = Norb*Norb;

[a_list,b_list] = ndgrid(1:Norb,1:Norb);
a_list = a_list(:); b_list = b_list(:);

Px = mesh.Px; Py = mesh.Py; Pz = mesh.Pz;

phase_plus  = zeros(Nky,Nkx,P,'like',1+1i);
phase_minus = zeros(Nky,Nkx,P,'like',1+1i);

parfor p = 1:P
  a = a_list(p); b = b_list(p);
  dt = tau(a,:) - tau(b,:);
  theta = Px*dt(1) + Py*dt(2) + Pz*dt(3);
  phase_plus(:,:,p)  = exp( 1i*theta );
  phase_minus(:,:,p) = exp(-1i*theta );
end

ph.P = P;
ph.Norb = Norb;
ph.a_list = a_list;
ph.b_list = b_list;
ph.phase_plus  = phase_plus;
ph.phase_minus = phase_minus;

end
