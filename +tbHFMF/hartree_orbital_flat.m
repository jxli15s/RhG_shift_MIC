function SigmaH = hartree_orbital_flat(V0_cont, rho_orb_flat, mesh, opts)
%HARTREE_ORBITAL_FLAT  Hartree term (q=0 only) from flat rho_orb.
%
% For translational invariant HF, Hartree uses q=0 component:
%   Σ_H,aa = Σ_b V_cont(q->0; z_a,z_b) * n_area(b)
%
% where areal density contributed by the patch:
%   n_area(b) = ∫_{patch} d^2k/(2π)^2 ρ_bb(k)
%            ≈ mesh.pref_patch * mean_k[ρ_bb(k)]
% Units:
%   V0_cont : eV·Å^2
%   n_area  : Å^-2
% => Σ_H    : eV

Norb = size(rho_orb_flat,1);

n_area = zeros(Norb,1);
for b = 1:Norb
  n_area(b) = mesh.pref_patch * real(mean(rho_orb_flat(b,b,:)));
end

% Optional: subtract uniform background to avoid large uniform Hartree shift
if isfield(opts,'hartree_subtract_mean') && opts.hartree_subtract_mean
  n_area = n_area - mean(n_area);
end

SigmaH = diag( real(V0_cont) * n_area );
SigmaH = (SigmaH + SigmaH')/2;

end
