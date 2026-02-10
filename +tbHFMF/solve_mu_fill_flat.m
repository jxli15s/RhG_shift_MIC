function mu = solve_mu_fill_flat(Ek_flat, mesh, opts)
%SOLVE_MU_FILL_FLAT  Solve chemical potential μ from flat eigenvalues Ek_flat.
%
% Ek_flat: Nact x Nk (energies in eV)
% Fermi:
%   f(E) = 1 / (1+exp((E-μ)/kT))
%
% Patch average:
%   <sum_n f>_k = mean over k of (sum over active bands)
%
% Two possible targets:
%   opts.n_target_cell : patch contribution to electrons per unit cell
%       n_cell(μ) = area_frac * <sum_n f>_k
%   opts.n_target_area : patch contribution to areal density (Å^-2)
%       n_area(μ) = pref_patch * <sum_n f>_k

Evec = Ek_flat(:);
emin = min(Evec); emax = max(Evec);

muL = emin - 50*opts.kT - 1;
muR = emax + 50*opts.kT + 1;

use_cell = isfield(opts,'n_target_cell') && ~isempty(opts.n_target_cell);
use_area = isfield(opts,'n_target_area') && ~isempty(opts.n_target_area);
if ~use_cell && ~use_area
  error('solve_mu_fill_flat: set opts.n_target_cell or opts.n_target_area.');
end

Nact = size(Ek_flat,1);

for it = 1:80
  muM = 0.5*(muL + muR);

  f = 1 ./ (1 + exp((Evec - muM)/opts.kT));
  occ_sum_per_k = mean(f) * Nact;   % mean over all states -> times Nact gives mean_k sum_n

  if use_cell
    nM = mesh.area_frac * occ_sum_per_k;
    target = opts.n_target_cell;
  else
    nM = mesh.pref_patch * occ_sum_per_k;
    target = opts.n_target_area;
  end

  if nM > target
    muR = muM;
  else
    muL = muM;
  end

  if abs(muR - muL) < 1e-12, break; end
end

mu = 0.5*(muL + muR);
end
