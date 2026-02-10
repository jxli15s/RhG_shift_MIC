function act = diag_hk_active(hk, act_idx)
%DIAG_HK_ACTIVE  Diagonalize h(k) for each k (PARFOR only), keep active bands.
%
% Input:
%   hk     : Norb x Norb x Nky x Nkx
%   act_idx: selected band indices after sorting eigenvalues ascending
%
% Output:
%   act.U_flat   : Norb x Nact x Nk    (Nk = Nky*Nkx)
%   act.eps_flat : Nact x Nk
%
% NOTE:
%   We keep k as a single "flat" dimension Nk for memory locality and easy parfor slicing.
%   This does NOT reduce raw memory (same #elements), but reduces temporary reshapes/copies.

[Norb,~,Nky,Nkx] = size(hk);
Nk   = Nky*Nkx;
Nact = numel(act_idx);

hk_flat = reshape(hk, Norb, Norb, Nk);

U_flat   = zeros(Norb, Nact, Nk, 'like', 1+1i);
eps_flat = zeros(Nact, Nk);

parfor ik = 1:Nk
  H = hk_flat(:,:,ik);
  H = (H + H')/2;        % enforce Hermitian numerically
  [V,D] = eig(H,'vector');
  [d,ord] = sort(real(D),'ascend');
  V = V(:,ord);
  U_flat(:,:,ik) = V(:,act_idx);
  eps_flat(:,ik) = d(act_idx);
end

act.Nact     = Nact;
act.idx      = act_idx;
act.U_flat   = U_flat;
act.eps_flat = eps_flat;

% Optional 4D views (ONLY for plotting/debug)
act.U   = reshape(U_flat,   Norb, Nact, Nky, Nkx);
act.eps = reshape(eps_flat, Nact, Nky, Nkx);

end
