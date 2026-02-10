function SigmaF_flat = fock_fft_atomgauge_flat(mesh, V, ph, rho_orb_flat)
%FOCK_FFT_ATOMGAUGE_FLAT  FFT-based Fock self-energy in atom gauge (flat-k storage).
%
% Target (scheme A, orbital basis):
%   Σ_F,ab(k) = - ∫_{patch} d^2k'/(2π)^2  V_cont(k-k'; z_a,z_b)
%               * exp[-i (k-k')·(τ_a-τ_b)] * ρ_ab(k')
%
% On a uniform patch:
%   ∫ d^2k'/(2π)^2 (...)  ≈  mesh.pref_patch * mean_k'(...)
% where mesh.pref_patch = area_frac/Acell (Å^-2).
%
% FFT trick:
%   Define ρ~_ab(p) = exp(+i p·dtau_ab) * ρ_ab(p)  (centered order).
%   Then:
%     Σ~_ab(p) = - mesh.pref_patch * [ V_cont(q)* ρ~ ](p)  (2D convolution)
%   Compute by FFT:
%     Σ~ = -mesh.pref_patch * IFFT( FFT(V) .* FFT(ρ~) )
%   Finally untwist:
%     Σ_ab(p) = exp(-i p·dtau_ab) * Σ~_ab(p)
%
% IMPORTANT about shifts:
%   Our arrays are stored in centered order (q=0 at center).
%   But fft2 expects FFT order (q=0 at (1,1)).
%   Therefore:
%     FFT(ρ~) uses fft2(mesh.to_fft(ρ~))   with mesh.to_fft = ifftshift
%     after ifft2, return to centered by mesh.from_fft = fftshift

Norb = ph.Norb;
Nkx  = mesh.Nkx;
Nky  = mesh.Nky;
Nk   = mesh.Nk;
P    = ph.P;

% --- Pack rho_orb_flat into pages rho_pages (Nky x Nkx x P), centered order ---
rho_pages = zeros(Nky,Nkx,P,'like',rho_orb_flat);
for p = 1:P
  a = ph.a_list(p);
  b = ph.b_list(p);
  rho_pages(:,:,p) = reshape(rho_orb_flat(a,b,:), Nky, Nkx);
end

% Twist: ρ~ = e^{+ip·dtau} ρ
rho_tilde = rho_pages .* ph.phase_plus;  % centered order

% FFT(ρ~) in FFT order
rho_r = fft2( mesh.to_fft(rho_tilde) );  % still FFT order in frequency domain

% Multiply by precomputed Vr = FFT(Vq) (also in FFT order)
Sigma_r = - rho_r .* V.Vr;

% Back to k-space (still FFT order), then return to centered order
Sigma_tilde = mesh.from_fft( ifft2(Sigma_r) );  % centered order

% Apply physical k-measure prefactor: convert V_cont (eV·Å^2) to Σ (eV)
Sigma_tilde = mesh.pref_patch * Sigma_tilde;

% Untwist: Σ = e^{-ip·dtau} Σ~
Sigma_pages = Sigma_tilde .* ph.phase_minus;

% --- Unpack pages back to flat SigmaF_flat ---
SigmaF_flat = zeros(Norb,Norb,Nk,'like',rho_orb_flat);
for p = 1:P
  a = ph.a_list(p);
  b = ph.b_list(p);
  SigmaF_flat(a,b,:) = reshape(Sigma_pages(:,:,p), 1,1,Nk);
end

% Hermitian cleanup per k
for ik = 1:Nk
  S = SigmaF_flat(:,:,ik);
  SigmaF_flat(:,:,ik) = (S + S')/2;
end

end
