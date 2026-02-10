function hk = build_hk_atomgauge(g, mesh)
%BUILD_HK_ATOMGAUGE  Bloch Hamiltonian in atom gauge:
%   h_{αβ}(k) = Σ_R t_{αβ}(R) * exp(+i k·R)
%
% Sign convention:
%   We use exp(+i k·R). Then in Fock we must be consistent with phases in dtau,
%   but the convolution itself is independent of this sign as long as you keep
%   everything consistent.
%
% Output:
%   hk: Norb x Norb x Nky x Nkx  (centered order in k around K0)

Norb = size(g.ham,1);
Nh   = size(g.ham,3);
Nky  = mesh.Nky;
Nkx  = mesh.Nkx;

% Convert hopping direction (fractional) to Cartesian displacement (Å)
Rcart = g.hopr * g.a;   % Nh x 3

Kx = mesh.Kx; Ky = mesh.Ky; Kz = mesh.Kz;

hk = zeros(Norb,Norb,Nky,Nkx,'like',1+1i);

for l = 1:Nh
  phase = exp( 1i*(Kx*Rcart(l,1) + Ky*Rcart(l,2) + Kz*Rcart(l,3)) ); % Nky x Nkx
  hk = hk + g.ham(:,:,l) .* reshape(phase,1,1,Nky,Nkx);
end

end
