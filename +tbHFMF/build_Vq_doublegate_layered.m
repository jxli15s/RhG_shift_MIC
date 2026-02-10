function V = build_Vq_doublegate_layered(g, mesh, pars)
%BUILD_VQ_DOUBLEGATE_LAYERED  Double-gate screened Coulomb kernel V_cont(q; z_a,z_b).
%
% Geometry:
%   Metallic gates at z = ±d_gate (Å).
%   Orbitals at z_a (Å), z_b (Å), |z| < d_gate.
%
% Kernel (one common closed form):
%   V_cont(q;z,z') = (2π e^2)/(ε q) * [ 2 sinh(q(d-z_>)) sinh(q(d+z_<)) ] / sinh(2qd)
% where z_> = max(z,z'), z_< = min(z,z').
%
% Units:
%   e2_over_eps : eV·Å
%   1/q         : Å
%   => V_cont   : eV·Å^2
%
% IMPORTANT:
%   Do NOT divide by Acell here. The k-measure prefactor mesh.pref_patch (Å^-2)
%   is applied later in Fock/Hartree, giving Σ in eV.

if ~isfield(pars,'q_small'),     pars.q_small = 1e-6; end
if ~isfield(pars,'include_q2'),  pars.include_q2 = false; end
if ~isfield(pars,'z0_mode'),     pars.z0_mode = 'center'; end

d = pars.d_gate;
q = mesh.qabs;                      % Nky x Nkx, centered order
Norb = size(g.wpos,1);
Nky  = mesh.Nky; Nkx = mesh.Nkx;

% --- set z origin (recommended: center of slab orbitals) ---
z = g.wpos(:,3);
if ischar(pars.z0_mode) && strcmpi(pars.z0_mode,'center')
  z = z - 0.5*(max(z)+min(z));
elseif isnumeric(pars.z0_mode)
  z = z - pars.z0_mode;
end
if any(abs(z) >= d)
  error('build_Vq_doublegate_layered: require |z|<d_gate.');
end

% --- V0(q->0) analytic limit (still eV·Å^2) ---
% Expand sinh(x)~x, sinh(2qd)~2qd:
%   V0 = (2π e^2/ε) * ((d-z_>)(d+z_<)/d)
V0 = zeros(Norb,Norb);
for a = 1:Norb
  for b = 1:Norb
    zgt = max(z(a), z(b));
    zlt = min(z(a), z(b));
    A = d - zgt;
    B = d + zlt;
    V0(a,b) = (2*pi*pars.e2_over_eps) * (A*B/d);
  end
end
V.V0 = V0;

% --- build pages Vq(:,:,p) for each (a,b) ---
[a_list,b_list] = ndgrid(1:Norb,1:Norb);
a_list = a_list(:); b_list = b_list(:);
P = Norb*Norb;

den   = sinh(2*q*d);
small = (q < pars.q_small);

Vq_pages = zeros(Nky,Nkx,P,'like',1+1i);

parfor p = 1:P
  a = a_list(p); b = b_list(p);
  zgt = max(z(a), z(b));
  zlt = min(z(a), z(b));
  A = d - zgt;
  B = d + zlt;

  num  = 2*sinh(q*A).*sinh(q*B);

  % regularize 1/q at q~0
  qreg = max(q, pars.q_small);
  Vq = (2*pi*pars.e2_over_eps) .* (num ./ den) ./ qreg;

  % small-q replacement
  Vq0 = V0(a,b);
  if pars.include_q2
    % optional q^2 correction (if you want)
    C2 = (A^2 + B^2)/6 - (2*d^2)/3;
    Vq(small) = Vq0 * (1 + (q(small).^2)*C2);
  else
    Vq(small) = Vq0;
  end

  Vq_pages(:,:,p) = Vq;
end

% --- FFT kernel in FFT order ---
% We store Vq_pages in centered order (q=0 at center). Before fft2:
%   centered -> FFT order via mesh.to_fft = ifftshift.
V.Vr = fft2( mesh.to_fft(Vq_pages) );

end
