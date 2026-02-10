function mesh = build_kmesh_patch(g, K0_frac, L_frac, Nk)
%BUILD_KMESH_PATCH  Uniform 2D patch mesh around K0 (in reduced coords).
%
% === What this mesh means ===
% We sample a small patch around a valley center K0 (e.g. (1/3,2/3)).
% Define local momentum p = k - K0. FFT convolution is done on this local p-grid.
%
% === Centered order (display/physics-friendly) ===
% We store arrays as centered order:
%   p=0 (equivalently k=K0) is at the CENTER index of the 2D array.
% This is convenient for:
%   - small-q regularization (q->0 at the center)
%   - symmetry checks and plotting
%
% === FFT order (fft2-friendly) ===
% MATLAB fft2 assumes the "zero lag / DC component" corresponds to index (1,1).
% Therefore before fft2 we must convert centered -> FFT order:
%   to_fft   = ifftshift  (move center to (1,1))
% and after ifft2 convert back FFT -> centered:
%   from_fft = fftshift   (move (1,1) back to center)
%
% We provide:
%   mesh.to_fft(A)   : centered -> FFT order  (ifftshift on dim 1&2)
%   mesh.from_fft(A) : FFT -> centered        (fftshift on dim 1&2)
%
% === Physical normalization ===
% In 2D continuum:  ∫ d^2k/(2π)^2 ...
% On a patch with uniform mesh:
%   ∫_{patch} d^2k/(2π)^2 (...)  ≈  pref_patch * mean_k(...)
% where
%   pref_patch = Apatch / (2π)^2 = area_frac / Acell   (units: Å^-2)
% This is CRITICAL because V_cont(q) from 2D FT has units eV·Å^2.
% Then pref_patch * V_cont -> eV.

Nkx = Nk(1);
Nky = Nk(2);
Lx  = L_frac(1);
Ly  = L_frac(2);

% --- Real-space unit-cell area Acell (Å^2) ---
a1 = g.a(1,:);
a2 = g.a(2,:);
Acell = norm(cross(a1,a2));   % Å^2
mesh.Acell = Acell;

% --- BZ area (Å^-2), ABZ = |b1×b2| = (2π)^2/Acell ---
b1 = g.b(1,:);
b2 = g.b(2,:);
ABZ = norm(cross(b1,b2));     % Å^-2
mesh.ABZ = ABZ;

% --- Patch fraction and physical prefactor ---
mesh.area_frac  = Lx * Ly;                % dimensionless fraction in reduced coords
mesh.Apatch     = mesh.area_frac * ABZ;   % Å^-2
mesh.pref_patch = mesh.area_frac / Acell; % Å^-2 = Apatch/(2π)^2

mesh.Nkx = Nkx;
mesh.Nky = Nky;
mesh.Nk  = Nkx * Nky;

mesh.K0_frac = K0_frac;
mesh.L_frac  = L_frac;

% --- Centered grid indices: ensures a "center index" representing p≈0 ---
ix = (0:Nkx-1) - floor(Nkx/2);
iy = (0:Nky-1) - floor(Nky/2);

px_frac = ix * (Lx / Nkx);
py_frac = iy * (Ly / Nky);
[PXf, PYf] = meshgrid(px_frac, py_frac);   % Nky x Nkx (MATLAB meshgrid convention)

KXf = K0_frac(1) + PXf;
KYf = K0_frac(2) + PYf;

% Reduced -> Cartesian using row-basis reciprocal vectors (rows of g.b)
Kx = KXf.*g.b(1,1) + KYf.*g.b(2,1);
Ky = KXf.*g.b(1,2) + KYf.*g.b(2,2);
Kz = KXf.*g.b(1,3) + KYf.*g.b(2,3);

% Local p in Cartesian (Å^-1): p = (PXf,PYf) mapped by g.b
Px = PXf.*g.b(1,1) + PYf.*g.b(2,1);
Py = PXf.*g.b(1,2) + PYf.*g.b(2,2);
Pz = PXf.*g.b(1,3) + PYf.*g.b(2,3);

mesh.PXf = PXf; mesh.PYf = PYf;
mesh.KXf = KXf; mesh.KYf = KYf;

mesh.Kx = Kx; mesh.Ky = Ky; mesh.Kz = Kz;
mesh.Px = Px; mesh.Py = Py; mesh.Pz = Pz;

mesh.qabs = sqrt(Px.^2 + Py.^2 + Pz.^2);  % |q| on the patch (Å^-1), centered order

% === FFT ordering helpers ===
% centered -> FFT order: q=0 at center -> q=0 at (1,1)
mesh.to_fft   = @(A) ifftshift(ifftshift(A,1),2);
% FFT order -> centered
mesh.from_fft = @(A) fftshift(fftshift(A,1),2);

end
