function [Unk, Enk, vxk, vyk, vzk] = get_bulk_plane_ham_velocity( ...
        Ham, Ham_kx, Ham_ky, Ham_kz)
% GET_BULK_PLANE_FROM_PRECOMPUTED_HAM
% Inputs:
%   Ham    : [nband, nband, nkx, nky]
%   Ham_kx : [nband, nband, nkx, nky]  dH/dkx  (same basis as Ham)
%   Ham_ky : [nband, nband, nkx, nky]  dH/dky
%   Ham_kz : [nband, nband, nkx, nky]  dH/dkz (optional, pass [] if not used)
%
% Outputs:
%   Unk : [nband, nband, nkx, nky]   eigenvectors (columns are bands)
%   Enk : [nkx, nky, nband]         eigenvalues (sorted ascending)
%   vxk : [nband, nband, nkx, nky]  band-basis v_x = U^† (dH/dkx) U
%   vyk : [nband, nband, nkx, nky]  band-basis v_y
%   vzk : [nband, nband, nkx, nky]  band-basis v_z (empty if Ham_kz not provided)

    % ---- sanity ----
    if ndims(Ham) ~= 4
        error('Ham must be 4D: [nband, nband, nkx, nky].');
    end
    [nb1,~, nkx, nky] = size(Ham);
    nband = nb1;
    itotal = nkx * nky;

    % ---- flatten to help parfor ----
    % store pages in linear index ll = sub2ind([nkx nky], i, j)
    unktem = zeros(nband, nband, itotal, 'like', Ham);
    enktem = zeros(itotal, nband, 'like', real(Ham(1)));

    vxtem  = zeros(nband, nband, itotal, 'like', Ham);
    vytem  = zeros(nband, nband, itotal, 'like', Ham);
    vztem = zeros(nband, nband, itotal, 'like', Ham);

    parfor ll = 1:itotal
        [i, j] = ind2sub([nkx, nky], ll);

        hk  = Ham(:,:,i,j);
        hkx = Ham_kx(:,:,i,j);
        hky = Ham_ky(:,:,i,j);
        hkz = Ham_kz(:,:,i,j);

        % enforce Hermitian (numerical safety)
        hk  = (hk  + hk')  / 2;
        hkx = (hkx + hkx') / 2;
        hky = (hky + hky') / 2;
        hkz = (hkz + hkz') / 2;

        % diagonalize
        [V, D] = eig(hk);
        [E, ind] = sort(real(diag(D)));   % energies should be real
        U = V(:, ind);

        % band-basis velocities
        Udag = U';
        vx_band = Udag * hkx * U;
        vy_band = Udag * hky * U;
        vz_band = Udag * hkz * U;

        enktem(ll,:)   = E;
        unktem(:,:,ll) = U;
        vxtem(:,:,ll)  = vx_band;
        vytem(:,:,ll)  = vy_band;
        vztem(:,:,ll) = vz_band;
    end

    % ---- reshape back ----
    Unk = reshape(unktem, [nband, nband, nkx, nky]);
    Enk = reshape(enktem, [nkx, nky, nband]);

    vxk = reshape(vxtem,  [nband, nband, nkx, nky]);
    vyk = reshape(vytem,  [nband, nband, nkx, nky]);
    vzk = reshape(vztem,  [nband, nband, nkx, nky]);
end
