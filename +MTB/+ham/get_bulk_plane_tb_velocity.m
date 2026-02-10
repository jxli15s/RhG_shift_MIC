function [Unk, Enk, vxk, vyk, vzk] = get_bulk_plane_tb_velocity(obj, Kx, Ky, Kz)
% GET_BULK_PLANE_TB_VELOCITY
%   Solve TB Hamiltonian on a 2D k-mesh and return eigenstates/energies and
%   band-basis velocity matrices.
%
% Outputs:
%   Unk : [dim_H, dim_H, knum, knum]   eigenvectors (columns are bands)
%   Enk : [knum, knum, dim_H]         eigenvalues (sorted ascending)
%   vxk : [dim_H, dim_H, knum, knum]  band-basis v_x = U^† (dH/dkx) U
%   vyk : [dim_H, dim_H, knum, knum]  band-basis v_y
%   vzk : [dim_H, dim_H, knum, knum]  band-basis v_z
%
% Notes:
%   - Here v_x etc. are computed from dH/dk (units: eV*Å). If you want
%     physical velocity (Å/s), divide by ħ (in eV*s): v_phys = (1/ħ) dH/dk.

    knum   = size(Kx, 1);
    itotal = knum^2;

    dim_H = size(obj.ham, 1);
    L     = size(obj.hopr, 1);

    % Flatten hopping matrices: [dim_H*dim_H, L]
    ham_flat = reshape(obj.ham, dim_H*dim_H, L);

    % Real-space displacement vectors for each hopping, in Cartesian (Å): [L,3]
    % (Your original code used (obj.hopr*obj.a)*k', so keep consistent.)
    Rcart = obj.hopr * obj.a;  % [L,3]

    % Allocate temporaries for parfor
    unktem = zeros(dim_H, dim_H, itotal);
    enktem = zeros(itotal, dim_H);

    vxtem  = zeros(dim_H, dim_H, itotal);
    vytem  = zeros(dim_H, dim_H, itotal);
    vztem  = zeros(dim_H, dim_H, itotal);

    parfor ll = 1:itotal
        [i, j] = ind2sub([knum, knum], ll);
        k = [Kx(i,j), Ky(i,j), Kz(i,j)];

        [Etem, Psik, vx_band, vy_band, vz_band] = solve_one_k_tb( ...
            k, dim_H, ham_flat, Rcart);

        enktem(ll,:)   = Etem;
        unktem(:,:,ll) = Psik;

        vxtem(:,:,ll)  = vx_band;
        vytem(:,:,ll)  = vy_band;
        vztem(:,:,ll)  = vz_band;
    end

    Unk = reshape(unktem, [dim_H, dim_H, knum, knum]);
    Enk = reshape(enktem, [knum, knum, dim_H]);

    vxk = reshape(vxtem, [dim_H, dim_H, knum, knum]);
    vyk = reshape(vytem, [dim_H, dim_H, knum, knum]);
    vzk = reshape(vztem, [dim_H, dim_H, knum, knum]);
end


function [E, Psik, vx_band, vy_band, vz_band] = solve_one_k_tb(k, dim_H, ham_flat, Rcart)
% Solve eigenproblem and compute band-basis velocity matrices at one k.

    % ----- phase factors e^{i k·R} -----
    phase = exp(1j * (Rcart * k(:)));  % [L,1]

    % ----- H(k) -----
    hk = ham_flat * phase;            % [dim_H^2, 1]
    hk = reshape(hk, dim_H, dim_H);
    hk = (hk + hk') / 2;              % enforce Hermitian

    % ----- dH/dk = sum_R t(R) * i R_alpha e^{ik·R} -----
    iphase = 1j * phase;              % i e^{ik·R}

    hkx = ham_flat * (iphase .* Rcart(:,1));  % [dim_H^2,1]
    hky = ham_flat * (iphase .* Rcart(:,2));
    hkz = ham_flat * (iphase .* Rcart(:,3));

    hkx = reshape(hkx, dim_H, dim_H); hkx = (hkx + hkx')/2;
    hky = reshape(hky, dim_H, dim_H); hky = (hky + hky')/2;
    hkz = reshape(hkz, dim_H, dim_H); hkz = (hkz + hkz')/2;

    % ----- diagonalize -----
    [V, D] = eig(hk);
    [E, ind] = sort(real(diag(D)));   % eigenvalues should be real
    Psik = V(:, ind);

    % ----- band-basis velocities: U^† (dH/dk) U -----
    Udag = Psik';
    vx_band = Udag * hkx * Psik;
    vy_band = Udag * hky * Psik;
    vz_band = Udag * hkz * Psik;
end
