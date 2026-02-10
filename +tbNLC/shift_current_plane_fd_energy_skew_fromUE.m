function [sigma_abc, out] = shift_current_plane_fd_energy_skew_fromUE( ...
    Kx, Ky, U, E, Eph_list, Ef, kT, eta, opts)
%SHIFT_CURRENT_PLANE_FD_ENERGY_SKEW_FROMUE
%==========================================================================
% 2D Shift Current 计算（输入为本征矢 U 与本征值 E）：
%
%   sigma^{abc}(Omega) = (pi * g_s * e^2 / hbar) *
%       ∫(d^2k/(2pi)^2) sum_{n!=m} f_nm *
%       Im[r^b_{mn} * r^c_{nm;a}] * delta_eV((E_m-E_n)-Omega)
%
% 其中
%   r^b_{mn}   = i <u_m | d_{k_b} u_n>
%   A^a_nn     = i <u_n | d_{k_a} u_n>
%   r^c_{mn;a} = d_{k_a} r^c_{mn} - i(A^a_mm - A^a_nn) r^c_{mn}
%
% 该函数将主脚本中的旧版局部实现迁移到 +tbNLC 包，接口保持兼容。
%==========================================================================

    arguments
        Kx double
        Ky double
        U  {mustBeNumeric}
        E  double
        Eph_list double
        Ef double
        kT double
        eta double
        opts struct = struct()
    end

    % -------------------- options --------------------
    periodicFD = tbNLC.get_opt(opts, 'periodicFD', false);
    trimBoundary = tbNLC.get_opt(opts, 'trimBoundary', true);
    symBC = tbNLC.get_opt(opts, 'symBC', true);
    verbose = tbNLC.get_opt(opts, 'verbose', true);
    g_s = tbNLC.get_opt(opts, 'g_s', 1);
    doGaugeFix = tbNLC.get_opt(opts, 'doGaugeFix', true);

    % -------------------- constants --------------------
    e_charge = 1.602176634e-19;
    hbar_Js = 1.054571817e-34;
    pref = pi * g_s * (e_charge^2) / hbar_Js;

    % -------------------- sizes & reshape --------------------
    [Nkx, Nky] = size(Kx);
    if ~isequal(size(Ky), [Nkx, Nky])
        error('Kx and Ky must have the same size Nkx x Nky.');
    end

    szU = size(U);
    if numel(szU) == 3
        nb = szU(1);
        nb_sel = szU(2);
        if szU(3) ~= Nkx*Nky
            error('If U is 3D, its 3rd dim must be Nkx*Nky.');
        end
        U = reshape(U, [nb, nb_sel, Nkx, Nky]);
    elseif numel(szU) == 4
        nb = szU(1);
        nb_sel = szU(2);
        if szU(3) ~= Nkx || szU(4) ~= Nky
            error('U must be nb x nb_sel x Nkx x Nky.');
        end
    else
        error('U must be 4D or 3D.');
    end

    if ~isequal(size(E), [Nkx, Nky, nb_sel])
        if isequal(size(E), [Nkx*Nky, nb_sel])
            E = reshape(E, [Nkx, Nky, nb_sel]);
        else
            error('E must be Nkx x Nky x nb_sel (or (Nkx*Nky) x nb_sel).');
        end
    end

    % 可选: 在响应计算阶段进一步筛选参与求和的能带
    band_list = tbNLC.get_opt(opts, 'band_list', []);
    if ~isempty(band_list)
        band_list = band_list(:).';
        if any(mod(band_list,1) ~= 0)
            error('opts.band_list must contain integer indices.');
        end
        if any(band_list < 1) || any(band_list > nb_sel)
            error('opts.band_list index out of range. Current nb_sel=%d.', nb_sel);
        end
        U = U(:, band_list, :, :);
        E = E(:, :, band_list);
        nb_sel = numel(band_list);
    end

    Eph_list = Eph_list(:).';
    Nw = numel(Eph_list);

    % -------------------- shared skew-grid ops --------------------
    ops = tbNLC.skewgrid_ops_fromK(Kx, Ky, periodicFD, trimBoundary);

    if verbose
        fprintf('[tbNLC.shift] Nkx=%d Nky=%d nb=%d nb_sel=%d\n', Nkx, Nky, nb, nb_sel);
        fprintf('[tbNLC.shift] dk1=(%.6e, %.6e), dk2=(%.6e, %.6e)\n', ...
            ops.dk1(1), ops.dk1(2), ops.dk2(1), ops.dk2(2));
        fprintf('[tbNLC.shift] det(J)=%.6e\n', ops.detJ);
        fprintf('[tbNLC.shift] periodicFD=%d, trimBoundary=%d, symBC=%d, g_s=%g, doGaugeFix=%d\n', ...
            periodicFD, trimBoundary, symBC, g_s, doGaugeFix);
    end

    % -------------------- gauge smoothing --------------------
    if doGaugeFix
        Ug = tbNLC.gauge_fix_parallel_transport(U);
    else
        Ug = U;
    end

    ix_list = ops.ix_list;
    iy_list = ops.iy_list;

    % =====================================================================
    % 1) eigenvector derivatives in Cartesian directions
    % =====================================================================
    [du_x, du_y] = ops.fd_du(Ug);

    % =====================================================================
    % 2) build A_diag and r_mn
    % =====================================================================
    A_diag = zeros(Nkx, Nky, nb_sel, 2);
    r_mn = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 'like', Ug);

    for ix = ix_list
        for iy = iy_list
            U0 = squeeze(Ug(:,:,ix,iy));
            dux = squeeze(du_x(:,:,ix,iy));
            duy = squeeze(du_y(:,:,ix,iy));

            Ax = 1i * diag(U0' * dux);
            Ay = 1i * diag(U0' * duy);
            A_diag(ix,iy,:,1) = real(Ax);
            A_diag(ix,iy,:,2) = real(Ay);

            r_mn(ix,iy,:,:,1) = 1i * (U0' * dux);
            r_mn(ix,iy,:,:,2) = 1i * (U0' * duy);
        end
    end

    off = ones(nb_sel) - eye(nb_sel);
    r_mn = r_mn .* reshape(off, [1,1,nb_sel,nb_sel,1]);

    % 与 MIC 保持一致的厄米投影策略
    r_mn = tbNLC.hermitianize_rmn(r_mn);
    r_mn = r_mn .* reshape(off, [1,1,nb_sel,nb_sel,1]);

    % =====================================================================
    % 3) derivatives of r_mn on skew grid
    % =====================================================================
    dr = ops.fd_rmn(r_mn);

    % =====================================================================
    % 4) covariant derivative r_cov = dr - i(Amm-Ann)*r
    % =====================================================================
    r_cov = zeros(size(dr), 'like', dr);
    for a = 1:2
        Amm = reshape(A_diag(:,:,:,a), [Nkx, Nky, nb_sel, 1]);
        Ann = reshape(A_diag(:,:,:,a), [Nkx, Nky, 1, nb_sel]);
        dA = Amm - Ann;
        for c = 1:2
            r_cov(:,:,:,:,c,a) = dr(:,:,:,:,c,a) - 1i * dA .* r_mn(:,:,:,:,c);
        end
    end

    % =====================================================================
    % 5) fermi occupations
    % =====================================================================
    f_n = tbNLC.fermi_dirac(E, Ef, kT);

    % =====================================================================
    % 6) main integration
    % =====================================================================
    sigma_abc = zeros(2,2,2,Nw);

    for ix = ix_list
        for iy = iy_list
            Ek = squeeze(E(ix,iy,:));
            fk = squeeze(f_n(ix,iy,:));

            for n = 1:nb_sel
                for m = 1:nb_sel
                    if m == n
                        continue;
                    end

                    dE = Ek(m) - Ek(n);
                    f_nm = fk(n) - fk(m);
                    if abs(f_nm) < 1e-14
                        continue;
                    end

                    delta_w = (1/pi) * eta ./ ((dE - Eph_list).^2 + eta^2);

                    for a = 1:2
                        for b = 1:2
                            for c = 1:2
                                I1 = imag(r_mn(ix,iy,m,n,b) * r_cov(ix,iy,n,m,c,a));
                                if symBC
                                    I2 = imag(r_mn(ix,iy,m,n,c) * r_cov(ix,iy,n,m,b,a));
                                    I = 0.5 * (I1 + I2);
                                else
                                    I = I1;
                                end
                                sigma_abc(a,b,c,:) = sigma_abc(a,b,c,:) + ...
                                    reshape(pref * f_nm * (I * ops.w_k) .* delta_w, [1,1,1,Nw]);
                            end
                        end
                    end
                end
            end
        end
    end

    % -------------------- outputs --------------------
    out = struct();
    out.E = E;
    out.Ug = Ug;
    out.A_diag = A_diag;
    out.r_mn = r_mn;
    out.dr = dr;
    out.r_cov = r_cov;
    out.dk1 = ops.dk1;
    out.dk2 = ops.dk2;
    out.J = ops.J;
    out.invJT = ops.invJT;
    out.w_k = ops.w_k;
    out.pref = pref;
    out.ix_list = ix_list;
    out.iy_list = iy_list;
    out.band_list = band_list;
end
