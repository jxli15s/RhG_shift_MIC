function [eta_abc, out] = mic_metric_plane_fromUE_mn( ...
    Kx, Ky, U, E, Eph_list, Ef, kT, eta, opts)
%MIC_METRIC_PLANE_FROMUE_MN
%==========================================================================
% Magnetic Injection Current (MIC, metric-type, linear polarization)
% 多带形式：sum_{m!=n}，在 2D k 网格上积分。
%
%   eta^{abc}(Eph) = -(pi * g_s * e^2)/(2*hbar) *
%       ∫ d^2k/(2*pi)^2 sum_{m!=n}
%       f_nm * Delta_v^a_mn * [2 g^{bc}_mn] *
%       delta_eV((E_m-E_n)-Eph)
%
% 其中:
%   2 g^{bc}_{mn} = r^b_{mn} r^c_{nm} + r^c_{mn} r^b_{nm}
%   r^b_{mn} = i <u_m | d_{k_b} u_n>
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
    doGaugeFix = tbNLC.get_opt(opts, 'doGaugeFix', true);
    g_s = tbNLC.get_opt(opts, 'g_s', 1);
    verbose = tbNLC.get_opt(opts, 'verbose', true);
    positiveDE = tbNLC.get_opt(opts, 'positiveDE', true);
    saveFullMN = tbNLC.get_opt(opts, 'saveFullMN', true);

    % -------------------- constants (SI) --------------------
    e_charge = 1.602176634e-19;
    hbar_Js = 1.054571817e-34;
    pref = -pi * g_s * (e_charge^2) / (2*hbar_Js);
    vel_pref = e_charge / hbar_Js;

    % -------------------- reshape & checks --------------------
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
        fprintf('[tbNLC.mic] Nkx=%d Nky=%d nb=%d nb_sel=%d\n', Nkx, Nky, nb, nb_sel);
        fprintf('[tbNLC.mic] dk1=(%.3e,%.3e), dk2=(%.3e,%.3e), det(J)=%.3e\n', ...
            ops.dk1(1), ops.dk1(2), ops.dk2(1), ops.dk2(2), ops.detJ);
        fprintf('[tbNLC.mic] periodicFD=%d, trimBoundary=%d, doGaugeFix=%d, g_s=%g, positiveDE=%d\n', ...
            periodicFD, trimBoundary, doGaugeFix, g_s, positiveDE);
    end

    % -------------------- gauge smoothing --------------------
    if doGaugeFix
        Ug = tbNLC.gauge_fix_parallel_transport(U);
    else
        Ug = U;
    end

    mask_k4 = ops.mask_k4;
    ix_list = ops.ix_list;
    iy_list = ops.iy_list;

    % -------------------- eigenvector derivatives --------------------
    [du_x, du_y] = ops.fd_du(Ug);

    % -------------------- r_mn(k) --------------------
    r_mn = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 'like', Ug);

    for ix = ix_list
        for iy = iy_list
            U0 = squeeze(Ug(:,:,ix,iy));
            dux = squeeze(du_x(:,:,ix,iy));
            duy = squeeze(du_y(:,:,ix,iy));
            r_mn(ix,iy,:,:,1) = 1i * (U0' * dux);
            r_mn(ix,iy,:,:,2) = 1i * (U0' * duy);
        end
    end

    off = ones(nb_sel) - eye(nb_sel);
    mask_off = reshape(off, [1,1,nb_sel,nb_sel,1]);

    % MIC 保持与原脚本一致：先厄米化再去对角
    r_mn_H = tbNLC.hermitianize_rmn(r_mn);
    r_mn_H = r_mn_H .* mask_off;

    r_nm_mn = permute(r_mn_H, [1 2 4 3 5]);

    % -------------------- metric kernel --------------------
    g2 = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 2, 'like', Ug);
    for b = 1:2
        for cc = 1:2
            g2(:,:,:,:,b,cc) = real( ...
                r_mn_H(:,:,:,:,b) .* r_nm_mn(:,:,:,:,cc) + ...
                r_mn_H(:,:,:,:,cc) .* r_nm_mn(:,:,:,:,b));
        end
    end

    % -------------------- energies: dE_mn --------------------
    Em = reshape(E, [Nkx, Nky, nb_sel, 1]);
    En = reshape(E, [Nkx, Nky, 1, nb_sel]);
    dE_mn = Em - En;

    if positiveDE
        mask_pos = double(dE_mn > 0);
    else
        mask_pos = ones(size(dE_mn));
    end

    % -------------------- fermi factor f_nm --------------------
    f = tbNLC.fermi_dirac(E, Ef, kT);
    fm = reshape(f, [Nkx, Nky, nb_sel, 1]);
    fn = reshape(f, [Nkx, Nky, 1, nb_sel]);
    f_nm = fn - fm;

    % -------------------- band velocities --------------------
    [dE_dx, dE_dy] = ops.fd_tensor(E);
    v_band = zeros(Nkx, Nky, nb_sel, 2);
    v_band(:,:,:,1) = vel_pref * dE_dx;
    v_band(:,:,:,2) = vel_pref * dE_dy;

    dv_mn = zeros(Nkx, Nky, nb_sel, nb_sel, 2);
    for a = 1:2
        Vm = reshape(v_band(:,:,:,a), [Nkx, Nky, nb_sel, 1]);
        Vn = reshape(v_band(:,:,:,a), [Nkx, Nky, 1, nb_sel]);
        dv_mn(:,:,:,:,a) = Vm - Vn;
    end

    % -------------------- total mask --------------------
    mask_total = mask_k4 .* reshape(off, [1,1,nb_sel,nb_sel]) .* mask_pos;

    % -------------------- integrate over photon energy --------------------
    eta_abc = zeros(2,2,2,Nw);

    for a = 1:2
        dvA = dv_mn(:,:,:,:,a);
        for b = 1:2
            for cc = 1:2
                gBC = g2(:,:,:,:,b,cc);
                Kabc = (f_nm .* dvA .* gBC) .* mask_total;
                for iw = 1:Nw
                    Eph = Eph_list(iw);
                    delta_w = (1/pi) * eta ./ ((dE_mn - Eph).^2 + eta^2);
                    eta_abc(a,b,cc,iw) = pref * ops.w_k * sum(Kabc .* delta_w, 'all');
                end
            end
        end
    end

    % -------------------- outputs --------------------
    out = struct();
    out.pref = pref;
    out.vel_pref = vel_pref;
    out.w_k = ops.w_k;
    out.J = ops.J;
    out.invJT = ops.invJT;
    out.dk1 = ops.dk1;
    out.dk2 = ops.dk2;
    out.mask_k = ops.mask_k;
    out.mask_off = off;
    out.positiveDE = positiveDE;
    out.Ug = Ug;
    out.band_list = band_list;

    if saveFullMN
        out.r_mn = r_mn_H;
        out.g2 = g2;
        out.dE_mn = dE_mn;
        out.f_nm = f_nm;
        out.v_band = v_band;
        out.dv_mn = dv_mn;
    end
end
