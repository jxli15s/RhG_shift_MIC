function [U, E] = diag_mesh_fromHk(Hk, opts)
%DIAG_MESH_FROMHK
%==========================================================================
% 在 2D k 网格上对哈密顿量逐点对角化（parfor-safe）。
%
% 输入:
%   Hk: nb x nb x Nkx x Nky, 每个 k 点的哈密顿量
%   opts.band_list (可选): 指定保留的能带索引（按升序本征值排序后）
%   opts.useParfor (可选): true/false, 是否在该函数内部使用 parfor (default=true)
%
% 输出:
%   U: nb x nb_sel x Nkx x Nky
%   E: Nkx x Nky x nb_sel
%==========================================================================

    arguments
        Hk {mustBeNumeric}
        opts struct = struct()
    end

    [nb, nb2, Nkx, Nky] = size(Hk);
    if nb ~= nb2
        error('Hk must be nb x nb x Nkx x Nky.');
    end

    if isfield(opts, 'band_list') && ~isempty(opts.band_list)
        band_list = opts.band_list(:).';
    else
        band_list = 1:nb;
    end
    if isfield(opts, 'useParfor') && ~isempty(opts.useParfor)
        useParfor = logical(opts.useParfor);
    else
        useParfor = true;
    end
    nb_sel = numel(band_list);

    itotal = Nkx * Nky;
    Utmp = zeros(nb, nb_sel, itotal, 'like', Hk);
    Etmp = zeros(itotal, nb_sel);

    if useParfor
        parfor ll = 1:itotal
            [ix, iy] = ind2sub([Nkx, Nky], ll);
            H = Hk(:,:,ix,iy);
            H = (H + H') / 2;

            [V, D] = eig(H);
            evals = real(diag(D));
            [evals, ind] = sort(evals, 'ascend');
            V = V(:, ind);

            Etmp(ll,:) = evals(band_list);
            Utmp(:,:,ll) = V(:, band_list);
        end
    else
        for ll = 1:itotal
            [ix, iy] = ind2sub([Nkx, Nky], ll);
            H = Hk(:,:,ix,iy);
            H = (H + H') / 2;

            [V, D] = eig(H);
            evals = real(diag(D));
            [evals, ind] = sort(evals, 'ascend');
            V = V(:, ind);

            Etmp(ll,:) = evals(band_list);
            Utmp(:,:,ll) = V(:, band_list);
        end
    end

    U = reshape(Utmp, [nb, nb_sel, Nkx, Nky]);
    E = reshape(Etmp, [Nkx, Nky, nb_sel]);
end
