clc;
clear;

%% 0) User setup
result_dir = fullfile(pwd, './');
% sigma_file = fullfile(result_dir, 'sigma_aggregate_121U.mat');
sigma_file = fullfile(result_dir, 'sigma_aggregate_121U_qvh.mat');
out_dir = fullfile(result_dir, 'symmetry_plots');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% symmetry setup
sym_group = 'C3v';   % 'C3' or 'C3v'
mirror_opt = 'Mx';   % 'Mx' / 'My' / numeric angle(rad)

% plotting setup
scale_factor = 1e6 / 10;
comp_main = 'xxx';   % 主分析分量（热图）
save_sym_mat = true;
do_blockwise_sym_compare = true; % 若 sigma_blocks_all 存在，则做 block-wise 对称化比较

%% 1) Load aggregated sigma
if ~isfile(sigma_file)
    error('sigma aggregate file not found: %s', sigma_file);
end
S = load(sigma_file);

if ~isfield(S, 'sigma_all') || ~isfield(S, 'U_list') || ~isfield(S, 'Eph_list')
    error('sigma_file missing required vars: sigma_all / U_list / Eph_list');
end

sigma_all = S.sigma_all;                    % 2x2x2xNw x NU
U_list = S.U_list(:).';
Eph_list = S.Eph_list(:);

[~, ~, ~, Nw, NU] = size(sigma_all);
if numel(U_list) ~= NU
    error('U_list length (%d) != sigma_all U dimension (%d)', numel(U_list), NU);
end
if numel(Eph_list) ~= Nw
    error('Eph_list length (%d) != sigma_all omega dimension (%d)', numel(Eph_list), Nw);
end

if isfield(S, 'success_mask') && numel(S.success_mask) == NU
    success_mask = logical(S.success_mask(:).');
else
    success_mask = true(1, NU);
end

if isfield(S, 'block_names_run')
    block_names_run = S.block_names_run;
else
    block_names_run = {'block1','block2','block3','block4'};
end

has_sigma_blocks = isfield(S, 'sigma_blocks_all');

fprintf('Loaded: Nw=%d, NU=%d, validU=%d\n', Nw, NU, nnz(success_mask));

%% 2) Symmetrize sigma(U) one-by-one
sigma_sym_all = nan(size(sigma_all));
rel_err_sum = nan(1, NU);

for iu = 1:NU
    if ~success_mask(iu)
        continue;
    end
    sig_u = sigma_all(:,:,:,:,iu);
    if all(isnan(sig_u), 'all')
        continue;
    end
    [sig_sym_u, info_u] = symmetrize_sigma_rank3(sig_u, sym_group, 'mirror', mirror_opt);
    sigma_sym_all(:,:,:,:,iu) = sig_sym_u;
    rel_err_sum(iu) = info_u.rel_error;
end

%% 3) Optional: block-wise symmetrize then sum
sigma_blocks_sym_all = [];
sigma_from_block_sym = [];
rel_err_blocks = [];

if do_blockwise_sym_compare && has_sigma_blocks
    sigma_blocks_all = S.sigma_blocks_all;  % 2x2x2xNw x nblock x NU
    nblock = size(sigma_blocks_all, 5);

    sigma_blocks_sym_all = nan(size(sigma_blocks_all));
    rel_err_blocks = nan(nblock, NU);

    for iu = 1:NU
        if ~success_mask(iu)
            continue;
        end
        for ib = 1:nblock
            sig_bu = sigma_blocks_all(:,:,:,:,ib,iu);
            if all(isnan(sig_bu), 'all')
                continue;
            end
            [sig_bu_sym, info_bu] = symmetrize_sigma_rank3(sig_bu, sym_group, 'mirror', mirror_opt);
            sigma_blocks_sym_all(:,:,:,:,ib,iu) = sig_bu_sym;
            rel_err_blocks(ib, iu) = info_bu.rel_error;
        end
    end

    sigma_from_block_sym = sum_dim5_omitnan(sigma_blocks_sym_all);
    fprintf('Block-wise symmetrization enabled. nblock=%d\n', nblock);
else
    fprintf('Block-wise symmetrization skipped (sigma_blocks_all not found or disabled).\n');
end

%% 4) Main component heatmaps
comp_main='xxx';
[a, b, c, comp_label] = component_index(comp_main);

comp_before = squeeze(real(sigma_all(a,b,c,:,:))) * scale_factor;      % Nw x NU
comp_after_sum = squeeze(real(sigma_sym_all(a,b,c,:,:))) * scale_factor;

if ~isempty(sigma_from_block_sym)
    comp_after_blk = squeeze(real(sigma_from_block_sym(a,b,c,:,:))) * scale_factor;
end

if ~isempty(sigma_from_block_sym)
    fig = figure('Color','w', 'Position', [100 100 1600 420]);
    tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

    nexttile;
    imagesc(U_list, Eph_list, comp_before);
    colormap(slanCM('RdBu'))%RdGy,RdBu,Blues,heat,thermal,seismic
    colormap(flipud(colormap));
    clim([-3e5,3e5])
    axis xy; colorbar;
    xlabel('U'); ylabel('\hbar\omega (eV)');
    title(sprintf('%s before sym', comp_label));

    nexttile;
    imagesc(U_list, Eph_list, comp_after_sum);
    clim([-2e5,2e5])
    axis xy; colorbar;
    xlabel('U'); ylabel('\hbar\omega (eV)');
    title(sprintf('%s sum-then-sym', comp_label));

    nexttile;
    imagesc(U_list, Eph_list, comp_after_blk);
    axis xy; colorbar;
    xlabel('U'); ylabel('\hbar\omega (eV)');
    title(sprintf('%s block-sym-then-sum', comp_label));
else
    fig = figure('Color','w', 'Position', [100 100 1100 420]);
    tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

    nexttile;
    imagesc(U_list, Eph_list, comp_before);
    axis xy; colorbar;
    xlabel('U'); ylabel('\hbar\omega (eV)');
    title(sprintf('%s before sym', comp_label));
    

    nexttile;
    imagesc(U_list, Eph_list, comp_after_sum);
    axis xy; colorbar;
    xlabel('U'); ylabel('\hbar\omega (eV)');
    title(sprintf('%s sum-then-sym', comp_label));
end
%%
% save("sigma_all_qvh.mat","U_list","Eph_list","sigma_all","sigma_sym_all","sigma_blocks_all","sigma_blocks_sym_all")
% save("sigma_all.mat","U_list","Eph_list","sigma_all","sigma_sym_all","sigma_blocks_all","sigma_blocks_sym_all")
%%
f1 = fullfile(out_dir, sprintf('heatmap_%s_%s_%s.png', comp_label, sym_group, string(mirror_opt)));
exportgraphics(fig, f1, 'Resolution', 300);
close(fig);

%% 5) Line plots for a few U
u_pick = unique([1, round(NU/2), NU]);
u_pick = u_pick(u_pick >= 1 & u_pick <= NU);

for ii = 1:numel(u_pick)
    iu = u_pick(ii);
    if ~success_mask(iu)
        continue;
    end

    fig = figure('Color','w', 'Position', [120 120 1200 520]);
    tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

    % left: all 8 components before/after
    nexttile; hold on;
    comp = ['x','y'];
    for ia = 1:2
        for ib = 1:2
            for ic = 1:2
                y_before = squeeze(real(sigma_all(ia,ib,ic,:,iu))) * scale_factor;
                y_after = squeeze(real(sigma_sym_all(ia,ib,ic,:,iu))) * scale_factor;
                tag = sprintf('%s%s%s', comp(ia), comp(ib), comp(ic));
                plot(Eph_list, y_before, ':', 'LineWidth', 0.9, 'DisplayName', [tag ' before']);
                plot(Eph_list, y_after, '-', 'LineWidth', 1.1, 'DisplayName', [tag ' sym']);
            end
        end
    end
    xlabel('\hbar\omega (eV)'); ylabel('scaled response');
    title(sprintf('U=%.3f : before vs sym', U_list(iu)));
    grid on; box on;

    % right: main component compare
    nexttile; hold on;
    y1 = squeeze(real(sigma_all(a,b,c,:,iu))) * scale_factor;
    y2 = squeeze(real(sigma_sym_all(a,b,c,:,iu))) * scale_factor;
    plot(Eph_list, y1, '--', 'LineWidth', 1.2, 'DisplayName', [comp_label ' before']);
    plot(Eph_list, y2, '-', 'LineWidth', 1.4, 'DisplayName', [comp_label ' sum-sym']);
    if ~isempty(sigma_from_block_sym)
        y3 = squeeze(real(sigma_from_block_sym(a,b,c,:,iu))) * scale_factor;
        plot(Eph_list, y3, '-.', 'LineWidth', 1.3, 'DisplayName', [comp_label ' block-sym']);
    end
    xlabel('\hbar\omega (eV)'); ylabel('scaled response');
    title(sprintf('U=%.3f : %s', U_list(iu), comp_label));
    legend('Location','best');
    grid on; box on;

    f2 = fullfile(out_dir, sprintf('line_U_%.3f_%s_%s_%s.png', U_list(iu), comp_label, sym_group, string(mirror_opt)));
    exportgraphics(fig, f2, 'Resolution', 300);
    close(fig);
end
%%
figure
hold on;
for ia = 1:1
    for ib = 1:1
        for ic = 1:1
            % y_before = squeeze(real(sigma_all(ia,ib,ic,:,50))) * scale_factor;
            y_after = squeeze(real(sigma_sym_all(ia,ib,ic,:,70))) * scale_factor;
            tag = sprintf('%s%s%s', comp(ia), comp(ib), comp(ic));
            plot(Eph_list, y_before, ':', 'LineWidth', 0.9, 'DisplayName', [tag ' before']);
            plot(Eph_list, y_after, '-', 'LineWidth', 1.1, 'DisplayName', [tag ' sym']);
        end
    end
end
%% 6) Symmetry error vs U
fig = figure('Color','w', 'Position', [150 150 1100 420]);
if ~isempty(rel_err_blocks)
    tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');
else
    tiledlayout(1,1, 'Padding','compact', 'TileSpacing','compact');
end

nexttile;
plot(U_list, rel_err_sum, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
xlabel('U'); ylabel('relative projection error');
title(sprintf('Sum tensor symmetry error (%s,%s)', sym_group, string(mirror_opt)));
grid on; box on;

if ~isempty(rel_err_blocks)
    nexttile; hold on;
    nblock = size(rel_err_blocks, 1);
    for ib = 1:nblock
        label = sprintf('block %d', ib);
        if ib <= numel(block_names_run)
            label = string(block_names_run{ib});
        end
        plot(U_list, rel_err_blocks(ib,:), '-', 'LineWidth', 1.1, 'DisplayName', label);
    end
    xlabel('U'); ylabel('relative projection error');
    title(sprintf('Block symmetry error (%s,%s)', sym_group, string(mirror_opt)));
    legend('Location','best');
    grid on; box on;
end

f3 = fullfile(out_dir, sprintf('sym_error_%s_%s.png', sym_group, string(mirror_opt)));
exportgraphics(fig, f3, 'Resolution', 300);
close(fig);

%% 7) Save symmetrized arrays (optional)
if save_sym_mat
    sym_mat_file = fullfile(out_dir, sprintf('sigma_symmetrized_%s_%s.mat', sym_group, string(mirror_opt)));
    if ~isempty(sigma_blocks_sym_all)
        save(sym_mat_file, ...
            'sigma_sym_all', 'sigma_from_block_sym', 'sigma_blocks_sym_all', ...
            'rel_err_sum', 'rel_err_blocks', ...
            'U_list', 'Eph_list', 'success_mask', 'block_names_run', ...
            'sym_group', 'mirror_opt', 'scale_factor', ...
            '-v7.3');
    else
        save(sym_mat_file, ...
            'sigma_sym_all', 'rel_err_sum', ...
            'U_list', 'Eph_list', 'success_mask', ...
            'sym_group', 'mirror_opt', 'scale_factor', ...
            '-v7.3');
    end
    fprintf('Saved symmetrized arrays: %s\n', sym_mat_file);
end

fprintf('\nDone. Figures saved in: %s\n', out_dir);

function Y = sum_dim5_omitnan(X)
% X: (..., dim5, ...), here fixed as dim5=sum dim
    X0 = X;
    X0(isnan(X0)) = 0;
    Y = sum(X0, 5);

    all_nan_mask = all(isnan(X), 5);
    Y(all_nan_mask) = nan;
end

function [a, b, c, tag] = component_index(tag)
    t = lower(char(tag));
    if numel(t) ~= 3
        error('component tag must be 3 chars, e.g. ''xxy''');
    end
    map = @(ch) 1 + (ch == 'y');
    if any(~ismember(t, ['x','y']))
        error('component tag must only contain x/y');
    end
    a = map(t(1));
    b = map(t(2));
    c = map(t(3));
    tag = t;
end

function [sigma_sym, info] = symmetrize_sigma_rank3(sigma_abc, group, varargin)
    p = inputParser;
    addOptional(p, 'mirror', 'Mx');
    parse(p, varargin{:});
    mirrorOpt = p.Results.mirror;

    R0 = eye(2);
    R1 = rot2(2*pi/3);
    R2 = rot2(4*pi/3);
    M = mirror2(mirrorOpt);

    switch lower(group)
        case 'c3'
            G = {R0, R1, R2};
        case 'c3v'
            G = {R0, R1, R2, M, M*R1, M*R2};
        otherwise
            error('Unknown group. Use ''C3'' or ''C3v''.');
    end

    sz = size(sigma_abc);
    if numel(sz) == 3
        sigma_abc = reshape(sigma_abc, 2,2,2,1);
    end
    Nw = size(sigma_abc,4);

    sigma_sym = zeros(size(sigma_abc));
    for iw = 1:Nw
        S = sigma_abc(:,:,:,iw);
        Sacc = zeros(2,2,2);
        for ig = 1:numel(G)
            Sacc = Sacc + apply_op_rank3(G{ig}, S);
        end
        sigma_sym(:,:,:,iw) = Sacc / numel(G);
    end

    diff = sigma_abc - sigma_sym;
    info.norm_sigma = norm(sigma_abc(:));
    info.norm_diff = norm(diff(:));
    info.rel_error = info.norm_diff / max(info.norm_sigma, 1e-30);
end

function S2 = apply_op_rank3(R, S)
    S2 = zeros(2,2,2);
    for a = 1:2
        for b = 1:2
            for c = 1:2
                tmp = 0;
                for ap = 1:2
                    for bp = 1:2
                        for cp = 1:2
                            tmp = tmp + R(a,ap) * R(b,bp) * R(c,cp) * S(ap,bp,cp);
                        end
                    end
                end
                S2(a,b,c) = tmp;
            end
        end
    end
end

function R = rot2(theta)
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
end

function M = mirror2(mirrorOpt)
    if isnumeric(mirrorOpt)
        phi = mirrorOpt;
        Mx = [1 0; 0 -1];
        M = rot2(phi) * Mx * rot2(-phi);
        return;
    end

    switch lower(string(mirrorOpt))
        case "mx"
            M = [1 0; 0 -1];
        case "my"
            M = [-1 0; 0 1];
        otherwise
            error('mirror must be ''Mx'', ''My'', or numeric angle.');
    end
end
