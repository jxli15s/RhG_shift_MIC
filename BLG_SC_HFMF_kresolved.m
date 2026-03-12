clc;
clear;

%% 0) User setup (fixed U)
U_fix = 0.0;
ne_tag = "0.0000e12";

% HFMF 数据目录（按需修改）
hfmf_data_dir = "/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data";
hfmf_file = fullfile(hfmf_data_dir, sprintf("ne=%s_U=%.3fdata.mat", ne_tag, U_fix));


% 计算设置
Eph_list = linspace(0.0, 0.1, 2000);
kT = 0.0;
eta = 1e-4;


% block 选择（默认四个全算）
block_names_run = {'K_up', 'Kp_up', 'K_dn', 'Kp_dn'};


%% 1) Output paths (single folder)
result_dir = fullfile(pwd, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

u_tag_mat = sprintf('U_%.3f', U_fix);
out_file = fullfile(result_dir, [u_tag_mat, '.mat']);

u_tag_fig = strrep(sprintf('U_%3.0f', U_fix), ' ', '0');
fig_png = fullfile(result_dir, [u_tag_fig, '.png']);

%% 2) Parallel pool (Threads)
p = gcp('nocreate');
if ~isempty(p)
    delete(p);
end
try
    parpool("Threads", 6);
catch ME
    warning('parpool not started (%s). Continue in serial mode.', ME.message);
end

%% 3) RhG geometry and K mesh
if ~isfile(hfmf_file)
    error('HFMF file not found: %s', hfmf_file);
end

S = load(hfmf_file, 'H_int', 'E_int');
if ~isfield(S, 'H_int') || ~isfield(S, 'E_int')
    error('HFMF file must contain H_int and E_int.');
end

H_int = S.H_int / 1000; % meV -> eV
E_int = S.E_int / 1000; % meV -> eV
Ef = tbNLC.calculate_ef(E_int(:), 0.5);

[Nkx_data, Nky_data, dimH1, dimH2] = size(H_int);


g = MTB.geometry("RhG");
a = 2.46; % Angstrom
g.a = [ ...
    1/2, -sqrt(3)/2, 0; ...
    1/2,  sqrt(3)/2, 0; ...
    0,    0,         1/a] * a;
g.b = inv(g.a') * 2*pi;

knum = Nkx_data - 1;
kxline = [-0.15, 0.15] / (2*pi);
kyline = [-0.15, 0.15] / (2*pi);
[Kx, Ky, ~] = g.get_Bulk2Dkmesh(kxline, kyline, knum);

if ~isequal(size(Kx), [Nkx_data, Nky_data])
    error('K-grid mismatch: size(Kx)=%s, H_int grid=%s', ...
        mat2str(size(Kx)), mat2str([Nkx_data, Nky_data]));
end

%% 4) Per-block SC + out blocks
H_blocks = tbNLC.extract_spinvalley_blocks(H_int, 1);

% 对角化参数
dopts = struct();
dopts.band_list = 1:2;

% Shift current 参数：这里明确开启中间量保存
sc_opts = struct();
sc_opts.band_list = 1:2;
sc_opts.periodicFD = false;
sc_opts.trimBoundary = true;
sc_opts.symBC = true;
sc_opts.verbose = false;
sc_opts.doGaugeFix = true;
sc_opts.g_s = 1;
sc_opts.saveIntermediates = true;   % 关键：保存 out 内中间量

nblock = numel(block_names_run);
Hk_blocks = cell(1, nblock);
for ib = 1:nblock
    bname = block_names_run{ib};
    if ~isfield(H_blocks, bname)
        error('Unknown block name: %s', bname);
    end
    Hk_blocks{ib} = permute(H_blocks.(bname), [3,4,1,2]);
end

U_blocks = cell(1, nblock);
E_blocks = cell(1, nblock);
sigma_blocks = cell(1, nblock);
out_sc_blocks = cell(1, nblock);

for ib = 1:nblock
    bname = block_names_run{ib};
    fprintf('SC block %d/%d: %s\n', ib, nblock, bname);

    Hk = Hk_blocks{ib};
    [U_blk, E_blk] = tbNLC.diag_mesh_fromHk(Hk, dopts);

    U_blocks{ib} = U_blk;
    E_blocks{ib} = E_blk;

    [sigma_blocks{ib}, out_sc_blocks{ib}] = tbNLC.shift_current_plane_fd_energy_skew_fromUE( ...
        Kx, Ky, U_blk, E_blk, Eph_list, Ef, kT, eta, sc_opts);
end

%% 5) Symmetrize each block, then sum
% sigma_blocks_sym = cell(1, nblock);
sym_info_blocks = cell(1, nblock);
for ib = 1:nblock
    [sigma_blocks_sym{ib}, sym_info_blocks{ib}] = symmetrize_sigma_rank3( ...
        sigma_blocks{ib}, 'C3v', 'mirror', 'Mx');
end

sigma_abc = sum_tensor_blocks(sigma_blocks);
sigma_abc_sym = sum_tensor_blocks(sigma_blocks_sym);
[sigma_C3v_from_sum, sym_info_from_sum] = symmetrize_sigma_rank3( ...
    sigma_abc, 'C3v', 'mirror', 'Mx');

%% 6) Plot
%%
fig = figure('Color', 'w');
hold on;
comp = ['x','y'];
for ia = 1:2
    for ib = 1:2
        for ic = 1:2
            % sig_curve = squeeze(sigma_abc(ia,ib,ic,:)) * 1e6 / 10;
            sig_curve = squeeze(sigma_blocks_sym{2}(ia,ib,ic,:)) * 1e6 / 10;
            plot(Eph_list, real(sig_curve), '--', ...
                'DisplayName', sprintf('%s%s%s', comp(ia), comp(ib), comp(ic)));
        end
    end
end
xlabel('\hbar\omega (eV)');
ylabel('\sigma_{abc} (scaled)');
title(sprintf('Shift Current (block-wise C3v sym), U=%.3f', U_fix));
legend('Location', 'best');
grid on;
box on;

exportgraphics(fig, fig_png, 'Resolution', 300);
close(fig);

%% 7) Save all results (including out_sc_blocks)
save(out_file, ...
    'hfmf_file', 'U_fix', 'ne_tag', ...
    'Kx', 'Ky', 'Eph_list', 'Ef', 'kT', 'eta', ...
    'block_names_run', 'dopts', 'sc_opts', ...
    'U_blocks', 'E_blocks', ...
    'sigma_blocks', 'out_sc_blocks', ...
    'sigma_blocks_sym', ...
    'sigma_abc', 'sigma_abc_sym', ...
    'sym_info_blocks', 'sigma_C3v_from_sum', 'sym_info_from_sum', ...
    '-v7.3');

fprintf('\nShift Current result saved:\n%s\n', out_file);
fprintf('Figure saved:\n%s\n', fig_png);

function tensor_sum = sum_tensor_blocks(tensor_blocks)
    idx = find(~cellfun(@isempty, tensor_blocks), 1, 'first');
    if isempty(idx)
        tensor_sum = [];
        return;
    end
    tensor_sum = zeros(size(tensor_blocks{idx}), 'like', tensor_blocks{idx});
    for ii = 1:numel(tensor_blocks)
        if ~isempty(tensor_blocks{ii})
            tensor_sum = tensor_sum + tensor_blocks{ii};
        end
    end
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
                            tmp = tmp + R(a,ap)*R(b,bp)*R(c,cp)*S(ap,bp,cp);
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
            error('mirror must be ''Mx'', ''My'', or numeric angle phi (rad).');
    end
end
