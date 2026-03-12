clc;
clear;

%% 1) User setup
% HFMF 输出目录（按需修改）
hfmf_data_dir = "/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data";
ne_tag = "0.0000e12";

% 多 U 扫描列表
% U_list = 1:31;
numU=31;
U_list=linspace(0,30,numU);
%%

% 计算设置
Eph_list = linspace(0.0, 0.1, 2000);
kT = 0.0;
eta = 5e-4;

% 只分析部分 valley 时可改这里，例如 {'K_up','K_dn'}
block_names_run = {'K_up', 'Kp_up', 'K_dn', 'Kp_dn'};

% 固定输出目录: 每个 U 存成单独 mat
result_dir = fullfile(pwd, 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

% 并行与保存配置
pool_profile = "Threads";     % "Threads" 或 "local"
pool_workers = 6;
save_mat_version = "-v7";   % 大文件建议 -v7.3

% Threads worker 不支持 save -v7.3，自动切换到 local 进程池
if strcmpi(pool_profile, "Threads") && strcmp(save_mat_version, "-v7.3")
    warning('Threads + -v7.3 save is unsupported. Switching pool_profile to \"local\".');
    pool_profile = "local";
end

%% 0) Parallel pool
p = gcp('nocreate');
if ~isempty(p)
    delete(p);
end
try
    parpool(pool_profile, pool_workers);
catch ME
    warning('parpool not started (%s). Continue in serial mode.', ME.message);
end

%% 2) RhG reciprocal basis (numeric, parfor friendly)
a = 2.46; % Angstrom
a_mat = [ ...
    1/2, -sqrt(3)/2, 0; ...
    1/2,  sqrt(3)/2, 0; ...
    0,    0,         1/a] * a;
b_mat = inv(a_mat') * 2*pi;

% K-mesh window in fractional coordinates
kxline = [-0.15, 0.15] / (2*pi);
kyline = [-0.15, 0.15] / (2*pi);

%% 3) NLC options
dopts = struct();
dopts.band_list = 1:2;
% 外层 U 已经并行，内部对角化禁用 parfor，避免嵌套并行冲突
dopts.useParfor = false;

sc_opts = struct();
sc_opts.band_list = 1:2;
sc_opts.periodicFD = false;
sc_opts.trimBoundary = true;
sc_opts.symBC = true;
sc_opts.verbose = false;
sc_opts.doGaugeFix = true;
sc_opts.g_s = 1;
sc_opts.saveIntermediates = false;

%% 4) Multi-U parallel run
NU = numel(U_list);
success_mask = false(1, NU);
status_msg = strings(1, NU);

parfor iu = 1:NU
    U_fix = U_list(iu);
    hfmf_file = fullfile(hfmf_data_dir, sprintf("ne=%s_U=%.3fdata.mat", ne_tag, U_fix));
    fprintf('[U %d/%d] Start U=%.3f, file=%s\n', iu, NU, U_fix, hfmf_file);

    try
        S = load(hfmf_file, 'H_int', 'E_int');
        H_int = S.H_int / 1000; % meV -> eV
        E_int = S.E_int / 1000; % meV -> eV
        Ef = tbNLC.calculate_ef(E_int(:), 0.5);

        [Nkx_data, Nky_data, ~, ~] = size(H_int);
        knum = Nkx_data - 1;
        [Kx, Ky] = build_bulk2d_kmesh(kxline, kyline, knum, b_mat);

        H_blocks = tbNLC.extract_spinvalley_blocks(H_int, 1);

        nblock = numel(block_names_run);
        U_blocks = cell(1, nblock);
        E_blocks = cell(1, nblock);
        sigma_blocks = cell(1, nblock);
        sigma_blocks_sym = cell(1, nblock);
        sym_info_blocks = cell(1, nblock);

        for ib = 1:nblock
            bname = block_names_run{ib};
            fprintf('[U %.3f] Block %d/%d: %s\n', U_fix, ib, nblock, bname);
            Hk = permute(H_blocks.(bname), [3,4,1,2]);

            [U_blk, E_blk] = tbNLC.diag_mesh_fromHk(Hk, dopts);
            U_blocks{ib} = U_blk;
            E_blocks{ib} = E_blk;

            sigma_blk = tbNLC.shift_current_plane_fd_energy_skew_fromUE( ...
                Kx, Ky, U_blk, E_blk, Eph_list, Ef, kT, eta, sc_opts);
            sigma_blocks{ib} = sigma_blk;

            [sigma_blocks_sym{ib}, sym_info_blocks{ib}] = symmetrize_sigma_rank3( ...
                sigma_blk, 'C3v', 'mirror', 'Mx');
        end

        sigma_abc = sum_tensor_blocks(sigma_blocks);
        sigma_abc_sym = sum_tensor_blocks(sigma_blocks_sym);
        [sigma_C3v_from_sum, sym_info_from_sum] = symmetrize_sigma_rank3( ...
            sigma_abc, 'C3v', 'mirror', 'Mx');

        u_tag_mat = sprintf('U_%.3f', U_fix);
        out_file = fullfile(result_dir, [u_tag_mat, '.mat']);
        out_payload = struct();
        out_payload.hfmf_file = hfmf_file;
        out_payload.U_fix = U_fix;
        out_payload.ne_tag = ne_tag;
        out_payload.Kx = Kx;
        out_payload.Ky = Ky;
        out_payload.Eph_list = Eph_list;
        out_payload.Ef = Ef;
        out_payload.kT = kT;
        out_payload.eta = eta;
        out_payload.block_names_run = block_names_run;
        out_payload.dopts = dopts;
        out_payload.sc_opts = sc_opts;
        out_payload.U_blocks = U_blocks;
        out_payload.E_blocks = E_blocks;
        out_payload.sigma_blocks = sigma_blocks;
        out_payload.sigma_blocks_sym = sigma_blocks_sym;
        out_payload.sigma_abc = sigma_abc;
        out_payload.sigma_abc_sym = sigma_abc_sym;
        out_payload.sym_info_blocks = sym_info_blocks;
        out_payload.sigma_C3v_from_sum = sigma_C3v_from_sum;
        out_payload.sym_info_from_sum = sym_info_from_sum;
        parsave_u_result(out_file, out_payload, save_mat_version);

        success_mask(iu) = true;
        status_msg(iu) = "OK: " + string(out_file);
        fprintf('[U %d/%d] Done U=%.3f -> %s\n', iu, NU, U_fix, out_file);
    catch ME
        success_mask(iu) = false;
        status_msg(iu) = "FAILED U=" + string(U_fix) + ": " + string(ME.message);
        fprintf('[U %d/%d] FAILED U=%.3f: %s\n', iu, NU, U_fix, ME.message);
    end
end

%% 5) Summary
for iu = 1:NU
    fprintf('%s\n', status_msg(iu));
end
fprintf('\nFinished: %d/%d succeeded.\n', nnz(success_mask), NU);

summary_file = fullfile(result_dir, 'Uscan_summary.mat');
save(summary_file, 'U_list', 'success_mask', 'status_msg', 'Eph_list', ...
    'block_names_run', 'dopts', 'sc_opts', 'hfmf_data_dir', 'ne_tag', ...
    'pool_profile', 'pool_workers', 'save_mat_version', '-v7.3');
fprintf('Summary saved: %s\n', summary_file);

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

function [Kx, Ky] = build_bulk2d_kmesh(kxline, kyline, knum, b_mat)
    kx_frac = linspace(kxline(1), kxline(2), knum + 1);
    ky_frac = linspace(kyline(1), kyline(2), knum + 1);
    [X, Y] = meshgrid(kx_frac, ky_frac);
    Kx = X * b_mat(1,1) + Y * b_mat(2,1);
    Ky = X * b_mat(1,2) + Y * b_mat(2,2);
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
        R = rot2(phi);
        M = R * Mx * R.';
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

function parsave_u_result(out_file, out_payload, save_mat_version)
    save(out_file, '-struct', 'out_payload', save_mat_version);
end
