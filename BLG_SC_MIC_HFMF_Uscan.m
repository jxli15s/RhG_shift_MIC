clc;
clear;

%% 0) Scan setup
% 目标: 扫描 U=1:31，逐 valley 计算 shift current + MIC，并仅保存轻量结果。
%
% 保存字段（默认）:
%   sigma_shift_u : 2x2x2xNw x NU
%   eta_mic_u     : 2x2x2xNw x NU
%   Kx, Ky, Eph_list, U_list, success_mask, Ef_list
%
% 不保存任何 K-resolved 中间量（不保存 out / r_mn / g2 等），降低存储占用。
% 环境变量控制:
%   SAVE_MIC_TENSOR=0  -> 仅保存 shift current
%   其他/未设置        -> 同时保存 MIC 张量
save_mic_tensor = ~strcmpi(getenv('SAVE_MIC_TENSOR'), '0');

% HFMF 数据目录与文件命名模板（按实际路径修改）
hfmf_data_dir = "/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data";
file_pattern = "ne=0.0000e12_U=%.3fdata.mat";

U_list = 1:31;

% 频率和展宽参数
Eph_list = linspace(0.0, 0.1, 5000); % eV
kT = 0.0;
eta = 5e-4;

% 输出路径
run_tag = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(pwd, 'outputs_hfmf', ['Uscan_', run_tag]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
out_file = fullfile(out_dir, 'scan_U1to31_sigma_mic.mat');

%% 1) 并行池（可选）
p = gcp('nocreate');
if ~isempty(p)
    delete(p);
end
try
    parpool('local', 6);
catch ME
    warning('parpool not started (%s). Continue in serial mode.', ME.message);
end

%% 2) RhG 几何（用于生成 Kx/Ky）
g = MTB.geometry("RhG");
a = 2.46; % Angstrom
g.a = [ ...
    1/2, -sqrt(3)/2, 0; ...
    1/2,  sqrt(3)/2, 0; ...
    0,    0,         1/a] * a;
g.b = inv(g.a') * 2*pi;

kxline = [-0.15, 0.15] / (2*pi);
kyline = [-0.15, 0.15] / (2*pi);

%% 3) tbNLC 参数（关闭大中间量输出）
dopts = struct();
dopts.band_list = 1:2;

sigma_opts = struct();
sigma_opts.band_list = 1:2;
sigma_opts.periodicFD = false;
sigma_opts.trimBoundary = true;
sigma_opts.symBC = true;
sigma_opts.verbose = false;
sigma_opts.doGaugeFix = true;
sigma_opts.g_s = 1;
sigma_opts.saveIntermediates = false; % 仅在请求第二输出时才有影响

mic_opts = struct();
mic_opts.band_list = 1:2;
mic_opts.periodicFD = false;
mic_opts.trimBoundary = true;
mic_opts.verbose = false;
mic_opts.doGaugeFix = true;
mic_opts.g_s = 1;
mic_opts.positiveDE = true;
mic_opts.saveIntermediates = false;
mic_opts.saveFullMN = false;

run_opts = struct();
run_opts.do_sigma = true;
run_opts.do_mic = true;
run_opts.sigma_opts = sigma_opts;
run_opts.mic_opts = mic_opts;

%% 4) 预分配（在首个有效文件上初始化）
NU = numel(U_list);
Nw = numel(Eph_list);

Kx = [];
Ky = [];
grid_initialized = false;

sigma_shift_u = [];
eta_mic_u = [];
success_mask = false(1, NU);
Ef_list = nan(1, NU);

block_names = {'K_up', 'Kp_up', 'K_dn', 'Kp_dn'};

%% 5) 主循环：U scan
for iu = 1:NU
    Uval = U_list(iu);
    fpath = fullfile(hfmf_data_dir, sprintf(file_pattern, Uval));

    if ~isfile(fpath)
        warning('Skip U=%g: file not found: %s', Uval, fpath);
        continue;
    end

    fprintf('\n[%d/%d] U=%g, file=%s\n', iu, NU, Uval, fpath);
    S = load(fpath, 'H_int', 'E_int');
    if ~isfield(S, 'H_int') || ~isfield(S, 'E_int')
        warning('Skip U=%g: H_int/E_int missing in file.', Uval);
        continue;
    end

    H_int = S.H_int / 1000; % meV -> eV
    E_int = S.E_int / 1000; % meV -> eV
    Ef = tbNLC.calculate_ef(E_int(:), 0.5);
    Ef_list(iu) = Ef;

    [Nkx_data, Nky_data, dimH1, dimH2] = size(H_int);
    if dimH1 ~= dimH2
        warning('Skip U=%g: invalid H_int shape.', Uval);
        continue;
    end

    if ~grid_initialized
        knum = Nkx_data - 1; % get_Bulk2Dkmesh returns (knum+1)x(knum+1)
        [Kx, Ky, ~] = g.get_Bulk2Dkmesh(kxline, kyline, knum);

        if ~isequal(size(Kx), [Nkx_data, Nky_data])
            error('K-grid mismatch on first valid file. size(Kx)=%s, H_int grid=%s', ...
                mat2str(size(Kx)), mat2str([Nkx_data, Nky_data]));
        end

        sigma_shift_u = nan(2,2,2,Nw,NU);
        eta_mic_u = nan(2,2,2,Nw,NU);
        grid_initialized = true;
    else
        if ~isequal(size(Kx), [Nkx_data, Nky_data])
            warning('Skip U=%g: K-grid size changed to [%d %d].', Uval, Nkx_data, Nky_data);
            continue;
        end
    end

    H_blocks = tbNLC.extract_spinvalley_blocks(H_int, 1);

    sigma_sum = zeros(2,2,2,Nw);
    eta_sum = zeros(2,2,2,Nw);

    for ib = 1:numel(block_names)
        name = block_names{ib};
        Hk = permute(H_blocks.(name), [3,4,1,2]);

        [U_blk, E_blk] = tbNLC.diag_mesh_fromHk(Hk, dopts);
        resp_blk = tbNLC.compute_nonlinear_conductivity_fromUE( ...
            Kx, Ky, U_blk, E_blk, Eph_list, Ef, kT, eta, run_opts);

        sigma_sum = sigma_sum + resp_blk.sigma_abc;
        eta_sum = eta_sum + resp_blk.eta_abc;
    end

    sigma_shift_u(:,:,:,:,iu) = sigma_sum;
    eta_mic_u(:,:,:,:,iu) = eta_sum;
    success_mask(iu) = true;
end

%% 6) 保存最小结果
if ~grid_initialized
    error('No valid U file was processed. Nothing to save.');
end

if save_mic_tensor
    save(out_file, ...
        'sigma_shift_u', 'eta_mic_u', ...
        'Kx', 'Ky', 'Eph_list', 'U_list', 'success_mask', 'Ef_list', ...
        '-v7.3');
else
    save(out_file, ...
        'sigma_shift_u', ...
        'Kx', 'Ky', 'Eph_list', 'U_list', 'success_mask', 'Ef_list', ...
        '-v7.3');
end

fprintf('\nSaved scan result to:\n%s\n', out_file);
fprintf('Valid U count: %d / %d\n', nnz(success_mask), NU);
