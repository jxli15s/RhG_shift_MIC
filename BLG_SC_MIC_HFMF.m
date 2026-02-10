clc;
clear;

%% 0) Runtime switches
% 集群/无界面运行建议:
% - is_cluster=true: 不弹图窗
% - save_figures=true: 自动保存图到输出目录
% - save_results=true: 保存 mat 结果
is_cluster = true;
save_figures = true;
save_results = true;
plot_each_block = true;      % true 时额外输出每个 valley block 的曲线
% 环境变量控制: USE_PARFOR_BLOCKS=1 -> 4 个 block 并行（内存充足时可开）
use_parfor_over_blocks = strcmp(getenv('USE_PARFOR_BLOCKS'), '1');

% HFMF 输出文件（按需修改）
hfmf_mat_file = "/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data/ne=0.0000e12_U=14.000data.mat";
% hfmf_mat_file = "/Volumes/T9/work/tb/matlab/data/RhG_HFMF/chiho/sequence/2001/ne=0.0000e12_U=0.000data.mat";

if is_cluster
    set(groot, 'defaultFigureVisible', 'off');
end

run_tag = datestr(now, 'yyyymmdd_HHMMSS');
output_root = fullfile(pwd, 'outputs_hfmf', run_tag);
fig_root = fullfile(output_root, 'figures');
if save_figures && ~exist(fig_root, 'dir')
    mkdir(fig_root);
end

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

%% 2) 构造 RhG 模型参数
g = MTB.geometry("RhG");
a = 2.46; % Angstrom
g.a = [ ...
    1/2, -sqrt(3)/2, 0; ...
    1/2,  sqrt(3)/2, 0; ...
    0,    0,         1/a] * a;
g.b = inv(g.a') * 2*pi;

Layer_N = 2;
pars = struct( ...
    'v0', 3.16, ...
    'gamma1', 0.46, ...
    'gamma2', -0.017, ...
    'gamma3', -0.30, ...
    'gamma4', -0.086, ...
    'uext', 0.015, ...
    'delta', -0.0011, ...
    'xi', 1, ...
    'N', Layer_N);

%% 3) 读取 HFMF 哈密顿量并建立 k 网格
if ~isfile(hfmf_mat_file)
    error('HFMF data file not found: %s', hfmf_mat_file);
end

S = load(hfmf_mat_file, 'H_int', 'E_int');
if ~isfield(S, 'H_int') || ~isfield(S, 'E_int')
    error('Input mat file must contain H_int and E_int.');
end

H_int = S.H_int / 1000; % meV -> eV
E_int = S.E_int / 1000; % meV -> eV

[Nkx_data, Nky_data, dimH1, dimH2] = size(H_int);
if dimH1 ~= dimH2
    error('H_int must have shape Nkx x Nky x Nh x Nh.');
end
if Nkx_data ~= Nky_data
    warning('Non-square k grid detected: Nkx=%d, Nky=%d.', Nkx_data, Nky_data);
end

% get_Bulk2Dkmesh 返回 (knum+1)x(knum+1)，因此这里按 H_int 反推 knum
knum = Nkx_data - 1;
kxline = [-0.15, 0.15] / (2*pi);
kyline = [-0.15, 0.15] / (2*pi);
[Kx, Ky, ~] = g.get_Bulk2Dkmesh(kxline, kyline, knum);
if ~isequal(size(Kx), [Nkx_data, Nky_data])
    error('K-grid size mismatch: size(Kx)=%s vs H_int grid=[%d %d].', ...
        mat2str(size(Kx)), Nkx_data, Nky_data);
end

%% 4) 按 valley block 对角化得到 U/E
dopts = struct();
dopts.band_list = 1:2;

H_blocks = tbNLC.extract_spinvalley_blocks(H_int, 1);
block_names = {'K_up', 'Kp_up', 'K_dn', 'Kp_dn'};
Hk_blocks = { ...
    permute(H_blocks.K_up,  [3,4,1,2]), ...
    permute(H_blocks.Kp_up, [3,4,1,2]), ...
    permute(H_blocks.K_dn,  [3,4,1,2]), ...
    permute(H_blocks.Kp_dn, [3,4,1,2])};

nblock = numel(Hk_blocks);
U_blocks = cell(1, nblock);
E_blocks = cell(1, nblock);

for ib = 1:nblock
    [U_blocks{ib}, E_blocks{ib}] = tbNLC.diag_mesh_fromHk(Hk_blocks{ib}, dopts);
end

all_eigs_cell = cellfun(@(x) x(:), E_blocks, 'UniformOutput', false);
Ef = tbNLC.calculate_ef(vertcat(all_eigs_cell{:}), 0.5);

%% 5) 频率网格与展宽
kT = 0.0;
eta = 5e-4;
Eph_list = linspace(0.0, 0.1, 5000);

%% 6) 配置 MIC / Shift current 参数
sigma_opts = struct();
sigma_opts.band_list = 1:2;
sigma_opts.periodicFD = false;
sigma_opts.trimBoundary = true;
sigma_opts.symBC = true;
sigma_opts.verbose = false;
sigma_opts.doGaugeFix = true;
sigma_opts.g_s = 1;
sigma_opts.saveIntermediates = false;

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
run_opts.do_sigma = true;   % 新增: 打开 shift current 计算
run_opts.do_mic = true;
run_opts.sigma_opts = sigma_opts;
run_opts.mic_opts = mic_opts;

%% 7) 对每个 valley block 分别计算，再求和
resp_blocks = cell(1, nblock);

tic;
if use_parfor_over_blocks
    parfor ib = 1:nblock
        resp_blocks{ib} = tbNLC.compute_nonlinear_conductivity_fromUE( ...
            Kx, Ky, U_blocks{ib}, E_blocks{ib}, Eph_list, Ef, kT, eta, run_opts);
    end
else
    for ib = 1:nblock
        fprintf('Compute block %d/%d: %s\n', ib, nblock, block_names{ib});
        resp_blocks{ib} = tbNLC.compute_nonlinear_conductivity_fromUE( ...
            Kx, Ky, U_blocks{ib}, E_blocks{ib}, Eph_list, Ef, kT, eta, run_opts);
    end
end
toc;

sigma_blocks = cell(1, nblock);
eta_blocks = cell(1, nblock);
for ib = 1:nblock
    if run_opts.do_sigma && isfield(resp_blocks{ib}, 'sigma_abc')
        sigma_blocks{ib} = resp_blocks{ib}.sigma_abc;
    end
    if run_opts.do_mic && isfield(resp_blocks{ib}, 'eta_abc')
        eta_blocks{ib} = resp_blocks{ib}.eta_abc;
    end
end

sigma_abc = [];
eta_abc = [];
sigma_K = [];
sigma_Kp = [];
eta_K = [];
eta_Kp = [];

if run_opts.do_sigma
    sigma_abc = sum_tensor_blocks(sigma_blocks);
    sigma_K = sum_tensor_blocks(sigma_blocks([1,3]));   % K valley (spin up/down)
    sigma_Kp = sum_tensor_blocks(sigma_blocks([2,4]));  % K' valley (spin up/down)
end

if run_opts.do_mic
    eta_abc = sum_tensor_blocks(eta_blocks);
    eta_K = sum_tensor_blocks(eta_blocks([1,3]));
    eta_Kp = sum_tensor_blocks(eta_blocks([2,4]));
end

%% 8) 作图（总和 + 可选每 block）
if run_opts.do_sigma
    plot_response_tensor(Eph_list, sigma_abc, 1e5, 'Shift Current \sigma_{abc} (sum all blocks)', ...
        fullfile(fig_root, 'sigma_total'), save_figures, is_cluster);
end
if run_opts.do_mic
    plot_response_tensor(Eph_list, eta_abc, 1e-8, 'MIC Metric \eta_{abc} (sum all blocks)', ...
        fullfile(fig_root, 'eta_total'), save_figures, is_cluster);
end

if plot_each_block
    for ib = 1:nblock
        if run_opts.do_sigma && ~isempty(sigma_blocks{ib})
            plot_response_tensor(Eph_list, sigma_blocks{ib}, 1e5, ...
                sprintf('Shift Current %s', block_names{ib}), ...
                fullfile(fig_root, sprintf('sigma_%s', block_names{ib})), ...
                save_figures, is_cluster);
        end
        if run_opts.do_mic && ~isempty(eta_blocks{ib})
            plot_response_tensor(Eph_list, eta_blocks{ib}, 1e-8, ...
                sprintf('MIC Metric %s', block_names{ib}), ...
                fullfile(fig_root, sprintf('eta_%s', block_names{ib})), ...
                save_figures, is_cluster);
        end
    end
end

%% 9) 保存结果
if save_results
    if ~exist(output_root, 'dir')
        mkdir(output_root);
    end
    save(fullfile(output_root, 'nonlinear_conductivity_hfmf.mat'), ...
        'hfmf_mat_file', 'block_names', ...
        'Eph_list', 'Ef', 'kT', 'eta', ...
        'dopts', 'sigma_opts', 'mic_opts', 'run_opts', ...
        'sigma_blocks', 'eta_blocks', ...
        'sigma_abc', 'eta_abc', 'sigma_K', 'sigma_Kp', 'eta_K', 'eta_Kp', ...
        '-v7.3');
end

function tensor_sum = sum_tensor_blocks(tensor_blocks)
%SUM_TENSOR_BLOCKS 对 cell 内 2x2x2xNw 张量按 block 求和。

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

function plot_response_tensor(Eph_list, tensor_abcw, scale_factor, fig_title, ...
    out_prefix, save_figure, close_after_save)
%PLOT_RESPONSE_TENSOR 将 2x2x2xNw 张量分量画在同一张图上并可选保存。

    if isempty(tensor_abcw)
        return;
    end

    fig = figure('Name', fig_title, 'Color', 'w');
    hold on;
    comp = ['x', 'y'];
    for a = 1:2
        for b = 1:2
            for c = 1:2
                y = squeeze(tensor_abcw(a,b,c,:)) * scale_factor;
                plot(Eph_list, real(y), 'LineWidth', 1.1, ...
                    'DisplayName', sprintf('%s%s%s', comp(a), comp(b), comp(c)));
            end
        end
    end
    xlabel('\hbar\omega (eV)');
    ylabel('Response (scaled)');
    title(fig_title);
    legend('Location', 'best');
    grid on;
    box on;

    if save_figure
        exportgraphics(fig, [out_prefix, '.png'], 'Resolution', 300);
        exportgraphics(fig, [out_prefix, '.pdf'], 'ContentType', 'vector');
        savefig(fig, [out_prefix, '.fig']);
    end

    if close_after_save
        close(fig);
    end
end
