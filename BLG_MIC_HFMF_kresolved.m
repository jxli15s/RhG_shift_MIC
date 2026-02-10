clc;
clear;

%% 0) User setup (fixed U)
% 固定 U 的 HFMF 文件参数
U_fix = 9.0;
U_fix = 11.0;
ne_tag = "0.0000e12";

% HFMF 输出目录（按需修改）
% hfmf_data_dir = "/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data";
% hfmf_data_dir = "/Volumes/T9/work/tb/matlab/data/RhG_HFMF/chiho/sequence/2001";
hfmf_data_dir = "/Volumes/T9/work/tb/matlab/data/RhG_HFMF/chiho/sequence/1001-1001";

env_u = getenv('U_FIX');
if ~isempty(env_u)
    U_try = str2double(env_u);
    if ~isnan(U_try)
        U_fix = U_try;
    end
end
hfmf_file = fullfile(hfmf_data_dir, sprintf("ne=%s_U=%.3fdata.mat", ne_tag, U_fix));
env_hfile = getenv('HFMF_FILE');
if ~isempty(env_hfile)
    hfmf_file = string(env_hfile);
end

% 计算设置
Eph_list = linspace(0.0, 0.2, 5000);
kT = 0.0;
eta = 5e-5;
env_nw = getenv('EPH_NW');
if ~isempty(env_nw)
    nw_try = round(str2double(env_nw));
    if ~isnan(nw_try) && nw_try > 1
        Eph_list = linspace(0.0, 0.1, nw_try);
    end
end

% 只分析部分 valley 时可改这里，例如 {'K_up'}
block_names_run = {'K_up', 'Kp_up', 'K_dn', 'Kp_dn'};
env_blocks = getenv('BLOCKS_RUN'); % 例如: K_up,K_dn
if ~isempty(env_blocks)
    block_names_run = cellstr(strtrim(split(string(env_blocks), ',')))';
    block_names_run = block_names_run(~cellfun(@isempty, block_names_run));
end

% 可选: 对 H_int 的 k 网格做步长抽样（用于快速测试）
k_stride = 1;
env_stride = getenv('K_STRIDE');
if ~isempty(env_stride)
    stride_try = round(str2double(env_stride));
    if ~isnan(stride_try) && stride_try >= 1
        k_stride = stride_try;
    end
end

% 输出目录
run_tag = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(pwd, 'outputs_hfmf', ['MIC_kresolved_U', num2str(U_fix), '_', run_tag]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
out_file = fullfile(out_dir, 'mic_kresolved.mat');

%% 1) Parallel pool (optional)
p = gcp('nocreate');
if ~isempty(p)
    delete(p);
end
try
    parpool('local', 6);
catch ME
    warning('parpool not started (%s). Continue in serial mode.', ME.message);
end

%% 2) RhG geometry for K mesh
g = MTB.geometry("RhG");
a = 2.46; % Angstrom
g.a = [ ...
    1/2, -sqrt(3)/2, 0; ...
    1/2,  sqrt(3)/2, 0; ...
    0,    0,         1/a] * a;
g.b = inv(g.a') * 2*pi;

%% 3) Load HFMF Hamiltonian
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

if k_stride > 1
    H_int = H_int(1:k_stride:end, 1:k_stride:end, :, :);
    fprintf('Apply k-stride downsampling: %d\n', k_stride);
end

[Nkx_data, Nky_data, dimH1, dimH2] = size(H_int);
if dimH1 ~= dimH2
    error('Invalid H_int shape: Nkx x Nky x Nh x Nh expected.');
end

% 与 HFMF 网格一致的 Kx/Ky
knum = Nkx_data - 1; % get_Bulk2Dkmesh -> (knum+1)x(knum+1)
kxline = [-0.15, 0.15] / (2*pi);
kyline = [-0.15, 0.15] / (2*pi);
[Kx, Ky, ~] = g.get_Bulk2Dkmesh(kxline, kyline, knum);
if ~isequal(size(Kx), [Nkx_data, Nky_data])
    error('K-grid mismatch: size(Kx)=%s, H_int grid=%s', ...
        mat2str(size(Kx)), mat2str([Nkx_data, Nky_data]));
end

%% 4) Valley blocks + eigensystem
H_blocks = tbNLC.extract_spinvalley_blocks(H_int, 1);

dopts = struct();
dopts.band_list = 1:2;

mic_opts = struct();
mic_opts.band_list = 1:2;
mic_opts.periodicFD = false;
mic_opts.trimBoundary = true;
mic_opts.verbose = false;
mic_opts.doGaugeFix = true;
mic_opts.g_s = 1;
mic_opts.positiveDE = true;

% 核心: 保留 k 分辨中间量（体积较大）
% mic_opts.saveIntermediates = true;
% mic_opts.saveFullMN = true;
mic_opts.saveIntermediates = false;
mic_opts.saveFullMN = false;

% 结果容器
nblock = numel(block_names_run);
eta_blocks = cell(1, nblock);
out_mic_blocks = cell(1, nblock);
E_blocks = cell(1, nblock);

for ib = 1:nblock
    bname = block_names_run{ib};
    if ~isfield(H_blocks, bname)
        error('Unknown block name: %s', bname);
    end
    fprintf('MIC block %d/%d: %s\n', ib, nblock, bname);

    Hk = permute(H_blocks.(bname), [3,4,1,2]);
    [U_blk, E_blk] = tbNLC.diag_mesh_fromHk(Hk, dopts);
    E_blocks{ib} = E_blk;
    tic;
    [eta_blocks{ib}, out_mic_blocks{ib}] = tbNLC.mic_metric_plane_fromUE_mn( ...
        Kx, Ky, U_blk, E_blk, Eph_list, Ef, kT, eta, mic_opts);
    toc;
end

%% 5) Sum over selected blocks
eta_sum = sum_tensor_blocks(eta_blocks);
eta_abc=sum_tensor_blocks(eta_blocks);
eta_abc1=eta_blocks{1};
%%
figure()
hold on
str=['x','y'];
for i=1:2
    for j=1:2
        for k=1:2
            sig_xxy = squeeze(eta_abc1(i,j,k,:))*10^-13*10^6/10;
            % sig_xxy = squeeze(sig_r(i,j,k,:))/10;
            % sig_xxy = squeeze(sigma_C3v(i,j,k,:))*10^6/10/10^20;
            plot(Eph_list, real(sig_xxy), '--','DisplayName',  sprintf('%s%s%s',str(i),str(j),str(k)));
            legend
        end
    end
end
%% 6) Save (k-resolved intermediates kept)
save(out_file, ...
    'hfmf_file', 'U_fix', 'ne_tag', ...
    'Kx', 'Ky', 'Eph_list', 'Ef', 'kT', 'eta', ...
    'block_names_run', 'dopts', 'mic_opts', ...
    'E_blocks', 'eta_blocks', 'eta_sum', 'out_mic_blocks', ...
    '-v7.3');

fprintf('\nMIC k-resolved result saved:\n%s\n', out_file);

function tensor_sum = sum_tensor_blocks(tensor_blocks)
% 对 cell 内 2x2x2xNw 张量求和。
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
