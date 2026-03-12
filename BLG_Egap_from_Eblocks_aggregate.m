clc;
clear;

%% 0) User setup
result_dir = fullfile(pwd, 'results');
eblocks_file = fullfile(result_dir, 'Eblocks_aggregate_121U.mat');
out_file = fullfile(result_dir, 'Egap_from_Eblocks_121U.mat');

% 你指定的 K 点索引
k_idx = [1001, 1001]; % [ikx, iky]

% 模式选择:
% 默认 true  -> 一次性读入 E_blocks_all 后向量化计算（最快，但吃内存）
% 设环境变量 EGAP_SLICE_MODE=1 可切到 false（省内存）
use_full_load = ~strcmpi(getenv('EGAP_SLICE_MODE'), '1');

%% 1) Basic checks and metadata
if ~isfile(eblocks_file)
    error('Eblocks aggregate file not found: %s', eblocks_file);
end

vars = whos('-file', eblocks_file);
var_names = string({vars.name});
if any(var_names == "E_blocks_all")
    evar = 'E_blocks_all';
elseif any(var_names == "Eb_blocks_all")
    evar = 'Eb_blocks_all';
else
    error('No E_blocks_all / Eb_blocks_all found in %s', eblocks_file);
end

m = matfile(eblocks_file);
sz = size(m, evar); % expected: Nkx x Nky x 2 x nblock x numU
if numel(sz) ~= 5
    error('%s must be 5D: (Nkx,Nky,2,nblock,numU). Got size=%s', evar, mat2str(sz));
end

Nkx = sz(1);
Nky = sz(2);
nb_sel = sz(3);
nblock = sz(4);
numU = sz(5);

if nb_sel < 2
    error('%s third dimension must be >=2 (Ev/Ec). Got %d', evar, nb_sel);
end
if k_idx(1) < 1 || k_idx(1) > Nkx || k_idx(2) < 1 || k_idx(2) > Nky
    error('k_idx=[%d,%d] out of range. Grid size is [%d,%d].', ...
        k_idx(1), k_idx(2), Nkx, Nky);
end

meta = load(eblocks_file, 'U_list', 'block_names_run', 'success_mask', 'status_msg');
if isfield(meta, 'U_list') && numel(meta.U_list) == numU
    U_list = meta.U_list(:).';
else
    U_list = 1:numU;
end
if isfield(meta, 'block_names_run')
    block_names_run = meta.block_names_run;
else
    block_names_run = arrayfun(@(x) sprintf('block%d', x), 1:nblock, 'UniformOutput', false);
end
if isfield(meta, 'success_mask') && numel(meta.success_mask) == numU
    success_mask = logical(meta.success_mask(:).');
else
    success_mask = true(1, numU);
end
if isfield(meta, 'status_msg')
    status_msg = meta.status_msg;
else
    status_msg = strings(1, numU);
end

%% 2) Compute Egap
% Egap_K(iu,ib)      = Ec(k_idx)-Ev(k_idx)
% Egap_global(iu,ib) = min(Ec(:,:)) - max(Ev(:,:))

t0 = tic;

if use_full_load
    fprintf('Mode: FULL LOAD (fastest). Loading %s ...\n', evar);
    L = load(eblocks_file, evar);
    Eall = L.(evar);

    [Egap_K, Egap_global] = calc_egap_from_full(Eall, k_idx);
else
    fprintf('Mode: PER-U SLICE (memory friendly).\n');
    Egap_K = nan(numU, nblock);
    Egap_global = nan(numU, nblock);

    for iu = 1:numU
        fprintf('[%d/%d] U=%g\n', iu, numU, U_list(iu));
        if strcmp(evar, 'E_blocks_all')
            Eu = m.E_blocks_all(:,:,:,:,iu); % Nkx x Nky x 2 x nblock
        else
            Eu = m.Eb_blocks_all(:,:,:,:,iu);
        end

        [eg_k_u, eg_g_u] = calc_egap_from_slice(Eu, k_idx);
        Egap_K(iu, :) = eg_k_u;
        Egap_global(iu, :) = eg_g_u;
    end
end

% 无效 U 置 NaN
invalid = ~success_mask;
if any(invalid)
    Egap_K(invalid, :) = nan;
    Egap_global(invalid, :) = nan;
end

t_elapsed = toc(t0);

%% 3) Save compact output
save(out_file, ...
    'Egap_K', 'Egap_global', ...
    'U_list', 'block_names_run', 'k_idx', ...
    'Nkx', 'Nky', 'nblock', 'numU', ...
    'success_mask', 'status_msg', ...
    'evar', 'use_full_load', 't_elapsed', ...
    '-v7.3');

fprintf('\nDone. Saved Egap file:\n%s\n', out_file);
fprintf('Elapsed time: %.2f s\n', t_elapsed);

function [Egap_K, Egap_global] = calc_egap_from_full(Eall, k_idx)
% Eall: Nkx x Nky x 2 x nblock x numU
    sz = size(Eall);
    Nkx = sz(1);
    Nky = sz(2);
    nblock = sz(4);
    numU = sz(5);

    Ev = Eall(:,:,1,:,:); % Nkx x Nky x 1 x nblock x numU
    Ec = Eall(:,:,2,:,:); % Nkx x Nky x 1 x nblock x numU

    % K-point gap
    tmpK = Ec(k_idx(1), k_idx(2), 1, :, :) - Ev(k_idx(1), k_idx(2), 1, :, :); % 1x1x1xnblockxnumU
    Egap_K = reshape(permute(tmpK, [5,4,1,2,3]), [numU, nblock]);

    % Global gap: min(Ec(:,:)) - max(Ev(:,:)) for each (block,U)
    Ev2 = reshape(Ev, Nkx * Nky, nblock, numU);
    Ec2 = reshape(Ec, Nkx * Nky, nblock, numU);

    tmpG = min(Ec2, [], 1) - max(Ev2, [], 1); % 1 x nblock x numU
    Egap_global = reshape(permute(tmpG, [3,2,1]), [numU, nblock]);
end

function [eg_k_u, eg_g_u] = calc_egap_from_slice(Eu, k_idx)
% Eu: Nkx x Nky x 2 x nblock (single U)
    sz = size(Eu);
    Nkx = sz(1);
    Nky = sz(2);
    nblock = sz(4);

    Ev = Eu(:,:,1,:); % Nkx x Nky x 1 x nblock
    Ec = Eu(:,:,2,:); % Nkx x Nky x 1 x nblock

    tmpK = Ec(k_idx(1), k_idx(2), 1, :) - Ev(k_idx(1), k_idx(2), 1, :); % 1x1x1xnblock
    eg_k_u = reshape(tmpK, [1, nblock]);

    Ev2 = reshape(Ev, Nkx * Nky, nblock);
    Ec2 = reshape(Ec, Nkx * Nky, nblock);
    eg_g_u = min(Ec2, [], 1) - max(Ev2, [], 1);
end
