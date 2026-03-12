clc;
clear;

%% 0) User setup
% 原始 U_*.mat 文件目录
result_dir = fullfile(pwd, 'results');

% 输出文件
out_file = fullfile(result_dir, 'Egap_from_U_files.mat');

% 指定 K 点索引（你定义的 K 点 gap）
k_idx = [2001, 2001]; % [ikx, iky]

% 是否强制使用固定 U_list；默认 false（自动扫描目录下 U_*.mat）
% 设环境变量 USE_FIXED_ULIST=1 可启用固定 U_list 模式
use_fixed_ulist = strcmpi(getenv('USE_FIXED_ULIST'), '1');

% 固定 U_list（仅在 use_fixed_ulist=true 时使用）
numU = 121;
U_list_fixed = linspace(0, 30, numU);

%% 1) Collect file list
if use_fixed_ulist
    U_list = U_list_fixed(:).';
    file_list = strings(1, numel(U_list));
    for iu = 1:numel(U_list)
        file_list(iu) = string(fullfile(result_dir, sprintf('U_%.3f.mat', U_list(iu))));
    end
else
    D = dir(fullfile(result_dir, 'U_*.mat'));
    if isempty(D)
        error('No U_*.mat found in %s', result_dir);
    end

    Uvals = nan(1, numel(D));
    files = strings(1, numel(D));
    for i = 1:numel(D)
        files(i) = string(fullfile(D(i).folder, D(i).name));
        tok = regexp(D(i).name, '^U_([-+]?\d*\.?\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            Uvals(i) = str2double(tok{1});
        end
    end

    valid = ~isnan(Uvals);
    if ~any(valid)
        error('No valid U_*.mat filename could be parsed in %s', result_dir);
    end

    Uvals = Uvals(valid);
    files = files(valid);
    [U_list, ord] = sort(Uvals, 'ascend');
    file_list = files(ord);
end

numU = numel(U_list);

%% 2) First valid file: infer dimensions
first_ok = false;
nblock = 0;
Nkx = 0;
Nky = 0;
nb_sel = 0;
block_names_run = {};

for iu = 1:numU
    fpath = char(file_list(iu));
    if ~isfile(fpath)
        continue;
    end

    S0 = load(fpath, 'E_blocks', 'block_names_run');
    if ~isfield(S0, 'E_blocks') || isempty(S0.E_blocks)
        continue;
    end

    nblock = numel(S0.E_blocks);
    sz0 = size(S0.E_blocks{1});
    if numel(sz0) ~= 3 || sz0(3) < 2
        error('Invalid E_blocks shape in %s, expected [Nkx,Nky,>=2]. got %s', fpath, mat2str(sz0));
    end

    Nkx = sz0(1);
    Nky = sz0(2);
    nb_sel = sz0(3);

    if isfield(S0, 'block_names_run')
        block_names_run = S0.block_names_run;
    else
        block_names_run = arrayfun(@(x) sprintf('block%d', x), 1:nblock, 'UniformOutput', false);
    end

    first_ok = true;
    break;
end

if ~first_ok
    error('No valid file containing E_blocks found.');
end

if k_idx(1) < 1 || k_idx(1) > Nkx || k_idx(2) < 1 || k_idx(2) > Nky
    error('k_idx=[%d,%d] out of range. Grid size=[%d,%d].', k_idx(1), k_idx(2), Nkx, Nky);
end

%% 3) Preallocate outputs
% Egap_K(iu,ib)      = Ec(k_idx)-Ev(k_idx)
% Egap_global(iu,ib) = min(Ec(:,:))-max(Ev(:,:))
Egap_K = nan(numU, nblock);
Egap_global = nan(numU, nblock);

success_mask = false(1, numU);
status_msg = strings(1, numU);

%% 4) Loop over U files
t0 = tic;
for iu = 1:numU
    fpath = char(file_list(iu));
    fprintf('[%d/%d] U=%.3f, file=%s\n', iu, numU, U_list(iu), fpath);

    if ~isfile(fpath)
        status_msg(iu) = "MISSING: " + string(fpath);
        continue;
    end

    try
        S = load(fpath, 'E_blocks');
        if ~isfield(S, 'E_blocks') || numel(S.E_blocks) ~= nblock
            status_msg(iu) = "SKIP(E_blocks missing/nblock mismatch): " + string(fpath);
            continue;
        end

        ok = true;
        for ib = 1:nblock
            Eb = S.E_blocks{ib};
            if ~isequal(size(Eb), [Nkx, Nky, nb_sel])
                ok = false;
                break;
            end

            Ev = Eb(:,:,1);
            Ec = Eb(:,:,2);

            Egap_K(iu, ib) = Ec(k_idx(1), k_idx(2)) - Ev(k_idx(1), k_idx(2));
            Egap_global(iu, ib) = min(Ec, [], 'all') - max(Ev, [], 'all');
        end

        if ~ok
            status_msg(iu) = "SKIP(E_blocks shape mismatch): " + string(fpath);
            continue;
        end

        success_mask(iu) = true;
        status_msg(iu) = "OK";
    catch ME
        status_msg(iu) = "FAILED: " + string(ME.message);
    end
end

t_elapsed = toc(t0);

%% 5) Save compact output
save(out_file, ...
    'Egap_K', 'Egap_global', ...
    'U_list', 'file_list', ...
    'block_names_run', 'k_idx', ...
    'Nkx', 'Nky', 'nb_sel', 'nblock', 'numU', ...
    'success_mask', 'status_msg', 't_elapsed', ...
    '-v7.3');

fprintf('\nDone. success=%d/%d\n', nnz(success_mask), numU);
fprintf('Saved: %s\n', out_file);
fprintf('Elapsed time: %.2f s\n', t_elapsed);
