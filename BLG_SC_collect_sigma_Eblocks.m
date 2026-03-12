clc;
clear;

%% 0) User setup
% 结果文件所在目录（就是 U_*.mat 所在目录）
result_dir = fullfile(pwd, 'results');

% U 列表
numU = 121;
U_list = linspace(0, 30, numU);
expected_Nw = 5000;

% 输出文件名
sigma_out_file = fullfile(result_dir, sprintf('sigma_aggregate_%dU.mat', numU));
eblocks_out_file = fullfile(result_dir, sprintf('Eblocks_aggregate_%dU.mat', numU));

%% 1) Pre-allocation handles
success_mask = false(1, numU);
status_msg = strings(1, numU);

sigma_all = [];
sigma_blocks_all = [];
Eph_list = [];
block_names_run = {};

E_blocks_all = [];

first_ok = false;
Nw = 0;
nblock = 0;

%% 2) Main loop: load each U_%.3f.mat and aggregate
for iu = 1:numU
    Uv = U_list(iu);
    fpath = fullfile(result_dir, sprintf('U_%.3f.mat', Uv));

    fprintf('[%d/%d] loading %s\n', iu, numU, fpath);

    if ~isfile(fpath)
        status_msg(iu) = "MISSING: " + string(fpath);
        continue;
    end

    try
        S = load(fpath, 'sigma_abc', 'sigma_blocks', 'E_blocks', 'Eph_list', 'block_names_run');

        if ~isfield(S, 'sigma_abc') || ~isfield(S, 'sigma_blocks')
            status_msg(iu) = "SKIP(no sigma vars): " + string(fpath);
            continue;
        end

        if ~first_ok
            Nw = size(S.sigma_abc, 4);
            if Nw ~= expected_Nw
                error('Expected Nw=%d, but got Nw=%d in first valid file.', expected_Nw, Nw);
            end
            nblock = numel(S.sigma_blocks);
            if nblock < 1
                error('sigma_blocks is empty in first valid file.');
            end
            if ~isfield(S, 'E_blocks') || numel(S.E_blocks) ~= nblock
                error('E_blocks missing or nblock mismatch in first valid file.');
            end
            e_shape = size(S.E_blocks{1});
            if numel(e_shape) ~= 3
                error('Expected E_blocks{ib} to be 3D (Nkx,Nky,nb_sel), got size=%s', mat2str(e_shape));
            end

            sigma_all = nan(2, 2, 2, expected_Nw, numU);
            sigma_blocks_all = nan(2, 2, 2, expected_Nw, nblock, numU);
            E_blocks_all = nan(e_shape(1), e_shape(2), e_shape(3), nblock, numU);

            if isfield(S, 'Eph_list')
                Eph_list = S.Eph_list;
            else
                Eph_list = [];
            end
            if isfield(S, 'block_names_run')
                block_names_run = S.block_names_run;
            else
                block_names_run = {};
            end

            first_ok = true;
            fprintf('Detected Nw=%d, nblock=%d from first valid file.\n', Nw, nblock);
        end

        % consistency checks
        if size(S.sigma_abc, 1) ~= 2 || size(S.sigma_abc, 2) ~= 2 || size(S.sigma_abc, 3) ~= 2
            status_msg(iu) = "SKIP(bad sigma_abc shape): " + string(fpath);
            continue;
        end
        if size(S.sigma_abc, 4) ~= expected_Nw
            status_msg(iu) = "SKIP(Nw mismatch): " + string(fpath);
            continue;
        end
        if numel(S.sigma_blocks) ~= nblock
            status_msg(iu) = "SKIP(nblock mismatch): " + string(fpath);
            continue;
        end

        % write sigma_abc
        sigma_all(:,:,:,:,iu) = S.sigma_abc;

        % write sigma_blocks
        block_ok = true;
        for ib = 1:nblock
            sb = S.sigma_blocks{ib};
            if ~isequal(size(sb), [2, 2, 2, Nw])
                block_ok = false;
                break;
            end
            sigma_blocks_all(:,:,:,:,ib,iu) = sb;
        end
        if ~block_ok
            status_msg(iu) = "SKIP(sigma_blocks shape mismatch): " + string(fpath);
            continue;
        end

        % write E_blocks (numeric array: Nkx x Nky x nb_sel x nblock x numU)
        if ~isfield(S, 'E_blocks') || numel(S.E_blocks) ~= nblock
            status_msg(iu) = "SKIP(E_blocks missing/nblock mismatch): " + string(fpath);
            continue;
        end
        e_ok = true;
        for ib = 1:nblock
            eb = S.E_blocks{ib};
            if ~isequal(size(eb), size(E_blocks_all(:,:,:,1,1)))
                e_ok = false;
                break;
            end
            E_blocks_all(:,:,:,ib,iu) = eb;
        end
        if ~e_ok
            status_msg(iu) = "SKIP(E_blocks shape mismatch): " + string(fpath);
            continue;
        end

        success_mask(iu) = true;
        status_msg(iu) = "OK";

    catch ME
        status_msg(iu) = "FAILED: " + string(ME.message);
        continue;
    end
end

if ~first_ok
    error('No valid U_*.mat file found in %s', result_dir);
end

%% 3) Save aggregated sigma
save(sigma_out_file, ...
    'sigma_all', 'sigma_blocks_all', ...
    'U_list', 'Eph_list', 'block_names_run', ...
    'success_mask', 'status_msg', 'Nw', 'nblock', ...
    '-v7.3');

%% 4) Save aggregated E_blocks
save(eblocks_out_file, ...
    'E_blocks_all', ...
    'U_list', 'block_names_run', ...
    'success_mask', 'status_msg', 'nblock', ...
    '-v7.3');

%% 5) Summary print
fprintf('\nDone. success = %d / %d\n', nnz(success_mask), numU);
fprintf('Sigma file  : %s\n', sigma_out_file);
fprintf('Eblocks file: %s\n', eblocks_out_file);
