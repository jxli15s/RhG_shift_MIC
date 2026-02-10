function H_blocks = extract_spinvalley_blocks(H, N)
%EXTRACT_SPINVALLEY_BLOCKS
%==========================================================================
% 从完整哈密顿量中提取 (K up, K' up, K down, K' down) 四个块。
%
% 输入:
%   H: nkx x nky x (8N) x (8N)
%   N: 层数
%
% 输出:
%   H_blocks.K_up, H_blocks.Kp_up, H_blocks.K_dn, H_blocks.Kp_dn
%   每个块尺寸: nkx x nky x (2N) x (2N)
%==========================================================================

    arguments
        H {mustBeNumeric}
        N (1,1) {mustBeInteger, mustBePositive}
    end

    [nkx, nky, ~, ~] = size(H);
    dim_layer = 2 * N;

    H_blocks = struct();
    H_blocks.K_up = zeros(nkx, nky, dim_layer, dim_layer, 'like', H);
    H_blocks.Kp_up = zeros(nkx, nky, dim_layer, dim_layer, 'like', H);
    H_blocks.K_dn = zeros(nkx, nky, dim_layer, dim_layer, 'like', H);
    H_blocks.Kp_dn = zeros(nkx, nky, dim_layer, dim_layer, 'like', H);

    idx_K_up = 1:dim_layer;
    idx_Kp_up = dim_layer + (1:dim_layer);
    idx_K_dn = 2*dim_layer + (1:dim_layer);
    idx_Kp_dn = 3*dim_layer + (1:dim_layer);

    H_blocks.K_up = H(:, :, idx_K_up, idx_K_up);
    H_blocks.Kp_up = H(:, :, idx_Kp_up, idx_Kp_up);
    H_blocks.K_dn = H(:, :, idx_K_dn, idx_K_dn);
    H_blocks.Kp_dn = H(:, :, idx_Kp_dn, idx_Kp_dn);
end

