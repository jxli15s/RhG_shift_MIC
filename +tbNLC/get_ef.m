function ef = get_ef(g, pars, model, nbands, knum)
%GET_EF
%==========================================================================
% 基于给定模型在指定 k 网格上计算本征值并估算费米能。
%
% 与原主脚本中的 get_ef 行为保持一致：
%   - 固定在 [-0.1,0.1] x [-0.1,0.1] 范围取样
%   - 使用半填充 u = 0.5
%==========================================================================

    arguments
        g
        pars struct
        model function_handle
        nbands (1,1) double {mustBeInteger, mustBePositive}
        knum (1,1) double {mustBeInteger, mustBePositive}
    end

    kxline = [-0.1, 0.1];
    kyline = [-0.1, 0.1];
    u = 0.5;

    [Kx, Ky, Kz] = g.get_Bulk2Dkmesh(kxline, kyline, knum);
    [~, Enk] = MTB.ham.get_bulk_plane_kp(pars, model, nbands, Kx, Ky, Kz);
    ef = tbNLC.calculate_ef(Enk(:), u);
end

