function [H_full, Hx_full, Hy_full, Hz_full] = get_ham_kp_RhG_valley(kpoint, pars)
% 构建包含K和K'能谷的完整哈密顿量
% 输出维度: 4N x 4N (层空间2N × 能谷空间2)
% 基矢排列: (A1_K, B1_K, ..., AN_K, BN_K, A1_K', B1_K', ..., AN_K', BN_K')

    % 获取K谷哈密顿量 (xi = +1)
    pars_K = pars;
    pars_K.xi = +1;
    [HK, HKx, HKy, HKz] = MTB.ham.get_ham_kp_RhG(kpoint, pars_K);
    
    % 获取K'谷哈密顿量 (xi = -1)
    pars_Kp = pars;
    pars_Kp.xi = -1;
    [HKp, HKpx, HKpy, HKpz] = MTB.ham.get_ham_kp_RhG(kpoint, pars_Kp);
    
    % 组装成块对角矩阵
    dim = 2 * pars.N;
    H_full  = blkdiag(HK, HKp);
    Hx_full = blkdiag(HKx, HKpx);
    Hy_full = blkdiag(HKy, HKpy);
    Hz_full = blkdiag(HKz, HKpz);  % 仍然是零矩阵
end