function [Ham_int, Kx_int, Ky_int, dHdkx_int, dHdky_int] = ...
         get_dHdk_2D(Ham, Kx, Ky)
%GET_H_AND_DHDK_FROM_MESH
%   从给定的 H(k)、Kx、Ky（可以是六角/斜BZ）数值上计算
%   直角坐标下的 dH/dkx, dH/dky，并且只保留内部 (N1-2)x(N2-2) 网格。
%
% 输入：
%   Ham : [n, n, N1, N2]  已经在整个k-mesh上的哈密顿量
%   Kx  : [N1, N2]        每个网格点的 kx 坐标
%   Ky  : [N1, N2]        每个网格点的 ky 坐标
%
% 输出（内部网格）：
%   Ham_int   : [n, n, N1-2, N2-2]   内部点的哈密顿量
%   Kx_int    : [N1-2, N2-2]         内部点的 kx
%   Ky_int    : [N1-2, N2-2]         内部点的 ky
%   dHdkx_int : [n, n, N1-2, N2-2]   直角坐标下的 dH/dkx
%   dHdky_int : [n, n, N1-2, N2-2]   直角坐标下的 dH/dky
%
% 说明：
%   - 网格可以是六角/斜的，只要在索引方向上是规则的即可。
%   - 程序自动从 Kx,Ky 的局部差分中提取两个方向的 k 步长向量，
%     然后通过一个 2x2 的线性变换得到直角坐标导数。
%   - 边界 i=1,N1 或 j=1,N2 不参与差分，直接丢弃。

    % ---------- 基本尺寸 ----------
    [n, ~, N1, N2] = size(Ham);

    if N1 < 3 || N2 < 3
        error('N1和N2都必须 >= 3，才能做中心差分并保留内部网格。');
    end

    % ---------- 1. 从 Kx,Ky 提取两个索引方向的步长向量 ----------
    % 选一个内部点作为参考点（一般 2,2 就可以）
    i0 = 2; 
    j0 = 2;

    % 索引方向 1：沿 i 方向前后各取一步
    dk1x = (Kx(i0+1,j0) - Kx(i0-1,j0)) / 2;
    dk1y = (Ky(i0+1,j0) - Ky(i0-1,j0)) / 2;

    % 索引方向 2：沿 j 方向前后各取一步
    dk2x = (Kx(i0,j0+1) - Kx(i0,j0-1)) / 2;
    dk2y = (Ky(i0,j0+1) - Ky(i0,j0-1)) / 2;

    % 组成 2x2 矩阵 C，并取逆：grad_k = C^{-1} * grad_kappa
    C = [dk1x, dk1y;
         dk2x, dk2y];

    if abs(det(C)) < 1e-12
        error('C 的行列式太小，说明两个索引方向在 kx-ky 平面里几乎共线，不能用于反演。');
    end

    invC = inv(C);
    c11 = invC(1,1); c12 = invC(1,2);
    c21 = invC(2,1); c22 = invC(2,2);

    % ---------- 2. 只保留内部 (2..N1-1, 2..N2-1) ----------
    N1_int = N1 - 2;
    N2_int = N2 - 2;

    Ham_int   = zeros(n, n, N1_int, N2_int);
    dHdkx_int = zeros(n, n, N1_int, N2_int);
    dHdky_int = zeros(n, n, N1_int, N2_int);

    Kx_int = zeros(N1_int, N2_int);
    Ky_int = zeros(N1_int, N2_int);

    % ---------- 3. 对内部点做索引方向差分 + 线性变换 ----------
    for I = 1:N1_int
        i = I + 1;           % 原网格索引：2..N1-1

        for J = 1:N2_int
            j = J + 1;       % 原网格索引：2..N2-1

            % 本点的 H 和 k 坐标
            Hc = Ham(:,:,i,j);
            Ham_int(:,:,I,J) = Hc;
            Kx_int(I,J) = Kx(i,j);
            Ky_int(I,J) = Ky(i,j);

            % --- (1) 在索引方向 κ1,κ2 上做中心差分 ---
            % κ1 对应 i 方向
            H_plus_i  = Ham(:,:,i+1,j);
            H_minus_i = Ham(:,:,i-1,j);
            dH_dkappa1 = (H_plus_i - H_minus_i) / 2;

            % κ2 对应 j 方向
            H_plus_j  = Ham(:,:,i,j+1);
            H_minus_j = Ham(:,:,i,j-1);
            dH_dkappa2 = (H_plus_j - H_minus_j) / 2;

            % 强制厄米（数值噪声会导致微小非厄米）
            dH_dkappa1 = 0.5*(dH_dkappa1 + dH_dkappa1');
            dH_dkappa2 = 0.5*(dH_dkappa2 + dH_dkappa2');

            % --- (2) 用 C^{-1} 把 grad_kappa 变成 grad_k ---
            dHdkx = c11*dH_dkappa1 + c12*dH_dkappa2;
            dHdky = c21*dH_dkappa1 + c22*dH_dkappa2;

            dHdkx_int(:,:,I,J) = dHdkx;
            dHdky_int(:,:,I,J) = dHdky;
        end
    end
end
