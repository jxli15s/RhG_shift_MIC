function [g, Q, F] = quantum_metric_general(obj, k, dk_list, Uk, Ek, band_index, delta)
% QUAN_METRIC  Compute quantum metric and Berry curvature tensor for a band.
%
%   [g, Q, F] = quan_metric(obj, k, dk_list, Uk, Ek, band_index, delta)
%
% 输入参数：
%   obj        : 拥有 get_hk(k) 方法的对象，返回 H(k) (nb x nb)
%   k          : 1 x D 的 k 点（行向量） 
%   dk_list    : Ndir x D 的矩阵，每一行是一个在 k 空间的微小位移方向
%                例如 2D 情况：dk_list = [b1/knum; b2/knum];
%   Uk         : nb x nb 的本征矢矩阵，Uk(:,n) = |u_n(k)>
%   Ek         : nb x 1 的本征值向量，对应 Uk
%   band_index : 目标能带编号（非简并）
%   delta      : 规避能隙太小的正数阈值（单位同 Ek^2）
%
% 输出参数：
%   g          : Ndir x Ndir 的量子度规张量 g_ij = Re(Q_ij)
%   Q          : Ndir x Ndir 的量子几何张量 Q_ij
%   F          : Ndir x Ndir 的 Berry 曲率张量 F_ij = 2 Im(Q_ij)
%
% 说明：
%   - 这里使用的是带间公式：
%       Q_ij = sum_{m != n} <n| H_i |m> <m| H_j |n> / (E_m - E_n)^2
%     其中 H_i = ∂H/∂k_i 用有限差分近似。
%   - dk_list 的每一行定义一个方向 i，对应 H_i。
%   - 对 2D 系统，通常 Ndir = 2，此时：
%       g_xx = g(1,1), g_yy = g(2,2), g_xy = g(1,2)
%       Ω_z  = F(1,2)   (因为 F_12 = -F_21 是唯一独立分量)

    % 保证 k 是行向量
    k = k(:).';

    nb = size(Uk, 1);             % 轨道/能带数
    Ndir = size(dk_list, 1);      % 定义的方向个数

    % 目标带的本征态与本征值
    u0 = Uk(:, band_index);       % |u_n>
    E0 = Ek(band_index);          % E_n

    % 去掉目标带，留下 m != n 的其它能带
    mask = true(nb, 1);
    mask(band_index) = false;
    Uk_others = Uk(:, mask);      % nb x (nb-1)
    Ek_others = Ek(mask);         % (nb-1) x 1

    % 能量差平方 (E_n - E_m)^2
    dE = (E0 - Ek_others.').^2;   % 1 x (nb-1)

    % 对太小的能隙做正则化，避免除以 0
    small_gap = (dE < delta);
    dE(small_gap) = dE(small_gap) + delta;

    % Uk_scaled(:,m) = |u_m> / (E_n - E_m)^2
    % 用 bsxfun 兼容老版本 MATLAB
    Uk_scaled = bsxfun(@rdivide, Uk_others, dE);   % nb x (nb-1)

    % H(k) 在当前点，只算一次
    H0 = obj.get_hk(k);

    % 准备存放 <n|H_i|m>/(ΔE^2) 和 <m|H_i|n>
    M = nb - 1;
    V = zeros(Ndir, M);   % V(i,m) = <n|H_i|m> / (ΔE)^2
    W = zeros(M, Ndir);   % W(m,j) = <m|H_j|n>

    for a = 1:Ndir
        dk = dk_list(a, :);
        % 有需要可以改成中心差分：
        % H_a = (obj.get_hk(k + dk) - obj.get_hk(k - dk)) / (2*norm(dk));
        H_a = (obj.get_hk(k + dk) - H0) ./ norm(dk);  % ∂H/∂k_a

        % <n|H_a|m> / ΔE^2  = u0' * H_a * (|m>/ΔE^2)
        V(a, :) = (u0' * H_a * Uk_scaled);

        % <m|H_a|n> = (Uk_others' * H_a * u0)(m)
        W(:, a) = (Uk_others' * H_a * u0);
    end

    % 量子几何张量 Q_ij = sum_m V(i,m) * W(m,j)
    Q = V * W;           % Ndir x Ndir 复数矩阵
    g = real(Q);         % 量子度规（对称）
    F = 2 * imag(Q);     % Berry 曲率张量（反对称）

end
