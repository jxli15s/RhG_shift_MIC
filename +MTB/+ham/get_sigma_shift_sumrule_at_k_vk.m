function nhc = get_sigma_shift_sumrule_at_k_vk(en,Dk, nbands, Egrid, egam, fermi)
% hamshift  —— 计算基于 sum-rule 的移位电流核谱 (3×3×3×Nw)
%
% 输入
%   hamall : [nb, nb, 4]  其中 hamall(:,:,1)=H(k)，2/3/4 分别是 Hx, Hy, Hz（单位 eV 与 eV·Å）
%   estart : 标谱能量起点 (eV)
%   enum   : 能量网格点数
%   estep  : 能量步进 (eV)
%   egam   : 高斯展宽 (eV)
%   fermi  : 费米能级 (eV)
%
% 输出
%   nhc : [3,3,3,enum]，每个能量点的核谱（与 Fortran 版本一致）
%
% 备注
% - 这是"能量规范"的实现，与你的 Fortran 代码逐项等价。
% - 这里的 gaussian(x,egam) 采用归一化 (1/(sqrt(pi)*egam))*exp(-(x/egam)^2)。

    enum = numel(Egrid);
    nb = nbands;
    evals=en;

    % 费米占据
    f = double(evals <= fermi);        % 0/1 占据

    % 初始化
    berry = zeros(3,3,3,nb,nb);        % 对应 Fortran 的 berry(i,j,k,n,m)
    nhc   = zeros(3,3,3,enum);
    gap   = evals.' - evals;           % gap(m,n) = E_m - E_n

    % —— 主循环：构造 berry(i,j,k,n,m)
    % 等价于你的 Fortran 双带循环 + 三个方向的求和
    deg_eps = 1e-5;                    % 简并判据（eV）

    for n = 1:nb
        for m = 1:nb
            if abs(f(n) - f(m)) > 0.05    % 跨带跃迁
                Enm  = evals(n) - evals(m);
                if abs(Enm) < deg_eps
                    continue;
                end
                invEnm   = 1.0 / Enm;
                invEnm2  = invEnm * invEnm;

                for i = 1:2
                    for j = 1:2
                        for k = 1:2

                            % Σ_p 部分
                            tmp = 0.0 + 0.0i;
                            for p = 1:nb
                                if (abs(evals(p)-evals(m)) > deg_eps) && ...
                                   (abs(evals(p)-evals(n)) > deg_eps)
                                    tmp = tmp ...
                                        +  Dk(n,p,k) * Dk(p,m,i) / (evals(p)-evals(m)) ...
                                        -  Dk(n,p,i) * Dk(p,m,k) / (evals(n)-evals(p));
                                end
                            end

                            % 直写 Fortran 里的那一行
                            termA = Dk(n,m,k) * (Dk(n,n,i) - Dk(m,m,i)) ...
                                  + Dk(n,m,i) * (Dk(n,n,k) - Dk(m,m,k));
                            factor = Dk(m,n,j) * ( invEnm2 * termA + invEnm * tmp );

                            berry(i,j,k,n,m) = berry(i,j,k,n,m) ...
                                - 1.0 * (f(n)-f(m)) * imag( invEnm * factor );
                        end
                    end
                end
            end
        end
    end

    % —— 频谱卷积（高斯展宽）
    for ei = 1:enum
        ep=Egrid(ei);
        % 只对 |gap-ep| < 5*egam 的 (m,n) 做加和（与 Fortran 一致）
        mask = (abs(gap - ep) < 5.0*egam) & (abs(f.' - f) > 0.05);
        [mm, nn] = find(mask);       % m,n 的索引对
        for t = 1:numel(mm)
            m = mm(t); n = nn(t);
            w = gaussian(gap(m,n) - ep, egam);
            nhc(:,:,:,ei) = nhc(:,:,:,ei) + w * berry(:,:,:,n,m);
        end
    end
end

% ---- 归一化高斯，与 Fortran 的 gaussian(gap-ep, egam) 对应
function g = gaussian(x, gamma)
    g = (1.0/(sqrt(pi)*gamma)) * exp(-(x./gamma).^2);
end