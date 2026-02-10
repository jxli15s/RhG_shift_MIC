function sigma_k = get_sigma_at_k(hamall, omega_eV, eta, fermi)
% GET_SIGMA_AT_K  计算单个 k 点处的线性光学电导 sigma_xx, sigma_yy(ω)
%
% 输入:
%   hamall : [nb, nb, 4]
%       hamall(:,:,1) = H(k)       —— 能量哈密顿量 (单位 eV)
%       hamall(:,:,2) = Hx(k)      —— 与 v_x 成正比的算符 (通常是 ∂H/∂k_x, 单位 eV·Å)
%       hamall(:,:,3) = Hy(k)      —— 与 v_y 成正比的算符
%       hamall(:,:,4) = Hz(k)      —— 与 v_z 成正比的算符
%
%   omega_eV : [1, Nw]
%       频率网格，单位 eV (这里是 ħω，而不是 ω 本身)
%
%   eta : 标量
%       洛伦兹展宽的宽度，单位 eV
%
%   fermi : 标量
%       费米能级，单位 eV
%
% 输出:
%   sigma_k : [2, Nw]
%       对应这个 k 点的频率依赖电导 (实部/虚部取决于 Llor 的定义和 H 的输入)
%       sigma_k(1, :) ≈ σ_xx(k, ω)
%       sigma_k(2, :) ≈ σ_yy(k, ω)
%
% 备注:
%   - 本例中 f(n) 使用 T→0 极限的 0/1 阶跃函数，而不是真正的 Fermi-Dirac。
%   - Llor 采用洛伦兹形式 δ(ħω-ΔE) ≈ (1/π)*η/[(ħω-ΔE)^2+η^2]，单位约为 1/eV。
%   - prefactor 目前写成 (π e^2 / ħ) * (f_n - f_m)/ΔE，与标准 Kubo 公式兼容。
%   - 这里未显式把 eV 转为 J（有可能需要整体单位换算，这里沿用你现有的约定）。
%

    % ========= 物理常数 =========
    eC   = 1.602176634e-19;   % 元电荷, 单位 C
    % h  = 6.62607015e-34;    % Planck 常数 (J·s) —— 当前没用到，可保留也可删掉
    hbar = 1.054571817e-34;   % 约化 Planck 常数 (J·s)

    % 如果将来想真正用 Fermi-Dirac，可以用下面这两行; 现在只用到 T 定义，没有真正用 beta
    kB   = 8.617333262e-5;    % Boltzmann 常数 (eV/K)
    T    = 0.001;             % 温度 (K)，几乎是 T→0
    beta = 1 / max(kB*T, 1e-12); %#ok<NASGU>  % 这里beta暂时没被使用

    A_to_m = 1e-10;           %#ok<NASGU>  % Å → m 的换算系数，这里也暂时没显式用
    w_k = 1;                  % 当前 k 点的权重（在外层积分时乘上 Brillouin zone 权重）

    % ========= 定义洛伦兹化的 δ 函数 =========
    % Lorentzian: δ(ħω-ΔE) ≈ (1/π)*η / [ (ω-ΔE)^2 + η^2 ], 其中 ω, ΔE 用 eV 记
    Llor = @(Eph, DE, etaL) (1/pi) * etaL ./ ((Eph - DE).^2 + etaL^2); % 输出 ~ 1/eV

    % ========= 取出哈密顿量和"速度算符"矩阵 =========
    H  = hamall(:,:,1);       % H(k)  (eV)
    Hx = hamall(:,:,2);       % 与 v_x 成正比的算符
    Hy = hamall(:,:,3);       % 与 v_y 成正比的算符
    Hz = hamall(:,:,4);       %#ok<NASGU> % 当前没用到 z 方向电导，可以以后扩展 σ_zz

    Nw = numel(omega_eV);     % 能量/频率网格点数
    nb = size(H, 1);          % 带数

    % ========= 本征化 H(k) =========
    % 为了避免数值误差，先强制 H 变成严格厄米
    H = (H + H') / 2;

    % H * U = U * D
    [U, D] = eig(H);
    evals = real(diag(D));    % 本征值 (eV)，取实部
    [evals, idx] = sort(evals(:));  % 排序
    U = U(:, idx);            % 对应排序后的本征矢

    % ========= 费米占据 (T→0 极限，用阶跃函数) =========
    % 若想用真正的 Fermi-Dirac: f = 1./(1 + exp((evals - fermi)*beta));
    f = double(evals <= fermi);    % evals <= μ: 占据态 (1), 其余为 0

    % ========= 把"速度算符"旋到本征基 =========
    % V_ij^α = <i|v_α|j> ≈ <i| (1/ħ) ∂H/∂k_α |j>
    Up = U';                          % U^\dagger
    V = zeros(nb, nb, 3);             % V(:,:,1/2/3) 分别对应 x,y,z 分量
    V(:,:,1) = Up * Hx * U;           % V^x
    V(:,:,2) = Up * Hy * U;           % V^y
    V(:,:,3) = Up * Hz * U;           % V^z (当前未用)

    Mxx_all = abs(V(:,:,1)).^2;   % nb x nb
    Myy_all = abs(V(:,:,2)).^2;

    % ========= 初始化 k 点电导 =========
    % sigma_k(1,:) -> σ_xx(k, ω)
    % sigma_k(2,:) -> σ_yy(k, ω)
    sigma_k = zeros(2, Nw);

    % ========= 双和: ∑_{n,m} (f_n - f_m) ... =========
    for n = 1:nb
        fn = f(n);
        for m = 1:nb
            fm = f(m);

            dE = evals(m) - evals(n);  % 能量差 ΔE = E_m - E_n  (eV)
            fv = fn - fm;              % Fermi 因子差 (f_n - f_m)

            % 跳过: (1) 能量差太小容易数值发散; (2) f_n ≈ f_m 时贡献极小/为 0
            if abs(dE) <= 1e-8 || abs(fv) < 1e-8
                continue;
            end

            % 速度矩阵元: <n|v_i|m>
            % Mx = V(n,m,1);
            % My = V(n,m,2);
            % 对角元: |<n|v_x|m>|^2, |<n|v_y|m>|^2
            % Mxx = Mx * conj(Mx);
            % Myy = My * conj(My);
            Mxx = Mxx_all(n,m);
            Myy = Myy_all(n,m);


            % 洛伦兹化的 δ(ħω - ΔE)
            % L(ω) ~ 1/eV, 大小为 [Nw, 1]
            L = (Llor(omega_eV, dE, eta)).';   % [Nw,1] → 转置成 [1,Nw]

            % prefactor: (π e^2 / ħ) * w_k * (f_n - f_m) / ΔE
            % 单位大致为: C^2 / (J·s) × 1/eV ≈ (S/m)*something，根据你整体单位约定
            pref = (pi * eC^2 / hbar) * w_k * (fv / dE);

            % 累加到 σ_xx(k,ω), σ_yy(k,ω)
            sigma_k(1,:) = sigma_k(1,:) + pref * Mxx .* L;
            sigma_k(2,:) = sigma_k(2,:) + pref * Myy .* L;
        end
    end
end