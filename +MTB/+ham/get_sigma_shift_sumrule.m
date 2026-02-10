function sigma = get_sigma_shift_sumrule(Kx,Ky,Kz, pars, get_hk, Egrid, etaE, Te, mu, t_eff)
% 与 cell 版等价，但用数值数组：
% V(:,:,a)        [nb,nb,3]
% Delta(:,:,a)    [nb,nb,3]
% r(:,:,a)        [nb,nb,3]
% r_cov(:,:,a,b)  [nb,nb,3,3]
% 6.62607015×10−34 J.s Planck's constants
% Inputs:
%   kmesh.klist [Nk,2], kmesh.weights [Nk,1] (包含 (2π)^-2 因子)
%   Egrid  [Nw,1]  —— 频谱能量网格 (eV)
%   etaE   —— 能量展宽 (eV)
%   Te(K), mu(eV), t_eff(Å) 可选
e = 1.602176634e-19; %ev to J
h = 6.582119569e-16; % ev.s Reduced Planck's constants
% fprintf('Test')
knum=size(Kx,1);
Nk = knum^2;
Nw = numel(Egrid);
wk=1/Nk;

% sigma_acc = zeros(3,3,3,Nw, 'like', 1i);  % 先累计无前因子的谱和
sigma_k = zeros(3,3,3,Nw,Nk, 'like', 1i);


parfor ik = 1:Nk
    % fprintf('Test\n')
    [i,j]=ind2sub([knum,knum],ik);
    % fprintf('Test',i,j)
    q =[Kx(i,j),Ky(i,j),Kz(i,j)];
    [E, U, Hx, Hy, Hz]=slovek(q,pars,get_hk);
    nb = numel(E);
    [Em,En] = meshgrid(E,E); dE = (Em-En);                 % [nb,nb]
    dE_vec = dE(:);                                              % [nb^2,1]
    % 
    % 旋到本征基：Ha_eig (eV·Å) —— 没有 ħ
    Ha = zeros(nb,nb,3);
    Ha(:,:,1) = U'*Hx*U;
    Ha(:,:,2) = U'*Hy*U;
    Ha(:,:,3) = U'*Hz*U;
    % fprintf('Test ik : %d',ik)
    % Δ~^a_{nm} = H^a_{nn} - H^a_{mm} （能量版差）
    Dlt = zeros(nb,nb,3);
    for a=1:3
        d = real(diag(Ha(:,:,a)));
        Dlt(:,:,a) = d.' - d;             % eV·Å
    end

    % % r^a_{nm} = H^a_{nm} / (i ΔE_{mn}) ; 置零对角
    r = zeros(nb,nb,3);
    for a=1:3
        E_floor = min(1e-6, 1e-3*etaE);        % 也可固定为 1e-6 eV
        den = 1i * dE;                          % 目标分母 iΔE
        mask = (abs(dE) < E_floor) & ~eye(nb);  % 仅限 m≠n 的近零项
        den(mask) = 1i * sign(dE(mask)) .* E_floor;
        Ra = Ha(:,:,a) ./ den;
        Ra(1:nb+1:end) = 0;
        Ra = 0.5*(Ra - Ra');                    % 可选：强制反厄米，清数值噪声
        r(:,:,a)=Ra
    end
    % 
    % % 线性 k·p → W^{ab}=⟨∂_b H_a⟩ = 0；构造能量版 r^{a;b}
    % % r^{a;b}_{nm} = (i/ΔE_{mn}) * [ H^a_{nm} Δ~^b_{nm} + H^b_{nm} Δ~^a_{nm}
    %                               + Σ_p ( H^a_{np} H^b_{pm}/ΔE_{pm} - H^b_{np} H^a_{pm}/ΔE_{np} ) ]
    r_cov = zeros(nb,nb,3,3);
    E_floor = min(1e-6, 1e-3*etaE);   % 很小的能量下限（eV），<< 展宽
    Ieye = ~eye(nb);
    % fprintf('Test ik: %d\n',ik)
    diagMask=eye(nb)~=0;
    offMask=~diagMask;

    for a = 1:3
        Ha_ab = Ha(:,:,a);
        ra_ab=r(:,:,a)
        for b = 1:3
            Hb_ab = Ha(:,:,b);
            rb_ab=r(:,:,b)
            % term1 = H^a_{nm} Δ~^b_{nm} + H^b_{nm} Δ~^a_{nm}
            
            % denom_nm = dE;                         % ε_nm
            % denom_nm(diagMask) = Inf;
            % small_nm = offMask & (abs(denom_nm) < E_floor);
            % denom_nm(small_nm) = sign(denom_nm(small_nm)).*E_floor;
            % term1 = Ha_ab .* Dlt(:,:,b) + Hb_ab .* Dlt(:,:,a);
            % term1 = term1./denom_nm;
            term1 = ra_ab .* Dlt(:,:,b) + rb_ab .* Dlt(:,:,a);
            S = zeros(nb,nb);

            for p = 1:nb
                denom_pm = dE(p,:);    % 1×nb,  ΔE_{pm}
                denom_np = dE(:,p);    % nb×1,  ΔE_{np}

                % 排除 m=p / n=p
                denom_pm(p) = Inf;     % 使 1./denom_pm(p)=0
                denom_np(p) = Inf;     % 使 1./denom_np(p)=0

                % 近简并 floor（同尺寸一维掩码，保持符号）
                idx_pm = abs(denom_pm) < E_floor;                 % 1×nb
                % denom_pm(idx_pm) = sign(denom_pm(idx_pm)) .* E_floor;
                denom_pm(idx_pm) = Inf;

                idx_np = abs(denom_np) < E_floor;                 % nb×1
                % denom_np(idx_np) = sign(denom_np(idx_np)) .* E_floor;
                denom_np(idx_np) = Inf;

                % 组装
                T1 = (Ha_ab(:,p) * Hb_ab(p,:)) ./ denom_pm;       % 隐式行广播
                T2 = (Hb_ab(:,p) * Ha_ab(p,:)) ./ denom_np;       % 隐式列广播
                S  = S + T1 - T2;
            end


            % rcov = (i/ΔE) * (term1 + S) ，对角设零，并做一次反厄米化消噪
            denom_nm = dE;
            % 对 n=m 的对角置 Inf，使得 1/denom=0
            denom_nm(Ieye) = denom_nm(Ieye);                    % off-diag 不动
            denom_nm(~Ieye) = Inf;                              % diag → Inf

            % 对近简并（但 n≠m）做 floor（保符号）
            mask_small = (abs(denom_nm) < E_floor) & Ieye;
            % denom_nm(mask_small) = sign(denom_nm(mask_small)) * E_floor;
            denom_nm(mask_small) = Inf;

            rcov = (1i ./ denom_nm) .* (term1 + S);
            rcov(~Ieye) = 0;                    % 明确置零对角
            rcov = 0.5*(rcov - rcov');          % 反厄米化（理论上应当反厄米）

            r_cov(:,:,a,b) = rcov;              % 单位：Å
        end
    end


    % ---- Fermi factors
    kB = 8.617333262e-5;          % eV/K
    if Te > 0
        f = 1./(1 + exp((E - mu)/(kB*Te)));   % E: [nb,1]
    else
        f = double(E <= mu);                  % T=0 极限
    end
    fnm = f.' - f;                            % [nb,nb]
    % fnm_vec = fnm(:);                         % [nb^2,1]

    % ---- Lorentzian broadening δ_η(E - ΔE)
    % Egrid: [Nw,1] provided by caller; dE_vec: [nb^2,1]
    % 新版 MATLAB（R2016b+）可直接写隐式扩展：
    % DeltaE = (etaE/pi) ./ ((Egrid(:).' - dE_vec).^2 + etaE^2);  % [Nw, nb^2]
    % Nw x nb^2
    % ---- 频谱卷积（每个 a,b,c）
% 先保证 DeltaE 是 Nw x nb^2
DeltaE = (etaE/pi) ./ ( bsxfun(@minus, Egrid(:), dE_vec.').^2 + etaE^2 );  % [Nw, nb^2]

sigma_loc = zeros(3,3,3,Nw, 'like', 1i);

for a=1:3
  for b=1:3
    for c=1:3
      % 注意用非共轭转置 .'
      M = fnm .* ( r(:,:,b).' .* r_cov(:,:,c,a) + r(:,:,c).' .* r_cov(:,:,b,a) );  % [nb,nb]
      Mm   = reshape(M, [], 1);                       % [nb^2,1]
      spec = DeltaE * Mm;                              % [Nw,1]
      sigma_loc(a,b,c,:) = wk * reshape(spec, 1,1,1,[]);  % 1×1×1×Nw
    end
  end
end

    sigma_k(:,:,:,:,ik) = sigma_loc;
end

% 循环外求和 + 乘全局前因子
sigma_acc = sum(sigma_k, 5);                 % [3,3,3,Nw]
% K = (-1i*pi*e) / (2*h);                    % A/V^2 per eV
K = e^3/h^2*10^6*pi; %A/V^2 per eV ~~~ e^3/h/ev
sigma = K * sigma_acc;                       % 2D 单位 A·m/V^2
% % if exist('t_eff','var') && ~isempty(t_eff), sigma = sigma / t_eff; end
end

function [E,Psik, Hx, Hy, Hz]=slovek(k,pars,get_hk)
    [hk,Hx,Hy,Hz]=get_hk(k,pars);
    [V,D]=eig(hk);
    [E,ind]=sort(diag(D));
    Psik=V(:,ind);
end
