function nhc_mic = get_sigma_mic_sumrule_at_k_vk_reQ(en, Dk, nbands, Egrid, egam, fermi, positiveDE)
%GET_SIGMA_MIC_SUMRULE_AT_K_VK_REQ
% Single-k MIC sum-rule kernel using ReQ -> g.
%
%   g_{mn}^{bc} = Re[Q_{mn}^{bc}],
%   Q_{mn}^{bc} = r_{mn}^b r_{nm}^c,
%   r_{mn}^b = -i * Dk(m,n,b) / (E_m - E_n)
%
% MIC kernel form at single k:
%   K^{abc}_{mn} ~ f_nm * Delta^{a}_{mn} * (2 g^{bc}_{mn})
% with
%   Delta^{a}_{mn} = Dk(m,m,a) - Dk(n,n,a)
%
% Inputs:
%   en, Dk, nbands, Egrid, egam, fermi: same style as shift sum-rule kernels
%   positiveDE (optional, default=true): keep only E_m-E_n > 0 transitions
%
% Output:
%   nhc_mic : [3,3,3,Nw]

    if nargin < 7
        positiveDE = true;
    end

    evals = en(:);
    nb = nbands;
    enum = numel(Egrid);

    f = double(evals <= fermi);

    berry = zeros(3,3,3,nb,nb);
    nhc_mic = zeros(3,3,3,enum);
    gap = evals.' - evals; % gap(m,n) = E_m - E_n

    deg_eps = 1e-8;

    for n = 1:nb
        for m = 1:nb
            if m == n
                continue;
            end

            dE = evals(m) - evals(n); % E_m - E_n
            if abs(dE) < deg_eps
                continue;
            end

            if positiveDE && (dE <= 0)
                continue;
            end

            f_nm = f(n) - f(m);
            if abs(f_nm) <= 0.05
                continue;
            end

            invdE2 = 1.0 / (dE * dE);

            for a = 1:3
                Delta_a = real(Dk(m,m,a) - Dk(n,n,a));
                for b = 1:3
                    for c = 1:3
                        % 2 g^{bc}_{mn} = Re[ r^b_mn r^c_nm + r^c_mn r^b_nm ]
                        % with r from Dk / dE, this becomes:
                        g2_bc = real( ...
                            Dk(m,n,b) * Dk(n,m,c) + ...
                            Dk(m,n,c) * Dk(n,m,b) ) * invdE2;

                        berry(a,b,c,n,m) = berry(a,b,c,n,m) + f_nm * Delta_a * g2_bc;
                    end
                end
            end
        end
    end

    for ei = 1:enum
        ep = Egrid(ei);
        mask = (abs(gap - ep) < 5.0*egam) & (abs(f.' - f) > 0.05);
        if positiveDE
            mask = mask & (gap > 0);
        end

        [mm, nn] = find(mask);
        for t = 1:numel(mm)
            m = mm(t);
            n = nn(t);
            w = gaussian(gap(m,n) - ep, egam);
            nhc_mic(:,:,:,ei) = nhc_mic(:,:,:,ei) + w * berry(:,:,:,n,m);
        end
    end
end

function g = gaussian(x, gamma)
    g = (1.0/(sqrt(pi)*gamma)) * exp(-(x./gamma).^2);
end
