function [nhc_nsc, nhc_msc] = get_sigma_shift_sumrule_at_k_vk_nsc_msc(en, Dk, nbands, Egrid, egam, fermi)
%GET_SIGMA_SHIFT_SUMRULE_AT_K_VK_NSC_MSC
% Single-k sum-rule kernel for
%   1) normal shift current (NSC, j<->k symmetric channel)
%   2) magnetic shift current (MSC, j<->k antisymmetric channel)
%
% Inputs:
%   en    : [nb,1] eigen energies at this k (eV)
%   Dk    : [nb,nb,3] matrix elements in eigen basis (typically vx,vy,vz)
%   nbands: scalar nb
%   Egrid : [Nw,1] or [1,Nw] photon energy grid (eV)
%   egam  : scalar Gaussian broadening (eV)
%   fermi : scalar Fermi level (eV)
%
% Outputs:
%   nhc_nsc, nhc_msc : [3,3,3,Nw]

    enum = numel(Egrid);
    nb = nbands;
    evals = en(:);

    f = double(evals <= fermi);

    berry_nsc = zeros(3,3,3,nb,nb);
    berry_msc = zeros(3,3,3,nb,nb);
    nhc_nsc = zeros(3,3,3,enum);
    nhc_msc = zeros(3,3,3,enum);
    gap = evals.' - evals;

    deg_eps = 1e-5;

    for n = 1:nb
        for m = 1:nb
            if abs(f(n) - f(m)) <= 0.05
                continue;
            end

            Enm = evals(n) - evals(m);
            if abs(Enm) < deg_eps
                continue;
            end

            invEnm = 1.0 / Enm;
            invEnm2 = invEnm * invEnm;

            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        tmp1 = 0.0 + 0.0i;
                        tmp2 = 0.0 + 0.0i;

                        for p = 1:nb
                            if (abs(evals(p)-evals(m)) > deg_eps) && (abs(evals(p)-evals(n)) > deg_eps)
                                tmp1 = tmp1 + Dk(n,p,k) * Dk(p,m,i) / (evals(p)-evals(m)) ...
                                            - Dk(n,p,i) * Dk(p,m,k) / (evals(n)-evals(p));
                                tmp2 = tmp2 + Dk(n,p,j) * Dk(p,m,i) / (evals(p)-evals(m)) ...
                                            - Dk(n,p,i) * Dk(p,m,j) / (evals(n)-evals(p));
                            end
                        end

                        termA1 = Dk(n,m,k) * (Dk(n,n,i) - Dk(m,m,i)) ...
                               + Dk(n,m,i) * (Dk(n,n,k) - Dk(m,m,k));
                        factor1 = Dk(m,n,j) * (invEnm2 * termA1 + invEnm * tmp1);

                        termA2 = Dk(n,m,j) * (Dk(n,n,i) - Dk(m,m,i)) ...
                               + Dk(n,m,i) * (Dk(n,n,j) - Dk(m,m,j));
                        factor2 = Dk(m,n,k) * (invEnm2 * termA2 + invEnm * tmp2);

                        raw1 = -1.0 * (f(n)-f(m)) * imag(invEnm * factor1);
                        raw2 = -1.0 * (f(n)-f(m)) * imag(invEnm * factor2);

                        berry_nsc(i,j,k,n,m) = berry_nsc(i,j,k,n,m) + 0.5 * (raw1 + raw2);
                        berry_msc(i,j,k,n,m) = berry_msc(i,j,k,n,m) + 0.5 * (raw1 - raw2);
                    end
                end
            end
        end
    end

    for ei = 1:enum
        ep = Egrid(ei);
        mask = (abs(gap - ep) < 5.0*egam) & (abs(f.' - f) > 0.05);
        [mm, nn] = find(mask);
        for t = 1:numel(mm)
            m = mm(t);
            n = nn(t);
            w = gaussian(gap(m,n) - ep, egam);
            nhc_nsc(:,:,:,ei) = nhc_nsc(:,:,:,ei) + w * berry_nsc(:,:,:,n,m);
            nhc_msc(:,:,:,ei) = nhc_msc(:,:,:,ei) + w * berry_msc(:,:,:,n,m);
        end
    end
end

function g = gaussian(x, gamma)
    g = (1.0/(sqrt(pi)*gamma)) * exp(-(x./gamma).^2);
end
