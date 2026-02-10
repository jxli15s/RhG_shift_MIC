function [sigma_nsc, sigma_msc, nhc_nsc_all, nhc_msc_all] = get_sigma_shift_sumrule_fast_nsc_msc( ...
    Enk, vxk, vyk, vzk, Egrid, egam, fermi, weight)
%GET_SIGMA_SHIFT_SUMRULE_FAST_NSC_MSC
% Fast k-mesh sum-rule evaluation for NSC and MSC using precomputed
% eigen energies and velocity-like matrix elements.
%
% Inputs:
%   Enk  : [Nkx,Nky,nb]
%   vxk  : [nb,nb,Nkx,Nky]
%   vyk  : [nb,nb,Nkx,Nky]
%   vzk  : [nb,nb,Nkx,Nky]
%   Egrid, egam, fermi, weight : same meaning as existing fast interface
%
% Outputs:
%   sigma_nsc, sigma_msc : [3,3,3,Nw]
%   nhc_nsc_all, nhc_msc_all (optional): [3,3,3,Nw,Nk]

    e = 1.602176634e-19;   % C
    h = 6.582119569e-16;   % eV*s
    Kconst = e / h * 1e6;

    [Nkx, Nky, nbands] = size(Enk);
    Nk = Nkx * Nky;
    Nw = numel(Egrid);

    sigma_nsc = zeros(3, 3, 3, Nw, 'like', 1i);
    sigma_msc = zeros(3, 3, 3, Nw, 'like', 1i);

    keep_kresolved = (nargout > 2);
    if keep_kresolved
        nhc_nsc_all = zeros(3, 3, 3, Nw, Nk, 'like', 1i);
        nhc_msc_all = zeros(3, 3, 3, Nw, Nk, 'like', 1i);
    else
        nhc_nsc_all = [];
        nhc_msc_all = [];
    end

    vxk = reshape(vxk, nbands, nbands, Nk);
    vyk = reshape(vyk, nbands, nbands, Nk);
    vzk = reshape(vzk, nbands, nbands, Nk);
    Enk = reshape(Enk, Nk, nbands);

    parfor ik = 1:Nk
        en = Enk(ik, :).';
        Dk = cat(3, vxk(:,:,ik), vyk(:,:,ik), vzk(:,:,ik));

        [nhc_nsc_k, nhc_msc_k] = MTB.ham.get_sigma_shift_sumrule_at_k_vk_nsc_msc( ...
            en, Dk, nbands, Egrid, egam, fermi);

        sigma_nsc = sigma_nsc + weight * nhc_nsc_k;
        sigma_msc = sigma_msc + weight * nhc_msc_k;

        if keep_kresolved
            nhc_nsc_all(:,:,:,:,ik) = weight * nhc_nsc_k;
            nhc_msc_all(:,:,:,:,ik) = weight * nhc_msc_k;
        end
    end

    sigma_nsc = Kconst * sigma_nsc;
    sigma_msc = Kconst * sigma_msc;
end
