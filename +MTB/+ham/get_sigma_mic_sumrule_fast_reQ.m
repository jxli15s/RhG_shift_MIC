function [eta_mic, nhc_all] = get_sigma_mic_sumrule_fast_reQ( ...
    Enk, vxk, vyk, vzk, Egrid, egam, fermi, weight, positiveDE, Kconst)
%GET_SIGMA_MIC_SUMRULE_FAST_REQ
% Fast k-mesh MIC sum-rule using ReQ->g with precomputed Enk and velocities.
%
% Inputs:
%   Enk  : [Nkx,Nky,nb]
%   vxk  : [nb,nb,Nkx,Nky]
%   vyk  : [nb,nb,Nkx,Nky]
%   vzk  : [nb,nb,Nkx,Nky]
%   Egrid, egam, fermi, weight : same as shift fast interface
%   positiveDE (optional, default=true)
%   Kconst     (optional, default=e/h*1e6)
%
% Output:
%   eta_mic : [3,3,3,Nw]
%   nhc_all (optional): [3,3,3,Nw,Nk]

    if nargin < 9 || isempty(positiveDE)
        positiveDE = true;
    end

    e = 1.602176634e-19;   % C
    h = 6.582119569e-16;   % eV*s
    if nargin < 10 || isempty(Kconst)
        Kconst = e / h * 1e6;
    end

    [Nkx, Nky, nbands] = size(Enk);
    Nk = Nkx * Nky;
    Nw = numel(Egrid);

    eta_mic = zeros(3, 3, 3, Nw, 'like', 1i);

    keep_kresolved = (nargout > 1);
    if keep_kresolved
        nhc_all = zeros(3, 3, 3, Nw, Nk, 'like', 1i);
    else
        nhc_all = [];
    end

    vxk = reshape(vxk, nbands, nbands, Nk);
    vyk = reshape(vyk, nbands, nbands, Nk);
    vzk = reshape(vzk, nbands, nbands, Nk);
    Enk = reshape(Enk, Nk, nbands);

    parfor ik = 1:Nk
        en = Enk(ik, :).';
        Dk = cat(3, vxk(:,:,ik), vyk(:,:,ik), vzk(:,:,ik));

        nhc_k = MTB.ham.get_sigma_mic_sumrule_at_k_vk_reQ( ...
            en, Dk, nbands, Egrid, egam, fermi, positiveDE);

        eta_mic = eta_mic + weight * nhc_k;

        if keep_kresolved
            nhc_all(:,:,:,:,ik) = weight * nhc_k;
        end
    end

    eta_mic = Kconst * eta_mic;
end
