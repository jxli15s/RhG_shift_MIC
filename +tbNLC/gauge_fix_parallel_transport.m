function Ug = gauge_fix_parallel_transport(U)
%GAUGE_FIX_PARALLEL_TRANSPORT
% 在 2D k 网格上做并行输运式相位平滑（逐带、非简并假设）。
%
% 输入:
%   U: nb x nb_sel x Nkx x Nky
%
% 输出:
%   Ug: 与 U 同尺寸，k 点间相位更平滑的本征态

    Ug = U;
    [~, nb_sel, Nkx, Nky] = size(U);

    % 先沿 kx 方向平滑
    for iy = 1:Nky
        for ix = 2:Nkx
            Uprev = squeeze(Ug(:,:,ix-1,iy));
            Ucur  = squeeze(Ug(:,:,ix,iy));
            for n = 1:nb_sel
                ov = Uprev(:,n)' * Ucur(:,n);
                ph = ov / max(abs(ov), 1e-30);
                Ucur(:,n) = Ucur(:,n) / ph;
            end
            Ug(:,:,ix,iy) = Ucur;
        end
    end

    % 再沿 ky 方向平滑
    for ix = 1:Nkx
        for iy = 2:Nky
            Uprev = squeeze(Ug(:,:,ix,iy-1));
            Ucur  = squeeze(Ug(:,:,ix,iy));
            for n = 1:nb_sel
                ov = Uprev(:,n)' * Ucur(:,n);
                ph = ov / max(abs(ov), 1e-30);
                Ucur(:,n) = Ucur(:,n) / ph;
            end
            Ug(:,:,ix,iy) = Ucur;
        end
    end
end

