function ops = skewgrid_ops_fromK(Kx, Ky, periodicFD, trimBoundary)
%SKEWGRID_OPS_FROMK Common skew-grid finite-difference utilities.
%
% 返回结构体字段:
%   ops.dk1, ops.dk2, ops.J, ops.invJT, ops.detJ, ops.w_k
%   ops.Nkx, ops.Nky, ops.ix_list, ops.iy_list, ops.mask_k, ops.mask_k4
%   ops.fd_du(U)       -> (du_x, du_y)
%   ops.fd_tensor(F)   -> (dF_dx, dF_dy)
%   ops.fd_rmn(r_mn)   -> dr
%
% 坐标变换:
%   [d/dkx; d/dky] = inv(J') * [d/dkappa1; d/dkappa2],  J=[dk1 dk2]

    arguments
        Kx double
        Ky double
        periodicFD logical = false
        trimBoundary logical = true
    end

    [Nkx, Nky] = size(Kx);
    if ~isequal(size(Ky), [Nkx, Nky])
        error('Kx and Ky must have the same size Nkx x Nky.');
    end

    [dk1, dk2] = grid_step_vectors(Kx, Ky);
    J = [dk1(:), dk2(:)];
    invJT = inv(J.');
    detJ = det(J);

    d2k = abs(detJ);
    w_k = d2k / (2*pi)^2;

    if trimBoundary && ~periodicFD
        ix_list = 2:(Nkx-1);
        iy_list = 2:(Nky-1);
        mask_k = zeros(Nkx, Nky);
        mask_k(ix_list, iy_list) = 1;
    else
        ix_list = 1:Nkx;
        iy_list = 1:Nky;
        mask_k = ones(Nkx, Nky);
    end

    ops = struct();
    ops.Nkx = Nkx;
    ops.Nky = Nky;
    ops.dk1 = dk1;
    ops.dk2 = dk2;
    ops.J = J;
    ops.invJT = invJT;
    ops.detJ = detJ;
    ops.d2k = d2k;
    ops.w_k = w_k;
    ops.ix_list = ix_list;
    ops.iy_list = iy_list;
    ops.mask_k = mask_k;
    ops.mask_k4 = reshape(mask_k, [Nkx, Nky, 1, 1]);

    ops.fd_du = @(U) fd_du_skew(U, invJT, periodicFD);
    ops.fd_tensor = @(F) fd_tensor_skew(F, invJT, periodicFD);
    ops.fd_rmn = @(r) fd_rmn_skew(r, invJT, periodicFD);
end

function [dk1, dk2] = grid_step_vectors(Kx, Ky)
%GRID_STEP_VECTORS Infer dk1/dk2 from K-grid arrays.

    [Nkx, Nky] = size(Kx);

    dKx1 = diff(Kx(:,1));
    dKy1 = diff(Ky(:,1));
    idx1 = find((abs(dKx1) + abs(dKy1)) > 0, 1, 'first');
    if isempty(idx1) || Nkx < 2
        error('Cannot infer dk1 from K-grid.');
    end
    dk1 = [dKx1(idx1); dKy1(idx1)];

    dKx2 = diff(Kx(1,:));
    dKy2 = diff(Ky(1,:));
    idx2 = find((abs(dKx2) + abs(dKy2)) > 0, 1, 'first');
    if isempty(idx2) || Nky < 2
        error('Cannot infer dk2 from K-grid.');
    end
    dk2 = [dKx2(idx2); dKy2(idx2)];
end

function [du_x, du_y] = fd_du_skew(U, invJT, periodicFD)
%FD_DU_SKEW Cartesian derivatives of eigenvectors on a skew grid.
% U: nb x nb_sel x Nkx x Nky
% du_x, du_y: same size (dU/dkx, dU/dky)

    du_x = zeros(size(U), 'like', U);
    du_y = zeros(size(U), 'like', U);

    [~, ~, Nkx, Nky] = size(U);

    if periodicFD
        idxp1 = [2:Nkx, 1];
        idxm1 = [Nkx, 1:Nkx-1];
        idxp2 = [2:Nky, 1];
        idxm2 = [Nky, 1:Nky-1];

        du_k1 = (U(:,:,idxp1,:) - U(:,:,idxm1,:)) / 2;
        du_k2 = (U(:,:,:,idxp2) - U(:,:,:,idxm2)) / 2;
    else
        du_k1 = zeros(size(U), 'like', U);
        du_k2 = zeros(size(U), 'like', U);

        ix = 2:(Nkx-1);
        iy = 2:(Nky-1);

        du_k1(:,:,ix,:) = (U(:,:,ix+1,:) - U(:,:,ix-1,:)) / 2;
        du_k2(:,:,:,iy) = (U(:,:,:,iy+1) - U(:,:,:,iy-1)) / 2;
    end

    du_x = invJT(1,1) * du_k1 + invJT(1,2) * du_k2;
    du_y = invJT(2,1) * du_k1 + invJT(2,2) * du_k2;
end

function [dF_dx, dF_dy] = fd_tensor_skew(F, invJT, periodicFD)
%FD_TENSOR_SKEW Central differences for tensor field F(ix,iy,...).

    sz = size(F);
    Nkx = sz(1);
    Nky = sz(2);
    rest = prod(sz(3:end));
    F3 = reshape(F, [Nkx, Nky, rest]);

    if periodicFD
        idxp1 = [2:Nkx, 1];
        idxm1 = [Nkx, 1:Nkx-1];
        idxp2 = [2:Nky, 1];
        idxm2 = [Nky, 1:Nky-1];

        d_k1 = (F3(idxp1,:,:) - F3(idxm1,:,:)) / 2;
        d_k2 = (F3(:,idxp2,:) - F3(:,idxm2,:)) / 2;
    else
        d_k1 = zeros(size(F3), 'like', F3);
        d_k2 = zeros(size(F3), 'like', F3);

        ix = 2:(Nkx-1);
        iy = 2:(Nky-1);

        d_k1(ix,:,:) = (F3(ix+1,:,:) - F3(ix-1,:,:)) / 2;
        d_k2(:,iy,:) = (F3(:,iy+1,:) - F3(:,iy-1,:)) / 2;
    end

    dF_dx3 = invJT(1,1) * d_k1 + invJT(1,2) * d_k2;
    dF_dy3 = invJT(2,1) * d_k1 + invJT(2,2) * d_k2;

    dF_dx = reshape(dF_dx3, sz);
    dF_dy = reshape(dF_dy3, sz);
end

function dr = fd_rmn_skew(r_mn, invJT, periodicFD)
%FD_RMN_SKEW
% r_mn: Nkx x Nky x nb x nb x 2(c)
% dr  : Nkx x Nky x nb x nb x 2(c) x 2(a)

    sz = size(r_mn);
    if numel(sz) ~= 5 || sz(5) ~= 2
        error('r_mn must be Nkx x Nky x nb x nb x 2.');
    end

    Nkx = sz(1);
    Nky = sz(2);
    nb = sz(3);
    nc = sz(5);

    r3 = reshape(r_mn, [Nkx, Nky, nb*nb*nc]);
    [dx3, dy3] = fd_tensor_skew(r3, invJT, periodicFD);
    dx = reshape(dx3, [Nkx, Nky, nb, nb, nc]);
    dy = reshape(dy3, [Nkx, Nky, nb, nb, nc]);

    dr = zeros(Nkx, Nky, nb, nb, nc, 2, 'like', r_mn);
    dr(:,:,:,:,:,1) = dx;
    dr(:,:,:,:,:,2) = dy;
end

