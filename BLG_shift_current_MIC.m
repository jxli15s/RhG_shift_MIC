clc;
clear;
delete(gcp('nocreate'))
% Nworkers  = 6;   % 外层 parfor worker 数
% Nthreads  = 1;    % 每个 worker 内部 BLAS 线程数
% pctRunOnAll maxNumCompThreads(Nthreads);
parpool('local',6);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            RhG Model      %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the structure
g = MTB.geometry("RhG");
a = 2.46;   % Ang
g.a =  [1/2,-sqrt(3)/2,0;...
        1/2,sqrt(3)/2,0;...
        0,0,1/a]*a;

g.b=inv(g.a')*2*pi;
Layer_N = 2;
pars = struct('v0',3.16, 'gamma1',0.46, 'gamma2',-0.017, ...   % eV*Ang (示例数)
              'gamma3',-0.30, 'gamma4',-0.086, ... % eV
              'uext',0.015, 'delta',-0.0011, 'xi',1, 'N',Layer_N);        % eV
% % 
% Layer_N = 3;
% pars = struct('v0',3.16, 'gamma1',0.435, 'gamma2',-0.0185, ...   % eV*Ang (示例数)
%               'gamma3',-0.322, 'gamma4',-0.0675, ... % eV
%               'uext',0.02, 'delta',-0.000147, 'xi',1, 'N',Layer_N);
% get the ham
%kpoint=[0.5,0.0,0.0];
%[h,hx,hy,hz] = MTB.ham.get_ham_kp_RhG(kpoint, pars);
model=@MTB.ham.get_ham_kp_RhG;
% model=@MTB.ham.get_ham_kp_RhG_valley;
g.dim_kp=Layer_N*2;
nbands=g.dim_kp;
knum=301;
tic;
efermi=tbNLC.get_ef(g,pars,model,nbands,knum);
toc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Get the 3D band in plane   %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=300;
kxline=[-0.03,0.03];
kyline=[-0.03,0.03];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
model=@MTB.ham.get_ham_kp_RhG;
nbands=g.dim_kp;
part=(kxline(2)-kxline(1))^2;
knum=size(Kx,1);
weights=1/knum^2/det(g.a)*part*pi*2;

tic;
% [Unk,Enk,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity(pars,model,nbands,Kx,Ky,Kz);
[Unk,Enk,ham,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity_withH(pars,model,nbands,Kx,Ky,Kz);
toc;

efermi=tbNLC.calculate_ef(Enk(:),0.5);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            RhG Model      %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the structure
g = MTB.geometry("RhG");
a = 2.46;   % Ang
g.a =  [1/2,-sqrt(3)/2,0;...
        1/2,sqrt(3)/2,0;...
        0,0,1/a]*a;

g.b=inv(g.a')*2*pi;
Layer_N = 2;
pars = struct('v0',3.16, 'gamma1',0.46, 'gamma2',-0.017, ...   % eV*Ang (示例数)
              'gamma3',-0.30, 'gamma4',-0.086, ... % eV
              'uext',0.015, 'delta',-0.0011, 'xi',1, 'N',Layer_N);        % eV
% % 
% Layer_N = 3;
% pars = struct('v0',3.16, 'gamma1',0.435, 'gamma2',-0.0185, ...   % eV*Ang (示例数)
%               'gamma3',-0.322, 'gamma4',-0.0675, ... % eV
%               'uext',0.02, 'delta',-0.000147, 'xi',1, 'N',Layer_N);
% get the ham
%kpoint=[0.5,0.0,0.0];
%[h,hx,hy,hz] = MTB.ham.get_ham_kp_RhG(kpoint, pars);
model=@MTB.ham.get_ham_kp_RhG;
% model=@MTB.ham.get_ham_kp_RhG_valley;
g.dim_kp=Layer_N*2;
nbands=g.dim_kp;
knum=501;
tic;
efermi=tbNLC.get_ef(g,pars,model,nbands,knum);
toc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Get the 3D band in plane   %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=200;
kxline=[-0.03,0.03];
kyline=[-0.03,0.03];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
model=@MTB.ham.get_ham_kp_RhG;
nbands=g.dim_kp;
part=(kxline(2)-kxline(1))^2;
knum=size(Kx,1);
weights=1/knum^2/det(g.a)*part*pi*2;

tic;
% [Unk,Enk,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity(pars,model,nbands,Kx,Ky,Kz);
[Unk,Enk,ham,vxk,vyk,vzk]=MTB.ham.get_bulk_plane_kp_velocity_withH(pars,model,nbands,Kx,Ky,Kz);
toc;

efermi=tbNLC.calculate_ef(Enk(:),0.5);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculate the shift current by Ham  %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== photon energy grid (eV) =====
Eph_list = linspace(0.0, 0.1, 2000);   % hbar*omega in eV
% ===== fermi + broadening =====
Ef  = efermi;      % eV
kT  = 0.0;       % eV (0 => step)
eta = 0.0001;      % eV (Lorentz broadening for delta(Em-En - Eph))

% ===== options =====
opts = struct();
opts.band_list   = 2:3;     % only sum within these bands (speedup)
opts.periodicFD  = false;   % 推荐先 false：丢边界做中心差分（最稳）
opts.trimBoundary= true;    % true -> 使用内部 (2..Nk-1)
opts.symBC       = true;    % symmetrize b<->c
opts.verbose     = true;
opts.useEmbedding = false;
opts.doGaugeFix = true;
opts.g_s = 1;
opts.saveIntermediates = false;

% opts.tau = [0,0,0];

% [sigma_abc, out] = shift_current_plane_fd_energy_skew( Kx, Ky, Kz, pars, model, Eph_list, Ef, kT, eta, opts);
[sigma_abc, out] = tbNLC.shift_current_plane_fd_energy_skew_fromUE(Kx, Ky, Unk, Enk, Eph_list, Ef, kT, eta, opts);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculate the shift current by Ham from HFMF  %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S=load("/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data/ne=0.0000e12_U=14.000data.mat");
S=load("/Volumes/T9/work/tb/matlab/data/RhG_HFMF/chiho/sequence/2001/ne=0.0000e12_U=0.000data.mat");
% S=load("/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/parallel/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data/ne=0.0000e12_U=6.000data.mat");
H_int=S.H_int/1000;
Ham_int=permute(H_int,[3,4,1,2]);
E_int=S.E_int/1000;
fermi=tbNLC.calculate_ef(E_int(:), 0.5);
knum=2000;
kxline=[-0.15,0.15]/2/pi;
kyline=[-0.15,0.15]/2/pi;
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
%%
dopts = struct();
dopts.band_list = 1:2; % or choose your relevant bands
[H_blocks] = tbNLC.extract_spinvalley_blocks(H_int,1);

[Unk1,Enk1] = tbNLC.diag_mesh_fromHk(permute(H_blocks.K_up,[3,4,1,2]), dopts);
% [Unk2,Enk2] = tbNLC.diag_mesh_fromHk(permute(H_blocks.Kp_up,[3,4,1,2]), dopts);
% [Unk3,Enk3] = tbNLC.diag_mesh_fromHk(permute(H_blocks.K_dn,[3,4,1,2]), dopts);
% [Unk4,Enk4] = tbNLC.diag_mesh_fromHk(permute(H_blocks.Kp_dn,[3,4,1,2]), dopts);
%%
% efermi=tbNLC.calculate_ef([Enk1(:);Enk2(:);Enk3(:);Enk4(:)],0.5);
efermi=tbNLC.calculate_ef(Enk1(:),0.5);
%%
% ===== photon energy grid (eV) =====
Eph_list = linspace(0.0, 0.1, 1000);   % hbar*omega in eV
%
% ===== fermi + broadening =====
Ef  = efermi;      % eV
kT  = 0.0;       % eV (0 => step)
eta = 0.001;      % eV (Lorentz broadening for delta(Em-En - Eph))

% ===== options =====
opts = struct();
opts.band_list   = 1:2;     % only sum within these bands (speedup)
opts.periodicFD  = false;   % 推荐先 false：丢边界做中心差分（最稳）
opts.trimBoundary= true;    % true -> 使用内部 (2..Nk-1)
opts.symBC       = true;    % symmetrize b<->c
opts.verbose     = true;
opts.useEmbedding = false;
opts.doGaugeFix = true;
opts.g_s = 1;
opts.saveIntermediates = false;


tic;
sigma_abc1 = tbNLC.shift_current_plane_fd_energy_skew_fromUE(Kx, Ky, Unk1, Enk1, Eph_list, Ef, kT, eta, opts); % spin up
% sigma_abc2 = tbNLC.shift_current_plane_fd_energy_skew_fromUE(Kx, Ky, Unk2, Enk2, Eph_list, Ef, kT, eta, opts); % spin up
% sigma_abc3 = tbNLC.shift_current_plane_fd_energy_skew_fromUE(Kx, Ky, Unk3, Enk3, Eph_list, Ef, kT, eta, opts); % spin dn
% sigma_abc4 = tbNLC.shift_current_plane_fd_energy_skew_fromUE(Kx, Ky, Unk4, Enk4, Eph_list, Ef, kT, eta, opts); % spin dn
toc;
%%
% sigma_abc=sigma_abc1+sigma_abc2+sigma_abc3+sigma_abc4;
sigma_abc=sigma_abc1
%%
figure()
hold on
str=['x','y'];
for i=1:2
    for j=1:2
        for k=1:2
            sig_xxy = squeeze(sigma_abc(i,j,k,:))*10^6/10;
            % sig_xxy = squeeze(sig_r(i,j,k,:))/10;
            % sig_xxy = squeeze(sigma_C3v(i,j,k,:))*10^6/10/10^20;
            plot(Eph_list, real(sig_xxy), '--','DisplayName',  sprintf('%s%s%s',str(i),str(j),str(k)));
            legend
        end
    end
end
%%
%%
% ===== photon energy grid (eV) =====
Eph_list = linspace(0.0, 0.2, 2000);   % hbar*omega in eV
% ===== fermi + broadening =====
Ef  = efermi;      % eV
kT  = 0.0;       % eV (0 => step)
eta = 0.001;      % eV (Lorentz broadening for delta(Em-En - Eph))

% ===== options =====
opts = struct();
opts.band_list   = 2:3;  
opts.periodicFD  = false;   % 推荐先 false：丢边界做中心差分（最稳）
opts.trimBoundary= true;    % true -> 使用内部 (2..Nk-1)
opts.symBC       = true;    % symmetrize b<->c
opts.verbose     = true;
opts.useEmbedding = false;
opts.doGaugeFix = true;
opts.g_s = 1;
opts.saveIntermediates = false;
opts.saveFullMN = false;
%%
eta_abc = tbNLC.mic_metric_plane_fromUE_mn(Kx, Ky, Unk, Enk, Eph_list, Ef, kT, eta, opts);
%% 
% [eta_abc, out] = mic_metric_plane_fromUE_mn_unified(Kx, Ky, Unk, Enk, Eph_list, Ef, kT, eta, opts);
%%
figure()
hold on
str=['x','y'];
for i=1:2
    for j=1:2
        for k=1:2
            sig_xxy = squeeze(eta_abc(i,j,k,:))*10^-13*10^6/10;
            % sig_xxy = squeeze(sigma_abc(i,j,k,:))*10^6/10;
            % sig_xxy = squeeze(sig_r(i,j,k,:))/10;
            % sig_xxy = squeeze(sigma_C3v(i,j,k,:))*10^6/10/10^20;
            plot(Eph_list, real(sig_xxy), '--','DisplayName',  sprintf('%s%s%s',str(i),str(j),str(k)));
            legend
        end
    end
end
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculate the magnetic injection current by Ham from HFMF  %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=load("/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data/ne=0.0000e12_U=8.000data.mat");
% S=load("/Volumes/T9/work/code/python/chiho_hfmf/HF_share/test/parallel/result/dsweep_251125/er=27.0_alp=0.030/SOC=0.00_fill=-1_LAFz_5.0_5.0_-5.0_-5.0_SOC_single_c3_spinless/matlab_data/ne=0.0000e12_U=6.000data.mat");
H_int=S.H_int/1000;
Ham_int=permute(H_int,[3,4,1,2]);
E_int=S.E_int/1000;
fermi=tbNLC.calculate_ef(E_int(:), 0.5);
knum=300;
kxline=[-0.15,0.15]/2/pi;
kyline=[-0.15,0.15]/2/pi;
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
%%
dopts = struct();
dopts.band_list = 1:2; % or choose your relevant bands
[H_blocks] = tbNLC.extract_spinvalley_blocks(H_int,1);

[Unk1,Enk1] = tbNLC.diag_mesh_fromHk(permute(H_blocks.K_up,[3,4,1,2]), dopts);
[Unk2,Enk2] = tbNLC.diag_mesh_fromHk(permute(H_blocks.Kp_up,[3,4,1,2]), dopts);
[Unk3,Enk3] = tbNLC.diag_mesh_fromHk(permute(H_blocks.K_dn,[3,4,1,2]), dopts);
[Unk4,Enk4] = tbNLC.diag_mesh_fromHk(permute(H_blocks.Kp_dn,[3,4,1,2]), dopts);
%%
efermi=tbNLC.calculate_ef([Enk1(:);Enk2(:);Enk3(:);Enk4(:)],0.5);
%%
% ===== photon energy grid (eV) =====
Eph_list = linspace(0.0, 0.2, 3000);   % hbar*omega in eV
% ===== fermi + broadening =====
Ef  = efermi;      % eV
kT  = 8.6173*1e-5*10;       % eV (0 => step)
eta = 0.001;      % eV (Lorentz broadening for delta(Em-En - Eph))

% ===== options =====
opts = struct();
opts.periodicFD  = false;   % 推荐先 false：丢边界做中心差分（最稳）
opts.trimBoundary= true;    % true -> 使用内部 (2..Nk-1)
opts.symBC       = true;    % symmetrize b<->c
opts.verbose     = true;
opts.useEmbedding = false;
opts.doGaugeFix = true;
opts.g_s = 1;
opts.saveIntermediates = false;
opts.saveFullMN = false;
%%
tic;
eta_abc1 = tbNLC.mic_metric_plane_fromUE_mn(Kx, Ky, Unk1, Enk1, Eph_list, Ef, kT, eta, opts);
eta_abc2 = tbNLC.mic_metric_plane_fromUE_mn(Kx, Ky, Unk2, Enk2, Eph_list, Ef, kT, eta, opts);
eta_abc3 = tbNLC.mic_metric_plane_fromUE_mn(Kx, Ky, Unk3, Enk3, Eph_list, Ef, kT, eta, opts);
eta_abc4 = tbNLC.mic_metric_plane_fromUE_mn(Kx, Ky, Unk4, Enk4, Eph_list, Ef, kT, eta, opts);
toc;

%%
eta_abc=eta_abc1+eta_abc2+eta_abc3+eta_abc4;
%%

figure()
hold on
str=['x','y'];
for i=1:2
    for j=1:2
        for k=1:2
            sig_xxy = squeeze(eta_abc4(i,j,k,:))*10^-13*10^6/10;
            % sig_xxy = squeeze(sigma_abc(i,j,k,:))*10^6/10;
            % sig_xxy = squeeze(sig_r(i,j,k,:))/10;
            % sig_xxy = squeeze(sigma_C3v(i,j,k,:))*10^6/10/10^20;
            plot(Eph_list, real(sig_xxy), '--','DisplayName',  sprintf('%s%s%s',str(i),str(j),str(k)));
            legend
        end
    end
end
%%

function ops = skewgrid_ops_fromK(Kx, Ky, periodicFD, trimBoundary)
%SKEWGRID_OPS_FROMK  Common skew-grid finite-difference utilities.
%
% Returns struct:
%   ops.dk1, ops.dk2, ops.J, ops.invJT, ops.detJ, ops.w_k
%   ops.Nkx, ops.Nky, ops.ix_list, ops.iy_list, ops.mask_k, ops.mask_k4
%   ops.fd_du(U)       -> (du_x, du_y)   for U: nb x nb_sel x Nkx x Nky
%   ops.fd_tensor(F)   -> (dF_dx,dF_dy)  for F: Nkx x Nky x ...
%   ops.fd_rmn(r_mn)   -> dr             for r_mn: Nkx x Nky x nb x nb x 2
%                         dr: Nkx x Nky x nb x nb x 2(c) x 2(a)
%
% Transform:
%   [∂/∂kx; ∂/∂ky] = inv(J') [∂/∂κ1; ∂/∂κ2],  J=[dk1 dk2]
%
% Notes:
% - We do NOT divide by physical step size because dk1,dk2 already carry that.
% - Central differences use /2 in index space, then inv(J') maps to Cartesian.

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

    [dk1, dk2] = grid_step_vectors(Kx, Ky);   % 2x1 each
    J     = [dk1(:), dk2(:)];                 % 2x2
    invJT = inv(J.');                         % J^{-T}
    detJ  = det(J);

    d2k = abs(detJ);
    w_k = d2k / (2*pi)^2;

    if trimBoundary && ~periodicFD
        ix_list = 2:(Nkx-1);
        iy_list = 2:(Nky-1);
        mask_k  = zeros(Nkx, Nky);
        mask_k(ix_list, iy_list) = 1;
    else
        ix_list = 1:Nkx;
        iy_list = 1:Nky;
        mask_k  = ones(Nkx, Nky);
    end
    mask_k4 = reshape(mask_k, [Nkx, Nky, 1, 1]);

    ops = struct();
    ops.Nkx = Nkx; ops.Nky = Nky;
    ops.dk1 = dk1; ops.dk2 = dk2;
    ops.J = J; ops.invJT = invJT; ops.detJ = detJ;
    ops.d2k = d2k; ops.w_k = w_k;
    ops.ix_list = ix_list; ops.iy_list = iy_list;
    ops.mask_k = mask_k; ops.mask_k4 = mask_k4;

    ops.fd_du     = @(U) fd_du_skew(U, invJT, periodicFD);
    ops.fd_tensor = @(F) fd_tensor_skew(F, invJT, periodicFD);
    ops.fd_rmn    = @(r) fd_rmn_skew(r, invJT, periodicFD);
end

% ======================================================================
% local helpers
% ======================================================================

function [dk1, dk2] = grid_step_vectors(Kx, Ky)
% infer dk1 (ix direction) and dk2 (iy direction) from Kx,Ky arrays
    [Nkx, Nky] = size(Kx);

    dKx1 = diff(Kx(:,1)); dKy1 = diff(Ky(:,1));
    idx1 = find((abs(dKx1)+abs(dKy1))>0, 1, 'first');
    if isempty(idx1) || Nkx < 2
        error('Cannot infer dk1 from K-grid.');
    end
    dk1 = [dKx1(idx1); dKy1(idx1)];

    dKx2 = diff(Kx(1,:)); dKy2 = diff(Ky(1,:));
    idx2 = find((abs(dKx2)+abs(dKy2))>0, 1, 'first');
    if isempty(idx2) || Nky < 2
        error('Cannot infer dk2 from K-grid.');
    end
    dk2 = [dKx2(idx2); dKy2(idx2)];
end

function [du_x, du_y] = fd_du_skew(U, invJT, periodicFD)
% Cartesian derivatives of eigenvectors on a skew grid
% U: nb x nb_sel x Nkx x Nky
% du_x, du_y: same size, representing ∂u/∂kx and ∂u/∂ky

    du_x = zeros(size(U), 'like', U);
    du_y = zeros(size(U), 'like', U);

    [~, ~, Nkx, Nky] = size(U);

    if periodicFD
        idxp1 = [2:Nkx, 1];
        idxm1 = [Nkx, 1:Nkx-1];
        idxp2 = [2:Nky, 1];
        idxm2 = [Nky, 1:Nky-1];

        du_k1 = (U(:,:,idxp1,:) - U(:,:,idxm1,:)) / 2; % along κ1
        du_k2 = (U(:,:,:,idxp2) - U(:,:,:,idxm2)) / 2; % along κ2
    else
        du_k1 = zeros(size(U), 'like', U);
        du_k2 = zeros(size(U), 'like', U);

        ix = 2:(Nkx-1);
        iy = 2:(Nky-1);

        du_k1(:,:,ix,:) = (U(:,:,ix+1,:) - U(:,:,ix-1,:)) / 2;
        du_k2(:,:,:,iy) = (U(:,:,:,iy+1) - U(:,:,:,iy-1)) / 2;
    end

    % [∂x; ∂y] = inv(J') [∂k1; ∂k2]
    du_x = invJT(1,1)*du_k1 + invJT(1,2)*du_k2;
    du_y = invJT(2,1)*du_k1 + invJT(2,2)*du_k2;
end

function [dF_dx, dF_dy] = fd_tensor_skew(F, invJT, periodicFD)
% Central differences for a tensor field F(ix,iy,...) on skew grid,
% then transform to Cartesian derivatives using inv(J').

    sz = size(F);
    Nkx = sz(1); Nky = sz(2);
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

    dF_dx3 = invJT(1,1)*d_k1 + invJT(1,2)*d_k2;
    dF_dy3 = invJT(2,1)*d_k1 + invJT(2,2)*d_k2;

    dF_dx = reshape(dF_dx3, sz);
    dF_dy = reshape(dF_dy3, sz);
end

function dr = fd_rmn_skew(r_mn, invJT, periodicFD)
% r_mn: Nkx x Nky x nb x nb x 2(c)
% dr  : Nkx x Nky x nb x nb x 2(c) x 2(a)

    sz = size(r_mn);
    if numel(sz) ~= 5 || sz(5) ~= 2
        error('r_mn must be Nkx x Nky x nb x nb x 2.');
    end
    Nkx = sz(1); Nky = sz(2);
    nb  = sz(3); nc = sz(5);

    r3 = reshape(r_mn, [Nkx, Nky, nb*nb*nc]);
    [dx3, dy3] = fd_tensor_skew(r3, invJT, periodicFD);
    dx = reshape(dx3, [Nkx, Nky, nb, nb, nc]);
    dy = reshape(dy3, [Nkx, Nky, nb, nb, nc]);

    dr = zeros(Nkx, Nky, nb, nb, nc, 2, 'like', r_mn);
    dr(:,:,:,:,:,1) = dx; % a=1 -> kx
    dr(:,:,:,:,:,2) = dy; % a=2 -> ky
end


function rH = hermitianize_rmn(r_mn)
%HERMITIANIZE_RMN  Enforce Hermiticity in band indices (m,n) for r_mn.
% rH = 0.5*(r + r^\dagger), where dagger means conj + swap (m<->n).
%
% This matches your MIC "r_mn_H" step exactly.

    rH = 0.5 * ( r_mn + permute(conj(r_mn), [1 2 4 3 5]) );
end


function [eta_abc, out] = mic_metric_plane_fromUE_mn( ...
    Kx, Ky, U, E, Eph_list, Ef, kT, eta, opts)
%MIC_METRIC_PLANE_FROMUE_MN
% =========================================================================
% Magnetic Injection Current (MIC) for LINEAR polarization (metric-type),
% multi-band form: SUM over all interband pairs (m,n), m≠n, on a 2D k-grid.
%
% Energy (eV) form with Lorentzian delta_eV (unit 1/eV):
%
%   eta^{abc}(Eph) = -(pi * g_s * e^2)/(2*hbar) *
%       ∫ d^2k/(2π)^2  Σ_{m≠n}  f_nm(k) * Δv^a_mn(k) * [2 g^{bc}_mn(k)] *
%       δ_eV( (E_m - E_n) - Eph )
%
% where (linear polarization => symmetric in b,c):
%   2 g^{bc}_{mn} = r^b_{mn} r^c_{nm} + r^c_{mn} r^b_{nm}
%   r^b_{mn} = i <u_m | ∂_{k_b} u_n>
%   v^a_n    = (1/ħ) ∂E_n/∂k_a  (E in eV => v = (e/ħ)∂E/∂k)
% =========================================================================

    arguments
        Kx double
        Ky double
        U  {mustBeNumeric}
        E  double
        Eph_list double
        Ef double
        kT double
        eta double
        opts struct = struct()
    end

    % -------------------- options --------------------
    periodicFD   = get_opt(opts,'periodicFD',false);
    trimBoundary = get_opt(opts,'trimBoundary',true);
    doGaugeFix   = get_opt(opts,'doGaugeFix',true);
    g_s          = get_opt(opts,'g_s',1);
    verbose      = get_opt(opts,'verbose',true);
    positiveDE   = get_opt(opts,'positiveDE',true);
    saveFullMN   = get_opt(opts,'saveFullMN',true);

    % -------------------- constants (SI) --------------------
    e_charge = 1.602176634e-19;   % C
    hbar_Js  = 1.054571817e-34;   % J*s

    pref     = -pi * g_s * (e_charge^2) / (2*hbar_Js);
    vel_pref =  e_charge / hbar_Js;

    % -------------------- reshape & checks --------------------
    [Nkx, Nky] = size(Kx);
    if ~isequal(size(Ky), [Nkx,Nky])
        error('Kx and Ky must have the same size Nkx x Nky.');
    end

    % reshape U to nb x nb_sel x Nkx x Nky
    szU = size(U);
    if numel(szU) == 3
        nb = szU(1); nb_sel = szU(2);
        if szU(3) ~= Nkx*Nky
            error('If U is 3D, its 3rd dim must be Nkx*Nky.');
        end
        U = reshape(U, [nb, nb_sel, Nkx, Nky]);
    elseif numel(szU) == 4
        nb = szU(1); nb_sel = szU(2);
        if szU(3)~=Nkx || szU(4)~=Nky
            error('U must be nb x nb_sel x Nkx x Nky.');
        end
    else
        error('U must be 4D or 3D.');
    end

    % reshape E to Nkx x Nky x nb_sel
    if ~isequal(size(E), [Nkx,Nky,nb_sel])
        if isequal(size(E), [Nkx*Nky, nb_sel])
            E = reshape(E, [Nkx, Nky, nb_sel]);
        else
            error('E must be Nkx x Nky x nb_sel (or (Nkx*Nky) x nb_sel).');
        end
    end

    Eph_list = Eph_list(:).'; % 1 x Nw
    Nw = numel(Eph_list);

    % -------------------- shared skew-grid ops --------------------
    ops = skewgrid_ops_fromK(Kx, Ky, periodicFD, trimBoundary);

    if verbose
        fprintf('[MIC metric mn] Nkx=%d Nky=%d nb=%d nb_sel=%d\n', Nkx, Nky, nb, nb_sel);
        fprintf('[MIC metric mn] dk1=(%.3e,%.3e), dk2=(%.3e,%.3e), det(J)=%.3e\n', ...
            ops.dk1(1),ops.dk1(2),ops.dk2(1),ops.dk2(2), ops.detJ);
        fprintf('[MIC metric mn] periodicFD=%d, trimBoundary=%d, doGaugeFix=%d, g_s=%g, positiveDE=%d\n', ...
            periodicFD, trimBoundary, doGaugeFix, g_s, positiveDE);
    end

    % -------------------- gauge smoothing --------------------
    if doGaugeFix
        Ug = gauge_fix_parallel_transport(U);
    else
        Ug = U;
    end

    % -------------------- masks / lists --------------------
    mask_k4 = ops.mask_k4;
    ix_list = ops.ix_list;
    iy_list = ops.iy_list;

    % -------------------- eigenvector derivatives --------------------
    [du_x, du_y] = ops.fd_du(Ug);

    % -------------------- r_mn(k) --------------------
    r_mn = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 'like', Ug);

    for ix = ix_list
        for iy = iy_list
            U0  = squeeze(Ug(:,:,ix,iy));    % nb x nb_sel
            dux = squeeze(du_x(:,:,ix,iy));  % nb x nb_sel
            duy = squeeze(du_y(:,:,ix,iy));  % nb x nb_sel

            r_mn(ix,iy,:,:,1) = 1i * (U0' * dux);
            r_mn(ix,iy,:,:,2) = 1i * (U0' * duy);
        end
    end

    % enforce off-diagonal only using a mask (branchless)
    off = ones(nb_sel) - eye(nb_sel);
    mask_off = reshape(off, [1,1,nb_sel,nb_sel,1]); % broadcast

    % -------------------- Hermitian projection (keep your original behavior) ----
    r_mn_H = hermitianize_rmn(r_mn);
    r_mn_H = r_mn_H .* mask_off;

    % For convenience: r_nm with (m,n) indexing again
    r_nm_mn = permute(r_mn_H, [1 2 4 3 5]);   % Nkx x Nky x m x n x 2

    % -------------------- metric kernel: g2^{bc}_{mn} --------------------
    g2 = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 2, 'like', Ug);
    for b = 1:2
        for cc = 1:2
            g2(:,:,:,:,b,cc) = real( ...
                r_mn_H(:,:,:,:,b)  .* r_nm_mn(:,:,:,:,cc) + ...
                r_mn_H(:,:,:,:,cc) .* r_nm_mn(:,:,:,:,b)  );
        end
    end

    % -------------------- energies: ΔE_mn = E_m - E_n --------------------
    Em = reshape(E, [Nkx, Nky, nb_sel, 1]); % ... m
    En = reshape(E, [Nkx, Nky, 1, nb_sel]); % ... n
    dE_mn = Em - En;                        % Nkx x Nky x m x n (eV)

    if positiveDE
        mask_pos = double(dE_mn > 0);
    else
        mask_pos = ones(size(dE_mn));
    end

    % -------------------- Fermi factor f_nm = f_n - f_m --------------------
    f = fermi_dirac(E, Ef, kT); % Nkx x Nky x nb_sel
    fm = reshape(f, [Nkx, Nky, nb_sel, 1]); % m
    fn = reshape(f, [Nkx, Nky, 1, nb_sel]); % n
    f_nm = fn - fm;

    % -------------------- velocities v_n = (e/ħ) ∂E/∂k --------------------
    [dE_dx, dE_dy] = ops.fd_tensor(E); % Nkx x Nky x nb_sel
    v_band = zeros(Nkx, Nky, nb_sel, 2);
    v_band(:,:,:,1) = vel_pref * dE_dx;  % vx
    v_band(:,:,:,2) = vel_pref * dE_dy;  % vy

    dv_mn = zeros(Nkx, Nky, nb_sel, nb_sel, 2);
    for a = 1:2
        Vm = reshape(v_band(:,:,:,a), [Nkx, Nky, nb_sel, 1]);
        Vn = reshape(v_band(:,:,:,a), [Nkx, Nky, 1, nb_sel]);
        dv_mn(:,:,:,:,a) = Vm - Vn;
    end

    % -------------------- total mask --------------------
    mask_total = mask_k4 .* reshape(off, [1,1,nb_sel,nb_sel]) .* mask_pos;

    % -------------------- integrate over photon energy --------------------
    eta_abc = zeros(2,2,2,Nw);

    for a = 1:2
        dvA = dv_mn(:,:,:,:,a);
        for b = 1:2
            for cc = 1:2
                gBC = g2(:,:,:,:,b,cc);
                Kabc = (f_nm .* dvA .* gBC) .* mask_total;

                for iw = 1:Nw
                    Eph = Eph_list(iw);
                    delta_w = (1/pi) * eta ./ ((dE_mn - Eph).^2 + eta^2); % 1/eV
                    eta_abc(a,b,cc,iw) = pref * ops.w_k * sum(Kabc .* delta_w, 'all');
                end
            end
        end
    end

    % -------------------- outputs --------------------
    out = struct();
    out.pref      = pref;
    out.vel_pref  = vel_pref;
    out.w_k       = ops.w_k;
    out.J         = ops.J;
    out.invJT     = ops.invJT;
    out.dk1       = ops.dk1;
    out.dk2       = ops.dk2;
    out.mask_k    = ops.mask_k;
    out.mask_off  = off;
    out.positiveDE = positiveDE;
    out.Ug        = Ug;

    if saveFullMN
        out.r_mn    = r_mn_H;
        out.g2      = g2;
        out.dE_mn   = dE_mn;
        out.f_nm    = f_nm;
        out.v_band  = v_band;
        out.dv_mn   = dv_mn;
    end
end

% ======================================================================
% helpers (same as your original style)
% ======================================================================

function val = get_opt(opts, name, default)
    if isstruct(opts) && isfield(opts, name) && ~isempty(opts.(name))
        val = opts.(name);
    else
        val = default;
    end
end

function Ug = gauge_fix_parallel_transport(U)
% band-wise phase smoothing (nondegenerate)
    Ug = U;
    [~, nb_sel, Nkx, Nky] = size(U);

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

function f = fermi_dirac(E, Ef, kT)
% energies in eV
    if kT <= 0
        f = double(E < Ef);
    else
        f = 1 ./ (1 + exp((E - Ef)/kT));
    end
end

function [sigma_abc, out] = shift_current_plane_fd_energy_skew_fromUE( ...
    Kx, Ky, U, E, Eph_list, Ef, kT, eta, opts)
%SHIFT_CURRENT_PLANE_FD_ENERGY_SKEW_FROMUE
% =========================================================================
% Shift current (2D) via eigenvector finite differences on a general 2D k-grid,
% taking eigenvectors U and eigenvalues E as INPUT.
%
% Energy-form:
%   sigma^{abc}(Omega) = (pi * g_s * e^2 / hbar) * ∫ (d^2k/(2pi)^2) Σ_{n≠m}
%        f_nm * Im[ r^b_{mn} * r^c_{nm;a} ] * δ_eV( (E_m-E_n) - Omega )
%   (+ b<->c sym if opts.symBC = true)
%
% Definitions:
%   r^b_{mn}   = i <u_m | ∂_{k_b} u_n>
%   A^a_nn     = i <u_n | ∂_{k_a} u_n>
%   r^c_{mn;a} = ∂_{k_a} r^c_{mn} - i (A^a_mm - A^a_nn) r^c_{mn}
% =========================================================================

    arguments
        Kx double
        Ky double
        U  {mustBeNumeric}
        E  double
        Eph_list double
        Ef double
        kT double
        eta double
        opts struct = struct()
    end

    % -------------------- options --------------------
    periodicFD   = get_opt(opts,'periodicFD',false);
    trimBoundary = get_opt(opts,'trimBoundary',true);
    symBC        = get_opt(opts,'symBC',true);
    verbose      = get_opt(opts,'verbose',true);
    g_s          = get_opt(opts,'g_s',1);
    doGaugeFix   = get_opt(opts,'doGaugeFix',true);

    % -------------------- constants --------------------
    e_charge = 1.602176634e-19;
    hbar_Js  = 1.054571817e-34;
    pref     = pi * g_s * (e_charge^2) / hbar_Js;

    % -------------------- sizes & reshape --------------------
    [Nkx, Nky] = size(Kx);
    if ~isequal(size(Ky), [Nkx,Nky])
        error('Kx and Ky must have the same size Nkx x Nky.');
    end

    szU = size(U);
    if numel(szU) == 3
        nb = szU(1); nb_sel = szU(2);
        if szU(3) ~= Nkx*Nky
            error('If U is 3D, its 3rd dim must be Nkx*Nky.');
        end
        U = reshape(U, [nb, nb_sel, Nkx, Nky]);
    elseif numel(szU) == 4
        nb = szU(1); nb_sel = szU(2);
        if szU(3)~=Nkx || szU(4)~=Nky
            error('U must be nb x nb_sel x Nkx x Nky.');
        end
    else
        error('U must be 4D or 3D.');
    end

    if ~isequal(size(E), [Nkx,Nky,nb_sel])
        if isequal(size(E), [Nkx*Nky, nb_sel])
            E = reshape(E, [Nkx, Nky, nb_sel]);
        else
            error('E must be Nkx x Nky x nb_sel (or (Nkx*Nky) x nb_sel).');
        end
    end

    Eph_list = Eph_list(:).';
    Nw = numel(Eph_list);

    % -------------------- shared skew-grid ops --------------------
    ops = skewgrid_ops_fromK(Kx, Ky, periodicFD, trimBoundary);

    if verbose
        fprintf('[shift_current_fromUE] Nkx=%d Nky=%d nb=%d nb_sel=%d\n', Nkx, Nky, nb, nb_sel);
        fprintf('[shift_current_fromUE] dk1=(%.6e, %.6e), dk2=(%.6e, %.6e)\n', ...
            ops.dk1(1),ops.dk1(2),ops.dk2(1),ops.dk2(2));
        fprintf('[shift_current_fromUE] det(J)=%.6e\n', ops.detJ);
        fprintf('[shift_current_fromUE] periodicFD=%d, trimBoundary=%d, symBC=%d, g_s=%g, doGaugeFix=%d\n', ...
            periodicFD, trimBoundary, symBC, g_s, doGaugeFix);
    end

    % -------------------- gauge smoothing --------------------
    if doGaugeFix
        Ug = gauge_fix_parallel_transport(U);
    else
        Ug = U;
    end

    ix_list = ops.ix_list;
    iy_list = ops.iy_list;

    % =====================================================================
    % 1) FD of eigenvectors in Cartesian directions
    % =====================================================================
    [du_x, du_y] = ops.fd_du(Ug);

    % =====================================================================
    % 2) Build A_diag and r_mn
    % =====================================================================
    A_diag = zeros(Nkx, Nky, nb_sel, 2);
    r_mn   = zeros(Nkx, Nky, nb_sel, nb_sel, 2, 'like', Ug);

    for ix = ix_list
        for iy = iy_list
            U0  = squeeze(Ug(:,:,ix,iy));
            dux = squeeze(du_x(:,:,ix,iy));
            duy = squeeze(du_y(:,:,ix,iy));

            Ax = 1i * diag(U0' * dux);
            Ay = 1i * diag(U0' * duy);
            A_diag(ix,iy,:,1) = real(Ax);
            A_diag(ix,iy,:,2) = real(Ay);

            r_mn(ix,iy,:,:,1) = 1i * (U0' * dux);
            r_mn(ix,iy,:,:,2) = 1i * (U0' * duy);
        end
    end

    % off-diagonal mask
    off = ones(nb_sel) - eye(nb_sel);
    r_mn = r_mn .* reshape(off, [1,1,nb_sel,nb_sel,1]);

    % IMPORTANT: make NSC use the SAME Hermitian projection convention as MIC
    r_mn = hermitianize_rmn(r_mn);
    r_mn = r_mn .* reshape(off, [1,1,nb_sel,nb_sel,1]);

    % =====================================================================
    % 3) FD of r_mn (Cartesian) using the same skew-grid ops
    % =====================================================================
    dr = ops.fd_rmn(r_mn);  % Nkx x Nky x m x n x 2(c) x 2(a)

    % =====================================================================
    % 4) Covariant derivative r_cov = dr - i(Amm-Ann)*r
    % =====================================================================
    r_cov = zeros(size(dr), 'like', dr);
    for a = 1:2
        Amm = reshape(A_diag(:,:,:,a), [Nkx, Nky, nb_sel, 1]);
        Ann = reshape(A_diag(:,:,:,a), [Nkx, Nky, 1, nb_sel]);
        dA  = Amm - Ann;

        for c = 1:2
            r_cov(:,:,:,:,c,a) = dr(:,:,:,:,c,a) - 1i * dA .* r_mn(:,:,:,:,c);
        end
    end

    % =====================================================================
    % 5) Fermi occupations
    % =====================================================================
    f_n = fermi_dirac(E, Ef, kT); % Nkx x Nky x nb_sel

    % =====================================================================
    % 6) Main integration
    % =====================================================================
    sigma_abc = zeros(2,2,2,Nw);

    for ix = ix_list
        for iy = iy_list
            Ek = squeeze(E(ix,iy,:));
            fk = squeeze(f_n(ix,iy,:));

            for n = 1:nb_sel
                for m = 1:nb_sel
                    if m==n, continue; end

                    dE   = Ek(m) - Ek(n);
                    f_nm = fk(n) - fk(m);
                    if abs(f_nm) < 1e-14, continue; end

                    delta_w = (1/pi) * eta ./ ((dE - Eph_list).^2 + eta^2);

                    for a = 1:2
                        for b = 1:2
                            for c = 1:2
                                I1 = imag( r_mn(ix,iy,m,n,b) * r_cov(ix,iy,n,m,c,a) );

                                if symBC
                                    I2 = imag( r_mn(ix,iy,m,n,c) * r_cov(ix,iy,n,m,b,a) );
                                    I  = 0.5*(I1 + I2);
                                else
                                    I  = I1;
                                end

                                sigma_abc(a,b,c,:) = sigma_abc(a,b,c,:) + ...
                                    reshape(pref * f_nm * (I * ops.w_k) .* delta_w, [1,1,1,Nw]);
                            end
                        end
                    end
                end
            end
        end
    end

    % -------------------- outputs --------------------
    out = struct();
    out.E = E;
    out.Ug = Ug;
    out.A_diag = A_diag;
    out.r_mn = r_mn;
    out.dr = dr;
    out.r_cov = r_cov;
    out.dk1 = ops.dk1;
    out.dk2 = ops.dk2;
    out.J = ops.J;
    out.invJT = ops.invJT;
    out.w_k = ops.w_k;
    out.pref = pref;
    out.ix_list = ix_list;
    out.iy_list = iy_list;
end


function [U, E] = diag_mesh_fromHk(Hk, opts)
%DIAG_MESH_FROMHK  (parfor-safe)
% Hk: nb x nb x Nkx x Nky  (Hermitian)
% U : nb x nb_sel x Nkx x Nky
% E : Nkx x Nky x nb_sel   (eV)
%
% opts.band_list (optional): choose subset of bands after sorting.


    [nb, nb2, Nkx, Nky] = size(Hk);
    if nb ~= nb2
        error('Hk must be nb x nb x Nkx x Nky.');
    end

    if isfield(opts,'band_list') && ~isempty(opts.band_list)
        band_list = opts.band_list(:).';
    else
        band_list = 1:nb;
    end
    nb_sel = numel(band_list);

    itotal = Nkx * Nky;

    % parfor-safe temporaries
    Utmp = zeros(nb, nb_sel, itotal, 'like', Hk);   % complex ok
    Etmp = zeros(itotal, nb_sel);                   % real energies

    parfor ll = 1:itotal
        % convert linear index -> (ix,iy)
        [ix, iy] = ind2sub([Nkx, Nky], ll);

        H = Hk(:,:,ix,iy);
        H = (H + H')/2;

        [V, D] = eig(H);
        evals = real(diag(D));
        [evals, ind] = sort(evals, 'ascend');
        V = V(:, ind);

        Etmp(ll,:) = evals(band_list);
        Utmp(:,:,ll) = V(:, band_list);
    end

    % reshape back
    U = reshape(Utmp, [nb, nb_sel, Nkx, Nky]);
    E = reshape(Etmp, [Nkx, Nky, nb_sel]);
end

function H_blocks = extract_spinvalley_blocks(H,N)
% 从k空间网格上的完整哈密顿量中提取各个自旋-能谷块
%
% 输入:
%   H:  完整哈密顿量, 维度 (nkx, nky, 8N, 8N)
%   Hx: ∂H/∂kx, 维度 (nkx, nky, 8N, 8N)
%   Hy: ∂H/∂ky, 维度 (nkx, nky, 8N, 8N)
%   Hz: ∂H/∂kz, 维度 (nkx, nky, 8N, 8N) (通常为0)
%   N:  层数
%
% 输出:
%   H_blocks: struct 包含各块
%     .K_up, .K_dn, .Kp_up, .Kp_dn
%     每个字段维度: (nkx, nky, 2N, 2N)
%   同样对于 Hx_blocks, Hy_blocks, Hz_blocks

arguments
    H 
    N
end

    % 获取k网格大小
    [nkx, nky, ~, ~] = size(H);
    dim_layer = 2*N;
    
    % ===== 基矢排列: (K↑, K'↑, K↓, K'↓) × layers =====
    
    % 预分配输出结构体
    H_blocks.K_up   = zeros(nkx, nky, dim_layer, dim_layer);
    H_blocks.Kp_up  = zeros(nkx, nky, dim_layer, dim_layer);
    H_blocks.K_dn   = zeros(nkx, nky, dim_layer, dim_layer);
    H_blocks.Kp_dn  = zeros(nkx, nky, dim_layer, dim_layer);
    
    % 定义块索引
    idx_K_up  = 1:dim_layer;
    idx_Kp_up = dim_layer + (1:dim_layer);
    idx_K_dn  = 2*dim_layer + (1:dim_layer);
    idx_Kp_dn = 3*dim_layer + (1:dim_layer);
    
    % ===== 提取块 =====
    % K谷, 自旋↑
    H_blocks.K_up   = H(:, :, idx_K_up, idx_K_up);   
    % K'谷, 自旋↑
    H_blocks.Kp_up  = H(:, :, idx_Kp_up, idx_Kp_up);
    % K谷, 自旋↓
    H_blocks.K_dn   = H(:, :, idx_K_dn, idx_K_dn);
    % K'谷, 自旋↓
    H_blocks.Kp_dn  = H(:, :, idx_Kp_dn, idx_Kp_dn);
end

function ef=get_ef(g,pars,model,nbands,knum)
    kxline=[-0.1,0.1];
    kyline=[-0.1,0.1];
    u=0.5;
    [Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
    [~,Enk]=MTB.ham.get_bulk_plane_kp(pars,model,nbands,Kx,Ky,Kz);
    ef=tbNLC.calculate_ef(Enk(:), u);
end

function efermi = calculate_ef(Enk, u)
    % 计算费米能级，基于填充因子 u
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)
    % u: 填充因子 (0 <= u <= 1)
    %
    % 输出:
    % efermi: 费米能级

    % 将能量值展平并取前 u*N 个最低能量值的最大值
    total_states = numel(Enk);                % 总的能量态数
    occupied_states = ceil(total_states * u); % 填充的态数
    efermi = max(mink(Enk(:), occupied_states));
end
