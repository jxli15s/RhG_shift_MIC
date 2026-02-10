function gs=get_supercell_wannier_3d_general(g,Umatrix,shift,mode)
% Slove the Hamiltonian on the K-mesh
% obj: The object of the continuum model system
% Kx,Ky: the k-mesh we will solve
% dim_H: num of the bands/orbits
% Unk: Norbits*Nands*knum*knum
% Enk: knum*knum*Nbands
gs = MTB.geometry("supercell");
gs.a=Umatrix*g.a;
gs.b=inv(gs.a')*2*pi;
pvolm=det(g.a); % primitive cell volum
svolm=det(gs.a); % supercell volum
c=fix(svolm/pvolm);
norbitalc = round(size(g.ham,1)); % number of orbitals in a cell
norbitals = round(size(g.ham,1)*svolm/pvolm); %total number of orbitals
natomc = size(g.atoms,1); % number of atoms in a cell
natoms = natomc*svolm/pvolm; %total number of atoms
gs.atoms=zeros(round(natoms),3);
gs.wpos=zeros(round(natoms),3);

if nargin < 4
    mode='fill';
elseif nargin < 3
    shift=[0,0,0];
end

if strcmp(mode, "brute")
    % mode
    rs3=replicate3d(g.atoms*g.a,g.a,c,c,c);
    rs1=brute_expand(rs3, gs.a(1,:), gs.a(2,:), gs.a(3,:), 2);
    rs1f=rs1*inv(gs.a)-round(rs1*inv(gs.a));
    gs.atoms=rs1f;
elseif strcmp(mode, "fill")
    % mode
    dim=3;
    [rs, sl] = generate_supercell_atoms(g, gs, c, dim);
    gs.atoms=rs*inv(gs.a);
    gs.sublattice=sl;
end
% 
% rs3=replicate3d(g.atoms*g.a,g.a,c,c,c);
% rs1=brute_expand(rs3, gs.a(1,:), gs.a(2,:), gs.a(3,:), 2);
% rs1f=rs1*inv(gs.a)-round(rs1*inv(gs.a));
% gs.atoms=rs1f;

for col_idx=1:3
    gs.atoms(gs.atoms(:,col_idx) > 0.999,col_idx)=gs.atoms(gs.atoms(:,col_idx) > 0.999,col_idx)-1;
    gs.atoms(gs.atoms(:,col_idx) < -0.0001,col_idx)=gs.atoms(gs.atoms(:,col_idx) < -0.0001,col_idx)+1;
end
gs.atoms=abs(gs.atoms);
[~, idx] = sortrows(gs.atoms, [1,2]);  % 先按 x 升序，再按 y 升序
gs.atoms = gs.atoms(idx, :);
gs.wpos=gs.atoms*gs.a;

% shift=[0.25,0,0];
gs.atoms=gs.atoms-shift;
for col_idx=1:3
    gs.atoms(gs.atoms(:,col_idx) > 0.999,col_idx)=gs.atoms(gs.atoms(:,col_idx) > 0.999,col_idx)-1;
    gs.atoms(gs.atoms(:,col_idx) < -0.0001,col_idx)=gs.atoms(gs.atoms(:,col_idx) < -0.0001,col_idx)+1;
end

gs.atoms=abs(gs.atoms);
[~, idx] = sortrows(gs.atoms, [1,2]);  % 先按 x 升序，再按 y 升序
gs.atoms = gs.atoms(idx, :);
gs.wpos=gs.atoms*gs.a;

g.atoms=g.atoms-shift*gs.a*inv(g.a);
g.wpos=g.atoms*g.a;
gs.wpos=gs.atoms*gs.a;

gs.sublattice=match_sublattice_cartesian(gs,g);
% gs.sublattice
% gs.suborbidx=zeros(size(gs.atoms,1),2);
end_idx=0;
for i=1:size(gs.atoms,1)
    orbnum=g.suborbidx(gs.sublattice(i),2)-g.suborbidx(gs.sublattice(i),1)+1;
    start_idx=end_idx+1;
    end_idx=start_idx+orbnum-1;
    gs.suborbidx(i,:)=[start_idx,end_idx];
    gs.orbnum_list=[gs.orbnum_list,orbnum];
end

plot_geometry_sub(gs, g)
plot_geometry(gs, g)

n1=max(Umatrix(:,1));
n2=max(Umatrix(:,2));
n3=max(Umatrix(:,3));
nx=max(g.hopr(:,1));
ny=max(g.hopr(:,2));
nz=max(g.hopr(:,3));
hop_dim=[ceil(nx/n1),ceil(ny/n2),ceil(nz/n3)];
nhoprs=(2*hop_dim(1)+1)*(2*hop_dim(2)+1)*(2*hop_dim(3)+1);
gs.ham=zeros(norbitals,norbitals,nhoprs);
mm=1;

for i=-hop_dim(1):hop_dim(1)
    for j=-hop_dim(2):hop_dim(2)
        for k=-hop_dim(3):hop_dim(3)
            gs.hopr(mm,:)=[i,j,k];
            gs.ham(:,:,mm)=get_t(gs,g,i,j,k,norbitals);
            mm=mm+1;
        end
    end
end

end



function [rs, sl] = generate_supercell_atoms(g, gs, c, dim)
% 扩展晶胞原子坐标，并筛选落在 supercell 内的原子
%
% 输入：
%   g      - 原始结构体（应包含字段 g.wpos, g.a, g.sublattice）
%   gs     - supercell结构体（应包含字段 gs.a）
%   c      - 扩展范围（邻居层数，例如 c=1 表示 -1:1）
%   dim    - 系统维度（1, 2 或 3）
%
% 输出：
%   rs     - 落在 supercell 内的原子坐标（笛卡尔）
%   sl     - 对应的 sublattice 标签（如有）

% 默认边界容差
d0 = -0.122132112;
d1 = 1.0 + d0;

% 坐标变换矩阵
L = inv(gs.a);

% 获取邻居胞的平移索引
inds = neighbor_cells(c, dim);

% 初始化输出
rs = [];
sl = [];

for idx = 1:size(inds,1)
    i = inds(idx,1);
    j = inds(idx,2);
    k = inds(idx,3);

    for ir = 1:size(g.wpos,1)
        ri = g.wpos(ir, :);                  % 原胞中原子坐标
        rj = ri + [i,j,k] * g.a;             % 平移后的位置
        rn = rj * L;                         % 转为 supercell 分数坐标

        n1 = rn(1); n2 = rn(2); n3 = rn(3);

        store = false;

        if dim == 3 && all([n1 n2 n3] > d0) && all([n1 n2 n3] < d1)
            store = true;
        elseif dim == 2 && all([n1 n2] > d0) && all([n1 n2] < d1)
            store = true;
        elseif dim == 1 && n1 > d0 && n1 < d1
            store = true;
        end

        if store
            rs(end+1,:) = rj;
            sl(end+1,1) = g.sublattice(ir);
        end
    end
end
end


function cells = neighbor_cells(num,dim)
% Return neighboring cell indices within `num` layers in `dim` dimensions.
% Sorted by distance from the origin [0 0 0]

if nargin < 2
    dim = 3; % default to 3D
end

cells = []; % initialize

if dim == 0
    return
elseif dim == 1
    for i = -num:num
        cells(end+1, :) = [i 0 0];
    end
elseif dim == 2
    for i = -num:num
        for j = -num:num
            cells(end+1, :) = [i j 0];
        end
    end
elseif dim == 3
    for i = -num:num
        for j = -num:num
            for k = -num:num
                cells(end+1, :) = [i j k];
            end
        end
    end
else
    error('Only supports dim = 0, 1, 2, or 3');
end

% Sort by squared Euclidean distance from [0 0 0]
d2 = sum(cells.^2, 2);
[~, idx] = sort(d2);
cells = cells(idx, :);  % sort from near to far
end



function ham=get_t(gs,g,i,j,k,norbitals)
ham=zeros(norbitals,norbitals);
natoms=size(gs.atoms,1);
for l=1:natoms
    for m=1:natoms
        r1=gs.wpos(l,:);
        r2=gs.wpos(m,:)+[i,j,k]*gs.a;
        r1_sub=g.wpos(gs.sublattice(l),:);
        r2_sub=g.wpos(gs.sublattice(m),:);
        % d_frac_r=(r2-r1-(r2_sub-r1_sub))*inv(g.a)
        d_frac_r=round((r2-r1-(r2_sub-r1_sub))*inv(g.a));
        temp=find(ismember(g.hopr,d_frac_r,'rows'));
        if temp
            ham(gs.suborbidx(l,1):gs.suborbidx(l,2),gs.suborbidx(m,1):gs.suborbidx(m,2))=...
                g.ham(g.suborbidx(gs.sublattice(l),1):g.suborbidx(gs.sublattice(l),2), ...
                g.suborbidx(gs.sublattice(m),1):g.suborbidx(gs.sublattice(m),2),temp);
        end
    end
end
end



function ro = replicate3d(rs, a, n1, n2, n3)
% replicate3d replicates a unit cell in 3D to build a supercell.
%
% Inputs:
%   rs  - nc x 3 matrix of atomic positions (cartesian)
%   a1, a2, a3 - 1x3 lattice vectors
%   n1, n2, n3 - number of repetitions in each direction
%
% Output:
%   ro  - (n1*n2*n3*nc) x 3 matrix of replicated positions
a1=a(1,:);
a2=a(2,:);
a3=a(3,:);
nc = size(rs, 1);  % number of atoms in the original cell
ro = zeros(n1 * n2 * n3 * nc, 3);  % allocate output array
ik = 1;  % MATLAB uses 1-based indexing

for i = 0:(n1-1)
    for j = 0:(n2-1)
        for l = 0:(n3-1)
            for k = 1:nc
                ro(ik, :) = i*a1 + j*a2 + l*a3 + rs(k, :);
                ik = ik + 1;
            end
        end
    end
end
end

function rs1 = brute_expand(rs3, a1, a2, a3, ncheck)
% BRUTE_EXPAND expand atomic positions rs3 by repeated lattice shifts
% until convergence (no new atoms added)
%
% Inputs:
%   rs3 - n0 x 3 array of initial atomic positions (cartesian)
%   a1, a2, a3 - 1x3 lattice vectors
%   ncheck - max number of lattice units to shift in each direction
%
% Output:
%   rs1 - expanded atomic positions

    rs1 = rs3;  % initial positions

    while true
        rs_old = rs1;  % store previous round
        n_before = size(rs1, 1);

        for i = -ncheck:ncheck
            for j = -ncheck:ncheck
                for k = -ncheck:ncheck
                    if i == 0 && j == 0 && k == 0
                        continue;  % skip (0,0,0)
                    end

                    % compute shifted positions
                    rs2 = rs1 + i*a1 + j*a2 + k*a3;

                    % find new atoms that are not in rs1
                    rs1 = return_unique(rs1, rs2);
                    size(rs1)
                end
            end
        end

        % if no new atoms added, stop
        if size(rs1, 1) == n_before
            break;
        end
    end
end

function rout = return_unique(rs1, rs2, tol)
% RETURN_UNIQUE returns only those positions in rs1 that do not appear in rs2
% Inputs:
%   rs1 - n1×3 原子位置
%   rs2 - n2×3 原子位置
%   tol - 容差 (默认为 sqrt(0.001))
%
% Output:
%   rout - 从 rs1 中选出的“唯一”点

    if nargin < 3
        tol = sqrt(0.001);  % default tolerance
    end

    rout = [];
    for i = 1:size(rs1, 1)
        ri = rs1(i, :);  % 当前 rs1 中的原子
        distsq = sum((rs2 - ri).^2, 2);  % 所有 rs2 点到 ri 的平方距离
        if min(distsq) > tol^2
            rout = [rout; ri];  % 只有与所有 rs2 点都不重合时才加入
        end
    end
end



function sublat = match_sublattice_cartesian(gs, g, tol)
% gs.wpos: N x 3 supercell positions (Cartesian)
% g.wpos:  n x 3 primitive positions (Cartesian)
% g.atoms: n x 1 sublattice label
% g.a: 3x3 lattice vectors (each row = a1, a2, a3)
% tol: tolerance, e.g., 1e-3

if nargin < 3
    tol = 1e-3;
end

% lattice vectors (3x3)
A = g.a;

% supercell atom count
N = size(gs.wpos,1);
n = size(g.wpos,1);
sublat = zeros(N,1);

% loop over all supercell atoms
for i = 1:N
    ri = gs.wpos(i,:);  % Cartesian
    frac = (A' \ ri')'; % fractional coordinate of this atom
    frac = frac - floor(frac); % map to [0,1)

    matched = false;
    for j = 1:n
        rj_frac = g.atoms(j,:);  % fractional
        x_mod=mod(frac-rj_frac,[1,1,1]);
        x_norm=norm(x_mod);
        if (abs(x_norm) < tol | abs(x_norm - 1) < tol | abs(x_norm - sqrt(2)) < tol | abs(x_norm - sqrt(3)) < tol)
            sublat(i) = j;  % or g.sublattice(j) if you have labels
            matched = true;
            break;
        end
    end

    if ~matched
        warning("Atom %d not matched!", i);
    end
end
end



function plot_geometry_sub(gs, g)
% 绘制 supercell 结构，并根据原胞 sublattice 编号区分颜色
% 输入：
%   gs.wpos: N×3 supercell 中原子的笛卡尔坐标
%   g.atoms: n×3 原胞原子的分数坐标
%   g.a: 3×3 晶格向量，每行一个 a1,a2,a3
%   g.sublattice: n×1 sublattice 标签（可选）
    % ===================== 绘图 =====================
    figure();
    hold on;

    % 不同子格不同颜色
    utypes = unique(gs.sublattice);
    cmap = lines(length(utypes));  % 颜色表

    for ii = 1:length(utypes)
        sel = (gs.sublattice == utypes(ii));
        scatter3(gs.wpos(sel,1), gs.wpos(sel,2), gs.wpos(sel,3), ...
                 60, cmap(ii,:), 'filled', 'DisplayName', sprintf('Sublattice %d', utypes(ii)));
    end

        % 1.1 给每个原子编号
    for i = 1:size(gs.wpos,1)
        pos = gs.wpos(i,:);
        text(pos(1), pos(2), pos(3), sprintf(' %d', i), 'FontSize', 10, 'Color', 'k');
    end

    % 2. 绘制晶格矢量
    origin = [0, 0, 0];
    quiver3(origin(1), origin(2), origin(3), gs.a(1,1), gs.a(1,2), gs.a(1,3), 0, 'r', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), gs.a(2,1), gs.a(2,2), gs.a(2,3), 0, 'g', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), gs.a(3,1), gs.a(3,2), gs.a(3,3), 0, 'b', 'LineWidth', 2);

    origin = [0, 0, 0];
    quiver3(origin(1), origin(2), origin(3), g.a(1,1), g.a(1,2), g.a(1,3), 0, 'r', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), g.a(2,1), g.a(2,2), g.a(2,3), 0, 'g', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), g.a(3,1), g.a(3,2), g.a(3,3), 0, 'b', 'LineWidth', 2);

    % scatter3(g.wpos(:,1), g.wpos(:,2), g.wpos(:,3), 60, 'filled','ro'); hold on;

    axis equal;
    grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    legend show;
    title('Supercell Structure with Sublattice Labels');
end


function plot_geometry(gs,g)
% 绘制原子位置和晶格矢量，并给原子编号
%
% gs.wpos: N x 3 原子笛卡尔坐标
% gs.a:    3 x 3 晶格矢量（每行一个向量）

    figure()
    
    % 1. 绘制原子坐标
    scatter3(gs.wpos(:,1), gs.wpos(:,2), gs.wpos(:,3), 60, 'filled'); hold on;

    % 1.1 给每个原子编号
    for i = 1:size(gs.wpos,1)
        pos = gs.wpos(i,:);
        text(pos(1), pos(2), pos(3), sprintf(' %d', i), 'FontSize', 10, 'Color', 'k');
    end

    % 2. 绘制晶格矢量
    origin = [0, 0, 0];
    quiver3(origin(1), origin(2), origin(3), gs.a(1,1), gs.a(1,2), gs.a(1,3), 0, 'r', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), gs.a(2,1), gs.a(2,2), gs.a(2,3), 0, 'g', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), gs.a(3,1), gs.a(3,2), gs.a(3,3), 0, 'b', 'LineWidth', 2);

    origin = [0, 0, 0];
    quiver3(origin(1), origin(2), origin(3), g.a(1,1), g.a(1,2), g.a(1,3), 0, 'k', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), g.a(2,1), g.a(2,2), g.a(2,3), 0, 'k', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), g.a(3,1), g.a(3,2), g.a(3,3), 0, 'k', 'LineWidth', 2);

        scatter3(g.wpos(:,1), g.wpos(:,2), g.wpos(:,3), 60, 'filled','ro'); hold on;

    % 3. 图形设置
    axis equal;
    grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Atomic Positions and Lattice Vectors');
    legend('Atoms', 'a_1', 'a_2', 'a_3');
end


% function hk=get_hk(obj,kpoint)
% %    kpoint=kpoint*obj.b
%     hk=zeros(size(obj.ham,1),size(obj.ham,1));
%     nrpts=size(obj.ham,3);
%     for j=1:nrpts
%         temp=obj.ham(:,:,j)*...
%         exp(1j*dot(kpoint,(obj.hopr(j,:)*obj.a))); 
%         hk=hk+temp;
%     end
% end