function gs=get_xy_open_wannier(g,n1,n2)
% Slove the Hamiltonian on the K-mesh
% obj: The object of the continuum model system
% Kx,Ky: the k-mesh we will solve
% dim_H: num of the bands/orbits
% Unk: Norbits*Nands*knum*knum
% Enk: knum*knum*Nbands
% kz=0;
gs = MTB.geometry("supercell");
gs.a=[n1,0,0;0,n2,0;0,0,1]*g.a;
gs.b=inv(gs.a')*2*pi;
norbitalc = size(g.wpos,1); % number of orbitals in a cell
norbitals = size(g.wpos,1)*n1*n2; %total number of orbitals
natomc = size(g.atoms,1); % number of atoms in a cell
natoms = natomc*n1*n2; %total number of atoms
% gs.atoms=zeros(natoms,3);
% gs.wpos=zeros(norbitals,3);


ncell = n1 * n2;
atoms_cell = cell(ncell, 1);
wpos_cell  = cell(ncell, 1);
cellindex = zeros(ncell, 3);

parfor idx = 1:ncell
    [i, j] = ind2sub([n1, n2], idx);
    shift_frac = [i - 1, j - 1, 0];
    shift_cart = shift_frac * g.a;

    atom_pos = (g.atoms + shift_frac) * g.a / gs.a;
    orb_pos  = g.wpos + shift_cart;
    % cellindex(idx,:) = shift_frac * g.a / gs.a;
    atoms_cell{idx} = atom_pos;
    wpos_cell{idx}  = orb_pos;
end

% 合并结果
gs.atoms = vertcat(atoms_cell{:});
gs.wpos  = vertcat(wpos_cell{:});

% 构造 hopping 哈希表：key = mat2str([dx,dy,dz]), value = index
hopr_map = containers.Map();
for h = 1:size(g.hopr,1)
    key = mat2str(g.hopr(h,:));
    hopr_map(key) = h;
end

% 并行构造所有 dz 的 hopping
nz = max(abs(g.hopr(:,3)));
dz_list = -nz : nz;
ndz = length(dz_list);
gs.hopr=zeros(ndz,3);
gs.ham =cell(1,ndz);
row_cell = cell(1, ndz);
col_cell = cell(1, ndz);
val_cell = cell(1, ndz);
parfor idx = 1:ndz
    dz = dz_list(idx);
    % rz = [0,0,dz] * sa;
    % [rr, cc, vv] = get_t_z_sparse_parallel2(gs, g, dz, cellindex, norbitalc, hopr_map);
    % [rr, cc, vv] = get_t_z_sparse_fast(gs, g, dz, cellindex, norbitalc, hopr_map);
    [rr,cc,vv] = get_t_z_sparse_fast_pruned_v2(gs, g, dz, n1, n2, norbitalc, hopr_map);
    % vv = vv * exp(1i * kz * rz');  % 统一加相因子
    row_cell{idx} = rr;
    col_cell{idx} = cc;
    val_cell{idx} = vv;
end

for idx=1:ndz
    dz = dz_list(idx);
    gs.hopr(idx,:)=[0,0,dz];
    gs.ham{idx} = sparse(row_cell{idx},col_cell{idx},val_cell{idx});
end

% 汇总稀疏矩阵所有元素
% row = vertcat(row_cell{:});
% col = vertcat(col_cell{:});
% val = vertcat(val_cell{:});
% assert(isequal(length(row), length(col), length(val)), 'Row, col, val length mismatch');
% gs.ham = sparse(row, col, val, norbitals, norbitals);
end


function [row, col, val] = get_t_z_sparse_parallel(gs, g, dz, cellindex, norbitalc, hopr_map)
% 返回 dz 方向的所有 hopping 项（三元组表示）

ncell = size(cellindex,1);
row = []; col = []; val = [];

parfor l = 1:ncell
    for m = 1:ncell
        r1 = cellindex(l,:);
        r2 = cellindex(m,:) + [0, 0, dz];
        d_frac_r = round((r2 - r1) * gs.a / g.a);  % 转为原胞中的分数坐标
        key = mat2str(d_frac_r);
        if hopr_map.isKey(key)
            temp = hopr_map(key);
            tblock = g.ham(:,:,temp);  % 原始 hopping block
            [Srow, Scol, Sval] = find(tblock);  % 稀疏三元组
            i0 = (l - 1) * norbitalc;
            j0 = (m - 1) * norbitalc;
            row = [row; i0 + Srow];
            col = [col; j0 + Scol];
            val = [val; Sval];
        end
    end
end
end

function [row, col, val] = get_t_z_sparse_fast_pruned_v2(gs, g, dz, n1, n2, norbitalc, hopr_map)
% 剪枝 + 邻接配对优化版：构造指定 dz 的 hopping 三元组 (row, col, val)

% 1. 获取有效 (l, m) 配对
[l_list, m_list] = build_adjacency(g, dz, n1, n2);
ncell = n1 * n2;

% 2. 初始化存储
row_cell = cell(length(l_list), 1);
col_cell = cell(length(l_list), 1);
val_cell = cell(length(l_list), 1);

% 3. 构造 hopping 跳跃项
parfor p = 1:length(l_list)
    l = l_list(p);
    m = m_list(p);

    [ix, iy] = ind2sub([n1, n2], l);
    [jx, jy] = ind2sub([n1, n2], m);

    r1 = [(ix - 1), (iy - 1), 0] * g.a / gs.a;
    r2 = [(jx - 1), (jy - 1), dz] * g.a / gs.a;
    d_frac_r = round((r2 - r1) * gs.a / g.a);
    key = mat2str(d_frac_r);

    if hopr_map.isKey(key)
        temp = hopr_map(key);
        tblock = g.ham(:,:,temp);
        [Srow, Scol, Sval] = find(tblock);
        i0 = (l - 1) * norbitalc;
        j0 = (m - 1) * norbitalc;
        row_cell{p} = i0 + Srow;
        col_cell{p} = j0 + Scol;
        val_cell{p} = Sval;
    else
        row_cell{p} = [];
        col_cell{p} = [];
        val_cell{p} = [];
    end
end

% 4. 汇总所有 hopping 三元组
row = vertcat(row_cell{:});
col = vertcat(col_cell{:});
val = vertcat(val_cell{:});
end

function [l_list, m_list] = build_adjacency(g, dz, n1, n2)
% 构建在给定 dz 方向下，xy 超胞内的 hopping 配对列表 (l, m)

ncell = n1 * n2;

% 每个 cell 的 xy 坐标
[i_grid, j_grid] = ndgrid(0:n1-1, 0:n2-1);
cell_pos = [i_grid(:), j_grid(:)];  % ncell × 2

% 所有 (l,m) 配对
[I1, I2] = ndgrid(1:ncell, 1:ncell);
lm_pair = [I1(:), I2(:)];

% 所有 dx,dy
delta = cell_pos(I2(:), :) - cell_pos(I1(:), :);  % 每对的 (dx, dy)

% 当前 dz 对应的 hopping 方向子集
hopr_xy = g.hopr(g.hopr(:,3)==dz, 1:2);  % 所有在此 dz 下的 dx dy
hopr_xy = unique(hopr_xy, 'rows');

% 判断哪些配对 (l,m) 符合 hopping 范围
is_valid = ismember(delta, hopr_xy, 'rows');

% 输出有效配对
l_list = lm_pair(is_valid, 1);
m_list = lm_pair(is_valid, 2);
end



% %
% % function ham=get_t_sparse(gs,g,k,cellindex,norbitalc,norbitals)
% %     ham=zeros(norbitals,norbitals);
% %     for l=1:size(cellindex,1)
% %         for m=1:size(cellindex,1)
% %             r1=cellindex(l,:);
% %             r2=cellindex(m,:)+[0,0,k];
% %             d_frac_r=round((r2-r1)*gs.a/g.a);
% %             temp=find(ismember(g.hopr,d_frac_r,'rows'));
% %            if temp
% %                ham((l-1)*norbitalc+1:l*norbitalc,(m-1)*norbitalc+1:m*norbitalc)=g.ham(:,:,temp);
% %            end
% %         end
% %     end
% % end

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