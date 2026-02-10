function rH = hermitianize_rmn(r_mn)
%HERMITIANIZE_RMN
% 在带指标 (m,n) 上做厄米化:
%   rH = 0.5 * (r + r^\dagger)
% 其中 dagger 表示共轭并交换 m<->n。

    rH = 0.5 * (r_mn + permute(conj(r_mn), [1 2 4 3 5]));
end

