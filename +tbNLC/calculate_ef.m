function efermi = calculate_ef(Enk, u)
%CALCULATE_EF
%==========================================================================
% 根据填充因子 u 估算费米能级：
%   从所有能量本征值中取前 ceil(numel(Enk)*u) 个最低态，
%   其最大值作为费米能级。
%==========================================================================

    arguments
        Enk {mustBeNumeric}
        u (1,1) double {mustBeGreaterThanOrEqual(u,0), mustBeLessThanOrEqual(u,1)}
    end

    total_states = numel(Enk);
    occupied_states = ceil(total_states * u);
    efermi = max(mink(Enk(:), occupied_states));
end

