clc;
clear;

% 可选：清理并重启并行池（若无并行工具箱，自动退化为串行）
p = gcp('nocreate');
if ~isempty(p)
    delete(p);
end
try
    parpool('local', 6);
catch ME
    warning('parpool not started (%s). Continue in serial mode.', ME.message);
end

%% 1) 构造 RhG 模型参数
% 这部分沿用原脚本设定，后续所有 nonlinear conductivity 计算都复用该模型。
g = MTB.geometry("RhG");
a = 2.46; % Angstrom
g.a = [ ...
    1/2, -sqrt(3)/2, 0; ...
    1/2,  sqrt(3)/2, 0; ...
    0,    0,         1/a] * a;
g.b = inv(g.a') * 2*pi;

Layer_N = 2;
pars = struct( ...
    'v0', 3.16, ...
    'gamma1', 0.46, ...
    'gamma2', -0.017, ...
    'gamma3', -0.30, ...
    'gamma4', -0.086, ...
    'uext', 0.015, ...
    'delta', -0.0011, ...
    'xi', 1, ...
    'N', Layer_N);

model = @MTB.ham.get_ham_kp_RhG;
g.dim_kp = Layer_N * 2;
nbands = g.dim_kp;

%% 2) 生成 2D k 网格并求解 U(k), E(k)
% 说明：
% - Kx/Ky 用于 2D 面内积分；
% - U/E 是后续 tbNLC 统一入口的直接输入。
knum = 800;
kxline = [-0.03, 0.03];
kyline = [-0.03, 0.03];
[Kx, Ky, Kz] = g.get_Bulk2Dkmesh(kxline, kyline, knum);

tic;
[Unk, Enk, ~, ~, ~, ~] = MTB.ham.get_bulk_plane_kp_velocity_withH( ...
    pars, model, nbands, Kx, Ky, Kz);
toc;

%% 3) 费米能、频率网格与展宽
Ef = tbNLC.calculate_ef(Enk(:), 0.5); % 半填充示例
kT = 0.0;          % eV, 0 表示 0K 阶跃分布
eta = 5e-4;        % eV, Lorentzian 展宽
Eph_list = linspace(0.0, 0.1, 5000); % eV

%% 4) 配置 sigma 与 MIC 的独立参数
% 你可以分别控制两个通道的细节（例如 symBC / positiveDE）。
sigma_opts = struct();
sigma_opts.band_list = 2:3;   % 给 shift current 用
sigma_opts.periodicFD = false;
sigma_opts.trimBoundary = true;
sigma_opts.symBC = true;
sigma_opts.verbose = true;
sigma_opts.doGaugeFix = true;
sigma_opts.g_s = 1;

mic_opts = struct();
mic_opts.band_list = 2:3;     % 给 MIC 用
mic_opts.periodicFD = false;
mic_opts.trimBoundary = true;
mic_opts.verbose = true;
mic_opts.doGaugeFix = true;
mic_opts.g_s = 1;
mic_opts.positiveDE = true;
mic_opts.saveFullMN = false; % 若只关心最终响应，建议 false 以节省内存

run_opts = struct();
run_opts.do_sigma = false;    % true: 计算 shift current
run_opts.do_mic = true;      % true: 计算 magnetic injection current
run_opts.sigma_opts = sigma_opts;
run_opts.mic_opts = mic_opts;

%% 5) 统一调用 tbNLC 工具箱
% 这是本次重构的核心：主文件不再直接维护大量局部函数，
% 只保留“输入准备 + 调用 + 后处理”。
tic;
[result, out_nlc] = tbNLC.compute_nonlinear_conductivity_fromUE( ...
    Kx, Ky, Unk, Enk, Eph_list, Ef, kT, eta, run_opts);
toc;

% 按需取出响应张量
% sigma_abc = result.sigma_abc; % 尺寸: 2 x 2 x 2 x Nw
eta_abc = result.eta_abc;     % 尺寸: 2 x 2 x 2 x Nw

%% 6) 作图示例（可按你的单位制继续调整缩放因子）
% plot_response_tensor(Eph_list, sigma_abc, 1e5, 'Shift Current \sigma_{abc}');
plot_response_tensor(Eph_list, eta_abc, 1e-8, 'MIC Metric \eta_{abc}');

%% 7) 保存结果（建议保留 out_nlc，便于后续排查数值细节）
% save('nonlinear_conductivity_result.mat', ...
%     'Eph_list', 'sigma_abc', 'eta_abc', 'Ef', 'kT', 'eta', ...
%     'run_opts', 'out_nlc', '-v7.3');

function plot_response_tensor(Eph_list, tensor_abcw, scale_factor, fig_title)
%PLOT_RESPONSE_TENSOR 将 2x2x2xNw 张量分量全部画在同一张图上。
%
% 输入:
%   Eph_list     : 频率点，1 x Nw
%   tensor_abcw  : 2 x 2 x 2 x Nw
%   scale_factor : 统一缩放因子
%   fig_title    : 图标题

    figure('Name', fig_title, 'Color', 'w');
    hold on;
    comp = ['x', 'y'];
    for a = 1:2
        for b = 1:2
            for c = 1:2
                y = squeeze(tensor_abcw(a,b,c,:)) * scale_factor;
                plot(Eph_list, real(y), 'LineWidth', 1.1, ...
                    'DisplayName', sprintf('%s%s%s', comp(a), comp(b), comp(c)));
            end
        end
    end
    xlabel('\hbar\omega (eV)');
    ylabel('Response (scaled)');
    title(fig_title);
    legend('Location', 'best');
    grid on;
    box on;
end
