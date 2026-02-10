function [resp, out] = compute_nonlinear_conductivity_fromUE( ...
    Kx, Ky, U, E, Eph_list, Ef, kT, eta, opts)
%COMPUTE_NONLINEAR_CONDUCTIVITY_FROMUE
%==========================================================================
% 统一入口：基于同一套 U/E 输入，一次性计算 nonlinear conductivity 响应。
%
% 支持的响应：
% 1) Shift current: sigma_abc
% 2) Magnetic injection current (metric form): eta_abc
%
% 输入参数与原始接口保持一致，便于直接替换主脚本中的调用。
%
% opts 字段：
%   do_sigma   (logical, default=true ) : 是否计算 sigma_abc
%   do_mic     (logical, default=true ) : 是否计算 eta_abc
%   sigma_opts (struct,  default=struct()) : 传给 shift current 子函数
%   mic_opts   (struct,  default=struct()) : 传给 MIC 子函数
%
% 返回：
%   resp.sigma_abc  (可选)
%   resp.eta_abc    (可选)
%   out.shift       shift current 细节输出
%   out.mic         MIC 细节输出
%   out.meta        运行配置摘要
%==========================================================================

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

    do_sigma   = tbNLC.get_opt(opts, 'do_sigma', true);
    do_mic     = tbNLC.get_opt(opts, 'do_mic', true);
    sigma_opts = tbNLC.get_opt(opts, 'sigma_opts', struct());
    mic_opts   = tbNLC.get_opt(opts, 'mic_opts', struct());

    if ~(do_sigma || do_mic)
        error('At least one of opts.do_sigma / opts.do_mic must be true.');
    end

    resp = struct();
    out = struct();

    if do_sigma
        [resp.sigma_abc, out.shift] = tbNLC.shift_current_plane_fd_energy_skew_fromUE( ...
            Kx, Ky, U, E, Eph_list, Ef, kT, eta, sigma_opts);
    end

    if do_mic
        [resp.eta_abc, out.mic] = tbNLC.mic_metric_plane_fromUE_mn( ...
            Kx, Ky, U, E, Eph_list, Ef, kT, eta, mic_opts);
    end

    out.meta = struct( ...
        'do_sigma', do_sigma, ...
        'do_mic', do_mic, ...
        'num_omega', numel(Eph_list), ...
        'kgrid_size', size(Kx));
end

