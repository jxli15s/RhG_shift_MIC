classdef HFMF_BandProj
%HFMF_BANDPROJ  Band-projected HF on uniform patch (FFT Fock), PARFOR-only.
%
% Storage strategy:
%   All k-dependent objects use flat index ik=1..Nk for speed and clarity:
%     rho_band_flat : Nact x Nact x Nk
%     rho_orb_flat  : Norb x Norb x Nk
%     Sigma_*_flat  : Norb x Norb x Nk
%     Sigma_band    : Nact x Nact x Nk
%     Ek_flat       : Nact x Nk
%     W_flat        : Nact x Nact x Nk
%
% FFT convolution is the only place requiring Nky x Nkx pages; we reshape per (a,b) page
% inside fock_fft_atomgauge_flat.

  properties
    g
    opts
    mesh
    hk
    act
    V
    ph
    rho_band_flat
    mu
    history
  end

  methods
    function obj = HFMF_BandProj(g, opts)
      obj.g = g;
      obj.opts = opts;

      obj.mesh = tbHFMF.build_kmesh_patch(g, opts.K0_frac, opts.L_frac, opts.Nk);
      obj.hk   = tbHFMF.build_hk_atomgauge(g, obj.mesh);

      obj.act  = tbHFMF.diag_hk_active(obj.hk, opts.act_idx);

      obj.V    = opts.Vbuilder(g, obj.mesh, opts.Vpars);
      obj.ph   = tbHFMF.precompute_dtau_phases(obj.mesh, g.wpos);

      obj = obj.init_rho_band_flat();

      obj.mu = 0;
      obj.history = struct();
    end

    function obj = init_rho_band_flat(obj)
      Nact = obj.act.Nact;
      Nk   = obj.mesh.Nk;

      mu0 = mean(obj.act.eps_flat(:));
      f0  = 1 ./ (1 + exp((obj.act.eps_flat - mu0)/obj.opts.kT));  % Nact x Nk

      rho = zeros(Nact,Nact,Nk,'like',1+1i);
      for ik = 1:Nk
        rho(:,:,ik) = diag(f0(:,ik));
      end

      if isfield(obj.opts,'seed_act') && ~isempty(obj.opts.seed_act)
        seed = obj.opts.seed_act;
        for ik = 1:Nk
          R = rho(:,:,ik) + seed;
          rho(:,:,ik) = (R + R')/2;
        end
      end

      obj.rho_band_flat = rho;
    end

    function out = run(obj)
      opts = obj.opts;
      mesh = obj.mesh;
      act  = obj.act;
      V    = obj.V;
      ph   = obj.ph;

      Norb = size(obj.g.ham,1);
      Nact = act.Nact;
      Nk   = mesh.Nk;

      Uflat   = act.U_flat;      % Norb x Nact x Nk
      epsflat = act.eps_flat;    % Nact x Nk

      rho_band_flat = obj.rho_band_flat;

      for it = 1:opts.max_iter

        % (1) rho_band -> rho_orb  (PARFOR)
        rho_orb_flat = zeros(Norb,Norb,Nk,'like',1+1i);
        parfor ik = 1:Nk
          U  = Uflat(:,:,ik);
          rb = rho_band_flat(:,:,ik);
          rho_orb_flat(:,:,ik) = U * rb * U';
        end

        % (2) Fock via FFT (centered<->FFT handled inside)
        SigmaF_orb_flat = tbHFMF.fock_fft_atomgauge_flat(mesh, V, ph, rho_orb_flat);

        % (3) Hartree (k-independent)
        if isfield(opts,'use_hartree') && opts.use_hartree
          SigmaH = tbHFMF.hartree_orbital_flat(V.V0, rho_orb_flat, mesh, opts);
        else
          SigmaH = zeros(Norb,Norb,'like',1+1i);
        end

        Sigma_orb_flat = SigmaF_orb_flat + SigmaH;  % implicit expansion

        % (4) Project to active-band space  (PARFOR)
        Sigma_band_flat = zeros(Nact,Nact,Nk,'like',1+1i);
        parfor ik = 1:Nk
          U = Uflat(:,:,ik);
          S = U' * Sigma_orb_flat(:,:,ik) * U;
          Sigma_band_flat(:,:,ik) = (S + S')/2;
        end

        % (5) Diagonalize Heff in active space (PARFOR)
        Ek_flat = zeros(Nact,Nk);
        W_flat  = zeros(Nact,Nact,Nk,'like',1+1i);
        parfor ik = 1:Nk
          Heff = diag(epsflat(:,ik)) + Sigma_band_flat(:,:,ik);
          Heff = (Heff + Heff')/2;

          [W,D] = eig(Heff,'vector');
          [E,ord] = sort(real(D),'ascend');
          W = W(:,ord);

          Ek_flat(:,ik)  = E;
          W_flat(:,:,ik) = W;
        end

        % (6) μ from filling constraint
        mu = tbHFMF.solve_mu_fill_flat(Ek_flat, mesh, opts);

        % (7) Update rho_band (PARFOR)
        rho_new_flat = zeros(Nact,Nact,Nk,'like',1+1i);
        parfor ik = 1:Nk
          W = W_flat(:,:,ik);
          f = 1 ./ (1 + exp((Ek_flat(:,ik) - mu)/opts.kT));
          R = W * diag(f) * W';
          rho_new_flat(:,:,ik) = (R + R')/2;
        end

        % mixing + convergence
        rho_mixed = (1-opts.mix)*rho_band_flat + opts.mix*rho_new_flat;
        err = max(abs(rho_mixed(:) - rho_band_flat(:)));
        rho_band_flat = rho_mixed;

        obj.history.it(it)  = it;
        obj.history.err(it) = err;
        obj.history.mu(it)  = mu;

        if isfield(opts,'verbose') && opts.verbose
          fprintf('it=%d  err=%.3e  mu=%.6f  pref_patch=%.3e\n', it, err, mu, mesh.pref_patch);
        end

        if err < opts.tol
          break;
        end
      end

      obj.rho_band_flat = rho_band_flat;
      obj.mu = mu;

      out = struct();
      out.rho_band_flat = rho_band_flat;
      out.mu = mu;
      out.history = obj.history;
    end

    function rho4 = rho_band_4d(obj)
      % convenience view for plotting
      Nact = obj.act.Nact;
      Nky  = obj.mesh.Nky;
      Nkx  = obj.mesh.Nkx;
      rho4 = reshape(obj.rho_band_flat, Nact, Nact, Nky, Nkx);
    end
  end
end
