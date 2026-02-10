function [H, Hx, Hy, Hz] = get_ham_kp_RhG(kpoint, pars)
% H_abcN_kp  ABC N-layer kp Hamiltonian for valley K_xi (xi=±1), N=2..5 (or larger)
% Basis: (A1,B1, A2,B2, ..., AN,BN)  size: 2N x 2N.
% Inputs:
%   kpoint = [qx; qy]  in 1/Ang
%   pars: struct with fields v0,gamma3,gamma4,gamma1,gamma2,uext,delta,xi,N
% Units:
%   v0, v3, v4 in eV*Ang (after rescale); gamma1,2, delta, uext in eV.

    arguments
        kpoint (3,1) double
        pars struct
    end

    a0   = 2.46;                % Ang
    % unpack
    v0    = pars.v0    * sqrt(3)*a0/2;   % eV*Ang
    v3    = pars.gamma3* sqrt(3)*a0/2;   % eV*Ang
    v4    = pars.gamma4* sqrt(3)*a0/2;   % eV*Ang
    g1    = pars.gamma1;                % eV
    g2    = pars.gamma2;                % eV
    uext  = pars.uext;                  % eV
    delta = pars.delta;                 % eV
    xi    = pars.xi;                    % +1 (K) or -1 (K')
    N     = pars.N;

    % complex momenta
    qx = kpoint(1);  qy = kpoint(2);
    pi  =  xi*qx + 1i*qy;
    pid =  xi*qx - 1i*qy;

    dim = 2*N;
    H   = zeros(dim, dim);
    Hx  = zeros(dim, dim);   % ∂H/∂q_x
    Hy  = zeros(dim, dim);   % ∂H/∂q_y
    Hz  = zeros(dim, dim);   % = 0 (no q_z dependence)

    % 2x2 blocks
    H0   = [ 0,     v0*pid; ...
             v0*pi, 0      ];
    dH0x = [ 0,     v0*xi; ...
             v0*xi, 0      ];
    dH0y = [ 0,    -1i*v0; ...
             1i*v0, 0      ];

    H1   = [ v4*pi,   g1; ...
             v3*pid,  v4*pi ];
    dH1x = [ v4*xi,   0; ...
             v3*xi,   v4*xi ];
    dH1y = [ 1i*v4,   0; ...
            -1i*v3,   1i*v4 ];

    H1d   = H1';
    dH1dx = dH1x';    % 导数块的厄米共轭
    dH1dy = dH1y';

    H2   = [ 0,     0; ...
             g2/2,  0 ];
    H2d  = H2';
    % dH2x = 0; dH2y = 0

    % ---- assemble block-tridiagonal H, Hx, Hy
    for l = 1:N
        idx = 2*(l-1)+(1:2);

        % intralayer
        H (idx,idx)  = H (idx,idx)  + H0;
        Hx(idx,idx)  = Hx(idx,idx)  + dH0x;
        Hy(idx,idx)  = Hy(idx,idx)  + dH0y;

        % nearest interlayer l <-> l+1
        if l < N
            idx2 = 2*l+(1:2);
            % (l+1,l): H1
            H (idx2,idx)  = H (idx2,idx)  + H1;
            Hx(idx2,idx)  = Hx(idx2,idx)  + dH1x;
            Hy(idx2,idx)  = Hy(idx2,idx)  + dH1y;
            % (l,l+1): H1^\dagger
            H (idx, idx2) = H (idx, idx2) + H1d;
            Hx(idx, idx2) = Hx(idx, idx2) + dH1dx;
            Hy(idx, idx2) = Hy(idx, idx2) + dH1dy;
        end

        % next-nearest interlayer l <-> l+2
        if l < N-1
            idx3 = 2*(l+1)+(1:2);
            H (idx3,idx)  = H (idx3,idx)  + H2;
            H (idx, idx3) = H (idx, idx3) + H2d;
            % derivatives are zero for H2
        end
    end

    % surface shift δ on A1 and BN (k-independent)
    H(1,1)       = H(1,1)       + delta;     % A1
    H(2*N,2*N)   = H(2*N,2*N)   + delta;     % BN

    % linear interlayer potential U (k-independent)
    if N == 1
        s = 0;
    else
        core = (N-3)/(N-1);
        s = [-1,  repmat(-core, 1, N-2),  1];   % length N
    end
    s2 = reshape([s; s], 1, []);            % [s1 s1 s2 s2 ... sN sN]
    Udiag = (uext/2) * s2;
    H = H + diag(Udiag);

    % (optional) enforce Hermiticity numerically
    H  = (H + H')/2;
    Hx = (Hx + Hx')/2;
    Hy = (Hy + Hy')/2;
    % Hz already zero
end

% function H = get_ham_kp_RhG(kpoint, pars)
% % H_abcN_kp  ABC N-layer kp Hamiltonian for valley K_xi (xi=±1), N=2..5 (or larger)
% % Basis: (A1,B1, A2,B2, ..., AN,BN)  -> size 2N x 2N.
% % Units: q in 1/Ang; v0,v3,v4 in eV*Ang; gamma1,gamma2, uext, delta in eV.
% 
%     arguments
%         kpoint (2,1) double
%         pars struct
%     end
%     % hbar = 6.582119569e-16;     % eV*s
%     % ang  = 1e-10;               % m
%     a0   = 2.46;
%     % --- unpack parameters
%     v0    = pars.v0*sqrt(3)*a0/2;         % eV*Ang
%     v3    = pars.gamma3*sqrt(3)*a0/2;         % eV*Ang
%     v4    = pars.gamma4*sqrt(3)*a0/2;         % eV*Ang
%     g1    = pars.gamma1;     % eV
%     g2    = pars.gamma2;     % eV
%     uext  = pars.uext;       % eV  (external interlayer potential scale)
%     delta = pars.delta;      % eV  (surface on-site shift)
%     xi    = pars.xi;         % 1 for k valley; -1 for k' valley
%     N     = pars.N;          % N number of layers
% 
%     % complex momenta
%     qx=kpoint(1);
%     qy=kpoint(2);
%     pi  =  xi*qx + 1i*qy;    % p_x + i p_y
%     pid =  xi*qx - 1i*qy;    % p_x - i p_y
% 
%     dim = 2*N;
%     H   = zeros(dim, dim);
% 
%     % ---- block helpers (2x2 each)
%     H0 = [ 0,     v0*pid; ...
%            v0*pi, 0      ];
% 
%     H1 = [ v4*pi,   g1; ...
%            v3*pid,  v4*pi ];
% 
%     H1d = H1';   % Hermitian conjugate
% 
%     H2  = [ 0,     0; ...
%             g2/2,  0 ];
% 
%     H2d = H2';   % conj transpose
% 
%     % ---- fill block-tridiagonal H_N
%     for l = 1:N
%         % intralayer block
%         idx = 2*(l-1)+(1:2);
%         H(idx, idx) = H(idx, idx) + H0;
% 
%         % nearest interlayer l <-> l+1 (H1 and H1^\dagger)
%         if l < N
%             idx2 = 2*l+(1:2);
%             H(idx2, idx) = H(idx2, idx) + H1;    % (l+1,l): H1
%             H(idx,  idx2)= H(idx,  idx2)+ H1d;   % (l,l+1): H1^\dagger
%         end
% 
%         % next-nearest interlayer l <-> l+2 (H2 and H2^\dagger)
%         if l < N-1
%             idx3 = 2*(l+1)+(1:2);
%             H(idx3, idx) = H(idx3, idx) + H2;    % (l+2,l): H2
%             H(idx,  idx3)= H(idx,  idx3)+ H2d;   % (l,l+2): H2^\dagger
%         end
%     end
% 
%     % ---- add surface shift H_delta = delta*diag(1,0,...,0,1)
%     H(1,1)       = H(1,1)       + delta;   % A1
%     H(2*N,2*N)   = H(2*N,2*N)   + delta;   % BN
% 
%     % ---- add linear layer potential U = (uext/2)*diag(...)
%     %     pattern per layer l=1..N: s_l repeated for A_l and B_l
%     if N == 1
%         s = 0;
%     else
%         core = (N-3)/(N-1);
%         s = [-1,  repmat(-core, 1, N-2),  1];   % length N
%     end
%     % expand to 2N by repeating each entry twice
%     s2 = reshape([s; s], 1, []);  % [s1 s1 s2 s2 ... sN sN]
%     Udiag = (uext/2) * s2;
%     H = H + diag(Udiag);
% end

        
            
            
            
            