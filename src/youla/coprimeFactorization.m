function cf = coprimeFactorization(G, K0, Fg, Fk0)
%COPRIMEFACTORIZATION  Doubly-coprime factorization of an LTI plant/controller pair.
%
%   CF = COPRIMEFACTORIZATION(G, K0, FG, FK0) computes a doubly-coprime
%   factorization of the plant G and the stabilizing controller K0 over
%   RH-infinity, returning all eight factors in the struct CF:
%
%       G  = N *inv(M)  = inv(Mt)*Nt
%       K0 = U *inv(V)  = inv(Vt)*Ut
%
%   Fields of CF: M, N, U, V, Mt, Nt, Ut, Vt  (the "t" prefix denotes the
%   left, or "tilde", factors), plus Fg and Fk0 as actually used.
%
%   Inputs:
%     G   - LTI plant (state-space)
%     K0  - LTI nominal/central stabilizing controller (state-space)
%     FG  - state-feedback gain stabilizing (A + B*Fg). Pass [] to fall back
%           on pole placement at the plant's own poles.
%     FK0 - state-feedback gain stabilizing (Ak0 + Bk0*Fk0). Pass [] to fall
%           back on pole placement at the controller's own poles.
%
%   The [] fallbacks keep the EXISTING poles, so they only produce factors in
%   RH-infinity when G and K0 are themselves stable. For an unstable plant or
%   controller, pass explicit stabilizing gains (e.g. from PLACE, LQR or the
%   LMIs in src/lmi).
%
%   Sign convention: positive feedback, u = K0*y, matching the "-r+G+n"
%   measurement of the synthesis interconnections (see docs/theory.md 11.2).
%
%   The construction is the standard observer/state-feedback one: Fg and Fk0
%   are the two design freedoms that shape the factorization and hence the
%   dynamics of the Youla parameter. Stability of the resulting closed loop
%   holds for *any* admissible pair, but the transient behaviour does depend
%   on the choice - this is flagged as open in the reference paper's
%   conclusion.
%
%   Every controller stabilizing G is then recovered as
%
%       K(Q) = (U + M*Q)*inv(V + N*Q)
%
%   for some stable Q, with Q = 0 giving back K0.
%
%   See also COPRIMEFACTORIZATIONLPV, YOULAPARAMETER, YKCONTROLLER.

arguments
    G
    K0
    Fg  = []
    Fk0 = []
end

Ac0 = K0.a;  Bc0 = K0.b;  Cc0 = K0.c;  Dc0 = K0.d;
At  = G.a;   Bt  = G.b;   Ct  = G.c;   Dt  = G.d;

if isempty(Fg)
    Fg = -place(At, Bt, pole(G).');
end
if isempty(Fk0)
    Fk0 = -place(Ac0, Bc0, pole(K0).');
end

% ---- Right coprime factors: [M; N] and [U; V] -------------------------
Ar = blkdiag(At + Bt*Fg, Ac0 + Bc0*Fk0);
Br = blkdiag(Bt, Bc0);
Cr = [Fg,         Cc0 + Dc0*Fk0;
      Ct + Dt*Fg, Fk0];
Dr = [eye(size(Dc0)), Dc0;
      Dt,             eye(size(Dt))];
sysR = ss(Ar, Br, Cr, Dr);

nU = size(Fg, 1);      % rows of the "upper" output block
nY = size(Ct, 1);      % rows of the "lower" output block
nG = size(Bt, 2);      % columns of the "left"  input block
nK = size(Bc0, 2);     % columns of the "right" input block

cf.M = sysR(1:nU,        1:nG);
cf.N = sysR(end-nY+1:end, 1:nG);
cf.U = sysR(1:nU,        end-nK+1:end);
cf.V = sysR(end-nY+1:end, end-nK+1:end);

% ---- Left coprime factors: [Vt -Ut] and [-Nt Mt] ----------------------
Y1 = inv(eye(size(Dc0)) - Dc0*Dt);
Z1 = inv(eye(size(Dt))  - Dt*Dc0);

Al = [At + Bt*Y1*Dc0*Ct,  Bt*Y1*Cc0;
      Bc0*Z1*Ct,          Ac0 + Bc0*Z1*Dt*Cc0];
Bl = [-Bt*Y1,             Bt*Y1*Dc0;
      -Bc0*Z1*Dt,         Bc0*Z1];
Cl = [Fg - Y1*Dc0*Ct,    -Y1*Cc0;
      Z1*Ct,             -(Fk0 - Z1*Dt*Cc0)];
Dl = [Y1,    -Y1*Dc0;
      Z1*Dt,  Z1];
sysL = ss(Al, Bl, Cl, Dl);

nY1 = size(Y1, 2);
nZ1 = size(Z1, 1);

cf.Vt =  sysL(1:nU,          1:nY1);
cf.Nt = -sysL(end-nZ1+1:end, 1:nY1);
cf.Ut = -sysL(1:nU,          end-nZ1+1:end);
cf.Mt =  sysL(end-nZ1+1:end, end-nZ1+1:end);

cf.Fg  = Fg;
cf.Fk0 = Fk0;

end
