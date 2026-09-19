function cf = coprimeFactorizationLpv(G, K0, Fg, Fk0)
%COPRIMEFACTORIZATIONLPV  Doubly-coprime factorization of an LPV plant/controller pair.
%
%   CF = COPRIMEFACTORIZATIONLPV(G, K0, FG, FK0) is the parameter-varying
%   counterpart of COPRIMEFACTORIZATION. G and K0 are LPVTools PSS objects
%   and FG, FK0 are parameter-varying state-feedback gains (PMAT or PSS)
%   satisfying LMIs (6) and (7) of the reference paper:
%
%       Fg(rho)  = V(rho) *inv(Xg(rho))     stabilizes  A(rho)  + B2(rho)*Fg(rho)
%       Fk0(rho) = W(rho) *inv(Xk0(rho))    stabilizes  AK0(rho) + Bk0(rho)*Fk0(rho)
%
%   The returned struct CF has the same eight fields as COPRIMEFACTORIZATION
%   (M, N, U, V, Mt, Nt, Ut, Vt), each a PSS defined over the same parameter
%   domain as G, together with Fg and Fk0.
%
%   Unlike the LTI version, no pole-placement fallback is provided: for an
%   LPV plant the stabilizing gains must come from the pLMI synthesis, since
%   pointwise pole placement carries no parameter-varying stability
%   guarantee.
%
%   See also COPRIMEFACTORIZATION, YOULAPARAMETER, LMISTATEFEEDBACKGRID.

arguments
    G
    K0
    Fg
    Fk0
end

Ac0 = K0.a;  Bc0 = K0.b;  Cc0 = K0.c;  Dc0 = K0.d;
At  = G.a;   Bt  = G.b;   Ct  = G.c;   Dt  = G.d;

% ---- Right coprime factors -------------------------------------------
Ar = [At + Bt*Fg,                             zeros(size(At,1),  size(Ac0,2));
      zeros(size(Ac0,1), size(At,2)),         Ac0 + Bc0*Fk0];
Br = [Bt,                            zeros(size(Bt,1),  size(Bc0,2));
      zeros(size(Bc0,1), size(Bt,2)), Bc0];
Cr = [Fg,         Cc0 + Dc0*Fk0;
      Ct + Dt*Fg, Fk0];
Dr = [eye(size(Dc0)), Dc0;
      Dt,             eye(size(Dt))];
sysR = ss(Ar, Br, Cr, Dr);

nU = size(Fg, 1);
nY = size(Ct, 1);
nG = size(Bt, 2);
nK = size(Bc0, 2);

cf.M = sysR(1:nU,                     1:nG);
cf.N = sysR(size(sysR,1)-nY+1:size(sysR,1), 1:nG);
cf.U = sysR(1:nU,                     size(sysR,2)-nK+1:size(sysR,2));
cf.V = sysR(size(sysR,1)-nY+1:size(sysR,1), size(sysR,2)-nK+1:size(sysR,2));

% ---- Left coprime factors --------------------------------------------
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

cf.Vt =  sysL(1:nU,                           1:nY1);
cf.Nt = -sysL(size(sysL,1)-nZ1+1:size(sysL,1), 1:nY1);
cf.Ut = -sysL(1:nU,                           size(sysL,2)-nZ1+1:size(sysL,2));
cf.Mt =  sysL(size(sysL,1)-nZ1+1:size(sysL,1), size(sysL,2)-nZ1+1:size(sysL,2));

cf.Fg  = Fg;
cf.Fk0 = Fk0;

end
