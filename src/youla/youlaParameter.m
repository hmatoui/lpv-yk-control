function Q = youlaParameter(G, K0, Ki, Fg, Fk0)
%YOULAPARAMETER  Youla-Kucera parameter Q mapping the nominal controller onto a local one.
%
%   Q = YOULAPARAMETER(G, K0, KI, FG, FK0) builds the Youla parameter Q that
%   encodes the local controller KI relative to the nominal controller K0,
%   such that
%
%       Ki = (U + M*Q)*inv(V + N*Q)
%
%   with (M,N,U,V) the coprime factors returned by COPRIMEFACTORIZATION.
%   This is eq. (10) of the reference paper in the grid-based case, and
%   eq. (32) evaluated at a vertex in the polytopic case.
%
%   Inputs:
%     G   - plant (LTI or PSS)
%     K0  - nominal / central controller
%     KI  - local controller to be encoded
%     FG  - plant state-feedback gain ([] for LTI pole-placement fallback)
%     FK0 - nominal-controller state-feedback gain ([] for fallback)
%
%   Structure of the result. Writing Ac/Bc/Cc/Dc for K0 and A1/B1/C1/D1 for
%   Ki, the realization is block upper-triangular:
%
%            | A + B*D1*C     B*C1      B*(D1-Dc)*Fk0 - B*Cc |
%       Aq = | B1*C           A1        B1*Fk0               |
%            | 0              0         Ac + Bc*Fk0          |
%
%   The two diagonal blocks are exactly
%     (a) the local closed loop  [A+B*D1*C, B*C1; B1*C, A1], stable by
%         assumption (A.2.1) - i.e. because Ki stabilizes G on its subregion;
%     (b) the stabilized nominal controller dynamics Ac + Bc*Fk0, stable by
%         LMI (7) / (29).
%
%   So Q is stable *by construction*, and Lemma A.1 (bounded off-diagonal
%   terms cannot destabilize a block-triangular system) does the rest. This
%   is precisely why the local designs never have to be co-designed: the only
%   thing the switching proof needs from Ki is that it stabilizes G locally.
%
%   Note also Dq = D1 - Dc: the Youla parameter measures the *difference*
%   between the local and the nominal controller, and Q = 0 iff Ki = K0.
%
%   See also COPRIMEFACTORIZATION, COPRIMEFACTORIZATIONLPV, YKCONTROLLER.

arguments
    G
    K0
    Ki
    Fg  = []
    Fk0 = []
end

A  = G.a;   B  = G.b;   C  = G.c;   D  = G.d;   %#ok<ASGLU>
Ac = K0.a;  Bc = K0.b;  Cc = K0.c;  Dc = K0.d;
A1 = Ki.a;  B1 = Ki.b;  C1 = Ki.c;  D1 = Ki.d;

if isempty(Fg)
    Fg = -place(A, B, pole(G).');
end
if isempty(Fk0)
    Fk0 = -place(Ac, Bc, pole(K0).');
end

Aq = [A + B*D1*C,                    B*C1,                       B*(D1-Dc)*Fk0 - B*Cc;
      B1*C,                          A1,                         B1*Fk0;
      zeros(size(Ac,1), size(C,2)),  zeros(size(Ac,1), size(A1,2)), Ac + Bc*Fk0];

Bq = [B*(D1 - Dc);
      B1;
      Bc];

Cq = [D1*C - Fg,  C1,  (D1 - Dc)*Fk0 - Cc];

Dq = D1 - Dc;

Q = ss(Aq, Bq, Cq, Dq);

end
