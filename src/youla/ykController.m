function [K, Kleft] = ykController(cf, Q)
%YKCONTROLLER  Realize a Youla-parameterized controller from coprime factors.
%
%   K = YKCONTROLLER(CF, Q) returns the right-factor realization
%
%       K = (U + M*Q)*inv(V + N*Q)
%
%   where CF is the struct produced by COPRIMEFACTORIZATION or
%   COPRIMEFACTORIZATIONLPV and Q is a stable Youla parameter.
%
%   [K, KLEFT] = YKCONTROLLER(CF, Q) additionally returns the dual
%   realization built from the left factors
%
%       Kleft = (Ut + Q*Mt)*inv(Vt + Q*Nt)
%
%   The two are mathematically identical; computing both is a cheap and
%   effective numerical check on the factorization, and the examples compare
%   them with VERIFYSTEPEQUALITY.
%
%   Key properties:
%     * K stabilizes the plant for ANY stable Q - stability is structural,
%       not the result of an optimization.
%     * Q = 0 recovers the nominal controller K0.
%     * The map from Q to the closed loop is affine, which is what makes
%       switching and interpolating over Q safe.
%
%   See also COPRIMEFACTORIZATION, YOULAPARAMETER, VERIFYSTEPEQUALITY.

arguments
    cf (1,1) struct
    Q
end

K = (cf.U + cf.M*Q) * inv(cf.V + cf.N*Q); %#ok<MINV>

if nargout > 1
    Kleft = (cf.Ut + Q*cf.Mt) * inv(cf.Vt + Q*cf.Nt); %#ok<MINV>
end

end
