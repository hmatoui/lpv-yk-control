function Pgen = stateFeedbackPlant(sys, ctrlWeight, asPss)
%STATEFEEDBACKPLANT  Build a generalized plant for state-feedback synthesis.
%
%   PGEN = STATEFEEDBACKPLANT(SYS, CTRLWEIGHT) wraps the system SYS into the
%   two-input / two-output generalized plant required by LPVSFSYN or
%   LMISTATEFEEDBACKGRID:
%
%       xdot = A*x + [0  B]*[w; u]
%        z   = [ 0 ] * x + [ 0   ctrlWeight ] * [w]
%              [ C ]       [ 0   0          ]   [u]
%
%   i.e. the controlled output penalizes the control effort (weighted by
%   CTRLWEIGHT) and the system output. Minimizing the induced L2 norm of this
%   channel yields a stabilizing state-feedback gain with a sensible
%   bandwidth, which is all the Youla construction needs.
%
%   PGEN = STATEFEEDBACKPLANT(SYS, CTRLWEIGHT, ASPSS) returns a PSS when
%   ASPSS is true (default: true when SYS is parameter varying).
%
%   Why this exists. Both Fg and Fk0 in the grid-based scheme are obtained by
%   solving the same kind of state-feedback problem - once for the plant
%   G(rho) (LMI 6) and once for the nominal controller dynamics K0(rho)
%   (LMI 7). Rather than duplicating the augmentation inline twice, it is
%   factored out here.
%
%   Reference values used in the paper's vehicle example:
%     plant G(rho)          : ctrlWeight = 1.8
%     nominal controller K0 : ctrlWeight = 1.0
%
%   See also LMISTATEFEEDBACKGRID, COPRIMEFACTORIZATIONLPV.

arguments
    sys
    ctrlWeight (1,1) double = 1
    asPss      (1,1) logical = true
end

nOut = size(sys.c, 1);
nIn  = size(sys.b, 2);

Bgen = [zeros(size(sys.b,1), 1), sys.b];
Cgen = [zeros(nOut, size(sys.a,2)); sys.c];
Dgen = [zeros(nOut, 1),  ctrlWeight*ones(nOut, nIn);
        zeros(nOut, 1),  zeros(nOut, nIn)];

if asPss
    Pgen = pss(sys.a, Bgen, Cgen, Dgen);
else
    Pgen = ss(sys.a, Bgen, Cgen, Dgen);
end

end
