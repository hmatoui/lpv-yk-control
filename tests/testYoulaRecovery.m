function tests = testYoulaRecovery()
%TESTYOULARECOVERY  Unit tests for the Youla-Kucera construction.
%
%   Run with:  runtests('tests')
%
%   These tests use a low-order LTI surrogate rather than the vehicle model,
%   so they execute in seconds and need only the Control System Toolbox.
%   They verify the properties the whole framework rests on:
%
%     1. the doubly-coprime factorization identities hold
%     2. Q = 0 recovers the nominal controller exactly
%     3. an arbitrary second stabilizing controller is recovered exactly
%        through its Youla parameter
%     4. the Youla parameter is stable by construction
%     5. ANY stable Q yields a stabilizing controller
%
%   Sign convention. The synthesis interconnection feeds the controller with
%   y - r (see DESIGNGRIDLPVCONTROLLER), so a controller K here closes the
%   loop as u = K*y, i.e. POSITIVE feedback. The coprime factorization uses
%   the same convention (well-posedness term I - Dk*Dg), and so do the
%   observer-based controllers built below and the FEEDBACK(...,+1) calls.
%
%   See also COPRIMEFACTORIZATION, YOULAPARAMETER, YKCONTROLLER.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setupPath();

% An unstable second-order plant and two observer-based stabilizing
% controllers. Note the controllers themselves are unstable (as H-infinity
% controllers for unstable plants often are), which is deliberate: it makes
% the tests exercise the stabilizing gains Fg and Fk0 for real.
s = tf('s');
G = ss((s + 1)/(s^2 + 0.2*s - 1));

A = G.a; B = G.b; C = G.c;

% Observer-based controller for u = K*y: state-feedback poles p, observer
% poles 2p. Realization (A - B*Kp - L*C, L, -Kp, 0).
buildK = @(p) ss(A - B*place(A, B, p) - place(A', C', 2*p)'*C, ...
                 place(A', C', 2*p)', ...
                 -place(A, B, p), 0);

K0 = buildK([-2 -3]);
K1 = buildK([-5 -6]);

% Stabilizing gains for the factorization: A + B*Fg and Ak0 + Bk0*Fk0 must
% be Hurwitz (LMIs (6)-(7) / (28)-(29) in the paper). The pole-placement
% fallback inside COPRIMEFACTORIZATION keeps the existing poles, which is
% only admissible when plant and controller are already stable - not here.
Fg  = -place(A, B, [-4 -5]);
Fk0 = -place(K0.a, K0.b, [-6 -7]);

testCase.TestData.G   = G;
testCase.TestData.K0  = K0;
testCase.TestData.K1  = K1;
testCase.TestData.Fg  = Fg;
testCase.TestData.Fk0 = Fk0;
testCase.TestData.w   = logspace(-2, 3, 200);   % comparison grid [rad/s]
testCase.TestData.tol = 1e-6;
end

function testSurrogateIsWellPosed(testCase)
% Guard the assumptions the other tests rely on.
d = testCase.TestData;
verifyGreaterThan(testCase, max(real(pole(d.G))), 0, 'Surrogate plant should be unstable.');
verifyLessThan(testCase, max(real(pole(closedLoop(d.G, d.K0)))), 0, 'K0 must stabilize G.');
verifyLessThan(testCase, max(real(pole(closedLoop(d.G, d.K1)))), 0, 'K1 must stabilize G.');
verifyLessThan(testCase, max(real(eig(d.G.a  + d.G.b *d.Fg))),  0, 'Fg must stabilize A + B*Fg.');
verifyLessThan(testCase, max(real(eig(d.K0.a + d.K0.b*d.Fk0))), 0, 'Fk0 must stabilize Ak0 + Bk0*Fk0.');
end

function testCoprimeIdentities(testCase)
d  = testCase.TestData;
cf = coprimeFactorization(d.G, d.K0, d.Fg, d.Fk0);

verifyLessThan(testCase, freqRespGap(d.G,  cf.N/cf.M,   d.w), d.tol, 'G  ~= N*inv(M)');
verifyLessThan(testCase, freqRespGap(d.G,  cf.Mt\cf.Nt, d.w), d.tol, 'G  ~= inv(Mt)*Nt');
verifyLessThan(testCase, freqRespGap(d.K0, cf.U/cf.V,   d.w), d.tol, 'K0 ~= U*inv(V)');
verifyLessThan(testCase, freqRespGap(d.K0, cf.Vt\cf.Ut, d.w), d.tol, 'K0 ~= inv(Vt)*Ut');

% All eight factors must live in RH-infinity
names = {'M','N','U','V','Mt','Nt','Ut','Vt'};
for k = 1:numel(names)
    verifyLessThan(testCase, max(real(pole(cf.(names{k})))), 0, ...
        sprintf('Factor %s is not stable.', names{k}));
end
end

function testNominalRecoveredAtZeroQ(testCase)
d  = testCase.TestData;
cf = coprimeFactorization(d.G, d.K0, d.Fg, d.Fk0);

Qzero = ss(zeros(size(d.K0)));
[Krec, KrecLeft] = ykController(cf, Qzero);

verifyLessThan(testCase, freqRespGap(d.K0, Krec,     d.w), d.tol, 'Q = 0 does not recover K0 (right factors).');
verifyLessThan(testCase, freqRespGap(d.K0, KrecLeft, d.w), d.tol, 'Q = 0 does not recover K0 (left factors).');
end

function testLocalControllerRecoveredExactly(testCase)
d  = testCase.TestData;
cf = coprimeFactorization(d.G, d.K0, d.Fg, d.Fk0);
Q  = youlaParameter(d.G, d.K0, d.K1, cf.Fg, cf.Fk0);

[Krec, KrecLeft] = ykController(cf, Q);

verifyLessThan(testCase, freqRespGap(d.K1, Krec, d.w), d.tol, ...
    'Right-factor YK realization does not reproduce K1.');
verifyLessThan(testCase, freqRespGap(d.K1, KrecLeft, d.w), d.tol, ...
    'Left-factor YK realization does not reproduce K1.');
end

function testYoulaParameterIsStable(testCase)
d  = testCase.TestData;
cf = coprimeFactorization(d.G, d.K0, d.Fg, d.Fk0);
Q  = youlaParameter(d.G, d.K0, d.K1, cf.Fg, cf.Fk0);

verifyLessThan(testCase, max(real(eig(Q.a))), 0, ...
    'Q must be stable by construction (Aq is block-triangular with stable diagonal).');

% Q measures the difference between the local and nominal controller
verifyEqual(testCase, Q.d, d.K1.d - d.K0.d, 'AbsTol', 1e-12);
end

function testClosedLoopStableForAnyStableQ(testCase)
% The defining YK property: ANY stable Q yields a stabilizing controller.
d  = testCase.TestData;
cf = coprimeFactorization(d.G, d.K0, d.Fg, d.Fk0);

rng(0);
for k = 1:5
    Qrand = rss(3, 1, 1);            % rss returns a stable model
    K  = ykController(cf, Qrand);
    CL = closedLoop(d.G, K);
    verifyLessThan(testCase, max(real(pole(CL))), 0, ...
        sprintf('Closed loop unstable for stable Q (trial %d).', k));
end
end

% ---- helpers (functiontests only treats local functions whose names -----
% ---- start with "test" as test points) ---------------------------------

function CL = closedLoop(G, K)
% Positive-feedback loop u = K*y, y = G*u - the convention of the framework.
CL = feedback(G*K, 1, +1);
end

function gap = freqRespGap(sys1, sys2, w)
% Largest frequency-response mismatch on the grid w. Unlike STEP this is
% well defined for unstable models, which the surrogate controllers are.
H1  = squeeze(freqresp(sys1, w));
H2  = squeeze(freqresp(sys2, w));
gap = max(abs(H1(:) - H2(:)));
end
