function sysOut = equalizeOrder(sysIn, padPole)
%EQUALIZEORDER  Pad an LTI array so every model has the same state dimension.
%
%   SYSOUT = EQUALIZEORDER(SYSIN) inflates every model in the LTI array SYSIN
%   to the largest state dimension found in the array, by appending fast,
%   decoupled modes:
%
%       A <- blkdiag(A, padPole*I),   B <- [B; 0],   C <- [C, 0],   D <- D
%
%   SYSOUT = EQUALIZEORDER(SYSIN, PADPOLE) uses PADPOLE (default -1e4) as the
%   eigenvalue of the appended modes.
%
%   Why this is needed. The grid-based LPV-YK controller stacks the N Youla
%   parameters into a single block-diagonal bank
%
%       Aq = diag(Aq_1, ..., Aq_N)
%
%   which a fixed-size real-time realization (Simulink / dSPACE) requires to
%   have a constant state dimension across the whole parameter grid. Local
%   controllers designed independently generally do not have identical
%   McMillan degree, so the shorter ones are padded here.
%
%   The padded modes are both uncontrollable and unobservable (zero rows in B
%   and zero columns in C) and are placed far in the left half plane, so the
%   input-output behaviour is unchanged and the added dynamics decay
%   essentially instantaneously. Verify with VERIFYSTEPEQUALITY if in doubt.
%
%   See also VERIFYSTEPEQUALITY, YOULAPARAMETER.

arguments
    sysIn
    padPole (1,1) double {mustBeNegative} = -1e4
end

sysOut  = sysIn;
nModels = nmodels(sysIn);   % NOT numel: an LTI array is a single object

sizeMax = 0;
for i = 1:nModels
    sizeMax = max(sizeMax, size(sysIn(:,:,i).a, 1));
end

for i = 1:nModels
    n = size(sysIn(:,:,i).a, 1);
    if n < sizeMax
        pad = sizeMax - n;
        A = blkdiag(sysIn(:,:,i).a, padPole*eye(pad));
        B = [sysIn(:,:,i).b; zeros(pad, size(sysIn(:,:,i).b, 2))];
        C = [sysIn(:,:,i).c, zeros(size(sysIn(:,:,i).c, 1), pad)];
        sysOut(:,:,i) = ss(A, B, C, sysIn(:,:,i).d);
    end
end

end
