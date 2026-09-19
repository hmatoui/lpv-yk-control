function [K, gopt] = lmiHinfStateFeedbackPolytope(listP, nstate, ncon, sp, percentage, solver)
%LMIHINFSTATEFEEDBACKPOLYTOPE  Polytopic LPV H-infinity state-feedback synthesis.
%
%   [K, GOPT] = LMIHINFSTATEFEEDBACKPOLYTOPE(LISTP, NSTATE, NCON, SP, PERCENTAGE, SOLVER)
%   computes vertex state-feedback gains for the polytopic generalized plant
%
%       xdot    A   Bw  Bu    x
%        z   =  Cz  Dzw Dzu   w
%
%   using a single constant Lyapunov matrix X shared across all vertices and
%   per-vertex gain variables Y{i}, so that K{i} = Y{i}*inv(X).
%
%   This implements LMIs (28)-(29) of the reference paper, which supply the
%   state-feedback gains Fg and Fk0 needed to build the coprime
%   factorization in the partitioned polytopic LPV-YK scheme.
%
%   Inputs:
%     LISTP      - 1-by-NV cell array of vertex generalized plants
%     NSTATE     - number of state variables (kept for signature symmetry)
%     NCON       - number of control inputs
%     SP         - kept for signature symmetry
%     PERCENTAGE - relaxation added to the optimal gamma, in percent
%     SOLVER     - YALMIP solver name
%
%   Outputs:
%     K    - cell array of vertex state-feedback gains
%     GOPT - achieved gamma
%
%   Original author: Olivier Sename (GIPSA-lab), Nov 2020.
%   Reproduced here with light cleanup. Requires YALMIP and an SDP solver.

sizeX = size(listP{1}.a, 1);
sizeZ = size(listP{1}, 1);
sizeW = size(listP{1}, 2) - ncon;
sizeU = ncon;
nV    = numel(listP);

epsi = 1e-6;

A = cell(1,nV); Bw = cell(1,nV); Bu = cell(1,nV);
Cz = cell(1,nV); Dzw = cell(1,nV); Dzu = cell(1,nV);
for i = 1:nV
    A{i}   = listP{i}.a(1:sizeX, 1:sizeX);
    Bw{i}  = listP{i}.b(1:sizeX, 1:sizeW);
    Bu{i}  = listP{i}.b(1:sizeX, sizeW+1:sizeW+sizeU);
    Cz{i}  = listP{i}.c(1:sizeZ, 1:sizeX);
    Dzw{i} = listP{i}.d(1:sizeZ, 1:sizeW);
    Dzu{i} = listP{i}.d(1:sizeZ, sizeW+1:sizeW+sizeU);
end

Y = cell(1,nV);
for i = 1:nV
    Y{i} = sdpvar(sizeU, sizeX);
end
X     = sdpvar(sizeX, sizeX, 'symmetric');
gamma = sdpvar(1, 1, 'full');

F = (X >= epsi);
for i = 1:nV
    M11 = A{i}*X + X*A{i}' + Bu{i}*Y{i} + (Bu{i}*Y{i})';
    M21 = Bw{i}';
    M31 = Cz{i}*X + Dzu{i}*Y{i};
    M32 = Dzw{i};
    H = [M11    M21'              M31';
         M21   -eye(sizeW)        M32';
         M31    M32              -gamma*eye(sizeZ)];
    F = [F, H <= -epsi]; %#ok<AGROW>
end

ops = sdpsettings('solver', solver, 'verbose', 0);
optimize(F, gamma, ops);

gopt = sqrt(value(gamma))*(1 + percentage/100);

Xv = value(X);
if any(eig(Xv) <= 0)
    warning('lmiHinfStateFeedbackPolytope:nonPositiveLyapunov', ...
        'X is not positive definite; the solution is unreliable.');
end

K = cell(1,nV);
for i = 1:nV
    K{i} = value(Y{i})/Xv;
end

end
