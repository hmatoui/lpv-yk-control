function [K, CL, gopt] = lmiStateFeedbackGrid(listP, rho, maxdrho, nmeas, ncon, sp, percentage, solver)
%LMISTATEFEEDBACKGRID  Grid-based LPV state-feedback synthesis with rate bounds.
%
%   [K, CL, GOPT] = LMISTATEFEEDBACKGRID(LISTP, RHO, MAXDRHO, NMEAS, NCON, SP, PERCENTAGE, SOLVER)
%   solves the rate-bounded parameter-dependent LMIs for an LPV state-feedback
%   gain on a grid, using the affine Lyapunov basis
%
%       X(rho) = X0 + rho*X1,     Y(rho) = Y0 + rho*Y1,     K(rho) = -Y(rho)*inv(X(rho))
%
%   The generalized plant is supplied as a cell array of models evaluated at
%   the grid points RHO:
%
%       xdot    A   Bw  Bu    x
%        z   =  Cz  Dzw Dzu   w
%
%   This implements LMIs (6)-(7) of the reference paper. Each grid point
%   contributes TWO inequalities, one for each extreme of the parameter rate
%   +/- MAXDRHO, which is the {nu_min, nu_max} notation used in the paper:
%   requiring the inequality at both vertices of the rate box implies it
%   throughout, by convexity in dX/drho.
%
%   Inputs:
%     LISTP      - 1-by-NR cell array of plants, one per grid point
%     RHO        - vector of grid points, same length as LISTP
%     MAXDRHO    - bound on |d(rho)/dt|
%     NMEAS      - number of trailing output rows of LISTP that are NOT
%                  controlled outputs; they are discarded, so the first
%                  size(LISTP{1},1) - NMEAS rows are taken as z. Pass 0 when
%                  every output of LISTP is a performance channel.
%     NCON       - number of control inputs
%     SP         - unused for state feedback, kept for signature symmetry
%     PERCENTAGE - relaxation added to the optimal gamma, in percent
%     SOLVER     - YALMIP solver name, e.g. 'sdpt3', 'sedumi'
%
%   Outputs:
%     K    - cell array of state-feedback gains, one per grid point
%     CL   - cell array of closed loops (u = -K*x)
%     GOPT - achieved gamma
%
%   Sign convention: K is returned such that the stabilized dynamics are
%   A - Bu*K. Callers that use the "+ B2*Fg" convention of the paper, as
%   COPRIMEFACTORIZATIONLPV does, must negate it.
%
%   Original author: Kazusa Yamada (GIPSA-lab), 2016.
%   Reproduced here with light cleanup. Requires YALMIP and an SDP solver.

%%% Size of the generalized plant
sizeX = size(listP{1}.a, 1);
sizeZ = size(listP{1}, 1) - nmeas;
sizeW = size(listP{1}, 2) - ncon;
sizeU = ncon;
sizeR = numel(rho);

assert(sizeR == numel(listP), 'lmiStateFeedbackGrid:mismatch', ...
    'Number of grid points (%d) and number of plants (%d) must match.', ...
    sizeR, numel(listP));

epsi = 1e-6;

%%% Partition the plants
A = cell(1,sizeR); Bw = cell(1,sizeR); Bu = cell(1,sizeR);
Cz = cell(1,sizeR); Dzw = cell(1,sizeR); Dzu = cell(1,sizeR);
for i = 1:sizeR
    A{i}   = listP{i}.a(1:sizeX, 1:sizeX);
    Bw{i}  = listP{i}.b(1:sizeX, 1:sizeW);
    Bu{i}  = listP{i}.b(1:sizeX, sizeW+1:sizeW+sizeU);
    Cz{i}  = listP{i}.c(1:sizeZ, 1:sizeX);
    Dzw{i} = listP{i}.d(1:sizeZ, 1:sizeW);
    Dzu{i} = listP{i}.d(1:sizeZ, sizeW+1:sizeW+sizeU);
end

%%% Affine Lyapunov / gain bases: X(rho) = X0 + rho*X1, Y(rho) = Y0 + rho*Y1
Y0 = sdpvar(sizeU, sizeX);
Y1 = sdpvar(sizeU, sizeX);
X0 = sdpvar(sizeX, sizeX, 'symmetric');
X1 = sdpvar(sizeX, sizeX, 'symmetric');

X = cell(1,sizeR); Y = cell(1,sizeR); dX = cell(1,sizeR);
for i = 1:sizeR
    X{i}  = X0 + rho(i)*X1;
    Y{i}  = Y0 + rho(i)*Y1;
    dX{i} = X1;                 % dX/drho for the affine basis
end

gamma = sdpvar(1, 1, 'full');

%%% LMI set
F = (gamma >= epsi);
for i = 1:sizeR
    F = [F, X{i} >= epsi]; %#ok<AGROW>
end

for i = 1:sizeR
    common = A{i}*X{i} + X{i}*A{i}' + Bu{i}*Y{i} + (Bu{i}*Y{i})';
    M21 = Bw{i}';
    M31 = Cz{i}*X{i} + Dzu{i}*Y{i};
    M32 = Dzw{i};

    % Both extremes of the parameter rate box
    for sgn = [+1, -1]
        M11 = common + sgn*dX{i}*maxdrho;
        H = [M11    M21'                M31';
             M21   -gamma*eye(sizeW)    M32';
             M31    M32                -gamma*eye(sizeZ)];
        F = [F, H <= -epsi]; %#ok<AGROW>
    end
end

%%% Minimize gamma
ops = sdpsettings('solver', solver, 'verbose', 0);
optimize(F, gamma, ops);
gopt = value(gamma)*(1 + percentage/100);

%%% Extract the gains
K  = cell(1,sizeR);
CL = cell(1,sizeR);
for i = 1:sizeR
    Xi = value(X{i});
    Yi = value(Y{i});
    if any(eig(Xi) <= 0)
        warning('lmiStateFeedbackGrid:nonPositiveLyapunov', ...
            'X(rho) is not positive definite at grid point %d.', i);
    end
    K{i}  = -Yi/Xi;
    CL{i} = ss(A{i} - Bu{i}*K{i}, Bw{i}, Cz{i} - Dzu{i}*K{i}, Dzw{i});
end

end
