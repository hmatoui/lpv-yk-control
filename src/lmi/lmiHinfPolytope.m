function [listK, listCL, gopt] = lmiHinfPolytope(listP, nmeas, ncon, sp, percentage, solver)
%LMIHINFPOLYTOPE  Polytopic LPV H-infinity output-feedback synthesis.
%
%   [LISTK, LISTCL, GOPT] = LMIHINFPOLYTOPE(LISTP, NMEAS, NCON, SP, PERCENTAGE, SOLVER)
%   solves the polytopic LPV H-infinity problem (Scherer & Weiland
%   formulation) for the generalized plant given as a cell array of vertex
%   models:
%
%       xdot     A   B1  B2     x
%        z    =  C1  D11 D12    w
%        y       C2  D21 D22    u
%
%   A single constant Lyapunov pair (X, Y) certifies all vertices, i.e.
%   quadratic stability - this is assumption (A.3.1) of the reference paper.
%
%   Inputs:
%     LISTP      - 1-by-NV cell array of vertex generalized plants
%     NMEAS      - number of measured outputs
%     NCON       - number of control inputs
%     SP         - 1 for a strictly proper controller (Dk = 0), 0 otherwise
%     PERCENTAGE - relaxation added to the optimal gamma before the second
%                  (conditioning) solve, in percent
%     SOLVER     - YALMIP solver name, e.g. 'sedumi', 'sdpt3'
%
%   Outputs:
%     LISTK  - cell array of vertex controllers
%     LISTCL - cell array of vertex closed loops
%     GOPT   - achieved gamma
%
%   The routine validates the polytopic structural requirements: B2, C2, D12
%   and D21 must be parameter independent and D22 must vanish. These are the
%   model assumptions stated in Sec. III-A of the reference paper.
%
%   Original author: Charles Poussot-Vassal.
%   Reproduced here with light cleanup (YALMIP syntax refresh, clearer
%   diagnostics). Requires YALMIP and an SDP solver.

%%% Size of the generalized plant
sizeX = size(listP{1}.a, 1);
sizeZ = size(listP{1}, 1) - nmeas;
sizeY = nmeas;
sizeW = size(listP{1}, 2) - ncon;
sizeU = ncon;

nV = numel(listP);

% Strict inequalities are not directly representable in YALMIP
epsi = 1e-6;

%%% Partition the vertex plants
A = cell(1,nV); B1 = cell(1,nV); B2 = cell(1,nV);
C1 = cell(1,nV); D11 = cell(1,nV); D12 = cell(1,nV);
C2 = cell(1,nV); D21 = cell(1,nV); D22 = cell(1,nV);

for i = 1:nV
    A{i}   = listP{i}.a(1:sizeX, 1:sizeX);
    B1{i}  = listP{i}.b(1:sizeX, 1:sizeW);
    B2{i}  = listP{i}.b(1:sizeX, sizeW+1:sizeW+sizeU);
    C1{i}  = listP{i}.c(1:sizeZ, 1:sizeX);
    D11{i} = listP{i}.d(1:sizeZ, 1:sizeW);
    D12{i} = listP{i}.d(1:sizeZ, sizeW+1:sizeW+sizeU);
    C2{i}  = listP{i}.c(sizeZ+1:sizeZ+sizeY, 1:sizeX);
    D21{i} = listP{i}.d(sizeZ+1:sizeZ+sizeY, 1:sizeW);
    D22{i} = listP{i}.d(sizeZ+1:sizeZ+sizeY, sizeW+1:sizeW+sizeU);
end

%%% Structural checks required by the polytopic formulation
for i = 1:nV
    assert(isequal(B2{1},  B2{i}),  'lmiHinfPolytope:B2Varying', ...
        'B2 must be parameter independent (vertex %d differs).', i);
    assert(isequal(D12{1}, D12{i}), 'lmiHinfPolytope:D12Varying', ...
        'D12 must be parameter independent (vertex %d differs).', i);
    assert(isequal(C2{1},  C2{i}),  'lmiHinfPolytope:C2Varying', ...
        'C2 must be parameter independent (vertex %d differs).', i);
    assert(isequal(D21{1}, D21{i}), 'lmiHinfPolytope:D21Varying', ...
        'D21 must be parameter independent (vertex %d differs).', i);
    assert(all(D22{i}(:) == 0),     'lmiHinfPolytope:D22NonZero', ...
        'D22 must be zero (vertex %d is non-zero).', i);
end

%%% Decision variables
At = cell(1,nV); Bt = cell(1,nV); Ct = cell(1,nV); Dt = cell(1,nV);
for i = 1:nV
    At{i} = sdpvar(sizeX, sizeX, 'full');
    Bt{i} = sdpvar(sizeX, sizeY, 'full');
    Ct{i} = sdpvar(sizeU, sizeX, 'full');
    if sp == 1
        Dt{i} = zeros(sizeU, sizeY);
    else
        Dt{i} = sdpvar(sizeU, sizeY, 'full');
    end
end
X     = sdpvar(sizeX, sizeX, 'symmetric');
Y     = sdpvar(sizeX, sizeX, 'symmetric');
gamma = sdpvar(1, 1, 'full');

%%% LMI set
H0 = [X eye(sizeX); eye(sizeX) Y];
F  = (H0 >= epsi);

M11 = cell(1,nV); M21 = cell(1,nV); M22 = cell(1,nV);
M31 = cell(1,nV); M32 = cell(1,nV);
M41 = cell(1,nV); M42 = cell(1,nV); M43 = cell(1,nV);

for i = 1:nV
    M11{i} = A{i}*X + X*A{i}' + B2{i}*Ct{i} + (B2{i}*Ct{i})';
    M21{i} = At{i} + (A{i} + B2{i}*Dt{i}*C2{i})';
    M22{i} = Y*A{i} + A{i}'*Y + Bt{i}*C2{i} + (Bt{i}*C2{i})';
    M31{i} = (B1{i} + B2{i}*Dt{i}*D21{i})';
    M32{i} = (Y*B1{i} + Bt{i}*D21{i})';
    M41{i} = C1{i}*X + D12{i}*Ct{i};
    M42{i} = C1{i} + D12{i}*Dt{i}*C2{i};
    M43{i} = D11{i} + D12{i}*Dt{i}*D21{i};

    H = [M11{i} M21{i}'  M31{i}'              M41{i}';
         M21{i} M22{i}   M32{i}'              M42{i}';
         M31{i} M32{i}  -gamma*eye(sizeW)     M43{i}';
         M41{i} M42{i}   M43{i}              -gamma*eye(sizeZ)];

    F = [F, H <= -epsi]; %#ok<AGROW>
end

%%% Minimize gamma
ops = sdpsettings('solver', solver, 'verbose', 0);
optimize(F, gamma, ops);
gopt = value(gamma)*(1 + percentage/100);

%%% Second pass: improve conditioning of the reconstruction at fixed gamma
alpha = sdpvar(1, 1, 'full');
H0 = [X -alpha*eye(sizeX); -alpha*eye(sizeX) Y];
F  = (H0 >= epsi);
for i = 1:nV
    H = [M11{i} M21{i}'  M31{i}'             M41{i}';
         M21{i} M22{i}   M32{i}'             M42{i}';
         M31{i} M32{i}  -gopt*eye(sizeW)     M43{i}';
         M41{i} M42{i}   M43{i}             -gopt*eye(sizeZ)];
    F = [F, H <= -epsi]; %#ok<AGROW>
end
optimize(F, alpha, ops);

%%% Extract the numerical solution
for i = 1:nV
    At{i} = value(At{i});
    Bt{i} = value(Bt{i});
    Ct{i} = value(Ct{i});
    Dt{i} = value(Dt{i});
end
X = value(X);
Y = value(Y);

if any(eig(X) <= 0) || any(eig(Y) <= 0)
    warning('lmiHinfPolytope:nonPositiveLyapunov', ...
        'Lyapunov matrices are not positive definite; the solution is unreliable.');
end

%%% Reconstruct the vertex controllers
[Usvd, S, Vsvd] = svd(eye(sizeX) - X*Y);
R = chol(S);
Mm = Usvd*R';
Nn = (R*Vsvd')';

listK  = cell(1,nV);
listCL = cell(1,nV);
for i = 1:nV
    Dc = Dt{i};
    Cc = (Ct{i} - Dc*C2{i}*X)/Mm';
    Bc = Nn\(Bt{i} - Y*B2{i}*Dc);
    Ac = Nn\(At{i} - Y*A{i}*X - Y*B2{i}*Dc*C2{i}*X ...
             - Nn*Bc*C2{i}*X - Y*B2{i}*Cc*Mm')/Mm';

    P = ss(listP{i}.a, listP{i}.b, listP{i}.c, listP{i}.d);
    listK{i}  = ss(Ac, Bc, Cc, Dc);
    listCL{i} = lft(P, listK{i});
end

end
