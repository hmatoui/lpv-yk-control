function [K, CL, gopt, P] = designPolytopicLpvController(Gv, We, Wu, opts)
%DESIGNPOLYTOPICLPVCONTROLLER  Polytopic LPV H-infinity output-feedback synthesis.
%
%   [K, CL, GOPT, P] = DESIGNPOLYTOPICLPVCONTROLLER(GV, WE, WU) synthesizes a
%   polytopic LPV controller for the vertex plants in the cell array GV,
%   using the same mixed-sensitivity interconnection as
%   DESIGNGRIDLPVCONTROLLER and solving the vertex LMIs with a single
%   constant Lyapunov pair (quadratic stability).
%
%   [K, CL, GOPT, P] = DESIGNPOLYTOPICLPVCONTROLLER(GV, WE, WU, OPTS) accepts
%   the name-value options:
%     StrictlyProper - enforce Dk = 0 (default true)
%     Percentage     - gamma relaxation before the conditioning solve, in
%                      percent (default 0)
%     Solver         - YALMIP solver name (default 'sedumi')
%
%   Outputs:
%     K    - cell array of vertex controllers, to be interpolated with
%            POLYTOPICINTERPOLATION at run time
%     CL   - cell array of vertex closed loops
%     GOPT - achieved gamma
%     P    - cell array of vertex generalized plants
%
%   This produces the nominal controller K0 of assumption (A.3.1): one LPV
%   controller quadratically stabilizing the plant over the *entire* convex
%   region P0. It is expected to be conservative - that is precisely why the
%   Youla parameters are then used to recover local performance.
%
%   Requires YALMIP, an SDP solver and the Robust Control Toolbox.
%
%   See also DESIGNGRIDLPVCONTROLLER, DESIGNLTICONTROLLER, LMIHINFPOLYTOPE.

arguments
    Gv (1,:) cell
    We
    Wu
    opts.StrictlyProper (1,1) logical = true
    opts.Percentage     (1,1) double  = 0
    opts.Solver         (1,:) char    = 'sedumi'
end

Wes = ss(We);
Wus = ss(Wu);

nV = numel(Gv);
P  = cell(1, nV);

for i = 1:nV
    G = Gv{i}; %#ok<NASGU>

    systemnames  = 'Wes Wus G'; %#ok<NASGU>
    inputvar     = '[ r(1); n(1); u(1)]'; %#ok<NASGU>
    outputvar    = '[ Wes; Wus; -r+G+n]'; %#ok<NASGU>
    input_to_Wes = '[-r+G+n]'; %#ok<NASGU>
    input_to_Wus = '[u]'; %#ok<NASGU>
    input_to_G   = '[u]'; %#ok<NASGU>
    sysoutname   = 'pp'; %#ok<NASGU>
    cleanupsysic = 'yes'; %#ok<NASGU>
    sysic;

    P{i} = pp;
end

nmeas = size(Gv{1}.c, 1);
ncon  = 1;

[K, CL, gopt] = lmiHinfPolytope(P, nmeas, ncon, ...
    double(opts.StrictlyProper), opts.Percentage, opts.Solver);

end
