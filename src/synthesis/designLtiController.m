function [K, CL, gam, P] = designLtiController(G, We, Wu)
%DESIGNLTICONTROLLER  LTI H-infinity synthesis at a single operating point.
%
%   [K, CL, GAM, P] = DESIGNLTICONTROLLER(G, WE, WU) synthesizes an LTI
%   H-infinity controller for the frozen plant G with the same mixed
%   sensitivity interconnection used by the LPV designs, so that the local
%   and global designs are directly comparable.
%
%   Inputs:
%     G  - LTI plant at one polytope vertex
%     WE - performance weight (see WEIGHTPERFORMANCE)
%     WU - control weight (see WEIGHTCONTROL)
%
%   Outputs:
%     K   - LTI controller
%     CL  - weighted closed loop
%     GAM - achieved gamma
%     P   - generalized plant
%
%   These are the vertex controllers Kij of assumption (A.3.2). Two points
%   from the reference paper matter when using this function:
%
%   1. Every Kij must be designed with the SAME weights and the SAME LMI
%      formulation, otherwise the controllers meeting at a shared vertex will
%      not coincide.
%   2. At a boundary shared by two adjacent subsets, the vertex controller
%      must be designed ONCE and reused - K_{i,3} == K_{i+1,2}. That identity
%      is what makes switching across the shared face a no-op, which is the
%      property validated experimentally by the hysteretic switching episode
%      in Sec. IV-C.
%
%   Requires the Robust Control Toolbox.
%
%   See also DESIGNPOLYTOPICLPVCONTROLLER, DESIGNGRIDLPVCONTROLLER.

arguments
    G
    We
    Wu
end

Wes = ss(We);
Wus = ss(Wu);
Glpv = G; %#ok<NASGU>   % sysic resolves this name from the local workspace

systemnames   = 'Wes Wus Glpv'; %#ok<NASGU>
inputvar      = '[ r(1); n(1); u(1)]'; %#ok<NASGU>
outputvar     = '[ Wes; Wus; -r+Glpv+n]'; %#ok<NASGU>
input_to_Wes  = '[-r+Glpv+n]'; %#ok<NASGU>
input_to_Wus  = '[u]'; %#ok<NASGU>
input_to_Glpv = '[u]'; %#ok<NASGU>
sysoutname    = 'P'; %#ok<NASGU>
cleanupsysic  = 'yes'; %#ok<NASGU>
sysic;

nmeas = size(G.c, 1);
ncon  = 1;

[K, ~, gam] = hinfsyn(P, nmeas, ncon);
CL = lft(P, K);

end
