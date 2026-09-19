function [K, CL, gam, P] = designGridLpvController(Glpv, We, Wu, Xb)
%DESIGNGRIDLPVCONTROLLER  Grid-based LPV H-infinity output-feedback synthesis.
%
%   [K, CL, GAM, P] = DESIGNGRIDLPVCONTROLLER(GLPV, WE, WU, XB) synthesizes a
%   grid-based LPV controller for the PSS plant GLPV, using the mixed
%   sensitivity interconnection
%
%                +----------------------+
%        r ----->|                      |----> z1 = We*(y + n - r)
%        n ----->|      P(rho)          |----> z2 = Wu*u
%        u ----->|                      |----> y  = Glpv*u + n - r
%                +----------------------+
%
%   and solves the rate-bounded parameter-dependent LMIs via LPVSYN with the
%   Lyapunov basis XB.
%
%   Inputs:
%     GLPV - LPV plant as an LPVTools PSS object
%     WE   - performance weight (see WEIGHTPERFORMANCE)
%     WU   - control weight (see WEIGHTCONTROL)
%     XB   - basis function vector for the parameter-dependent Lyapunov
%            matrix, e.g. [basis(1,0); basis(rho,1)] for the affine form
%            X(rho) = X0 + X1*rho used in the reference paper
%
%   Outputs:
%     K   - synthesized LPV controller (PSS, rate-independent realization)
%     CL  - weighted closed loop
%     GAM - achieved gamma performance level
%     P   - the generalized plant that was built
%
%   This single function replaces the four near-identical design scripts of
%   the original codebase. The only thing that differed between them was the
%   choice of WE and WU, so those are now arguments:
%
%     nominal K0 over the full region:
%       We = weightPerformance(1e-3, 2, 1);  Wu = weightControl(1e-2, 1, 5);
%     local Ki over a subregion:
%       We = weightPerformance(1e-3, 2, 1);  Wu = weightControl(1e-2, 1, 10);
%
%   Requires LPVTools and the Robust Control Toolbox.
%
%   See also DESIGNPOLYTOPICLPVCONTROLLER, DESIGNLTICONTROLLER,
%            WEIGHTPERFORMANCE, WEIGHTCONTROL.

arguments
    Glpv
    We
    Wu
    Xb
end

% sysic resolves the block names below from this workspace, so the local
% variable names Wes, Wus and Glpv are load bearing - do not rename them.
Wes = ss(We);
Wus = ss(Wu);

systemnames  = 'Wes Wus Glpv'; %#ok<NASGU>
inputvar     = '[ r(1); n(1); u(1)]'; %#ok<NASGU>
outputvar    = '[ Wes; Wus; -r+Glpv+n]'; %#ok<NASGU>
input_to_Wes = '[-r+Glpv+n]'; %#ok<NASGU>
input_to_Wus = '[u]'; %#ok<NASGU>
input_to_Glpv = '[u]'; %#ok<NASGU>
sysoutname   = 'P'; %#ok<NASGU>
cleanupsysic = 'yes'; %#ok<NASGU>
sysic;

nmeas = size(Glpv.Data.c(:,:,1), 1);
ncon  = 1;

[Kfull, gam] = lpvsyn(P, nmeas, ncon, Xb, Xb);

% Freeze the rate variable and drop it from the domain, giving a controller
% that is scheduled on rho alone and therefore implementable.
Krate = lpvinterp(Kfull, 'rhoDot', 0);
K     = lpvelimiv(Krate);

CL = lft(P, K);

end
