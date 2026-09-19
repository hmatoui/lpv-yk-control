%% EXAMPLE 01 - Grid-based LPV-YK switching controller
%
% Reproduces the grid-based design of Section II / Section IV-B1 of:
%
%   H. Atoui, O. Sename, V. Milanes, J.-J. Martinez,
%   "Advanced Design of Multiple LPV Controllers with Youla-Kucera
%    Parametrization: Real-world Validation on Autonomous Vehicles".
%
% Pipeline:
%   1. Partition the speed range [5, 30] m/s into 5 subsets
%   2. Design a nominal LPV controller K0 over the full region
%   3. Design local LPV controllers Ki over each subset, INDEPENDENTLY
%   4. Solve LMIs (6)-(7) for the state-feedback gains Fg and Fk0
%   5. Build the doubly-coprime factorization -> the fixed block J(rho)
%   6. Build one Youla parameter Qi per subset
%   7. Verify each Ki is recovered EXACTLY by the YK architecture
%   8. Equalize orders and discretize for real-time implementation
%
% The key structural point: steps 2 and 3 are completely decoupled from step
% 4. Nothing in the LMIs mentions the number of subsets, the local
% controllers, or the switching signal - which is why a sixth controller
% could be added later without redoing any of this.
%
% Requirements: LPVTools, Robust Control Toolbox, YALMIP + SDP solver.

clear; clc; close all;
setupPath;

%% 1 - Configuration --------------------------------------------------
speedRange = [5 30];      % longitudinal speed range [m/s]
nSubsets   = 5;           % number of parameter subsets
rateBound  = 5;           % |d(vx)/dt| bound [m/s^2]
Ts         = 0.01;        % implementation sample time [s]
nGridPerSubset = 6;       % grid points per subset

subsets = partitionInterval(speedRange, nSubsets);

fprintf('Parameter partition (m/s):\n');
disp(subsets);

%% 2 - Build the gridded LPV plant ------------------------------------
veh = vehicleParameters();
act = steeringActuator();

% Build a global grid whose points include every subset boundary
gridPoints = [];
subGrids   = zeros(nSubsets, nGridPerSubset);
for p = 1:nSubsets
    subGrids(p,:) = linspace(subsets(p,1), subsets(p,2), nGridPerSubset);
    gridPoints    = [gridPoints, subGrids(p,1:end-1)]; %#ok<AGROW>
end
gridPoints = [gridPoints, subGrids(end,end)];

rho = pgrid('rho', gridPoints, [-rateBound rateBound]);

% Affine Lyapunov basis: X(rho) = X0 + X1*rho
Xb = [basis(1, 0); basis(rho, 1)];

[A, B, C, D] = bicycleModelGrid(veh, act, rho);
Glpv = pss(A, B, C, D);

fprintf('LPV plant: %d states, %d grid points.\n', ...
        size(Glpv.a,1), numel(gridPoints));

%% 3 - Nominal controller K0 over the full region ---------------------
% Detuned on purpose: K0 only has to be a globally valid safety net, local
% performance is recovered later through the Youla parameters.
We0 = weightPerformance(1e-3, 2, 1);
Wu0 = weightControl(1e-2, 1, 5);

fprintf('\nDesigning nominal controller K0 over [%g %g] m/s ...\n', speedRange);
[K0, CL0, gam0] = designGridLpvController(Glpv, We0, Wu0, Xb);
fprintf('  gamma_inf,0 = %.4f\n', gam0);

%% 4 - Local controllers Ki, designed independently -------------------
% More control authority (higher Wu corner) because each only covers 5 m/s.
We = weightPerformance(1e-3, 2, 1);
Wu = weightControl(1e-2, 1, 10);

Gsub = cell(1, nSubsets);
Ksub = cell(1, nSubsets);
gamLocal = zeros(1, nSubsets);

fprintf('\nDesigning local controllers (Table I):\n');
for p = 1:nSubsets
    Gsub{p} = lpvsplit(Glpv, 'rho', [subsets(p,1) subsets(p,2)]);
    [Ksub{p}, ~, gamLocal(p)] = designGridLpvController(Gsub{p}, We, Wu, Xb);
    fprintf('  P%d = [%2g %2g] m/s :  gamma_inf,%d = %.4f\n', ...
            p, subsets(p,1), subsets(p,2), p, gamLocal(p));
end

fprintf('\n  max_i gamma_inf,i = %.4f   <-- guaranteed global level\n', max(gamLocal));
fprintf('  (a single LPV controller over the full region with the same\n');
fprintf('   weights achieves gamma_inf = 1.42 in the reference paper)\n');

%% 5 - State-feedback gains Fg and Fk0 (LMIs 6-7) ---------------------
% These are the ONLY LMIs in the whole scheme, and neither involves the
% local controllers or the switching signal.
fprintf('\nSolving LMIs (6)-(7) for Fg and Fk0 ...\n');

GgenPlant = stateFeedbackPlant(Glpv, 1.8);
Fg = lpvsfsyn(GgenPlant, 1, 'L2');

GgenCtrl = stateFeedbackPlant(K0, 1.0);
Fk0 = lpvsfsyn(GgenCtrl, 1, 'L2');

fprintf('  Fg  : %d-by-%d parameter-varying gain\n', size(Fg,1), size(Fg,2));
fprintf('  Fk0 : %d-by-%d parameter-varying gain\n', size(Fk0,1), size(Fk0,2));

%% 6 - Doubly-coprime factorization -> the fixed block J(rho) ---------
cf = coprimeFactorizationLpv(Glpv, K0, Fg, Fk0);

% Verify the factorization identities at one grid point
k = 2;
fprintf('\nVerifying coprime factorization at rho = %.2f m/s:\n', gridPoints(k));
verifyStepEquality( ...
    {Glpv.Data(:,:,k), ...
     cf.N.Data(:,:,k)/cf.M.Data(:,:,k), ...
     cf.Mt.Data(:,:,k)\cf.Nt.Data(:,:,k)}, ...
    {'G', 'N*inv(M)', 'inv(Mt)*Nt'});
verifyStepEquality( ...
    {K0.Data(:,:,k), ...
     cf.U.Data(:,:,k)/cf.V.Data(:,:,k), ...
     cf.Vt.Data(:,:,k)\cf.Ut.Data(:,:,k)}, ...
    {'K0', 'U*inv(V)', 'inv(Vt)*Ut'});

%% 7 - Youla parameter Qi per subset ----------------------------------
fprintf('\nBuilding the Youla parameters Qi and verifying exact recovery:\n');

Q     = cell(1, nSubsets);
Qimpl = cell(1, nSubsets);

for p = 1:nSubsets
    rhoRange = [subsets(p,1) subsets(p,2)];

    Q{p} = youlaParameter( ...
        lpvsplit(Glpv, 'rho', rhoRange), ...
        lpvsplit(K0,   'rho', rhoRange), ...
        Ksub{p}, ...
        lpvsplit(Fg,   'rho', rhoRange), ...
        lpvsplit(Fk0,  'rho', rhoRange));

    % Exact-recovery check: Fl(J, Qi) must equal Ki, not approximate it
    cfSub.M = lpvsplit(cf.M, 'rho', rhoRange);
    cfSub.N = lpvsplit(cf.N, 'rho', rhoRange);
    cfSub.U = lpvsplit(cf.U, 'rho', rhoRange);
    cfSub.V = lpvsplit(cf.V, 'rho', rhoRange);
    Krec = ykController(cfSub, Q{p});

    % Compare the underlying ss arrays (one model per grid point) rather than
    % the PSS objects, so the check runs on plain Control System Toolbox code.
    fprintf('  P%d : ', p);
    verifyStepEquality({Ksub{p}.Data, Krec.Data}, {sprintf('K%d',p), 'Fl(J,Q)'});

    % The implemented parameter is Q*inv(V0) - see docs/theory.md
    Qimpl{p} = Q{p} * inv(lpvsplit(cf.V, 'rho', rhoRange)); %#ok<MINV>
end

%% 8 - Equalize orders and discretize ---------------------------------
% The switched bank Aq = diag(Aq_1, ..., Aq_N) needs a constant state
% dimension for a fixed-size real-time realization.
fprintf('\nPreparing the implementation realizations ...\n');

nq = cellfun(@(q) size(minreal(q.Data(:,:,1)).a, 1), Qimpl);
fprintf('  Youla parameter orders before padding: %s\n', mat2str(nq));

Qd = cell(1, nSubsets);
for p = 1:nSubsets
    Qr    = equalizeOrder(minreal(Qimpl{p}.Data), -1e4);
    Qd{p} = c2d(Qr, Ts, 'Tustin');
end

K0d = c2d(K0.Data, Ts, 'Tustin');
Md  = c2d(equalizeOrder(minreal(cf.M.Data)),  Ts, 'Tustin');
Nd  = c2d(equalizeOrder(minreal(cf.N.Data)),  Ts, 'Tustin');
Ud  = c2d(equalizeOrder(minreal(cf.U.Data)),  Ts, 'Tustin');
Vid = c2d(equalizeOrder(minreal(inv(cf.V.Data))), Ts, 'Tustin');

%% 9 - Package the design ---------------------------------------------
design = struct();
design.method      = 'grid-based LPV-YK';
design.subsets     = subsets;
design.gridPoints  = gridPoints;
design.Ts          = Ts;
design.gammaNominal = gam0;
design.gammaLocal   = gamLocal;
design.gammaGlobal  = max(gamLocal);
design.K0d         = K0d;
design.Md          = Md;
design.Nd          = Nd;
design.Ud          = Ud;
design.Vinvd       = Vid;
design.Qd          = Qd;

outFile = fullfile(fileparts(mfilename('fullpath')), 'gridLpvYkDesign.mat');
save(outFile, '-struct', 'design');

fprintf('\nDesign complete. Saved to:\n  %s\n', outFile);
fprintf('\nSummary\n');
fprintf('  nominal gamma           : %.4f\n', gam0);
fprintf('  worst local gamma       : %.4f\n', max(gamLocal));
fprintf('  switching signal        : sigma_i in {0,1}, toggled at boundaries\n');
fprintf('  stability guarantee     : exponential, for ANY bounded sigma\n');
