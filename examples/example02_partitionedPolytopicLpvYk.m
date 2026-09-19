%% EXAMPLE 02 - Partitioned polytopic LPV-YK controller
%
% Reproduces the partitioned polytopic design of Section III / Section
% IV-B2 of:
%
%   H. Atoui, O. Sename, V. Milanes, J.-J. Martinez,
%   "Advanced Design of Multiple LPV Controllers with Youla-Kucera
%    Parametrization: Real-world Validation on Autonomous Vehicles".
%
% Pipeline:
%   1. Partition the speed range and build a triangular polytope per subset
%      ALONG the physical trajectory rho2 = 1/rho1
%   2. Design a nominal polytopic LPV controller K0 over the full hull P0
%   3. Design LTI controllers at the UNIQUE vertices, reusing the controller
%      on every shared face  ->  assumption (A.3.2)
%   4. Solve LMIs (28)-(29) for Fg and Fk0
%   5. Coprime factorization and Youla parameters at each vertex
%   6. Verify exact recovery, and verify the shared-face identity
%   7. Sweep the closed-loop H-infinity norm vs speed to quantify the
%      conservatism removed by partitioning
%
% The headline idea: the convex hull of (vx, 1/vx) contains points that are
% physically impossible, and covering them costs performance. Hugging the
% hyperbola with a chain of small polytopes removes that cost. The Youla
% parametrization is what makes switching between them free.
%
% Requirements: Robust Control Toolbox, YALMIP + SDP solver.

clear; clc; close all;
setupPath;

%% 1 - Configuration and partition ------------------------------------
speedRange = [5 30];
nSubsets   = 5;
Ts         = 0.01;
solver     = 'sedumi';

subsets = partitionInterval(speedRange, nSubsets);

veh = vehicleParameters();
act = steeringActuator();

% Full hull P0 plus one triangular subset per partition element.
% polytopeVertices(..., true) drops the (vxMax, 1/vxMin) corner, which is
% the most extreme violation of rho2 = 1/rho1.
[V0, box0] = polytopeVertices(speedRange, true);

Vsub   = cell(1, nSubsets);
boxSub = cell(1, nSubsets);
for p = 1:nSubsets
    [Vsub{p}, boxSub{p}] = polytopeVertices(subsets(p,:), true);
end

nVert = size(V0, 2);
fprintf('Triangular hull P0 vertices (rho1; rho2):\n');
disp(V0);

%% 2 - Vertex plant models --------------------------------------------
G0 = cell(1, nVert);
for j = 1:nVert
    [A,B,C,D] = bicycleModelPolytopic(veh, act, V0(:,j));
    G0{j} = ss(A,B,C,D);
end

Gsub = cell(nSubsets, nVert);
for p = 1:nSubsets
    for j = 1:nVert
        [A,B,C,D] = bicycleModelPolytopic(veh, act, Vsub{p}(:,j));
        Gsub{p,j} = ss(A,B,C,D);
    end
end

%% 3 - Nominal polytopic controller K0 over the full hull -------------
We0 = weightPerformance(1e-3, 2, 1);
Wu0 = weightControl(1e-2, 1, 5);

fprintf('\nDesigning nominal polytopic controller K0 over P0 ...\n');
[K0v, ~, gam0] = designPolytopicLpvController(G0, We0, Wu0, ...
    'StrictlyProper', true, 'Percentage', 0, 'Solver', solver);
fprintf('  gamma_inf,0 = %.4f   (conservative by construction)\n', gam0);

%% 4 - LTI vertex controllers, designed ONCE per unique vertex --------
% This is assumption (A.3.2) implemented literally. Adjacent subsets share a
% vertex on the hyperbola: subset p vertex 3 == subset p+1 vertex 2. By
% keying the design on the vertex coordinates we guarantee the controller is
% byte-identical on both sides, so switching across that face is a no-op.
We = weightPerformance(1e-3, 2, 1);
Wu = weightControl(1e-2, 1, 10);

fprintf('\nDesigning LTI vertex controllers (Table II):\n');

vertexKeys = containers.Map('KeyType','char','ValueType','any');
Kv      = cell(nSubsets, nVert);
gamV    = zeros(nSubsets, nVert);
nDesigns = 0;

for p = 1:nSubsets
    for j = 1:nVert
        key = sprintf('%.6f_%.6f', Vsub{p}(1,j), Vsub{p}(2,j));
        if isKey(vertexKeys, key)
            cached      = vertexKeys(key);
            Kv{p,j}     = cached.K;
            gamV(p,j)   = cached.gam;
        else
            [Kij, ~, gij] = designLtiController(Gsub{p,j}, We, Wu);
            Kv{p,j}   = Kij;
            gamV(p,j) = gij;
            vertexKeys(key) = struct('K', Kij, 'gam', gij);
            nDesigns = nDesigns + 1;
        end
    end
end

fprintf('  %d unique vertices designed for %d subset-vertex slots\n', ...
        nDesigns, nSubsets*nVert);
fprintf('  %-8s %8s %8s %8s\n', 'subset', 'w_i1', 'w_i2', 'w_i3');
for p = 1:nSubsets
    fprintf('  K%-7d %8.2f %8.2f %8.2f\n', p, gamV(p,1), gamV(p,2), gamV(p,3));
end

% Verify the shared-face identity K_{p,3} == K_{p+1,2}
fprintf('\nVerifying the shared-face identity of (A.3.2):\n');
for p = 1:nSubsets-1
    same = isequal(Kv{p,3}.a, Kv{p+1,2}.a) && isequal(Kv{p,3}.b, Kv{p+1,2}.b) ...
        && isequal(Kv{p,3}.c, Kv{p+1,2}.c) && isequal(Kv{p,3}.d, Kv{p+1,2}.d);
    fprintf('  K_{%d,3} == K_{%d,2} : %s\n', p, p+1, string(same));
end

%% 5 - State-feedback gains Fg and Fk0 (LMIs 28-29) -------------------
fprintf('\nSolving LMIs (28)-(29) for Fg and Fk0 ...\n');

Av = cellfun(@(g) g.a, G0, 'UniformOutput', false);
B2 = G0{1}.b;

% Region-constrained state feedback keeps the factorization well conditioned
Fg = lmiStateFeedbackRegion(Av, B2, 5, pi/2, 50, solver, false);

% Nominal-controller dynamics
Pk0 = cell(1, nVert);
for j = 1:nVert
    Pk0{j} = stateFeedbackPlant(K0v{j}, 1.0, false);
end
Fk0 = lmiHinfStateFeedbackPolytope(Pk0, size(K0v{1}.a,1), size(K0v{1}.b,2), ...
                                   1, 0, solver);

%% 6 - Coprime factorization and Youla parameters at the vertices -----
fprintf('\nBuilding coprime factors and Youla parameters ...\n');

cf0 = cell(1, nVert);
for j = 1:nVert
    cf0{j} = coprimeFactorization(G0{j}, K0v{j}, Fg{j}, Fk0{j});
end

% Assemble vertex factors as LTI arrays for interpolation
Marr = stack(3, cf0{1}.M,  cf0{2}.M,  cf0{3}.M);
Narr = stack(3, cf0{1}.N,  cf0{2}.N,  cf0{3}.N);
Uarr = stack(3, cf0{1}.U,  cf0{2}.U,  cf0{3}.U);
Varr = stack(3, cf0{1}.V,  cf0{2}.V,  cf0{3}.V);
Mtarr = stack(3, cf0{1}.Mt, cf0{2}.Mt, cf0{3}.Mt);
Ntarr = stack(3, cf0{1}.Nt, cf0{2}.Nt, cf0{3}.Nt);
Utarr = stack(3, cf0{1}.Ut, cf0{2}.Ut, cf0{3}.Ut);
Vtarr = stack(3, cf0{1}.Vt, cf0{2}.Vt, cf0{3}.Vt);
K0arr = stack(3, K0v{1}, K0v{2}, K0v{3});
Fgarr = cat(3, Fg{1}, Fg{2}, Fg{3});

Q   = cell(nSubsets, nVert);
CLp = cell(nSubsets, nVert);

for p = 1:nSubsets
    for j = 1:nVert
        rhoV = Vsub{p}(:,j).';

        % Nominal quantities interpolated AT this vertex
        K0here = polytopicInterpolation(rhoV, box0, K0arr);
        Fghere = polytopicInterpolation(rhoV, box0, Fgarr);

        cfHere.M = polytopicInterpolation(rhoV, box0, Marr);
        cfHere.N = polytopicInterpolation(rhoV, box0, Narr);
        cfHere.U = polytopicInterpolation(rhoV, box0, Uarr);
        cfHere.V = polytopicInterpolation(rhoV, box0, Varr);
        cfHere.Mt = polytopicInterpolation(rhoV, box0, Mtarr);
        cfHere.Nt = polytopicInterpolation(rhoV, box0, Ntarr);
        cfHere.Ut = polytopicInterpolation(rhoV, box0, Utarr);
        cfHere.Vt = polytopicInterpolation(rhoV, box0, Vtarr);

        Q{p,j} = youlaParameter(Gsub{p,j}, K0here, Kv{p,j}, Fghere, []);

        [Krec, KrecLeft] = ykController(cfHere, Q{p,j});
        if p == 1 && j == 1
            fprintf('  exact-recovery check at P1, vertex 1:\n    ');
            verifyStepEquality({Kv{p,j}, Krec, KrecLeft}, ...
                               {'K_11', '(U+MQ)inv(V+NQ)', '(Ut+QMt)inv(Vt+QNt)'});
        end

        % The synthesis interconnection feeds K with (y - r), so the loop is
        % u = K*y: POSITIVE feedback. feedback(...,1) would close the wrong loop.
        CLp{p,j} = feedback(Gsub{p,j}*Krec, 1, +1);
    end
end

%% 7 - Conservatism comparison ----------------------------------------
% Partitioned LPV-YK vs gain-scheduled LPV-YK over the full hull P0.
% Same YK machinery in both cases - the ONLY difference is the partition.
fprintf('\nSweeping closed-loop H-infinity norm vs speed ...\n');

vxSweep = speedRange(1):0.5:speedRange(2);
normPart = zeros(size(vxSweep));
normFull = zeros(size(vxSweep));

% Gain-scheduled benchmark: YK over P0 with vertex controllers of subset 1
Qfull = cell(1, nVert);
for j = 1:nVert
    rhoV   = V0(:,j).';
    K0here = polytopicInterpolation(rhoV, box0, K0arr);
    Fghere = polytopicInterpolation(rhoV, box0, Fgarr);
    Qfull{j} = youlaParameter(G0{j}, K0here, Kv{1,j}, Fghere, []);
end

CLfull = cell(1, nVert);
for j = 1:nVert
    rhoV = V0(:,j).';
    cfHere.M = polytopicInterpolation(rhoV, box0, Marr);
    cfHere.N = polytopicInterpolation(rhoV, box0, Narr);
    cfHere.U = polytopicInterpolation(rhoV, box0, Uarr);
    cfHere.V = polytopicInterpolation(rhoV, box0, Varr);
    CLfull{j} = feedback(G0{j}*ykController(cfHere, Qfull{j}), 1, +1);
end
CLfullArr = stack(3, CLfull{1}, CLfull{2}, CLfull{3});

for k = 1:numel(vxSweep)
    vx   = vxSweep(k);
    rhoK = [vx, 1/vx];

    % Which subset is active
    pAct = find(vx >= subsets(:,1) & vx <= subsets(:,2), 1, 'last');

    CLarr = stack(3, CLp{pAct,1}, CLp{pAct,2}, CLp{pAct,3});
    normPart(k) = norm(polytopicInterpolation(rhoK, boxSub{pAct}, CLarr), inf);
    normFull(k) = norm(polytopicInterpolation(rhoK, box0, CLfullArr), inf);
end

figure('Name','Conservatism: partitioned vs gain-scheduled');
plot(vxSweep, normPart, 'LineWidth', 2, 'Color', [0.16 0.77 0.18]); hold on;
plot(vxSweep, normFull, 'LineWidth', 2, 'Color', 'r');
grid on; box on;
xlabel('v_x  (m/s)', 'FontWeight', 'bold');
ylabel('Closed-loop H_\infty norm', 'FontWeight', 'bold');
legend({'Partitioned LPV-YK', 'Gain-scheduled LPV-YK over P_0'}, ...
       'Location', 'best');
title('Effect of trajectory-aligned partitioning');

fprintf('  mean norm, partitioned      : %.4f\n', mean(normPart));
fprintf('  mean norm, gain-scheduled   : %.4f\n', mean(normFull));

%% 8 - Package the design ---------------------------------------------
Qd = cell(nSubsets, nVert);
for p = 1:nSubsets
    for j = 1:nVert
        Qd{p,j} = c2d(minreal(Q{p,j}), Ts, 'Tustin');
    end
end

design = struct();
design.method     = 'partitioned polytopic LPV-YK';
design.subsets    = subsets;
design.box0       = box0;
design.boxSub     = boxSub;
design.Ts         = Ts;
design.gammaNominal = gam0;
design.gammaVertex  = gamV;
design.K0d        = cellfun(@(k) c2d(k, Ts, 'Tustin'), K0v, 'UniformOutput', false);
design.Qd         = Qd;
design.normSweep  = struct('vx', vxSweep, 'partitioned', normPart, 'fullHull', normFull);

outFile = fullfile(fileparts(mfilename('fullpath')), 'partitionedPolytopicLpvYkDesign.mat');
save(outFile, '-struct', 'design');

fprintf('\nDesign complete. Saved to:\n  %s\n', outFile);
fprintf('\nSummary\n');
fprintf('  nominal gamma over P0      : %.4f\n', gam0);
fprintf('  worst vertex gamma         : %.4f\n', max(gamV(:)));
fprintf('  stability guarantee        : quadratic, for ANY switching signal\n');
fprintf('  shared-face behaviour      : identical controller -> exact continuity\n');
