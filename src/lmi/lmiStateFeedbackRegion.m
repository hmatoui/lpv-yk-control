function [F, P] = lmiStateFeedbackRegion(A, B, alpha, theta, r, solver, showPlot)
%LMISTATEFEEDBACKREGION  Polytopic state feedback with LMI pole-region constraints.
%
%   [F, P] = LMISTATEFEEDBACKREGION(A, B, ALPHA, THETA, R, SOLVER) computes
%   vertex state-feedback gains that place the closed-loop eigenvalues of
%   every vertex inside an LMI region defined by the intersection of:
%
%     * a half plane   Re(s) <= -ALPHA      (minimum decay rate)
%     * a disc         |s|   <=  R          (maximum natural frequency)
%     * a conic sector of half-angle THETA  (minimum damping)
%
%   A is a cell array of vertex state matrices sharing the common input
%   matrix B. A single constant Lyapunov matrix P is used for all vertices,
%   so the result is a quadratically stabilizing polytopic gain.
%
%   [F, P] = LMISTATEFEEDBACKREGION(..., SHOWPLOT) draws the region and the
%   resulting closed-loop eigenvalues when SHOWPLOT is true (default false).
%
%   Sign convention: F is returned for POSITIVE feedback, i.e. the stabilized
%   dynamics are A{i} + B*F{i}. This matches the "+ B2*Fg" convention used
%   throughout the Youla-Kucera construction, so the gains can be passed
%   directly to COPRIMEFACTORIZATION.
%
%   Why a region constraint rather than plain stabilization: the choice of Fg
%   shapes the dynamics of the Youla parameter Q, and hence the closed-loop
%   transient during switching. Confining the factorization poles to a
%   well-conditioned region keeps Q numerically well behaved and avoids
%   very fast modes that would stiffen the real-time discretization.
%
%   Requires YALMIP and an SDP solver.

arguments
    A        (1,:) cell
    B        double
    alpha    (1,1) double
    theta    (1,1) double
    r        (1,1) double
    solver   (1,:) char
    showPlot (1,1) logical = false
end

epsi = 1e-10;
[nx, nu] = size(B);
nV = numel(A);

P = sdpvar(nx, nx);
Y   = cell(1,nV);
AtP = cell(1,nV);
PA  = cell(1,nV);
for i = 1:nV
    Y{i}   = sdpvar(nu, nx, 'full');
    AtP{i} = P*A{i}' - Y{i}'*B';
    PA{i}  = A{i}*P - B*Y{i};
end

LMIs = (P >= epsi);
for i = 1:nV
    % Decay rate
    LMIs = [LMIs, P*A{i}' + A{i}*P - Y{i}'*B' - B*Y{i} + 2*alpha*P <= -epsi]; %#ok<AGROW>
    % Disc of radius r
    LMIs = [LMIs, [-r*P, PA{i}; AtP{i}, -r*P] <= -epsi]; %#ok<AGROW>
    % Conic sector of half-angle theta
    LMIs = [LMIs, [sin(theta)*(PA{i} + AtP{i}), cos(theta)*(PA{i} - AtP{i}); ...
                   cos(theta)*(AtP{i} - PA{i}), sin(theta)*(PA{i} + AtP{i})] <= -epsi]; %#ok<AGROW>
end

ops = sdpsettings('solver', solver, 'verbose', 0);
optimize(LMIs, [], ops);

Pv = value(P);
if any(eig(Pv) <= 0)
    warning('lmiStateFeedbackRegion:nonPositiveLyapunov', ...
        'P is not positive definite; the solution is unreliable.');
elseif any(eig(Pv) <= 1e-8)
    warning('lmiStateFeedbackRegion:illConditioned', ...
        'P has very small eigenvalues; the solution may be ill conditioned.');
end

F  = cell(1,nV);
VP = cell(1,nV);
for i = 1:nV
    F{i}  = -(value(Y{i})/Pv);       % positive-feedback convention
    VP{i} = eig(A{i} + B*F{i});
end

P = Pv;

if showPlot
    maxAbs = r + 2;
    figure; hold on; grid on;
    plot([-maxAbs, 2], [0 0], 'k--');
    plot([0 0], [-maxAbs, maxAbs], 'k--');
    plot([-alpha -alpha], [-maxAbs maxAbs], 'LineWidth', 1);
    plot([-maxAbs 0], [-maxAbs 0]*tan(theta),  'LineWidth', 1);
    plot([-maxAbs 0], [-maxAbs 0]*tan(-theta), 'LineWidth', 1);
    plot(r*cos(0:0.01:2*pi), r*sin(0:0.01:2*pi), 'LineWidth', 1);
    for i = 1:nV
        plot(real(VP{i}), imag(VP{i}), '+', 'LineWidth', 2);
    end
    xlim([-maxAbs, 2]); ylim([-maxAbs, maxAbs]);
    xlabel('Re'); ylabel('Im'); title('Closed-loop eigenvalues and LMI region');
end

end
