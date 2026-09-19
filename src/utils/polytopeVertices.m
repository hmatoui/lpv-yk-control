function [V, box] = polytopeVertices(speedRange, reduced)
%POLYTOPEVERTICES  Vertices of the (rho1, rho2) = (vx, 1/vx) polytope.
%
%   [V, BOX] = POLYTOPEVERTICES(SPEEDRANGE) returns the vertices of the
%   rectangular convex hull spanned by rho1 = vx in SPEEDRANGE and
%   rho2 = 1/vx. V is a 2-by-NV matrix whose columns are the vertices and
%   BOX is the 2-by-2 matrix [rho1Min rho1Max; rho2Min rho2Max].
%
%   [V, BOX] = POLYTOPEVERTICES(SPEEDRANGE, TRUE) drops the vertex
%   (rho1Max, rho2Max) = (vxMax, 1/vxMin), producing the *triangular* hull
%   used in the reference paper. That corner is the most extreme violation of
%   the physical constraint rho2 = 1/rho1 - it asserts simultaneously the
%   highest speed and the lowest speed - so removing it eliminates a large
%   part of the over-bounding at zero cost in coverage.
%
%   Vertex ordering (columns of V):
%     1: (rho1Min, rho2Min)
%     2: (rho1Min, rho2Max)
%     3: (rho1Max, rho2Min)
%     4: (rho1Max, rho2Max)   -- omitted when REDUCED is true
%
%   This ordering must stay consistent with POLYTOPICINTERPOLATION, which
%   assumes the reduced 3-vertex simplex.
%
%   See also POLYTOPICINTERPOLATION, BICYCLEMODELPOLYTOPIC.

arguments
    speedRange (1,2) double {mustBePositive}
    reduced    (1,1) logical = true
end

rho1 = [min(speedRange), max(speedRange)];
rho2 = [1/max(speedRange), 1/min(speedRange)];
box  = [rho1; rho2];

V = [rho1(1) rho1(1) rho1(2) rho1(2);
     rho2(1) rho2(2) rho2(1) rho2(2)];

if reduced
    V = V(:, 1:3);
end

end
