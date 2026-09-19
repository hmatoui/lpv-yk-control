function G = polytopicInterpolation(rho, box, Gv)
%POLYTOPICINTERPOLATION  Convex interpolation over the reduced 3-vertex simplex.
%
%   G = POLYTOPICINTERPOLATION(RHO, BOX, GV) evaluates the polytopic
%   combination
%
%       G(rho) = sum_j alpha_j(rho) * Gv_j,   sum_j alpha_j = 1, alpha_j >= 0
%
%   at the operating point RHO = [rho1 rho2] inside the triangular parameter
%   region whose bounding box is BOX = [rho1Min rho1Max; rho2Min rho2Max].
%
%   GV holds the three vertex objects, ordered as produced by
%   POLYTOPEVERTICES with REDUCED = true:
%       1: (rho1Min, rho2Min)   2: (rho1Min, rho2Max)   3: (rho1Max, rho2Min)
%
%   GV may be either an LTI array (indexed GV(:,:,j)) or any object array for
%   which scalar multiplication and addition are defined. Both paths are
%   supported because the calling code interpolates coprime factors, state
%   feedback gains and closed loops alike.
%
%   The barycentric coordinates follow directly from the simplex geometry:
%       alpha_2 = |rho2 - rho2Min| / (rho2Max - rho2Min)
%       alpha_3 = |rho1 - rho1Min| / (rho1Max - rho1Min)
%       alpha_1 = 1 - alpha_2 - alpha_3
%
%   See also POLYTOPEVERTICES.

arguments
    rho (1,2) double
    box (2,2) double
    Gv
end

a2 = abs(rho(2) - box(2,1)) / (box(2,2) - box(2,1));
a3 = abs(rho(1) - box(1,1)) / (box(1,2) - box(1,1));
a1 = 1 - (a2 + a3);

if isnumeric(Gv) || size(Gv,1) >= 3 || size(Gv,2) >= 3
    G = a1*Gv(:,:,1) + a2*Gv(:,:,2) + a3*Gv(:,:,3);
else
    A = a1*Gv.a(:,:,1) + a2*Gv.a(:,:,2) + a3*Gv.a(:,:,3);
    B = a1*Gv.b(:,:,1) + a2*Gv.b(:,:,2) + a3*Gv.b(:,:,3);
    C = a1*Gv.c(:,:,1) + a2*Gv.c(:,:,2) + a3*Gv.c(:,:,3);
    D = a1*Gv.d(:,:,1) + a2*Gv.d(:,:,2) + a3*Gv.d(:,:,3);
    G = ss(A, B, C, D);
end

end
