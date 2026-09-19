function [A, B, C, D] = bicycleModelPolytopic(p, act, rho)
%BICYCLEMODELPOLYTOPIC  Lateral bicycle model in affine (polytopic) LPV form.
%
%   [A,B,C,D] = BICYCLEMODELPOLYTOPIC(P, ACT, RHO) builds the state-space
%   data of the lateral vehicle dynamics augmented with the steering
%   actuator, using the *affine* parameter dependence of eq. (42) in the
%   reference paper, with
%
%       rho = [rho1; rho2],   rho1 = vx,   rho2 = 1/vx
%
%   Treating vx and 1/vx as two independent coordinates is exactly what makes
%   the model affine - and exactly what introduces over-bounding, since the
%   convex hull of the two coordinates contains points off the physical
%   hyperbola rho2 = 1/rho1. Partitioning the hull along that hyperbola is
%   the mechanism by which the partitioned polytopic LPV-YK scheme recovers
%   the lost performance.
%
%   The state, input and output definitions match BICYCLEMODELGRID.
%
%   Inputs:
%     P   - vehicle parameter struct from VEHICLEPARAMETERS
%     ACT - actuator struct from STEERINGACTUATOR
%     RHO - 2-element vector [vx; 1/vx]
%
%   See also BICYCLEMODELGRID, POLYTOPEVERTICES, VEHICLEPARAMETERS.

arguments
    p   (1,1) struct
    act (1,1) struct
    rho (2,1) double
end

m  = p.m;   I  = p.Iz;
Cf = p.Cf;  Cr = p.Cr;
lf = p.lf;  lr = p.lr;

Kp = act.Kp;  Tw = act.Tw;  zeta = act.zeta;
num = act.numDelay;  den = act.denDelay;

rho1 = rho(1);   % vx
rho2 = rho(2);   % 1/vx

% Affine in (rho1, rho2)
A22 = -rho2*(Cr + Cf)/m;
A23 =  rho2*(Cr*lr - Cf*lf)/m - rho1;
A32 =  rho2*(Cr*lr - Cf*lf)/I;
A33 =  rho2*(-Cr*lr^2 - Cf*lf^2)/I;

B21 = Cf/m;
B31 = Cf*lf/I;

B51 = -den(2)*num(1)/den(1) + num(2);
B61 =  num(3) - den(3)*num(1)/den(1);

A = [A22  A23  B21        0                 0            0;
     A32  A33  B31        0                 0            0;
     0    0    0          1                 0            0;
     0    0   -1/(Tw^2)  -2*zeta/Tw   Kp/(Tw^2)/den(1)   0;
     0    0    0          0          -den(2)/den(1)      1;
     0    0    0          0          -den(3)/den(1)      0];

B = [0; 0; 0; Kp/(Tw^2); B51; B61];
C = [0 1 0 0 0 0];
D = zeros(size(C,1), size(B,2));

end
