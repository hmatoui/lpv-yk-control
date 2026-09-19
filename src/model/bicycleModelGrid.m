function [A, B, C, D] = bicycleModelGrid(p, act, rho)
%BICYCLEMODELGRID  Lateral bicycle model in gridded (non-affine) LPV form.
%
%   [A,B,C,D] = BICYCLEMODELGRID(P, ACT, RHO) builds the state-space data of
%   the lateral vehicle dynamics augmented with the steering actuator, using
%   the *non-affine* parameter dependence of eq. (41) in the reference paper.
%   The scheduling parameter is the longitudinal speed RHO = vx, which may be
%   an LPVTools PGRID object (yielding parameter-varying matrices) or a
%   numeric scalar.
%
%   Because the grid-based synthesis evaluates the pLMIs pointwise, the
%   1/vx terms need no affine rewriting - this is the structural advantage of
%   the gridded formulation over the polytopic one.
%
%   States:  x = [vy; w; delta; delta_dot; xd1; xd2]
%     vy         - lateral velocity in the body frame [m/s]
%     w          - yaw rate [rad/s]
%     delta      - front steering angle [rad]
%     delta_dot  - steering rate [rad/s]
%     xd1, xd2   - Pade delay states
%   Input:   u = steering angle command [rad]
%   Output:  y = yaw rate w [rad/s]
%
%   Inputs:
%     P   - vehicle parameter struct from VEHICLEPARAMETERS
%     ACT - actuator struct from STEERINGACTUATOR
%     RHO - longitudinal speed vx (PGRID or scalar) [m/s]
%
%   See also BICYCLEMODELPOLYTOPIC, VEHICLEPARAMETERS, STEERINGACTUATOR.

m  = p.m;   I  = p.Iz;
Cf = p.Cf;  Cr = p.Cr;
lf = p.lf;  lr = p.lr;

Kp = act.Kp;  Tw = act.Tw;  zeta = act.zeta;
num = act.numDelay;  den = act.denDelay;

vx = rho;

% Lateral / yaw dynamics - note the explicit 1/vx dependence
A22 = -(Cr + Cf)/(m*vx);
A23 = (Cr*lr - Cf*lf)/(m*vx) - vx;
A32 = (Cr*lr - Cf*lf)/(I*vx);
A33 = (-Cr*lr^2 - Cf*lf^2)/(I*vx);

B21 = Cf/m;
B31 = Cf*lf/I;

% Pade delay realization
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
