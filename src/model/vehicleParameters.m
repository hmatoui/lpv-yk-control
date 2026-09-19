function p = vehicleParameters()
%VEHICLEPARAMETERS  Physical parameters of the Renault ZOE test vehicle.
%
%   P = VEHICLEPARAMETERS() returns a struct holding the lateral-dynamics
%   parameters of the robotized electric Renault ZOE used in the Satory
%   experimental campaign.
%
%   Fields:
%     lf, lr  - distance from CoG to front/rear axle [m]
%     L       - wheelbase, lf + lr [m]
%     m       - vehicle mass [kg]
%     Iz      - yaw inertia [kg.m^2]
%     Cf, Cr  - front/rear cornering stiffness [N/rad]
%     mu      - road friction coefficient [-]
%     g       - gravitational acceleration [m/s^2]
%     Kus     - understeer gradient, derived [s^2/m]
%     vChar   - characteristic speed, derived [m/s]
%
%   See also BICYCLEMODELGRID, BICYCLEMODELPOLYTOPIC.

p.lf = 1.139;
p.lr = 1.449;
p.L  = p.lf + p.lr;

p.m  = 1664;
p.Ix = 380;
p.Iy = 1529;
p.Iz = 2198;

p.mu = 1;
p.Cf = 75000/2;
p.Cr = 120000/2;
p.g  = 9.81;

% Derived quantities
p.Kus   = (p.m/p.L) * (p.lr/(2*p.Cf) - p.lf/(2*p.Cr));
p.vChar = sqrt(p.L/p.Kus);

end
