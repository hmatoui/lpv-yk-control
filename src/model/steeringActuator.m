function a = steeringActuator()
%STEERINGACTUATOR  Identified steering actuator dynamics of the test vehicle.
%
%   A = STEERINGACTUATOR() returns a struct describing the identified
%   steering column dynamics: a second-order low-pass in series with a pure
%   input delay approximated by a 2nd-order Pade expansion.
%
%       Sigma_act(s) = Kp / (1 + 2*zeta*Tw*s + (Tw*s)^2) * exp(-Td*s)
%
%   Fields:
%     Kp        - static gain [-]
%     Tw        - time constant [s]
%     zeta      - damping ratio [-]
%     Td        - transport delay [s]
%     numDelay  - numerator of the 2nd-order Pade approximation
%     denDelay  - denominator of the 2nd-order Pade approximation
%     tf        - the full continuous-time transfer function
%
%   The Pade coefficients are returned explicitly because the bicycle model
%   builders embed the delay states directly in the state-space realization
%   rather than composing transfer functions.
%
%   See also BICYCLEMODELGRID, BICYCLEMODELPOLYTOPIC, PADE.

a.Kp   = 1.0012;
a.Tw   = 0.09842;
a.zeta = 0.75013;
a.Td   = 0.060756;

[a.numDelay, a.denDelay] = pade(a.Td, 2);

s = tf('s');
a.tf = (a.Kp/(1 + 2*a.zeta*a.Tw*s + (a.Tw*s)^2)) * tf(a.numDelay, a.denDelay);

end
