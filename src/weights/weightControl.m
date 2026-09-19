function Wu = weightControl(eps, Mu, wu)
%WEIGHTCONTROL  Control-effort weighting filter Wu(s).
%
%   WU = WEIGHTCONTROL(EPS, MU, WU_) returns the first-order control weight
%   that limits actuator activity and enforces roll-off:
%
%                 s + wu/Mu
%       Wu(s) = --------------
%                 eps*s + wu
%
%   Interpretation of the three knobs:
%     EPS - sets the high-frequency gain 1/eps, i.e. how strongly control
%           action is penalized beyond the actuator bandwidth.
%     MU  - low-frequency control authority bound.
%     WU_ - actuator roll-off corner frequency [rad/s].
%
%   Typical values from the reference paper (Sec. IV-B). Note the local
%   controllers are given a *higher* corner frequency, i.e. more authority,
%   because each only has to cover a narrow speed band:
%     nominal controller K0 : weightControl(1e-2, 1, 5)
%     local controllers  Ki : weightControl(1e-2, 1, 10)
%
%   See also WEIGHTPERFORMANCE.

arguments
    eps (1,1) double {mustBePositive}
    Mu  (1,1) double {mustBePositive}
    wu  (1,1) double {mustBePositive}
end

s  = tf('s');
Wu = (s + wu/Mu)/(eps*s + wu);

end
