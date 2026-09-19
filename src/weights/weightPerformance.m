function We = weightPerformance(eps, Ms, wb)
%WEIGHTPERFORMANCE  Tracking-performance weighting filter We(s).
%
%   WE = WEIGHTPERFORMANCE(EPS, MS, WB) returns the first-order performance
%   weight used to shape the sensitivity function in the H-infinity designs:
%
%                 wb + s/Ms
%       We(s) = --------------
%                 s + wb*eps
%
%   Interpretation of the three knobs:
%     EPS - low-frequency sensitivity bound; 1/eps is the DC gain of We, so
%           small eps enforces tight steady-state tracking.
%     MS  - high-frequency bound on |We|, i.e. the maximum sensitivity peak
%           allowed (typically 2).
%     WB  - desired closed-loop bandwidth [rad/s].
%
%   Typical values from the reference paper (Sec. IV-B):
%     nominal controller K0 : weightPerformance(1e-3, 2, 1)
%     local controllers  Ki : weightPerformance(1e-3, 2, 1)
%
%   See also WEIGHTCONTROL.

arguments
    eps (1,1) double {mustBePositive}
    Ms  (1,1) double {mustBePositive}
    wb  (1,1) double {mustBePositive}
end

s  = tf('s');
We = (wb + s/Ms)/(s + wb*eps);

end
