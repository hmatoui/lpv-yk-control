function subsets = partitionInterval(range, nSubsets)
%PARTITIONINTERVAL  Split a scalar parameter range into contiguous subsets.
%
%   SUBSETS = PARTITIONINTERVAL(RANGE, NSUBSETS) divides the interval
%   RANGE = [rhoMin rhoMax] into NSUBSETS equal closed subintervals that
%   share their endpoints. The result is an NSUBSETS-by-2 matrix whose i-th
%   row is [lower_i upper_i], with upper_i == lower_{i+1}.
%
%   Adjacent subsets deliberately share a boundary point: that shared point
%   is the switching surface on which the LPV-YK switching signal sigma
%   toggles. Overlapping subsets are not required by the theory.
%
%   Example - the partition used in the reference paper, eq. (43):
%       subsets = partitionInterval([5 30], 5)
%       %  5   10
%       % 10   15
%       % 15   20
%       % 20   25
%       % 25   30
%
%   See also POLYTOPEVERTICES.

arguments
    range    (1,2) double
    nSubsets (1,1) double {mustBeInteger, mustBePositive}
end

if range(2) <= range(1)
    error('partitionInterval:badRange', ...
          'RANGE must be increasing, got [%g %g].', range(1), range(2));
end

edges   = linspace(range(1), range(2), nSubsets+1);
subsets = [edges(1:end-1).', edges(2:end).'];

end
