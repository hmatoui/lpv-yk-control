function tests = testUtilities()
%TESTUTILITIES  Unit tests for the partitioning and interpolation helpers.
%
%   Run with:  runtests('tests')
%
%   See also PARTITIONINTERVAL, POLYTOPEVERTICES, POLYTOPICINTERPOLATION,
%            EQUALIZEORDER.

tests = functiontests(localfunctions);
end

function setupOnce(~)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
setupPath();
end

function testPartitionMatchesPaper(testCase)
subsets = partitionInterval([5 30], 5);
expected = [5 10; 10 15; 15 20; 20 25; 25 30];
verifyEqual(testCase, subsets, expected, 'AbsTol', 1e-12);
end

function testPartitionSubsetsShareBoundaries(testCase)
subsets = partitionInterval([5 30], 5);
verifyEqual(testCase, subsets(1:end-1,2), subsets(2:end,1), 'AbsTol', 1e-12, ...
    'Adjacent subsets must share their boundary - that is the switching surface.');
end

function testReducedPolytopeDropsInfeasibleCorner(testCase)
[Vfull, ~] = polytopeVertices([5 30], false);
[Vred,  ~] = polytopeVertices([5 30], true);

verifySize(testCase, Vfull, [2 4]);
verifySize(testCase, Vred,  [2 3]);

% The dropped corner asserts max speed and min speed simultaneously
dropped = Vfull(:,4);
verifyEqual(testCase, dropped, [30; 1/5], 'AbsTol', 1e-12);
end

function testAdjacentSubsetsShareAHyperbolaVertex(testCase)
% Subset p vertex 3 must coincide with subset p+1 vertex 2 - this is what
% makes assumption (A.3.2) implementable.
subsets = partitionInterval([5 30], 5);
V = cell(1,5);
for p = 1:5
    V{p} = polytopeVertices(subsets(p,:), true);
end
for p = 1:4
    verifyEqual(testCase, V{p}(:,3), V{p+1}(:,2), 'AbsTol', 1e-12, ...
        sprintf('Subsets %d and %d do not share a vertex.', p, p+1));
end
end

function testInterpolationRecoversVertices(testCase)
[V, box] = polytopeVertices([5 30], true);
Gv = cat(3, 1, 2, 3);

for j = 1:3
    val = polytopicInterpolation(V(:,j).', box, Gv);
    verifyEqual(testCase, val, Gv(:,:,j), 'AbsTol', 1e-10, ...
        sprintf('Interpolation does not reproduce vertex %d.', j));
end
end

function testEqualizeOrderPreservesResponse(testCase)
sys1 = ss(tf(1, [1 1]));          % order 1
sys2 = ss(tf(1, [1 2 1]));        % order 2
arr  = stack(3, sys1, sys2);      % non-uniform ss array

padded = equalizeOrder(arr, -1e4);

verifyEqual(testCase, nmodels(padded), 2, 'Array length must be preserved.');
verifyEqual(testCase, size(padded(:,:,1).a,1), size(padded(:,:,2).a,1), ...
    'Padded models must share a state dimension.');
verifyEqual(testCase, size(padded(:,:,1).a,1), 2, 'Padding must reach the largest order.');

t  = 0:0.01:10;
y0 = step(sys1, t);
y1 = step(padded(:,:,1), t);
verifyLessThan(testCase, sqrt(mean((y0 - y1).^2)), 1e-6, ...
    'Padding must not change the input-output behaviour.');
end
