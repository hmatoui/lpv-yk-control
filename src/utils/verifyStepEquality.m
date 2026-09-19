function [ok, err] = verifyStepEquality(sysList, names, tol, tFinal)
%VERIFYSTEPEQUALITY  Numerically check that several models realize one system.
%
%   [OK, ERR] = VERIFYSTEPEQUALITY(SYSLIST) compares the step responses of
%   the models in the cell array SYSLIST pairwise and returns OK = true when
%   every pairwise RMS difference is below the tolerance. ERR holds those RMS
%   differences.
%
%   [OK, ERR] = VERIFYSTEPEQUALITY(SYSLIST, NAMES, TOL, TFINAL) additionally
%   labels the models in the printed report, sets the tolerance (default
%   1e-3) and the simulation horizon (default 50 s).
%
%   This is the sanity check used throughout the examples to confirm the two
%   Youla-Kucera identities hold numerically:
%
%       G  = N*inv(M)        = inv(Mt)*Nt
%       K0 = U*inv(V)        = inv(Vt)*Ut
%       Ki = (U + M*Q)*inv(V + N*Q) = (Ut + Q*Mt)*inv(Vt + Q*Nt)
%
%   The last identity is the important one: it confirms the YK architecture
%   reproduces each local controller *exactly* rather than approximating it,
%   which is what the performance-recovery argument relies on.
%
%   See also COPRIMEFACTORIZATION, YOULAPARAMETER, YKCONTROLLER.

arguments
    sysList (1,:) cell
    names   (1,:) cell   = {}
    tol     (1,1) double = 1e-3
    tFinal  (1,1) double = 50
end

n = numel(sysList);
if n < 2
    error('verifyStepEquality:tooFew', 'Need at least two models to compare.');
end

if isempty(names)
    names = arrayfun(@(k) sprintf('sys%d', k), 1:n, 'UniformOutput', false);
end

t = 0:0.01:tFinal;
y = cell(1, n);
for k = 1:n
    y{k} = step(sysList{k}, t);
end

err = [];
ok  = true;
for a = 1:n-1
    for b = a+1:n
        e = sqrt(mean((y{a}(:) - y{b}(:)).^2));
        err(end+1) = e; %#ok<AGROW>
        if ~(abs(e) < tol)
            ok = false;
            fprintf('  [FAIL] %s vs %s : RMS = %.3e (tol %.1e)\n', ...
                    names{a}, names{b}, e, tol);
        end
    end
end

if ok
    fprintf('  [ OK ] %s agree (max RMS = %.3e)\n', strjoin(names, ', '), max(err));
end

end
