function ok = checkDependencies()
%CHECKDEPENDENCIES  Report which required toolboxes and packages are present.
%
%   OK = CHECKDEPENDENCIES() prints a report of the toolboxes and third-party
%   packages this repository needs and returns true when everything required
%   by at least one of the two design pipelines is available.
%
%   Dependency matrix:
%
%     Component                     Grid-based   Partitioned polytopic
%     ---------------------------   ----------   ---------------------
%     Control System Toolbox        required     required
%     Robust Control Toolbox        required     required
%     LPVTools                      required     not used
%     YALMIP                        optional*    required
%     SDP solver (SeDuMi/SDPT3)     optional*    required
%
%     * the grid pipeline can obtain Fg and Fk0 either from LPVTools'
%       LPVSFSYN or from the bundled LMISTATEFEEDBACKGRID, which needs YALMIP.
%
%   See also SETUPPATH.

fprintf('LPV-YK dependency check\n');
fprintf('=======================\n\n');

items = {
    'Control System Toolbox', @() ~isempty(ver('control')),    'required'
    'Robust Control Toolbox', @() ~isempty(ver('robust')),     'required'
    'LPVTools (pgrid/pss)',   @() exist('pgrid','file') == 2 || exist('pgrid','file') == 6, 'grid pipeline'
    'YALMIP (sdpvar)',        @() exist('sdpvar','file') == 2, 'polytopic pipeline'
    'SeDuMi',                 @() exist('sedumi','file') == 2, 'SDP solver'
    'SDPT3',                  @() exist('sdpt3','file') == 2,  'SDP solver'
};

status = false(size(items,1), 1);
for i = 1:size(items,1)
    try
        status(i) = items{i,2}();
    catch
        status(i) = false;
    end
    if status(i)
        mark = 'found   ';
    else
        mark = 'MISSING ';
    end
    fprintf('  [%s] %-24s  (%s)\n', mark, items{i,1}, items{i,3});
end

hasCore    = status(1) && status(2);
hasGrid    = hasCore && status(3);
hasPolytop = hasCore && status(4) && (status(5) || status(6));

fprintf('\n');
fprintf('  Grid-based pipeline      : %s\n', ternary(hasGrid,    'ready', 'unavailable'));
fprintf('  Partitioned polytopic    : %s\n', ternary(hasPolytop, 'ready', 'unavailable'));

ok = hasGrid || hasPolytop;

if ~ok
    fprintf('\n  Neither pipeline can run. See README.md for installation links.\n');
end

end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
