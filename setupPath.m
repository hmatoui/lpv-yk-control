function setupPath()
%SETUPPATH  Add the LPV-YK source tree to the MATLAB path.
%
%   SETUPPATH() adds every folder under src/ to the MATLAB path for the
%   current session. Run it once before using any function in this
%   repository, or call it at the top of your own script as the bundled
%   examples do.
%
%   The path is not saved. To make it permanent, run SETUPPATH followed by
%   SAVEPATH, or add a call to SETUPPATH in your startup.m.
%
%   See also CHECKDEPENDENCIES.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'src')));

end
