function BuildVerificationTest(varargin)
% BuildVerificationTest  Compile VVT, PVT, or PCT
%
%   Generates an executable for either the Verasonics Verification
%   Test (VVT), Production Verification Test (PVT), or Power Cycling
%   Test (PCT).
%
%   BuildVerificationTest creates an executable VVT.exe.
%   BuildVerificationTest(appName) creates an executable for the specified
%   app, either VVT, PVT, or PCT.
%
%   Example:
%   BuildWithMcc('PVT') creates PVT\Deploy\PVT.exe
%   BuildWithMcc('PCT') creates PCT\Deploy\PCT.exe
%
% VVT can be built from a standard Verasonics distribution.
% PVT must be built from a Verasonics manufacturing test distribution.
% PCT must be built from a Verasonics manufacturing test distribution.

% Copyright (C) 2001-2024, Verasonics, Inc.
% All worldwide rights and remedies under all intellectual
% property laws and industrial property laws are reserved.

activate();
vpfRoot = fileparts(which('activate'));

if ~ispc
    error('%s only supports Windows.', mfilename)
end

if ~license('test', 'compiler')
    error('MATLAB Compiler must be installed to run %s', mfilename)
end

app = 'VVT';
simFile = false;
if nargin > 0
    app = varargin{1};
end
if nargin > 1
    simFile = varargin{2};
end
if nargin > 2
    error('Invalid number of input arguments');
end
isVVT = strcmp(app, 'VVT');
isPVT = strcmp(app, 'PVT');
isPCT = strcmp(app, 'PCT');
if ~isPVT && ~isVVT && ~isPCT
    error('The app name must be VVT, PVT, or PCT');
end
% The entry-point script for the application.
if isPCT
    appScript = 'PowerCyclingTest(varargin{:})';
else
    appScript = app;
end

% ManufacturingTest needed for PVT and PCT.
if ~isVVT && ~isfolder('ManufacturingTest')
    error(['ManufacturingTest directory not found. Manufacturing or ' ...
           'source installation is required for building PCT and PVT']);
end

% License not needed for PCT.
if ~isPCT && ~isfile('System/licenseMgr.p')
    error(['Required file licenseMgr.p not found. Please copy ' ...
           'license.enc to the root directory, restart MATLAB, ' ...
           'and run activate.'])
end

% Warn about removing output directory.
projectName = app;
projectDir = ['./' projectName];
if isfolder(projectDir)
    warningState = warning();
    warning('off', 'backtrace')
    warning('Contents of the output directory "%s" will be deleted', ...
            projectName);
    warning(warningState);
end

% Create directories for project files and compiler output.
% Note: Forward slashes are used for all mcc arguments, because it
% leads to more consistently formatted debug statements.
outputDir = [projectDir '/Deploy'];
if isfolder(projectDir)
    projectDirFullPath = fullfile(vpfRoot, projectDir);
    if isOnPath(projectDirFullPath)
        rmpath(projectDirFullPath);
    end
    outputDirFullPath = fullfile(vpfRoot, outputDir);
    if isOnPath(outputDirFullPath)
        rmpath(outputDirFullPath);
    end
    rmdir(projectDir, 's');
end
mkdir(projectDir);
mkdir(outputDir);

function writeline(fstring, varargin)
    fprintf(fid, sprintf('%s\n', fstring), varargin{:});
end

% Create launch script.
launchScript = [projectDir '/launch.m'];
fid = fopen(launchScript, 'w');
writeline("function launch(varargin)");
writeline("activate();");
if simFile
    writeline("%s;", simFile);
end
writeline("try");
writeline("    %s;", appScript);
writeline("catch ex");
writeline("    disp(getReport(ex));");
writeline("end");
fclose(fid);

if isPCT
    getExecutablePath = [projectDir '/getPctExePath.m'];
    fid = fopen(getExecutablePath, 'w');
    writeline("function filePath = getPctExePath()");
    writeline("[~, result] = system('set PATH');");
    writeline("filePath = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));");
    fclose(fid);
end

additionalAdds = {};
if isPVT || isPCT
    additionalAdds = [additionalAdds {'-a', './ManufacturingTest'}];
end
if isPCT
    additionalAdds = [additionalAdds {'-a', getExecutablePath}];
end
if simFile
    simPath = which(simFile);
    if isempty(simPath)
        error("Simulation script not found '%s'", simFile)
    end
    additionalAdds = [additionalAdds {'-a', simPath}];
end
if isPVT || isVVT
    additionalAdds = [additionalAdds { ...
        '-a', './System/licenseMgr.p', ...
        '-a', './VSX.m', ...
    }];
end

fprintf(1, "Creating executable '%s'.\n", [projectName '.exe']);
mcc('-m', launchScript, ...
    '-R', ['-logfile,' projectName '.log'], ...
    '-o', projectName, ...
    '-d', outputDir, ...
    '-a', './activate.p', ...
    '-a', './HardwareTest', ...
    '-a', './Hal/lib', ... % TODO: Remove if VTS-668 is implemented
    '-a', './HwDiag/lib', ... % TODO: Remove if VTS-667 is implemented
    '-a', './System/*', ... % Exclude System/wd-install and System/app-install
    '-a', './System/+vsv', ...
    '-a', './System/P_Files', ...
    '-a', './System/Resource', ...
    '-a', './Utilities', ...
    '-a', fullfile(toolboxdir('compiler'), 'mcrversion.m'), ... % runAcq has this dependency for deployed application
    '-a', fullfile(toolboxdir('signal'), 'signal', 'sgolayfilt.m'), ... % Dependency analyzer was not including sgolayfilt
    additionalAdds{:}, ...
    '-a', projectDir)

close all;  % Close figures created during execution of mcc
end

function onPath = isOnPath(dir)
    pathCell = regexp(path(), pathsep, 'split');
    if ispc  % Windows is not case-sensitive
        onPath = any(strcmpi(dir, pathCell));
    else
        onPath = any(strcmp(dir, pathCell));
    end
end
