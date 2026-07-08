function BuildReprogramHardware()
% BuildReprogramHardware  Compile reprogramHardware

% Copyright (C) 2018-2024, Verasonics, Inc.
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

if ~isfile('System/licenseMgr.p')
    error(['Required file licenseMgr.p not found. Please copy ' ...
           'license.enc to the root directory, restart MATLAB, ' ...
           'and run activate.'])
end

% Warn about removing output directory.
projectName = 'reprogramHardware';
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
writeline("try");
writeline("    reprogramHardware(varargin{:});");
writeline("catch ex");
writeline("    disp(getReport(ex));");
writeline("end");
writeline("com.verasonics.hal.hardware.Hardware.closeHardware();");
fclose(fid);

fprintf(1, "Creating executable '%s'.\n", [projectName '.exe']);
mcc('-m', launchScript, ...
    '-R', ['-logfile,' projectName '.log'], ...
    '-o', projectName, ...
    '-d', outputDir, ...
    '-a', './activate.p', ...
    '-a', './reprogramHardware.p', ...
    '-a', './Hal/lib/fpga/', ...
    '-a', './HwDiag/lib/CommandLineSubstitutions.txt', ...
    '-a', './HwDiag/lib/CommandLineSubstitutionsFilename.txt', ...
    '-a', './System/+vsv/+multi/+util/isSecondarySystem.m', ... % dependency of activate
    '-a', './System/matlab-verasonics-loader-0.1.0.jar', ...
    '-a', './System/libverasonicshal.dll', ...
    '-a', './System/verasonics-common-0.1.0.jar', ...
    '-a', './System/verasonics-jhal-0.1.2.jar', ...
    '-a', './System/libVerasonicsHal-c.dll', ...
    '-a', './System/libVerasonicsJniHal.dll', ...
    '-a', './System/libVerasonicsCommon.dll', ...
    '-a', './System/libVerasonicsCommon-Jni.dll', ...
    '-a', './System/libwindriver-hwadapter.dll', ...
    '-a', './System/libwinpthread-1.dll', ...
    '-a', './System/licenseMgr.p', ...
    '-a', './System/HwDiag.exe', ...
    '-a', './Utilities/hwdiag.m', ...
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
