function gpudiag()
% Call the gpudiag executable.

% Copyright 2001-2023 Verasonics, Inc.  All world-wide rights and
% remedies under all intellectual property laws and industrial
% property laws are reserved.  Verasonics Registered U.S. Patent and
% Trademark Office.

    if ispc()
        callGpuDiag('.exe')
    elseif isunix() && ~ismac() % islinux
        callGpuDiag('.linux64')
    elseif ismac()
        warning('gpudiag:unsupportedOS', 'GPU Toolkit is not supported for macOS');
    else
        warning('gpudiag:unknownOS', 'Unable to detect OS');
    end
end

function callGpuDiag(appSuffix)
    appPath = fullfile(vsv.file.getVSXDir(), 'System', ['gpuDiag' appSuffix]);
    [status, result] = system(appPath);
    if status == 0 && ~isempty(result)
        disp(result)
    else
        warning('gpuDiag command failed.  Try restarting Matlab to fix the problem.\nError reported is: \n\n%s', result);
    end
end
