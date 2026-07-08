% Compile L11-5vCDIExternalProc script.

% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2018-2024 Verasonics, Inc.

setUpScriptName = 'SetUpL11_5vCDIExternalProc';
matFileName = 'L11-5vCDIExternalProc';
projectDir = './Example_Scripts/Specialty_Applications/ColorDopplerImagingExternal/';

BuildWithMcc(setUpScriptName, matFileName, ...
    {'-a', fullfile(projectDir, 'cdiprocessing.m'), ...
     '-a', fullfile(projectDir, 'vsCDIAdaptiveThreshold.m')});
