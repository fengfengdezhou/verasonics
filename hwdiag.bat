@echo off
rem
rem  Copyright 2001-2018 Verasonics, Inc.  All world-wide rights and
rem  remedies under all intellectual property laws and industrial
rem  property laws are reserved.  Verasonics Registered U.S. Patent and
rem  Trademark Office.
rem
rem  hwdiag.bat - Wrapper script to set the environment for running
rem               the HwDiag binary that exists in the same directory
rem               tree as this wrapper script.
rem

pushd "%~dp0"

System\HwDiag.exe %*

popd

exit /B %ERRORLEVEL%
