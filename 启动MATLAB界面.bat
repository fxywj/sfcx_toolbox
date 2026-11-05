@echo off
rem ============================================================
rem  SFCX MATLAB app launcher (double-click me)
rem  Opens MATLAB, adds the toolbox to the path, runs sfcx.app
rem ============================================================
setlocal

rem -- edit this line to point at YOUR matlab.exe -------------------------
set "MATLAB_EXE=C:\Path\To\Your\MATLAB\bin\matlab.exe"
rem ----------------------------------------------------------------------

if not exist "%MATLAB_EXE%" (
    echo [ERROR] MATLAB not found at:
    echo   %MATLAB_EXE%
    echo Edit the MATLAB_EXE line inside this file, then run it again.
    pause
    exit /b 1
)

set "TOOLBOX=%~dp0matlab"

echo ============================================================
echo   SFCX - structure-function coupling explorer (MATLAB)
echo   Toolbox : %TOOLBOX%
echo   MATLAB is starting, the app window follows shortly.
echo ============================================================
start "" "%MATLAB_EXE%" -r "addpath('%TOOLBOX%'); sfcx.app"
