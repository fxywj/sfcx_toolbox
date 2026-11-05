@echo off
rem ============================================================
rem  SFCX web GUI launcher (double-click me)
rem  Starts a local server and opens the browser automatically.
rem ============================================================
setlocal
cd /d "%~dp0"
set "PYTHONPATH=%~dp0"

rem -- pick the first Python that actually has the dependencies -------
rem -- add your own interpreter to the list if needed ------------------
set "PY="
for %%C in ("python" "py") do (
  if not defined PY (
    %%~C -c "import numpy,scipy,matplotlib" >nul 2>nul && set "PY=%%~C"
  )
)

if not defined PY (
    echo [ERROR] No Python with numpy/scipy/matplotlib found.
    echo Install Python 3.9+ with those packages, or edit this file.
    pause
    exit /b 1
)

set "PORT=8765"
set "SESSION=%~dp0sfcx_web_session"

echo ============================================================
echo   SFCX - structure-function coupling explorer
echo   Python  : %PY%
echo   Web UI  : http://127.0.0.1:%PORT%/
echo   Session : %SESSION%
echo   Close this window to stop the server.
echo ============================================================
"%PY%" -m sfcx.webapp --port %PORT% --open --root "%SESSION%"
pause
