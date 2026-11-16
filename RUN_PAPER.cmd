@echo off
setlocal
cd /d "%~dp0"
set "PYTHONUTF8=1"
if not exist ".venv\Scripts\python.exe" (
  echo Install the Python environment first. See the setup instructions.
  exit /b 1
)
".venv\Scripts\python.exe" -X utf8 -m paper_reproduction.workflow --only nominal --output new_runs\paper %*
set "RUN_STATUS=%ERRORLEVEL%"
if "%~1"=="" pause
exit /b %RUN_STATUS%
