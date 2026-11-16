@echo off
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Please create a Python environment and install requirements.txt first.
  pause
  exit /b 1
)
".venv\Scripts\python.exe" -X utf8 run_cases.py --case all
if errorlevel 1 (
  echo Case execution failed. Review the message above.
  pause
  exit /b 1
)
start "" "results\report.html"
pause
