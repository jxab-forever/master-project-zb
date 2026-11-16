@echo off
cd /d "%~dp0"
".venv\Scripts\python.exe" -X utf8 run_cases.py --plots-only
if errorlevel 1 (
  pause
  exit /b 1
)
start "" "results\report.html"
