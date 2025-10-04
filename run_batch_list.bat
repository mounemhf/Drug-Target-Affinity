@echo off
cd /d %~dp0
if not exist .venv (python -m venv .venv)
call .\.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install --upgrade -r requirements.txt
if "%~1"=="" (
  echo Usage: run_batch_list.bat path\to\batch_list.csv
  pause
  exit /b 2
)
python run_batch_from_csv.py "%~1"
echo.
echo Batch finished. Check the "results" folder.
pause
