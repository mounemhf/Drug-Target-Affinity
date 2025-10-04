@echo off
cd /d %~dp0
if not exist .venv (python -m venv .venv)
call .\.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install --upgrade -r requirements.txt

set "DRUG_NAME="
set "SMILES="
set "KEEPALL=--pubchem-keep-all"

echo.
echo Run one drug against ALL FASTAs in example_inputs\cancer_targets
echo (leave Drug name empty if you will provide SMILES)
set /p DRUG_NAME=Drug name (optional if SMILES provided): 
set /p SMILES=SMILES (optional; overrides drug name if provided): 

if not "%SMILES%"=="" (
  python run_one_drug_all_targets.py --smiles "%SMILES%" --drug-name "%DRUG_NAME%" %KEEPALL%
) else (
  python run_one_drug_all_targets.py --drug-name "%DRUG_NAME%" %KEEPALL%
)

echo.
echo Done. Check per-target result folders under .\results\
pause
