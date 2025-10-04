@echo off
cd /d %~dp0
if not exist .venv (python -m venv .venv)
call .\.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install --upgrade -r requirements.txt
python binding_fetch_online.py --drug-name Gefitinib --protein example_inputs\cancer_targets\egfr.fasta --outdir results --pubchem-keep-all
python make_per_source_reports.py
echo Done. Open results\report_online.md
pause
