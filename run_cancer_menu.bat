@echo off
cd /d %~dp0
if not exist .venv (python -m venv .venv)
call .\.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install --upgrade -r requirements.txt
:menu
cls
echo ======================================
echo   Cancer Binding Runner (Solid + Heme + IO)
echo ======================================
echo -- EGFR / HER2 / ALK / BRAF --
echo  1) Gefitinib         vs EGFR (P00533)
echo  2) Lapatinib         vs HER2/ERBB2 (P04626)
echo  3) Crizotinib        vs ALK (Q9UM73)
echo  4) Vemurafenib       vs BRAF (P15056)
echo.
echo -- DNA damage / PARP --
echo  5) Olaparib          vs PARP1 (P09874)
echo  6) Talazoparib       vs PARP1 (P09874)
echo.
echo -- Microtubules (Taxanes) --
echo  7) Docetaxel         vs TUBB (P07437)
echo  8) Cabazitaxel       vs TUBB (P07437)
echo  9) Paclitaxel        vs TUBB (P07437)
echo.
echo -- PI3K/AKT pathway & others --
echo 10) Ipatasertib       vs AKT1 (P31749)
echo 11) Alpelisib         vs PIK3CA (P42336)
echo 12) Cabozantinib      vs MET (P08581)
echo.
echo -- Hematology (TKIs) --
echo 13) Ibrutinib         vs BTK (Q06187)
echo 14) Acalabrutinib     vs BTK (Q06187)
echo 15) Zanubrutinib      vs BTK (Q06187)
echo 16) Imatinib          vs ABL1 (P00519)
echo 17) Dasatinib         vs ABL1 (P00519)
echo 18) Midostaurin       vs FLT3 (P36888)
echo 19) Gilteritinib      vs FLT3 (P36888)
echo.
echo -- Apoptosis / Cell cycle --
echo 20) Venetoclax        vs BCL2 (P10415)
echo 21) Palbociclib       vs CDK4 (P11802)
echo 22) Ribociclib        vs CDK6 (Q00534)
echo.
echo -- Immuno-oncology (small molecules) --
echo 23) BMS-202 (PD-L1 SM) vs PD-L1/CD274 (Q9NZQ7)
echo     Note: Antibodies (e.g., Pembrolizumab) may have sparse numeric data in these sources.
echo.
echo -- Metabolism / Fusions --
echo 24) Ivosidenib        vs IDH1 (O75874)
echo 25) Selpercatinib     vs RET (P07949)
echo 26) Entrectinib       vs ROS1 (P08922)
echo.
echo 27) Custom (drug or SMILES + FASTA path)
echo  0) Exit
echo.
set /p choice=Select option: 
if "%choice%"=="1"  set DRUG=Gefitinib&      set FASTA=example_inputs\cancer_targets\egfr.fasta& goto runpair
if "%choice%"=="2"  set DRUG=Lapatinib&      set FASTA=example_inputs\cancer_targets\erb2_her2.fasta& goto runpair
if "%choice%"=="3"  set DRUG=Crizotinib&     set FASTA=example_inputs\cancer_targets\alk.fasta& goto runpair
if "%choice%"=="4"  set DRUG=Vemurafenib&    set FASTA=example_inputs\cancer_targets\braf.fasta& goto runpair
if "%choice%"=="5"  set DRUG=Olaparib&       set FASTA=example_inputs\cancer_targets\parp1.fasta& goto runpair
if "%choice%"=="6"  set DRUG=Talazoparib&    set FASTA=example_inputs\cancer_targets\parp1.fasta& goto runpair
if "%choice%"=="7"  set DRUG=Docetaxel&      set FASTA=example_inputs\cancer_targets\tubb.fasta& goto runpair
if "%choice%"=="8"  set DRUG=Cabazitaxel&    set FASTA=example_inputs\cancer_targets\tubb.fasta& goto runpair
if "%choice%"=="9"  set DRUG=Paclitaxel&     set FASTA=example_inputs\cancer_targets\tubb.fasta& goto runpair
if "%choice%"=="10" set DRUG=Ipatasertib&    set FASTA=example_inputs\cancer_targets\akt1.fasta& goto runpair
if "%choice%"=="11" set DRUG=Alpelisib&      set FASTA=example_inputs\cancer_targets\pik3ca.fasta& goto runpair
if "%choice%"=="12" set DRUG=Cabozantinib&   set FASTA=example_inputs\cancer_targets\met.fasta& goto runpair
if "%choice%"=="13" set DRUG=Ibrutinib&      set FASTA=example_inputs\cancer_targets\btk.fasta& goto runpair
if "%choice%"=="14" set DRUG=Acalabrutinib&  set FASTA=example_inputs\cancer_targets\btk.fasta& goto runpair
if "%choice%"=="15" set DRUG=Zanubrutinib&   set FASTA=example_inputs\cancer_targets\btk.fasta& goto runpair
if "%choice%"=="16" set DRUG=Imatinib&       set FASTA=example_inputs\cancer_targets\abl1.fasta& goto runpair
if "%choice%"=="17" set DRUG=Dasatinib&      set FASTA=example_inputs\cancer_targets\abl1.fasta& goto runpair
if "%choice%"=="18" set DRUG=Midostaurin&    set FASTA=example_inputs\cancer_targets\flt3.fasta& goto runpair
if "%choice%"=="19" set DRUG=Gilteritinib&   set FASTA=example_inputs\cancer_targets\flt3.fasta& goto runpair
if "%choice%"=="20" set DRUG=Venetoclax&     set FASTA=example_inputs\cancer_targets\bcl2.fasta& goto runpair
if "%choice%"=="21" set DRUG=Palbociclib&    set FASTA=example_inputs\cancer_targets\cdk4.fasta& goto runpair
if "%choice%"=="22" set DRUG=Ribociclib&     set FASTA=example_inputs\cancer_targets\cdk6.fasta& goto runpair
if "%choice%"=="23" set DRUG=BMS-202&        set FASTA=example_inputs\cancer_targets\cd274_pdl1.fasta& goto runpair
if "%choice%"=="24" set DRUG=Ivosidenib&     set FASTA=example_inputs\cancer_targets\idh1.fasta& goto runpair
if "%choice%"=="25" set DRUG=Selpercatinib&  set FASTA=example_inputs\cancer_targets\ret.fasta& goto runpair
if "%choice%"=="26" set DRUG=Entrectinib&    set FASTA=example_inputs\cancer_targets\ros1.fasta& goto runpair
if "%choice%"=="27" goto custom
if "%choice%"=="0" goto end
goto menu
:runpair
echo.
echo === Running %DRUG% vs %FASTA% ===
python binding_fetch_online.py --drug-name "%DRUG%" --protein "%FASTA%" --outdir results --pubchem-keep-all
python make_per_source_reports.py
echo.
pause
goto menu
:custom
set "DRUG_NAME="
set "SMILES="
set "PROT_PATH="
set /p DRUG_NAME=Drug name (leave empty to use SMILES): 
set /p SMILES=SMILES (optional; overrides drug name if provided): 
echo Tip: drag your FASTA file here to auto-fill the path.
set /p PROT_PATH=Protein FASTA path: 
if not "%SMILES%"=="" (
  python binding_fetch_online.py --smiles "%SMILES%" --drug-name "%DRUG_NAME%" --protein "%PROT_PATH%" --outdir results --pubchem-keep-all
) else (
  python binding_fetch_online.py --drug-name "%DRUG_NAME%" --protein "%PROT_PATH%" --outdir results --pubchem-keep-all
)
python make_per_source_reports.py
pause
goto menu
:end
