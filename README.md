📄 HOW_TO_USE.txt
🔬 Drug–Target Affinity Fetcher

This script allows you to query ChEMBL, PubChem, IUPHAR, and BindingDB to retrieve experimental assay results for a drug–protein pair.

⚙️ Requirements

Python ≥ 3.9

Packages listed in requirements.txt

Internet connection (for online fetch)

Install dependencies (inside project folder):

python -m venv .venv
.\.venv\Scripts\activate   # Windows
pip install --upgrade pip
pip install -r requirements.txt

▶️ Usage
1. Run with a drug name
python binding_fetch_online.py --drug-name Cabazitaxel --protein example_inputs/cancer_targets/tubb.fasta --outdir results --pubchem-keep-all

2. Run with a SMILES string (preferred for accuracy)
python binding_fetch_online.py --smiles "CC[C@H]1C(=O)O..." --protein example_inputs/cancer_targets/egfr.fasta --outdir results

3. Run one drug against all targets in a folder
python run_one_drug_all_targets.py --drug-name Gefitinib --targets-dir example_inputs/cancer_targets --outroot results

4. Batch mode (CSV list of drugs and proteins)
run_batch_list.bat batch_list.csv

📂 Output

Results are saved under the results/ folder:

report_online.md → Global summary & interpretation

report_chembl.md → ChEMBL-specific

report_pubchem.md → PubChem-specific

report_iuphar.md → IUPHAR-specific

report_bindingdb.md → BindingDB-specific

summary.json + CSV files with raw records

⚠️ Notes

If no records are found, it means there are no curated assays for that drug–target pair.

The tool aggregates experimental data only. For novel molecules without assays, you will need predictive models (e.g., DeepDTA).
