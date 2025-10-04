📄 HOW_TO_USE.txt
🔬 Drug–Target Affinity Fetcher

This script allows you to query ChEMBL, PubChem, IUPHAR, and BindingDB to retrieve experimental assay results for a drug–protein pair.

⚙️ Requirements

Python ≥ 3.9

Packages listed in requirements.txt

Internet connection (for online fetch)

⚙️ How to use

1️⃣ Install dependencies with:

pip install -r requirements.txt


2️⃣ Run with a drug name:

python binding_fetch_online.py --drug-name Cabazitaxel --protein example_inputs/cancer_targets/tubb.fasta --outdir results --pubchem-keep-all


3️⃣ Or with a SMILES string (preferred for novel molecules):

python binding_fetch_online.py --smiles "CC[C@H]1C(=O)O..." --protein example_inputs/cancer_targets/egfr.fasta --outdir results


4️⃣ Results are saved under results/:

report_online.md → global summary

report_chembl.md, report_pubchem.md, report_iuphar.md, report_bindingdb.md

⚠️ Note: this tool aggregates existing experimental data. For completely new molecules with no assays, the next step is to integrate deep learning predictors (e.g., DeepDTA, GraphDTA) for computational forecasts before lab validation.
