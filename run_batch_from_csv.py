# -*- coding: utf-8 -*-
import csv, sys, os, hashlib, shlex, subprocess
from pathlib import Path

PY = sys.executable  # current python

def slugify(s: str) -> str:
    s = ''.join(c if c.isalnum() or c in ('-','_') else '_' for c in s.strip())
    while '__' in s: s = s.replace('__','_')
    return s.strip('_') or 'item'

def outdir_for(drug_name, smiles, fasta_path):
    tgt = Path(fasta_path).stem
    if smiles:
        h = hashlib.sha1(smiles.encode('utf-8')).hexdigest()[:8]
        return f'./results/{tgt}__smiles_{h}'
    if drug_name:
        return f'./results/{tgt}__{slugify(drug_name)}'
    return f'./results/{tgt}__item'

def run_one(drug_name, smiles, fasta_path, outdir):
    cmd = [PY, 'binding_fetch_online.py', '--protein', fasta_path, '--outdir', outdir, '--pubchem-keep-all']
    if smiles:
        cmd += ['--smiles', smiles]
        if drug_name:
            cmd += ['--drug-name', drug_name]
    else:
        if not drug_name:
            print('[SKIP] Need at least drug_name or smiles:', fasta_path)
            return 0
        cmd += ['--drug-name', drug_name]
    print('>>', ' '.join(shlex.quote(c) for c in cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print(f'[WARN] fetch failed (code {r.returncode}) for: drug={drug_name} smiles={bool(smiles)} fasta={fasta_path}')
        return r.returncode
    r2 = subprocess.run([PY, 'make_per_source_reports.py'])
    if r2.returncode != 0:
        print(f'[WARN] per-source reports failed (code {r2.returncode}) for outdir={outdir}')
    return 0

def main():
    if len(sys.argv) < 2:
        print('Usage: python run_batch_from_csv.py input.csv')
        return 2
    csv_path = Path(sys.argv[1])
    if not csv_path.exists():
        print('CSV not found:', csv_path)
        return 2
    with csv_path.open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader, start=1):
            drug_name = (row.get('drug_name') or '').strip()
            smiles = (row.get('smiles') or '').strip()
            fasta_path = (row.get('fasta_path') or '').strip().replace('\','/')
            outdir = (row.get('outdir') or '').strip()
            if not outdir:
                outdir = outdir_for(drug_name, smiles, fasta_path)
            Path(outdir).mkdir(parents=True, exist_ok=True)
            print(f'\n=== [{i}] {drug_name or "(SMILES)"} vs {fasta_path} -> {outdir} ===')
            run_one(drug_name, smiles, fasta_path, outdir)
    print('\n[Done] Batch complete.')
    return 0

if __name__ == '__main__':
    sys.exit(main())
