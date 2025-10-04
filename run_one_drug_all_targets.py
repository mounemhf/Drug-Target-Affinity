# -*- coding: utf-8 -*-
import argparse, hashlib, subprocess, shlex, sys
from pathlib import Path

PY = sys.executable

def slugify(s: str) -> str:
    s = ''.join(c if c.isalnum() or c in ('-','_') else '_' for c in s.strip())
    while '__' in s: s = s.replace('__','_')
    return s.strip('_') or 'item'

def smiles_hash(smiles: str) -> str:
    import hashlib
    return hashlib.sha1(smiles.encode('utf-8')).hexdigest()[:8]

def main():
    ap = argparse.ArgumentParser(description="Run one drug against ALL FASTAs in a folder.")
    ap.add_argument("--drug-name", type=str, help="Ligand name (optional if --smiles provided)")
    ap.add_argument("--smiles", type=str, help="Ligand SMILES (overrides --drug-name if provided)")
    ap.add_argument("--targets-dir", default="example_inputs/cancer_targets", help="Folder containing .fasta targets")
    ap.add_argument("--outroot", default="results", help="Root results folder")
    ap.add_argument("--pubchem-keep-all", action="store_true", help="Keep all PubChem assays (no gene/target filter)")
    args = ap.parse_args()

    tdir = Path(args.targets_dir)
    if not tdir.exists():
        print("Targets folder not found:", tdir); return 2
    fastas = sorted([p for p in tdir.glob("*.fasta") if p.is_file()])
    if not fastas:
        print("No FASTA files in:", tdir); return 2

    if not (args.smiles or args.drug_name):
        print("Provide --drug-name or --smiles"); return 2

    if args.smiles:
        label = f"smiles_{smiles_hash(args.smiles)}"
    else:
        label = slugify(args.drug_name)

    for i, fa in enumerate(fastas, start=1):
        target = fa.stem
        outdir = Path(args.outroot) / f"{target}__{label}"
        outdir.mkdir(parents=True, exist_ok=True)

        cmd = [PY, "binding_fetch_online.py", "--protein", str(fa), "--outdir", str(outdir)]
        if args.pubchem_keep_all:
            cmd += ["--pubchem-keep-all"]
        if args.smiles:
            cmd += ["--smiles", args.smiles]
            if args.drug_name:
                cmd += ["--drug-name", args.drug_name]
        else:
            cmd += ["--drug-name", args.drug_name]

        print(f"\n=== [{i}/{len(fastas)}] {args.drug_name or '(SMILES)'} vs {fa.name} → {outdir} ===")
        print(">>", " ".join(shlex.quote(c) for c in cmd))
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print(f"[WARN] fetch failed ({r.returncode}) for {fa.name}")
            continue

        r2 = subprocess.run([PY, "make_per_source_reports.py", "--outdir", str(outdir)])
        if r2.returncode != 0:
            print(f"[WARN] reports failed ({r2.returncode}) for {outdir}")

    print("\n[Done] All targets processed.")
    print(f"See per-target folders under: {args.outroot}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
