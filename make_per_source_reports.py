# -*- coding: utf-8 -*-
import json, sys, math, argparse
from pathlib import Path
import pandas as pd
import numpy as np
from pandas.errors import EmptyDataError

RUBRIC=[("Very high",0,1.0),("High",1.0,10.0),("Strong",10.0,100.0),("Moderate",100.0,1000.0),("Weak",1000.0,10000.0),("Very weak/None",10000.0,float("inf"))]

def to_nm(value, units):
    if value in (None,"","NA"): 
        return None
    try:
        v=float(value)
    except Exception:
        return None
    u="nM"
    if units is not None:
        try:
            if isinstance(units,float) and math.isnan(units):
                u="nM"
            else:
                u=str(units)
        except Exception:
            u="nM"
    u=u.lower().replace("µ","u")
    if u.startswith("nm"): return v
    if u.startswith("um"): return v*1000.0
    if u.startswith("mm"): return v*1_000_000.0
    if u.startswith("pm"): return v*0.001
    return v

def classify(nm):
    for label,lo,hi in RUBRIC:
        if lo<=nm<hi: return label
    return "Unknown"

def summarize_numeric(vals):
    vals=[x for x in vals if x is not None and not pd.isna(x)]
    if not vals: return {}
    arr=np.array(vals,dtype=float)
    return {"n":int(arr.size),"median_nM":float(np.median(arr)),"min_nM":float(np.min(arr)),"max_nM":float(np.max(arr))}

def best_available_key(summaries):
    for k in ["Ki","Kd","IC50","EC50"]:
        if k in summaries: return k
    return None

def potency_sentence(s):
    med_class = classify(s["median_nM"])
    best_class = classify(s["min_nM"])
    return f"Median potency is **{med_class.lower()}** (median {s['median_nM']:.3g} nM); best reported case is **{best_class.lower()}** (min {s['min_nM']:.3g} nM)."

def variability_sentence(s):
    span = s["max_nM"]/max(s["min_nM"], 1e-12)
    if span>=1e5: return "Assay results vary **over ≥5 orders of magnitude**, indicating strong context dependence (assay types, cell systems, readouts)."
    if span>=1e3: return "Assay results vary **over ≥3 orders of magnitude**, suggesting meaningful context or protocol effects."
    if span>=1e2: return "Assay results vary **over ~2 orders of magnitude**, typical across heterogeneous literature assays."
    return "Assay spread is **modest**, suggesting reasonably consistent results across experiments."

def evidence_sentence(n):
    if n>=1000: return f"Evidence base is **very strong** (n={n} measurements)."
    if n>=200:  return f"Evidence base is **strong** (n={n})."
    if n>=50:   return f"Evidence base is **moderate** (n={n})."
    if n>=10:   return f"Evidence base is **limited** (n={n})."
    return f"Evidence base is **sparse** (n={n})."

def caveats_block(source_name):
    items = [
        "Assays differ in format and conditions (biochemical vs cellular; readouts like inhibition vs viability).",
        "Units and conversions are normalized to nM; residual heterogeneity may remain.",
        "Outliers can reflect specific constructs, mutations, or co-factors.",
        "Cross-database curation standards vary; prefer original publications for critical decisions."
    ]
    out=["## Interpretation (detailed)","",
         f"This section interprets **{source_name}** results in plain language."]
    out += [f"- {x}" for x in items]
    return "\n".join(out)

def make_interpretation(title, summaries, n_records):
    lines=[]
    lines.append("## Biological meaning")
    key = best_available_key(summaries)
    if not key:
        lines.append("- No quantitative values parsed for this source; consider alternative names/SMILES or different target identifiers.")
        return "\n".join(lines)

    s = summaries[key]
    lines.append(f"- Using **{key}** as the primary metric:")
    lines.append(f"  - {potency_sentence(s)}")
    lines.append(f"  - {variability_sentence(s)}")
    lines.append(f"  - {evidence_sentence(s['n'])}")
    lines.append("")
    lines.append("## Practical take")
    med_cls = classify(s["median_nM"]).lower()
    if med_cls in ("very high","high","strong"):
        lines.append("- Overall activity profile is **compatible with drug-like potency**; multiple assays support pharmacological relevance.")
    elif med_cls in ("moderate",):
        lines.append("- Overall activity is **moderate**; potency may be context-dependent or require optimization/dose selection.")
    else:
        lines.append("- Overall activity appears **weak** in this source; consider orthogonal evidence or alternative targets/chemotypes.")
    lines.append("")
    lines.append("## Suggested next steps")
    nexts = [
        "Cross-check other sources (ChEMBL, PubChem, IUPHAR, BindingDB) for concordance.",
        "Inspect per-assay metadata (journal/year, assay type) to explain outliers.",
        "If applicable, align target identifiers (UniProt/Gene) and drug synonyms/SMILES."
    ]
    lines += [f"- {x}" for x in nexts]
    return "\n".join(lines)

def report_lines(title, meta, summaries, source_n):
    lines=[]
    lines.append(f"# {title}")
    lines.append("")
    lines.append(f"**Ligand**: `{meta.get('drug_name','')}` | **SMILES**: `{meta.get('smiles','')}` | **CID**: `{meta.get('cid','')}`")
    lines.append(f"**Target**: `{meta.get('protein_name','')}` | **UniProt**: `{meta.get('uniprot','')}` | **Gene**: `{meta.get('gene','')}`")
    lines.append("")
    lines.append("## Summary (normalized to nM)")
    if not summaries:
        lines.append("_No quantitative affinities found._")
    else:
        for k in ["Ki","Kd","IC50","EC50"]:
            s=summaries.get(k)
            if not s: continue
            best=classify(s["min_nM"])
            lines.append(f"- **{k}**: n={s['n']} | median={s['median_nM']:.3g} nM | min={s['min_nM']:.3g} nM | max={s['max_nM']:.3g} nM → **{best}** (best)")
    lines.append("")
    lines.append(caveats_block(title.split()[0]))
    lines.append("")
    lines.append(make_interpretation(title, summaries, source_n))
    return "\n".join(lines)

def safe_read_csv(path: Path)->pd.DataFrame:
    if not path.exists(): return pd.DataFrame()
    try:
        if path.stat().st_size==0: return pd.DataFrame()
    except Exception: return pd.DataFrame()
    try:
        return pd.read_csv(path, dtype=str)
    except EmptyDataError:
        return pd.DataFrame()
    except Exception:
        return pd.DataFrame()

def write_chembl_report(outdir, meta):
    df=safe_read_csv(outdir/"chembl_records.csv")
    summaries={}
    if not df.empty and "standard_type" in df.columns:
        for t in ["Ki","Kd","IC50","EC50"]:
            sub=df[df["standard_type"].astype(str).str.upper()==t]
            if sub.empty: continue
            arr=[to_nm(v,u) for v,u in zip(sub.get("standard_value",[]), sub.get("standard_units",[]))]
            s=summarize_numeric(arr)
            if s: summaries[t]=s
    n = int(df.shape[0]) if not df.empty else 0
    (outdir/"report_chembl.md").write_text(report_lines("ChEMBL Report", meta, summaries, n), encoding="utf-8")

def write_pubchem_report(outdir, meta):
    df=safe_read_csv(outdir/"pubchem_records.csv")
    summaries={}
    if not df.empty:
        for t in ["Ki","Kd","IC50","EC50"]:
            if t in df.columns:
                arr=[to_nm(x,"nM") for x in df[t].tolist()]
                s=summarize_numeric(arr)
                if s: summaries[t]=s
    n = int(df.shape[0]) if not df.empty else 0
    (outdir/"report_pubchem.md").write_text(report_lines("PubChem Report", meta, summaries, n), encoding="utf-8")

def write_iuphar_report(outdir, meta):
    df=safe_read_csv(outdir/"iuphar_records.csv")
    summaries={}
    if not df.empty and "type" in df.columns:
        for t in ["Ki","Kd","IC50","EC50"]:
            sub=df[df["type"].astype(str).str.upper()==t]
            if sub.empty: continue
            values=[to_nm(v,u) for v,u in zip(sub.get("value",[]), sub.get("units",[]))]
            s=summarize_numeric(values)
            if s: summaries[t]=s
    n = int(df.shape[0]) if not df.empty else 0
    (outdir/"report_iuphar.md").write_text(report_lines("IUPHAR Report", meta, summaries, n), encoding="utf-8")

def write_bindingdb_note(outdir, meta):
    p=outdir/"bindingdb_online_raw.csv"
    title="BindingDB Online Report"
    if not p.exists() or p.stat().st_size==0:
        text=[f"# {title}","","No BindingDB online rows saved."]
    else:
        text=[f"# {title}","",f"**Ligand**: `{meta.get('drug_name','')}` | **Target**: `{meta.get('protein_name','')}`","",
              "Raw HTML rows saved (not fully numeric-parsed in this helper).",
              f"- File: `{p.name}`"]
    (outdir/"report_bindingdb.md").write_text("\n".join(text), encoding="utf-8")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results", help="Folder that contains summary.json and CSV outputs")
    args=ap.parse_args()
    outdir=Path(args.outdir)
    meta_path=outdir/"summary.json"
    if not meta_path.exists():
        print(f"ERROR: {meta_path} not found. Run binding_fetch_online.py first."); return 2
    meta=json.loads(meta_path.read_text(encoding="utf-8")).get("meta", {})
    write_chembl_report(outdir, meta)
    write_pubchem_report(outdir, meta)
    write_iuphar_report(outdir, meta)
    write_bindingdb_note(outdir, meta)
    print(f"Per-source reports written to {outdir}:")
    print(" - report_chembl.md")
    print(" - report_pubchem.md")
    print(" - report_iuphar.md")
    print(" - report_bindingdb.md")
    return 0

if __name__=="__main__":
    sys.exit(main())
