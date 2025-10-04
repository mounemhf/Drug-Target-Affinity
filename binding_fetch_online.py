#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, json, re, sys, time, math
from pathlib import Path
from typing import Optional, Dict, List
import pandas as pd
import numpy as np

def eprint(*a, **k): print(*a, file=sys.stderr, **k)
def read_text(p: Path)->str: return Path(p).read_text(encoding='utf-8')

def parse_fasta_header_and_seq(path: Path):
    header=''; seq=[]
    for i, line in enumerate(read_text(path).splitlines()):
        if line.startswith('>'):
            if i==0: header=line[1:].strip()
            continue
        seq.append(line.strip())
    return header, ''.join(seq)

def extract_uniprot_from_header(h: str)->Optional[str]:
    m=re.search(r'\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])\b', h)
    return m.group(1) if m else None

def extract_gene_from_header(h:str)->Optional[str]:
    m=re.search(r'\bGN=([A-Za-z0-9_-]+)', h); return m.group(1) if m else None

def extract_protein_name_from_header(h:str)->Optional[str]:
    if '|' in h:
        parts=h.split('|')
        if len(parts)>=3: return parts[2].split()[0]
    return h.split()[0] if h else None

def http_get_json(url, params=None, timeout=20):
    import requests, time
    for _ in range(3):
        try:
            r=requests.get(url, params=params, timeout=timeout, headers={'User-Agent':'DTA-OnlineFetcher/1.0'})
            if r.ok: return r.json()
        except Exception:
            pass
        time.sleep(1.2)
    return None

def resolve_pubchem_by_name(name:str)->Dict[str,Optional[str]]:
    import requests
    out={}
    try:
        url=f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES,InChIKey/JSON"
        data=http_get_json(url)
        if data and 'PropertyTable' in data:
            props=data['PropertyTable']['Properties'][0]
            out['smiles']=props.get('IsomericSMILES'); out['inchikey']=props.get('InChIKey')
        url2=f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/cids/JSON"
        data2=http_get_json(url2)
        if data2 and 'IdentifierList' in data2 and data2['IdentifierList'].get('CID'):
            out['cid']=str(data2['IdentifierList']['CID'][0])
    except Exception: pass
    return out

def resolve_pubchem_by_smiles(smiles:str)->Dict[str,Optional[str]]:
    out={'smiles':smiles}
    try:
        import requests as rq
        url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/IsomericSMILES,InChIKey/JSON"
        data=http_get_json(url, params={'smiles':smiles})
        if data and 'PropertyTable' in data:
            props=data['PropertyTable']['Properties'][0]
            out['smiles']=props.get('IsomericSMILES', smiles)
            out['inchikey']=props.get('InChIKey')
        url2="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"
        r=rq.post(url2, data={'smiles':smiles}, timeout=20)
        if r.ok:
            d=r.json()
            if 'IdentifierList' in d and d['IdentifierList'].get('CID'):
                out['cid']=str(d['IdentifierList']['CID'][0])
    except Exception: pass
    return out

def pubchem_assay_summary(cid:str)->pd.DataFrame:
    data=http_get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON")
    rows=[]
    try:
        for a in data.get('AssaySummaries',{}).get('AssaySummary',[]):
            rows.append({'source':'pubchem','AID':a.get('AID'),'TargetName':a.get('TargetName'),
                         'GeneSymbol':a.get('GeneSymbol'),'ActivityOutcome':a.get('ActivityOutcome'),
                         'AC50':a.get('AC50'),'IC50':a.get('IC50'),'EC50':a.get('EC50'),'Ki':a.get('Ki'),'Kd':a.get('Kd'),
                         'PMID':a.get('PMID')})
    except Exception: pass
    return pd.DataFrame(rows)

def chembl_targets_by_uniprot(uniprot:str)->List[str]:
    base="https://www.ebi.ac.uk/chembl/api/data/target.json"
    data=http_get_json(base, params={'target_components__accession':uniprot,'limit':1000})
    return [t.get('target_chembl_id') for t in (data or {}).get('targets',[]) if t.get('target_chembl_id')]

def chembl_molecule_ids_by_name(name:str)->List[str]:
    base="https://www.ebi.ac.uk/chembl/api/data/molecule.json"
    ids=set()
    d=http_get_json(base, params={'molecule_synonyms__icontains':name,'limit':100})
    for m in (d or {}).get('molecules',[]): 
        mid=m.get('molecule_chembl_id'); 
        if mid: ids.add(mid)
    d2=http_get_json(base, params={'pref_name__iexact':name,'limit':50})
    for m in (d2 or {}).get('molecules',[]): 
        mid=m.get('molecule_chembl_id'); 
        if mid: ids.add(mid)
    return list(ids)

def chembl_activities(target_ids:List[str], molecule_ids:List[str])->pd.DataFrame:
    import requests
    rows=[]; base="https://www.ebi.ac.uk/chembl/api/data/activity.json"
    std_types={'Ki','Kd','IC50','EC50'}
    for tid in target_ids:
        off=0
        while True:
            r=requests.get(base, params={'target_chembl_id':tid,'limit':200,'offset':off},
                           timeout=25, headers={'User-Agent':'DTA-OnlineFetcher/1.0'})
            if not r.ok: break
            acts=r.json().get('activities',[])
            if not acts: break
            for a in acts:
                st=(a.get('standard_type') or '').upper()
                if st in std_types:
                    if molecule_ids and a.get('molecule_chembl_id') not in set(molecule_ids): 
                        continue
                    rows.append({
                        'source':'chembl','target_chembl_id':tid,'molecule_chembl_id':a.get('molecule_chembl_id'),
                        'standard_type':a.get('standard_type'),'standard_value':a.get('standard_value'),
                        'standard_units':a.get('standard_units'),'relation':a.get('standard_relation'),
                        'ligand_name':a.get('molecule_pref_name'),'PMID':a.get('pmid'),'DOI':a.get('doi'),
                        'Journal':a.get('journal'),'Year':a.get('year')
                    })
            off+=200
            if off>6000: break
    return pd.DataFrame(rows)

def iuphar_ligand_ids_by_name(name:str)->List[int]:
    data=http_get_json("https://www.guidetopharmacology.org/services/ligands", params={'name':name})
    return [int(d['ligandId']) for d in (data or []) if 'ligandId' in d]

def iuphar_affinities(ligand_ids:List[int], uniprot:Optional[str])->pd.DataFrame:
    rows=[]
    for lid in ligand_ids:
        data=http_get_json(f"https://www.guidetopharmacology.org/services/ligands/{lid}/interactions")
        if not isinstance(data,list): continue
        for it in data:
            tgt=it.get('target',{}) if isinstance(it.get('target'),dict) else {}
            up=tgt.get('uniprotId')
            if uniprot and up and up.upper()!=uniprot.upper(): continue
            aff=it.get('affinity')
            atype=relation=units=None; value=None
            if isinstance(aff,dict):
                atype=aff.get('type'); relation=aff.get('relation'); value=aff.get('value'); units=aff.get('units')
            elif isinstance(aff,list) and aff:
                a0=aff[0]
                if isinstance(a0,dict):
                    atype=a0.get('type'); relation=a0.get('relation'); value=a0.get('value'); units=a0.get('units')
                else:
                    value=str(a0)
            elif isinstance(aff,str):
                value=aff
            ref=it.get('reference'); pmid=None
            if isinstance(ref,dict): pmid=ref.get('pubmedId')
            rows.append({'source':'iuphar','ligandId':lid,'target_name':tgt.get('name'),
                         'uniprot':up,'type':atype,'relation':relation,
                         'value':value,'units':units,'PMID':pmid})
    return pd.DataFrame(rows)

def bindingdb_online(drug_name:str, protein:str)->pd.DataFrame:
    try:
        import requests
        from bs4 import BeautifulSoup
    except Exception:
        eprint('[WARN] BeautifulSoup not installed; skip BindingDB online')
        return pd.DataFrame()
    base="https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SummaryBindingPage.jsp"
    r=requests.get(base, params={'LigandSearch':drug_name,'target':protein}, timeout=20, headers={'User-Agent':'DTA-OnlineFetcher/1.0'})
    if not r.ok: return pd.DataFrame()
    soup=BeautifulSoup(r.text,'html.parser')
    rows=[]
    for tbl in soup.find_all('table'):
        headers=[th.get_text(strip=True) for th in tbl.find_all('th')]
        if not headers or not any(h.lower().startswith(('ki','kd','ic50','ec50')) for h in headers): 
            continue
        for tr in tbl.find_all('tr'):
            tds=[td.get_text(' ', strip=True) for td in tr.find_all('td')]
            if len(tds)>=3: rows.append({'raw':' | '.join(tds)})
    df=pd.DataFrame(rows); df['source']='bindingdb-online'; return df

RUBRIC=[('Very high',0,1.0),('High',1.0,10.0),('Strong',10.0,100.0),('Moderate',100.0,1000.0),('Weak',1000.0,10000.0),('Very weak/None',10000.0,float('inf'))]

def to_nm(val, units)->Optional[float]:
    if val in (None,'','NA'): return None
    try: v=float(val)
    except Exception: return None
    u='nM'
    if units is not None:
        try:
            if isinstance(units,float) and math.isnan(units):
                u='nM'
            else:
                u=str(units)
        except Exception:
            u='nM'
    u=u.lower().replace('µ','u')
    if u.startswith('nm'): return v
    if u.startswith('um'): return v*1000.0
    if u.startswith('mm'): return v*1_000_000.0
    if u.startswith('pm'): return v*0.001
    return v

def summarize_numeric(vals)->dict:
    arr=[float(x) for x in vals if x is not None and not pd.isna(x)]
    if not arr: return {}
    return {'n':len(arr),'median_nM':float(np.median(arr)),'min_nM':float(np.min(arr)),'max_nM':float(np.max(arr))}

def classify(x:float)->str:
    for label,lo,hi in RUBRIC:
        if lo<=x<hi: return label
    return 'Unknown'

def variability_sentence(s):
    span = s['max_nM']/max(s['min_nM'], 1e-12)
    if span>=1e5: return 'Assay results vary over ≥5 orders of magnitude, indicating strong context dependence.'
    if span>=1e3: return 'Assay results vary over ≥3 orders of magnitude, suggesting assay/context effects.'
    if span>=1e2: return 'Assay results vary over ~2 orders of magnitude, typical across heterogeneous assays.'
    return 'Assay spread is modest, indicating reasonably consistent results.'

def evidence_sentence(n):
    if n>=1000: return f'Evidence base is very strong (n={n} measurements).'
    if n>=200:  return f'Evidence base is strong (n={n}).'
    if n>=50:   return f'Evidence base is moderate (n={n}).'
    if n>=10:   return f'Evidence base is limited (n={n}).'
    return f'Evidence base is sparse (n={n}).'

def best_key(summaries):
    for k in ['Ki','Kd','IC50','EC50']:
        if k in summaries: return k
    return None

def render_report(meta, summaries, out_path:Path):
    lines=[]
    lines.append('# Online Binding Report'); lines.append('')
    lines.append(f"**Ligand**: `{meta.get('drug_name','')}` | **SMILES**: `{meta.get('smiles','')}` | **CID**: `{meta.get('cid','')}`")
    lines.append(f"**Target**: `{meta.get('protein_name','')}` | **UniProt**: `{meta.get('uniprot','')}` | **Gene**: `{meta.get('gene','')}`")
    lines.append(''); lines.append('## Summary (normalized to nM)')
    if not summaries:
        lines.append('_No quantitative affinities found online._')
    else:
        for key in ['Ki','Kd','IC50','EC50']:
            s=summaries.get(key); 
            if not s: continue
            lines.append(f"- **{key}**: n={s['n']}, median={s['median_nM']:.3g} nM, min={s['min_nM']:.3g} nM, max={s['max_nM']:.3g} nM → **{classify(s['min_nM'])}** (best)")
    lines.append('')
    lines.append('## Interpretation (detailed)')
    lines.append('- Assays across sources (ChEMBL/PubChem/IUPHAR/BindingDB) can differ in format and conditions; values are normalized to nM, but heterogeneity remains.')
    key = best_key(summaries)
    if key:
        s = summaries[key]
        med_cls = classify(s['median_nM']).lower()
        best_cls = classify(s['min_nM']).lower()
        lines.append('')
        lines.append(f"**Primary metric: {key}**")
        lines.append(f"- Median potency is **{med_cls}** (median {s['median_nM']:.3g} nM); best case is **{best_cls}** (min {s['min_nM']:.3g} nM)." )
        lines.append(f"- {variability_sentence(s)}")
        lines.append(f"- {evidence_sentence(s['n'])}")
        lines.append('')
        lines.append('### Practical take')
        if med_cls in ('very high','high','strong'):
            lines.append('- Overall activity profile is **compatible with drug-like potency**; multiple assays support pharmacological relevance.')
        elif med_cls in ('moderate',):
            lines.append('- Overall activity is **moderate**; potency may be context-dependent or require optimization/dose selection.')
        else:
            lines.append('- Overall activity appears **weak** in this aggregate; consider orthogonal evidence or alternative targets/chemotypes.')
        lines.append('')
        lines.append('### Suggested next steps')
        lines.append('- Cross-check per-source reports for concordance and outliers.')
        lines.append('- Inspect assay metadata (journal/year, assay type) to explain extremes.')
        lines.append('- If needed, try alternative drug synonyms or run with SMILES to ensure correct compound mapping.')
    else:
        lines.append('')
        lines.append('- No quantitative values parsed; try other names/SMILES or check UniProt mapping in the FASTA header.')
    out_path.write_text('\n'.join(lines), encoding='utf-8')

def main(argv=None):
    ap=argparse.ArgumentParser(description='Online DTA fetcher (ChEMBL, PubChem, IUPHAR, BindingDB)')
    ap.add_argument('--drug-name', type=str, help='Ligand name (e.g., Lapatinib)')
    ap.add_argument('--smiles', type=str, help='Ligand SMILES (overrides --drug-name)')
    ap.add_argument('--protein', type=str, required=True, help='Protein FASTA path')
    ap.add_argument('--outdir', type=str, default='./results')
    ap.add_argument('--pubchem-keep-all', action='store_true',
                    help='Do not filter PubChem assays by gene/target name.')
    args=ap.parse_args(argv)

    outdir=Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    header, seq=parse_fasta_header_and_seq(Path(args.protein))
    uniprot=extract_uniprot_from_header(header) or ''
    gene=extract_gene_from_header(header) or ''
    pname=extract_protein_name_from_header(header) or (header.split()[0] if header else '')

    dname=''; smiles=''; cid=None
    if args.smiles:
        dname=args.drug_name or ''
        res=resolve_pubchem_by_smiles(args.smiles)
        smiles=res.get('smiles') or args.smiles
        cid=res.get('cid')
    elif args.drug_name:
        dname=args.drug_name
        res=resolve_pubchem_by_name(dname)
        smiles=res.get('smiles') or ''
        cid=res.get('cid')
    if not (args.smiles or args.drug_name):
        eprint('ERROR: provide --drug-name or --smiles'); return 2

    chembl_t=chembl_targets_by_uniprot(uniprot) if uniprot else []
    chembl_m=chembl_molecule_ids_by_name(dname) if dname else []
    df_chembl=chembl_activities(chembl_t, chembl_m) if chembl_t else pd.DataFrame()

    df_pubchem=pd.DataFrame()
    if dname and not cid:
        res=resolve_pubchem_by_name(dname); cid=res.get('cid') or cid
    if cid:
        df_pubchem=pubchem_assay_summary(cid)
        if (not args.pubchem_keep_all) and (gene or pname) and not df_pubchem.empty:
            patt=(gene or pname or '').lower()
            df_pubchem = df_pubchem.assign(
                _t=df_pubchem['GeneSymbol'].astype(str).str.lower() + ' ' +
                   df_pubchem['TargetName'].astype(str).str.lower()
            )
            df_pubchem = df_pubchem.loc[df_pubchem['_t'].str.contains(patt, na=False)]
            df_pubchem = df_pubchem.drop(columns=['_t'])

    df_iuphar=pd.DataFrame()
    if dname:
        lids=iuphar_ligand_ids_by_name(dname)
        if lids: df_iuphar=iuphar_affinities(lids, uniprot if uniprot else None)

    df_bdb=pd.DataFrame()
    if dname and (uniprot or pname): df_bdb=bindingdb_online(dname, uniprot or pname)

    summaries={}
    if not df_chembl.empty:
        for t in ['Ki','Kd','IC50','EC50']:
            vals=df_chembl[df_chembl['standard_type'].astype(str).str.upper()==t][['standard_value','standard_units']]
            arr=[to_nm(v,u) for v,u in zip(vals['standard_value'], vals['standard_units'])]
            s=summarize_numeric(arr)
            if s: summaries[t]=s
    if not df_pubchem.empty:
        for t in ['Ki','Kd','IC50','EC50']:
            if t in df_pubchem.columns:
                arr=[to_nm(x,'nM') for x in df_pubchem[t].tolist()]
                s=summarize_numeric(arr)
                if s:
                    if t in summaries:
                        cur=summaries[t]
                        merged = summarize_numeric([cur['median_nM'], cur['min_nM'], cur['max_nM']] + arr)
                        summaries[t]=merged
                    else:
                        summaries[t]=s

    df_chembl.to_csv(outdir/'chembl_records.csv', index=False, encoding='utf-8')
    df_pubchem.to_csv(outdir/'pubchem_records.csv', index=False, encoding='utf-8')
    df_iuphar.to_csv(outdir/'iuphar_records.csv', index=False, encoding='utf-8')
    df_bdb.to_csv(outdir/'bindingdb_online_raw.csv', index=False, encoding='utf-8')

    meta={'drug_name':dname,'smiles':smiles,'cid':cid or '',
          'uniprot':uniprot,'gene':gene,'protein_name':pname}
    (outdir/'summary.json').write_text(json.dumps({'meta':meta, 'summaries':summaries}, ensure_ascii=False, indent=2), encoding='utf-8')
    render_report(meta, summaries, outdir/'report_online.md')

    print('[OK] Online fetch complete.')
    print(f" ChEMBL rows: {len(df_chembl)} | PubChem assays: {len(df_pubchem)} | IUPHAR rows: {len(df_iuphar)} | BindingDB rows: {len(df_bdb)}")
    if summaries:
        parts=[]
        for k in ['Ki','Kd','IC50','EC50']:
            if k in summaries:
                s=summaries[k]
                parts.append(f"{k}: n={s['n']} med={s['median_nM']:.3g} best={s['min_nM']:.3g} nM")
        print(' Summary:', ' | '.join(parts))
    else:
        print(' Summary: No quantitative values parsed (try other names/SMILES or check UniProt mapping).')
    return 0

if __name__=='__main__':
    sys.exit(main())
