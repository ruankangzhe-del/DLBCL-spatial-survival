"""Portable deterministic replay of frozen clinical LP and expression modules.
No outcomes are supplied to the scoring core. Missing nonzero genes fail closed.
GSE87371 keeps its historical clinical LP incompatibility and is never upgraded.
"""
from pathlib import Path
import argparse,csv,hashlib,json,sys
import numpy as np
import pandas as pd
from frozen_scoring_core import frozen_lp,module_scores

def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest().upper()
def safe(root,rel):
    p=(root/rel).resolve()
    if Path(rel).is_absolute() or not p.is_relative_to(root):raise ValueError('CONFIG_PATH_ESCAPES_REPOSITORY')
    return p
def verify(root):
    with (root/'provenance/source_access_manifest.tsv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    for r in rows:
        if digest(safe(root,r['repository_path']))!=r['sha256']:raise ValueError('INPUT_HASH_MISMATCH:'+r['repository_path'])
    return len(rows)
def read(p,**kw):return pd.read_csv(p,sep='\t',float_precision='round_trip',**kw)
def main(root,out):
    root=root.resolve();out=out.resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('OUTPUT_NOT_EMPTY')
    out.mkdir(parents=True,exist_ok=True)
    count=verify(root)
    cfg=json.loads((root/'config/scoring.json').read_text())
    state=json.loads(safe(root,cfg['state']).read_text());universe=read(safe(root,cfg['universe'])).gene.tolist();modules=read(safe(root,cfg['modules']))
    if len(universe)!=2099 or len(state['candidate_genes_order'])!=100 or len(state['nonzero_genes_order'])!=10:raise ValueError('FROZEN_MODEL_DIMENSION_CHANGED')
    checks=[];combined=[]
    for cohort,setting in cfg['cohorts'].items():
        expression=read(safe(root,setting['expression']),index_col=0).sort_index()
        if expression.index.duplicated().any() or len(expression)!=setting['expected_n']:raise ValueError('PATIENT_ID_OR_COUNT_MISMATCH')
        lp,raw,contrib,ranks=frozen_lp(expression,universe,state)
        m,coverage,standard=module_scores(expression,modules,cohort)
        # The expected reference is read only after deterministic expression-only scoring.
        expected=read(safe(root,cfg['expected_scores']))
        expected=expected[expected.cohort.eq(cohort)].set_index('patient_id').loc[expression.index]
        if setting['clinical_status']=='READY':
            delta=float(np.max(np.abs(lp-expected.clinical_LP_per_development_SD)))
            if delta>=cfg['strict_absolute_tolerance']:raise ValueError('FROZEN_LP_DIFFERENCE:'+cohort)
            checks.append(dict(cohort=cohort,check='frozen_LP',n=len(lp),max_abs_difference=delta,status='PASS'))
        else:
            historical=read(root/'data/internal_restricted/clinical/40_locked_patient_scores.tsv')
            historical=historical[historical.cohort.eq(cohort)&historical.model.eq('elastic_net')].set_index('patient_id')
            candidate=pd.Series(lp,index=expression.index)
            delta=float(np.max(np.abs(candidate.loc[historical.index]-historical.risk_score_per_development_sd)))
            if delta<cfg['strict_absolute_tolerance']:raise ValueError('HISTORICAL_COMPATIBILITY_GATE_CHANGED')
            recorded=read(safe(root,cfg['compatibility_reference']));rec=float(recorded[recorded.cohort.eq(cohort)].max_absolute_score_delta.iloc[0])
            if abs(delta-rec)>1e-10:raise ValueError('COMPATIBILITY_GAP_NOT_REPRODUCED')
            checks.append(dict(cohort=cohort,check='clinical_LP_compatibility',n=len(lp),max_abs_difference=delta,status='EXPECTED_BLOCKED_FROZEN_LP_COMPATIBILITY'))
            lp[:]=np.nan;raw[:]=np.nan
        for column in m.columns:
            delta=float(np.max(np.abs(m[column].to_numpy()-expected[column].to_numpy())))
            if not np.isfinite(delta) or delta>1e-10:raise ValueError('MODULE_DIFFERENCE:'+cohort+':'+column)
            checks.append(dict(cohort=cohort,check=column,n=len(m),max_abs_difference=delta,status='PASS'))
        result=m.copy();result.insert(0,'clinical_LP_per_development_SD',lp);result.insert(0,'clinical_LP_raw',raw);result.insert(0,'cohort',cohort);result['clinical_status']=setting['clinical_status'];combined.append(result.reset_index())
        pd.DataFrame(coverage).to_csv(out/(cohort+'_module_coverage.tsv'),sep='\t',index=False)
    pd.concat(combined,ignore_index=True).to_csv(out/'frozen_scores.tsv.gz',sep='\t',index=False,float_format='%.17g',na_rep='NA')
    pd.DataFrame(checks).to_csv(out/'score_comparison.tsv',sep='\t',index=False)
    (out/'validation.json').write_text(json.dumps({'status':'PASS_SCOPED_FROZEN_SCORING_WITH_EXPECTED_GSE87371_BLOCK','source_hash_checks':count,'patient_expression_records':sum(x['expected_n'] for x in cfg['cohorts'].values()),'clinical_LP_reproduced_patients':1481,'module_checks':18,'model_training':False,'outcomes_used_in_scoring':False,'external_null_evaluated':False,'checks':checks},indent=2)+'\n')
    print('PASS:1481 frozen clinical LPs;1704 expression records;18 module checks;GSE87371 compatibility block preserved')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args();main(args.root,args.out)
