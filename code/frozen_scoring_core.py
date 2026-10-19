"""Exact frozen scoring functions extracted without numerical changes; source in provenance."""
import numpy as np
import pandas as pd

def frozen_lp(expression,universe,state):
    """Expression-only function: patient x original-unit genes, no clinical fields."""
    prohibited={'os_time_years','os_event','pfs_time_years','pfs_event','event','treatment','coo','ipi','response','EBV'}
    if set(expression.columns)&prohibited:raise ValueError('NONEXPRESSION_FIELD')
    nonzero=state['nonzero_genes_order'];missing=set(nonzero)-set(expression.columns)
    if missing:raise ValueError('MISSING_NONZERO:'+','.join(sorted(missing)))
    background=[g for g in universe if g in expression]
    x=expression[background]
    if not np.isfinite(x.to_numpy()).all():raise ValueError('NONFINITE_CLINICAL_BACKGROUND')
    ranks=x.rank(axis=1,method='average',pct=True)
    columns=state['candidate_genes_order'];mu=np.array(state['scaler_mean']);sd=np.array(state['scaler_scale']);beta=np.array(state['coefficients'])
    rx=pd.DataFrame(np.broadcast_to(mu,(len(x),len(mu))).copy(),index=x.index,columns=columns)
    measured=[g for g in columns if g in ranks];rx.loc[:,measured]=ranks[measured]
    contrib=(rx.to_numpy()-mu)/sd*beta
    raw=contrib.sum(axis=1)-state['coxnet_offset'][0]
    scaled=(raw-state['development_center'])/state['development_sd']
    contributions=pd.DataFrame(contrib/state['development_sd'],index=x.index,columns=columns)[nonzero]
    return scaled,raw,contributions,ranks


def module_scores(expression,modules,cohort):
    """No outcome or covariate parameters; retain whole-cohort z means/SDs."""
    result={};coverage=[];standard=[]
    for name,m in modules.groupby('module',sort=False):
        available=[];weights=[];zvalues=[];missing=[]
        for row in m.itertuples():
            status='AVAILABLE';mean=sd=np.nan
            if row.gene not in expression:status='NOT_MEASURED_OR_IDENTITY_UNRESOLVED'
            else:
                x=expression[row.gene].to_numpy(float)
                if not np.isfinite(x).all():status='NONFINITE_ANY_PATIENT_EXCLUDES_GENE_COHORT'
                else:
                    mean=float(x.mean());sd=float(x.std(ddof=0))
                    if sd<=0:status='CONSTANT_GENE'
                    else:available.append(row.gene);weights.append(row.weight);zvalues.append((x-mean)/sd)
            if status!='AVAILABLE':missing.append(row.gene+':'+status)
            standard.append({'cohort':cohort,'module':name,'gene':row.gene,'mean':mean,'population_sd':sd,'ddof':0,'weight':row.weight,'status':status})
        c=len(available)/len(m);den=float(np.abs(weights).sum());full=float(m.weight.abs().sum())
        status='READY' if c>=.5 and den>0 else 'BLOCKED_LOW_COVERAGE'
        y=np.column_stack(zvalues)@np.array(weights)/den if status=='READY' else np.full(len(expression),np.nan)
        ys=float(np.std(y,ddof=0)) if np.isfinite(y).all() else np.nan
        if status=='READY' and ys<=0:status='BLOCKED_CONSTANT_MODULE';y[:]=np.nan
        coverage.append({'cohort':cohort,'module':name,'n_patients':len(expression),'defined_gene_n':len(m),'available_gene_n':len(available),'gene_coverage':c,'absolute_weight_coverage':den/full,'available_weight_denominator':den,'module_score_population_sd':ys,'status':status,'missing_reasons':';'.join(missing)})
        result[name]=y
    return pd.DataFrame(result,index=expression.index),coverage,standard


