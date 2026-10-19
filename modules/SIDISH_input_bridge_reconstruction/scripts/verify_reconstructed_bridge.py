from pathlib import Path
import json,hashlib,csv,shutil,sys,argparse
import numpy as np,pandas as pd,anndata as ad
parser=argparse.ArgumentParser()
parser.add_argument('--input-root',type=Path,required=True)
parser.add_argument('--output-root',type=Path,required=True)
args=parser.parse_args()
W=args.output_root;A=W/'audit';OUT=W/'outputs'
I=args.input_root
f=ad.read_h5ad(I/'processed_adata_sidish.h5ad');r=ad.read_h5ad(OUT/'reconstructed_processed_adata_sidish.h5ad')
fx=f.X.tocsr();rx=r.X.tocsr();fx.sort_indices();rx.sort_indices()
checks={'cell_order':np.array_equal(f.obs_names,r.obs_names),'gene_order':np.array_equal(f.var_names,r.var_names),'CSR_indptr':np.array_equal(fx.indptr,rx.indptr),'CSR_indices':np.array_equal(fx.indices,rx.indices),'CSR_values':np.array_equal(fx.data,rx.data)}
for c in ['n_genes','n_genes_by_counts','total_counts','total_counts_mt','pct_counts_mt']:
 checks['QC_'+c]=np.array_equal(f.obs[c],r.obs[c])
fb=pd.read_csv(I/'processed_bulk_sidish.csv',index_col=0);rb=pd.read_csv(OUT/'reconstructed_processed_bulk_sidish.csv',index_col=0)
checks['bulk_row_order']=fb.index.equals(rb.index);checks['bulk_column_order']=fb.columns.equals(rb.columns)
checks['bulk_reopened_exact']=np.array_equal(fb.to_numpy(),rb.to_numpy())
checks['bulk_reopened_abs_tolerance_1e_14']=np.allclose(fb.to_numpy(),rb.to_numpy(),rtol=0,atol=1e-14)
qc=pd.read_csv(OUT/'cell_selection_audit.tsv',sep='\t');qc['sample']=qc.cell_id.str.split('_').str[0]
flow=[]
for sample,v in qc.groupby('sample'):
 ng=v.n_genes_after_gene_filter>=5000;mt=v.mt_percent_after_gene_filter>=10
 flow.append(dict(sample=sample,input_cells=len(v),excluded_gene_only=int((ng&~mt).sum()),excluded_mt_only=int((mt&~ng).sum()),excluded_both=int((mt&ng).sum()),retained=int(v.keep.sum()),id_membership_mismatch=int((v.keep!=v.in_frozen).sum())))
pd.DataFrame(flow).to_csv(A/'sample_cell_flow.tsv',sep='\t',index=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest().upper()
inputs=list(I.glob('*'));sources=[]
for p in inputs:
 if p.name not in ['sc_cells.tsv','sc_genes.tsv','sc_metadata.csv','sc_counts_genes_by_cells.mtx','bulk_survival.csv','bulk_expression_samples_by_gene.csv','processed_adata_sidish.h5ad','processed_bulk_sidish.csv']:continue
 sources.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
expected={'processed_adata_sidish.h5ad':'415C78C140CDF1F6835D1926E16FD0762092630D17501520E768EF16825307AE','processed_bulk_sidish.csv':'42473A8A58EC81CCDE33400BBBCC30987750FEC52E39EEE8985596AD4450AF9E'}
checks['frozen_inputs_hash_unchanged']=all(x['sha256']==expected[Path(x['path']).name] for x in sources if Path(x['path']).name in expected)
pd.DataFrame(sources).to_csv(A/'input_sha256.tsv',sep='\t',index=False)
result={'checks':{k:bool(v) for k,v in checks.items()},'cell_flow':flow,'scope':'Independent reopen of serialized reconstructed outputs; CSR components compared exactly; original archived hashes verified','original_code_recovered':False,'training_repeated':False}
(A/'independent_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
assert all(v for k,v in checks.items() if k!='bulk_reopened_exact'),result
print(json.dumps(result,ensure_ascii=False,indent=2))
