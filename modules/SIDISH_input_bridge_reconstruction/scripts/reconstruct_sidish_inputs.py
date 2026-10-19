"""Independent reconstruction; never overwrite historical inputs or checkpoints."""
from pathlib import Path
import json,hashlib,datetime,gc,os,sys,argparse
O=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input-root',type=Path,required=True)
parser.add_argument('--output-root',type=Path,default=O/'input_bridge_reconstruction')
args=parser.parse_args(); W=args.output_root; A=W/'audit'; OUT=W/'outputs'
assert W.resolve()!=args.input_root.resolve() and args.input_root.resolve() not in W.resolve().parents, 'Output must not be inside historical input directory'
A.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
os.environ['MPLCONFIGDIR']=str(W/'runtime/matplotlib');os.environ['NUMBA_CACHE_DIR']=str(W/'runtime/numba')
import numpy as np,pandas as pd,scipy,scipy.sparse as sp
from scipy.io import mmread
import anndata as ad,scanpy as sc,h5py
I=args.input_root
def save(n,x):
 (A/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda z:z.item() if isinstance(z,np.generic) else str(z)),encoding='utf-8')
def log(x):print(x,flush=True)
def matrix_compare(a,b):
 if a.shape!=b.shape:return {'shape_equal':False,'left':a.shape,'right':b.shape}
 a=sp.csr_matrix(a);b=sp.csr_matrix(b);d=a.astype(np.float64)-b.astype(np.float64);d.eliminate_zeros()
 return {'shape_equal':True,'shape':a.shape,'different_values':d.nnz,'max_abs_difference':float(np.max(np.abs(d.data))) if d.nnz else 0.,'exact_equal':d.nnz==0,'a_dtype':str(a.dtype),'b_dtype':str(b.dtype)}
log('Loading frozen input and raw exported counts')
frozen=ad.read_h5ad(I/'processed_adata_sidish.h5ad')
save('frozen_object_structure.json',{'shape':frozen.shape,'obs':list(frozen.obs),'var':list(frozen.var),'uns':list(frozen.uns),'layers':list(frozen.layers),'X_dtype':str(frozen.X.dtype),'raw_shape':frozen.raw.shape if frozen.raw else None})
genes=pd.read_csv(I/'sc_genes.tsv',sep='\t').iloc[:,0].astype(str).to_numpy();cells=pd.read_csv(I/'sc_cells.tsv',sep='\t').iloc[:,0].astype(str).to_numpy()
X=mmread(I/'sc_counts_genes_by_cells.mtx').T.tocsr();assert X.shape==(len(cells),len(genes))
meta=pd.read_csv(I/'sc_metadata.csv',index_col=0);meta.index=meta.index.astype(str);assert np.array_equal(meta.index,cells)
adata=ad.AnnData(X,obs=meta.copy(),var=pd.DataFrame(index=genes));del X;gc.collect()
save('environment.json',{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'scanpy':sc.__version__,'anndata':ad.__version__,'h5py':h5py.__version__})
log('Replaying published raw-input preprocessing rules')
adata.var_names_make_unique();sc.pp.filter_cells(adata,min_genes=3);sc.pp.filter_genes(adata,min_cells=400)
adata.var['mt']=adata.var_names.str.startswith(('MT-','mt-'))
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
gene_filtered_shape=adata.shape
keep=(adata.obs.n_genes_by_counts<5000)&(adata.obs.pct_counts_mt<10)
cell_cmp={'raw_cells':len(cells),'after_gene_filter_shape':gene_filtered_shape,'reconstructed_cells':int(keep.sum()),'frozen_cells':frozen.n_obs,'same_cell_order':np.array_equal(adata.obs_names[keep],frozen.obs_names),'missing_frozen_ids':len(set(frozen.obs_names)-set(adata.obs_names[keep])),'extra_reconstructed_ids':len(set(adata.obs_names[keep])-set(frozen.obs_names))}
save('cell_identity.json',cell_cmp);log(cell_cmp)
pd.DataFrame({'cell_id':adata.obs_names,'n_genes_after_gene_filter':adata.obs.n_genes_by_counts,'mt_percent_after_gene_filter':adata.obs.pct_counts_mt,'keep':keep,'in_frozen':adata.obs_names.isin(frozen.obs_names)}).to_csv(OUT/'cell_selection_audit.tsv',sep='\t',index=False)
adata=adata[keep,:].copy();raw=adata.copy()
qc_compare={}
if cell_cmp['same_cell_order']:
 for col in ['n_genes','n_genes_by_counts','total_counts','total_counts_mt','pct_counts_mt']:
  if col in frozen.obs and col in adata.obs:
   a=np.asarray(adata.obs[col],float);b=np.asarray(frozen.obs[col],float);qc_compare[col]={'exact_equal':bool(np.array_equal(a,b)),'max_abs_difference':float(np.max(np.abs(a-b)))}
save('recomputed_QC_comparison.json',qc_compare)
sc.pp.normalize_total(adata,target_sum=1e4);sc.pp.log1p(adata);sc.pp.highly_variable_genes(adata,flavor='seurat')
bulk=pd.read_csv(I/'bulk_expression_samples_by_gene.csv',index_col=0)
bulk.index=bulk.index.astype(str)
hvg=adata.var.highly_variable.to_numpy();overlap=adata.var_names.isin(bulk.columns)
selected=hvg&overlap;pred_genes=adata.var_names[selected]
gene_cmp={'after_gene_filter_genes':gene_filtered_shape[1],'HVG_count':int(hvg.sum()),'HVG_bulk_intersection':int(selected.sum()),'frozen_genes':frozen.n_vars,'same_gene_order':np.array_equal(pred_genes,frozen.var_names),'missing_frozen_genes':sorted(set(frozen.var_names)-set(pred_genes)),'extra_reconstructed_genes':sorted(set(pred_genes)-set(frozen.var_names))}
save('gene_identity.json',gene_cmp);log(gene_cmp)
adata.var.assign(bulk_overlap=overlap,selected_for_model=selected).to_csv(OUT/'gene_selection_audit.tsv',sep='\t')
rebuilt=raw[:,selected].copy()
if 'cells' in rebuilt.obs:rebuilt.obs_names=rebuilt.obs['cells'].astype(str)
if 'cellType' in rebuilt.obs and 'celltype_major' not in rebuilt.obs:rebuilt.obs.rename(columns={'cellType':'celltype_major'},inplace=True)
comparison={'raw_count_matrix':matrix_compare(rebuilt.X,frozen.X),'normalized_matrix':matrix_compare(adata[:,selected].X,frozen.X),'same_cell_order':cell_cmp['same_cell_order'],'same_gene_order':gene_cmp['same_gene_order']}
save('expression_comparison.json',comparison);log(comparison)
fbulk=pd.read_csv(I/'processed_bulk_sidish.csv',index_col=0);fbulk.index=fbulk.index.astype(str)
surv=pd.read_csv(I/'bulk_survival.csv',index_col=0);surv.index=surv.index.astype(str)
assert list(surv.columns)==['futime','fustat']
reb_bulk=pd.concat([surv.rename(columns={'futime':'duration','fustat':'event'}),bulk.loc[:,pred_genes]],axis=1)
bc={'same_patient_order':np.array_equal(reb_bulk.index,fbulk.index),'same_columns':np.array_equal(reb_bulk.columns,fbulk.columns),'shape':reb_bulk.shape}
if bc['same_patient_order'] and bc['same_columns']:
 diff=reb_bulk.to_numpy(dtype=float)-fbulk.to_numpy(dtype=float)
 bc.update(exact_equal=bool(np.array_equal(reb_bulk.to_numpy(),fbulk.to_numpy())),max_abs_difference=float(np.nanmax(np.abs(diff))),different_values=int(np.count_nonzero(diff)))
save('bulk_comparison.json',bc);log(bc)
rebuilt.write_h5ad(OUT/'reconstructed_processed_adata_sidish.h5ad',compression='gzip')
reb_bulk.to_csv(OUT/'reconstructed_processed_bulk_sidish.csv')
passed=cell_cmp['same_cell_order'] and gene_cmp['same_gene_order'] and comparison['raw_count_matrix'].get('exact_equal',False) and bc.get('exact_equal',False)
save('reconstruction_result.json',{'status':'EXACT_NUMERIC_RECONSTRUCTION' if passed else 'PARTIAL_OR_MISMATCH','cell':cell_cmp,'gene':gene_cmp,'expression':comparison,'bulk':bc,'scope':'Reconstructed implementation from saved exported counts and bulk; not recovered historical entrypoint; no model training','parameters':{'min_genes_cell':3,'min_cells_gene':400,'n_genes_by_counts_lt':5000,'pct_counts_mt_lt':10,'normalization_target_sum_for_HVG':10000,'HVG_flavor':'seurat','batch_correction':False,'returned_expression':'counts before normalize/log, restricted to HVGs intersecting bulk'}})
log('Completed: '+str(passed))
