"""Hash-gated execution of the portable fixed-score clinical R analysis."""
from pathlib import Path
import argparse,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--rscript',default='Rscript');a=p.parse_args()
root=a.root.resolve();out=a.out.resolve()
if out==root or (out.exists() and any(out.iterdir())):raise ValueError('OUTPUT_MUST_BE_NEW_EMPTY_DIRECTORY')
subprocess.run([sys.executable,str(root/'code/check_inputs.py'),'--root',str(root)],check=True)
subprocess.run([a.rscript,str(root/'code/clinical_cohort_recheck.R'),str(root),str(out)],check=True)
