"""Verify the public file inventory without reading participant data."""
from pathlib import Path
import csv,hashlib,json
root=Path(__file__).resolve().parents[1]
manifest=root/'checksums/public_sha256.tsv'
with manifest.open(encoding='utf-8',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
seen=set()
for r in rows:
 p=(root/r['path']).resolve()
 if not p.is_relative_to(root) or r['path'] in seen: raise ValueError('INVALID_MANIFEST_PATH')
 seen.add(r['path'])
 if not p.is_file() or p.stat().st_size!=int(r['bytes']) or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']: raise ValueError('FILE_MISMATCH:'+r['path'])
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and '.git' not in p.relative_to(root).parts and '__pycache__' not in p.relative_to(root).parts and p.suffix not in {'.pyc','.pyo'}}
if actual!=seen|{'checksums/public_sha256.tsv'}: raise ValueError('UNLISTED_OR_MISSING_FILES')
print(json.dumps({'status':'PASS_PUBLIC_RELEASE_INTEGRITY','files':len(rows),'patient_data_included':False}))
