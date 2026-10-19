"""Read-only integrity validation of frozen copied sources."""
from pathlib import Path
import csv,hashlib,json,argparse
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root.resolve()
with (root/'provenance/source_access_manifest.tsv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
 x=(root/r['repository_path']).resolve()
 if not x.is_relative_to(root):raise ValueError('MANIFEST_ESCAPE')
 with x.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest().upper()
 if h!=r['sha256']:raise ValueError('INPUT_HASH_MISMATCH:'+r['repository_path'])
print(json.dumps({'status':'PASS_COPIED_SOURCE_HASHES','files':len(rows)}))
