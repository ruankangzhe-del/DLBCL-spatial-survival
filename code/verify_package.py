"""Verify the delivered payload without importing scientific dependencies."""
from pathlib import Path
import argparse, csv, hashlib, json

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
a = p.parse_args(); root = a.root.resolve()
with (root/'checksums/payload_sha256.tsv').open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
seen = set()
for r in rows:
    path = (root/r['path']).resolve()
    if Path(r['path']).is_absolute() or not path.is_relative_to(root) or r['path'] in seen:
        raise ValueError('INVALID_MANIFEST_PATH:'+r['path'])
    seen.add(r['path'])
    with path.open('rb') as f:
        sha = hashlib.file_digest(f,'sha256').hexdigest().upper()
    if path.stat().st_size != int(r['bytes']) or sha != r['sha256'].upper():
        raise ValueError('PAYLOAD_MISMATCH:'+r['path'])
actual = {x.relative_to(root).as_posix() for x in root.rglob('*') if x.is_file()}
if actual != seen | {'checksums/payload_sha256.tsv'}:
    raise ValueError('UNLISTED_OR_MISSING_PAYLOAD_FILES')
print(json.dumps(dict(status='PASS_DELIVERY_PAYLOAD', files=len(rows), manifest_self_excluded=True)))
