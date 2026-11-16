"""Read-only portable release integrity check; Python standard library only."""
from pathlib import Path
import csv,hashlib,json
R=Path(__file__).resolve().parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
m=json.loads((R/'SHA256SUMS.json').read_text())
for n,h in m.items():assert (R/n).is_file() and sha(R/n)==h,n
actual={p.relative_to(R).as_posix() for p in R.rglob('*') if p.is_file() and p.name!='SHA256SUMS.json' and not any(x in p.parts for x in ['.git','__pycache__']) and not p.name.startswith('.DS_Store')}
# Older nested version manifests remain covered by the root manifest.
actual.update(n for n in m if Path(n).name=='SHA256SUMS.json')
assert actual==set(m),('Unlisted/missing files',sorted(actual-set(m)),sorted(set(m)-actual))
V=R/'revisions/revision_v25_2_1'
rows=list(csv.DictReader((V/'Publication_table_source_index.csv').open()))
assert len(rows)==82 and len({r['data_block'] for r in rows})==82
for r in rows:assert sha(V/r['release_relative_path'])==r['sha256'],r['data_block']
ledger=json.loads((V/'metadata/Preserved_source_archive_members.json').read_text())
for r in ledger['members']:assert sha(R/r['current_path'])==r['sha256'],r['original_member']
print(json.dumps(dict(status='PASS',release_payload_files=len(m),publication_source_blocks=len(rows),previous_archive_members_preserved=len(ledger['members']),scope='Integrity and source links only; no experimental or upstream-reproduction claim.'),indent=2))
