"""Read-only public RC2 verifier; distinguishes intentional private FCS exclusions."""
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
actual={p.relative_to(R).as_posix() for p in R.rglob('*') if p.is_file() and not any(x in p.parts for x in ['.git','__pycache__']) and p.name!='.DS_Store'}-{'SHA256SUMS.json'}
retained=json.loads((R/'EXISTING_GITHUB_FILES_RETAINED.json').read_text())['files'];extras=actual-set(m);assert extras<=set(retained),('Unlisted files',sorted(extras-set(retained)))
for n in extras:
 b=(R/n).read_bytes();assert len(b)==retained[n]['bytes'] and hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==retained[n]['git_blob_sha1'],n
assert set(m)<=actual
V=R/'revisions/revision_v25_2_2';rows=list(csv.DictReader((V/'Publication_table_source_index.csv').open()));assert len(rows)==82
for row in rows:assert sha(V/row['release_relative_path'])==row['sha256']
ledger=json.loads((V/'metadata/RC2_prior_candidate_preservation_ledger.json').read_text());private=0
for row in ledger['members']:
 if row['disposition']=='WITHHELD_PRIVATE_HUMAN_FCS':
  private+=1;assert not (R/row['RC1_path']).exists();assert row['RC1_sha256'] not in m.values()
 else:assert sha(R/row['current_path'])==row['RC1_sha256']
assert private==34
fcs=list(R.rglob('*.fcs'));assert len(fcs)==3
print(json.dumps(dict(status='PASS',release_payload_files=len(m),source_blocks=len(rows),RC1_members_accounted=len(ledger['members']),private_human_FCS_excluded=private,retained_existing_GitHub_files_verified=len(extras),public_sample_free_buffer_FCS=3,raw_human_event_reanalysis='UNAVAILABLE_IN_PUBLIC_ARCHIVE',reproduction_scope=json.loads((R/'REPRODUCIBILITY_SCOPE.json').read_text()),scope='Public integrity and source links. Figure replay from saved bins is not raw human pilot replay or experimental validation.'),indent=2))
