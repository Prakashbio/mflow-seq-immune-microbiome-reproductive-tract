"""Verify the reviewer metadata snapshots without changing scientific data.
Usage: python verify_reviewer_metadata_v24_6.py /path/to/reviewer_metadata_v24_6
"""
from pathlib import Path
import csv,json,hashlib,sys
folder=Path(sys.argv[1]);plan=json.loads((folder/'Workbook_patch_plan.json').read_text());ledger=json.loads((folder/'Submission_language_cleanup.json').read_text())['metadata_description_changes'];allowed={(x['file'],x['csv_row'],x['column']) for x in ledger}
checks=[]
for b in plan['blocks']:
 p=folder/b['file'];old=folder/'historical_originals_NOT_FOR_SUBMISSION'/b['file'];a=list(csv.reader(old.open(newline='')));z=list(csv.reader(p.open(newline='')))
 assert hashlib.sha256(old.read_bytes()).hexdigest()==b['source_sha256']
 assert len(a)==len(z)==b['rows']+1 and all(len(x)==len(y)==b['columns'] for x,y in zip(a,z))
 changed=[]
 for i,(aa,zz) in enumerate(zip(a,z),1):
  for j,(x,y) in enumerate(zip(aa,zz),1):
   if x!=y:assert (b['file'],i,j) in allowed;changed.append((i,j))
 assert z==[[('' if v is None else str(v)) for v in row] for row in b['values']]
 checks.append(dict(file=b['file'],rows=b['rows'],columns=b['columns'],original_sha256=b['source_sha256'],reviewer_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),presentation_only_changes=len(changed)))
print(json.dumps(dict(status='PASS',checks=checks,numeric_data_and_missing_values_unchanged=True),indent=2))
