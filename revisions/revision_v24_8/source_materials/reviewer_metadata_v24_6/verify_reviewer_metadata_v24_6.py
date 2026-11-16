"""Verify public metadata hashes and correspondence with canonical analysis inputs."""
from pathlib import Path
import csv,hashlib,json,sys
folder=Path(sys.argv[1])
hashes=json.loads((folder/'Current_metadata_hashes.json').read_text())
for n,h in hashes.items():assert hashlib.sha256((folder/n).read_bytes()).hexdigest()==h,n
def rows(p):return list(csv.DictReader(p.open(newline='')))
a=rows(folder/'Analysis_metadata_149.csv');b=rows(folder.parent/'baseline_v20_1_analysis_release/bundle/source_inputs/source_metadata_original_codes.csv')
assert len(a)==len(b)==149
by={r['SourceCode']:r for r in b}
assert {r['Analysis library ID'] for r in a}==set(by)
for r in a:
 ref=by[r['Analysis library ID']]
 for x,y in [('Participant code','Participant'),('Anatomical site','Location'),('FACS fraction','Fraction'),('Preparation key','SpecimenKey')]:assert r[x]==ref[y],(r['Analysis library ID'],x)
assert len({r['Participant code'] for r in a})==15
assert len({r['Preparation key'] for r in a})==40
assert len(rows(folder/'Separate_archived_200_report_ledger.csv'))==200
assert len(rows(folder/'Seven_historical_PCR_annotations.csv'))==7
print(json.dumps(dict(status='PASS',unchanged_public_metadata_files=len(hashes),libraries=149,women=15,preparations=40,archived_reports=200,historical_annotations=7),indent=2))
