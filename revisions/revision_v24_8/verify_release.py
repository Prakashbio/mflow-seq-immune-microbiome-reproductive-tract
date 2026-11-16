from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parent
h=json.loads((r/'SHA256SUMS.json').read_text())
for name,expected in h.items():
 p=r/name; a=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):a.update(b)
 assert a.hexdigest()==expected,name
print('PASS:',len(h),'release files match SHA256; no raw FASTQ or laboratory validation is inferred.')
