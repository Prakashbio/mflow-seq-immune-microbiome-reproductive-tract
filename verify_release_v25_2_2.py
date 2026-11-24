"""Read-only public-release integrity, source-index and included-FCS checks."""
from pathlib import Path
import csv,hashlib,json
from openpyxl import load_workbook
from zipfile import ZipFile
import xml.etree.ElementTree as ET
R=Path(__file__).resolve().parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
m=json.loads((R/'SHA256SUMS.json').read_text())
for n,h in m.items():assert (R/n).is_file() and sha(R/n)==h,n
actual={p.relative_to(R).as_posix() for p in R.rglob('*') if p.is_file() and not any(x in p.parts for x in ['.git','__pycache__']) and p.name!='.DS_Store'}-{'SHA256SUMS.json'}
retained=json.loads((R/'EXISTING_GITHUB_FILES_RETAINED.json').read_text())['files']
extras=actual-set(m)
assert extras<=set(retained),('Unlisted files',sorted(extras-set(retained)))
for n in extras:
 b=(R/n).read_bytes()
 assert len(b)==retained[n]['bytes'] and hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==retained[n]['git_blob_sha1'],n
assert set(m)<=actual
V=R/'revisions/revision_v25_2_2'
rows=list(csv.DictReader((V/'Publication_table_source_index.csv').open()))
assert len(rows)==80
assert not any(row['worksheet']=='S2C.4' for row in rows)
for row in rows:assert sha(V/row['release_relative_path'])==row['sha256'],row['data_block']
workbook=load_workbook(V/'selected_tables/Supplementary_Tables_CLEAN.xlsx',read_only=True,data_only=False)
assert len(workbook.sheetnames)==38 and 'S2C.4' not in workbook.sheetnames
for sheet in workbook:
 for cells in sheet:
  for cell in cells:
   assert 'S2C.4' not in str(cell.value),(sheet.title,cell.coordinate)
workbook.close()
withdrawn_csv_hashes={
 'a41b0b57af0ff146a43f8952fa8f5b26e9c698781c7003e6b4b0900e4bc64612',
 'a6538af14cb75d4217162d9b61d4ddd4bf96e68312e700fee3db1d7dadc85302',
}
for name in actual:
 path=R/name
 assert 'Table_S2C_4' not in path.name,('Withdrawn table export',name)
 if path.suffix.lower()=='.csv':
  assert sha(path) not in withdrawn_csv_hashes,('Duplicate withdrawn CSV',name)
 if path.suffix.lower()=='.xlsx':
  with ZipFile(path) as archive:
   xml=ET.fromstring(archive.read('xl/workbook.xml'))
   sheets=xml.findall('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}sheets/{http://schemas.openxmlformats.org/spreadsheetml/2006/main}sheet')
   assert all(sheet.get('name')!='S2C.4' for sheet in sheets),('Withdrawn table in older workbook',name)
access=json.loads((R/'SOURCE_ACCESS.json').read_text())
assert access['withheld_human_sample_FCS_count']==34
fcs=[p for p in R.rglob('*') if p.is_file() and p.suffix.lower()=='.fcs']
assert len(fcs)==3 and {p.name for p in fcs}=={'HUT093_BSA control.fcs','HUT093_PBS1 control.fcs','HUT093_PBS2 control.fcs'}
print(json.dumps(dict(status='PASS',release_payload_files=len(m),source_blocks=len(rows),retained_existing_GitHub_files_verified=len(extras),public_sample_free_buffer_FCS=3,withheld_human_FCS_count_from_access_statement=34,scope='Public manifest/source agreement and included-FCS inventory. Private source identities are not published or required for this check.'),indent=2))
