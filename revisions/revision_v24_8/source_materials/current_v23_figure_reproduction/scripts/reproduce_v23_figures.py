"""Rebuild the integrated v23 figures from the included frozen plotting inputs.

Usage: python3 scripts/reproduce_v23_figures.py INPUT_FOLDER NEW_OUTPUT_FOLDER
INPUT_FOLDER contains source_data/, scripts/, captions/ and figures/.
The output must not exist. No input file or analytical result is modified.
"""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys

src=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve()
assert not out.exists(),'Use a new output folder; existing outputs are never overwritten.'
out.mkdir(parents=True)
for name in ['source_data','scripts']:shutil.copytree(src/name,out/name)
for name in ['figures','captions','audit']:(out/name).mkdir()
for p in (src/'captions').glob('*.txt'):shutil.copy2(p,out/'captions'/p.name)
for name in ['Main_Figure_1.png','Supplementary_Figure_S1.emf','Supplementary_Figure_S2.emf','Supplementary_Figure_S3.emf']:
 shutil.copy2(src/'figures'/name,out/'figures'/name)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((src/'audit/Immutable_source_inputs.json').read_text())
assert len(manifest)==376
assert all(sha(out/x['file'])==x['sha256'] for x in manifest)
(out/'audit/Immutable_source_inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')
with (out/'audit/Build_log.txt').open('w') as log:
 subprocess.run(['Rscript','--vanilla',str(out/'scripts/redesign_v22_figures.R'),str(out)],check=True,stdout=log,stderr=subprocess.STDOUT)
subprocess.run([sys.executable,str(out/'scripts/audit_v22_figure_design.py'),str(out)],check=True)
# The own-PRE plot is now part of Main 5B. S13 contains the four-index display.
for ext in ['png','pdf','svg']:
 old=out/'figures'/f'Supplementary_Figure_S13.{ext}'
 old.unlink()
 (out/'figures'/f'OPTIONAL_S18_Alpha_metrics_by_site_and_fraction.{ext}').rename(old)
expected=json.loads((src/'audit/Current_v23_PNG_EMF_hashes.json').read_text())
assert len(expected)==23
for rel,h in expected.items():assert sha(out/rel)==h,rel
assert all(sha(src/x['file'])==x['sha256'] for x in manifest)
receipt=dict(status='PASS',reproduced_current_figure_assets=23,identical_PNGs=20,
 unchanged_native_EMFs=3,CSV_RDS_inputs=376,new_statistical_tests=0,
 scope='Pixel-identical current figure rebuild from stored results. Does not reproduce or authenticate the unresolved upstream FASTQ pipeline or validate FACS gates.')
(out/'audit/Current_v23_reproduction_result.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
