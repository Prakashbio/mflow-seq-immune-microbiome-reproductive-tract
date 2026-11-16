"""Rebuild current v23.2 figures from frozen v23 inputs and original FACS evidence.

Usage: python3 reproduce_v23_2_figures.py PREPARED_RELEASE_ROOT NEW_OUTPUT_FOLDER
No analytical models are refitted and no original input is overwritten.
"""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys
src=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve()
assert not out.exists(),'Use a new output folder.'
out.mkdir(parents=True);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
old=src/'current_v23_figure_reproduction'
subprocess.run([sys.executable,str(old/'scripts/reproduce_v23_figures.py'),str(old),str(out/'frozen_v23_rebuild')],check=True)
dst=out/'figures';dst.mkdir()
mapping={**{i:i for i in range(1,8)},12:8,13:9,14:10,15:11,16:12,17:13,8:14,9:15,10:16,11:17}
import re
for p in (out/'frozen_v23_rebuild/figures').iterdir():
 if p.suffix not in ['.png','.emf','.pdf','.svg']:continue
 mm=re.fullmatch(r'Supplementary_Figure_S(\d+)',p.stem)
 if mm:name=f'Supplementary_Figure_S{mapping[int(mm[1])]}'+p.suffix
 elif p.stem.startswith('Main_Figure_'):name=p.name
 else:continue
 shutil.copy2(p,dst/name)
evidence=src/'facs_control_pilot_display_v23_2'
subprocess.run([sys.executable,str(evidence/'scripts/plot_v23_2_facs_evidence.py'),str(evidence),str(out/'FACS_display_rebuild')],check=True)
for p in (out/'FACS_display_rebuild/figures').iterdir():shutil.copy2(p,dst/p.name)
current=src/'current_v23_2_figure_reproduction'
shutil.copytree(current/'captions',out/'captions')
expected=json.loads((current/'Current_v23_2_graphic_hashes.json').read_text())
assert len(expected)==29
for name,h in expected.items():assert sha(dst/name)==h,name
report=dict(status='PASS',current_selected_image_assets=29,unchanged_existing_assets=23,new_FACS_plates=6,statistical_models_refitted=0,
 scope='Reproduces current figure graphics and display binning, not upstream FASTQ processing, biological specificity or quantitative DAPI validation.')
(out/'Current_v23_2_reproduction_result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
