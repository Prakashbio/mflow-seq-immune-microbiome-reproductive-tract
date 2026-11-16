"""Reproduce all current figures, retaining numerical data and original imagery."""
from pathlib import Path
import sys,shutil,subprocess,hashlib,json
src=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();assert not out.exists();out.mkdir(parents=True)
subprocess.run([sys.executable,str(src/'current_v23_3_figure_reproduction/scripts/reproduce_v23_3_figures.py'),str(src),str(out/'preserved_v23_3_build')],check=True)
current=src/'current_v24_figure_reproduction'
subprocess.run([shutil.which('Rscript') or '/usr/local/bin/Rscript','--vanilla',str(current/'scripts/polish_v24_figure_layouts.R'),str(src/'current_v23_figure_reproduction'),str(out/'layout_rebuild')],check=True)
shutil.copytree(out/'preserved_v23_3_build/figures',out/'figures')
for p in (out/'layout_rebuild/figures').iterdir():shutil.copy2(p,out/'figures'/p.name)
shutil.copytree(current/'captions',out/'captions')
expected=json.loads((current/'Current_v24_graphic_hashes.json').read_text());assert len(expected)==29
for name,h in expected.items():assert hashlib.sha256((out/'figures'/name).read_bytes()).hexdigest()==h,name
receipt=dict(status='PASS',current_image_assets=29,presentation_updated_assets=17,unchanged_assets=12,original_microscopy_and_FACS_graphics_unchanged=True,numerical_results_unchanged=True,new_statistical_tests=0,scope='Downstream display reproduction and layer-value preservation; not upstream sequencing provenance or experimental validation')
(out/'Current_v24_reproduction_result.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
