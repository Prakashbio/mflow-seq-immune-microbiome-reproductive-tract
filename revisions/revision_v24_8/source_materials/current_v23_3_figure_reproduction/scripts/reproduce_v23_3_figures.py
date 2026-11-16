"""Reproduce v23.3 graphics from the preserved baseline and Figure 2 layout script."""
from pathlib import Path
import sys,subprocess,shutil,hashlib,json
src=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();assert not out.exists();out.mkdir(parents=True)
subprocess.run([sys.executable,str(src/'current_v23_2_figure_reproduction/scripts/reproduce_v23_2_figures.py'),str(src),str(out/'preserved_v23_2_build')],check=True)
current=src/'current_v23_3_figure_reproduction'
subprocess.run([shutil.which('Rscript') or '/usr/local/bin/Rscript','--vanilla',str(current/'scripts/plot_figure2_observed_PRE_v23_3.R'),str(src/'current_v23_figure_reproduction'),str(out/'Figure2_observed_profiles_build')],check=True)
subprocess.run([shutil.which('Rscript') or '/usr/local/bin/Rscript','--vanilla',str(current/'scripts/plot_figure5_panelA_legends_v23_3.R'),str(src/'current_v23_figure_reproduction'),str(out/'Figure5A_horizontal_legend_build')],check=True)
shutil.copytree(out/'preserved_v23_2_build/figures',out/'figures')
for p in (out/'Figure2_observed_profiles_build/figures').iterdir():shutil.copy2(p,out/'figures'/p.name)
for p in (out/'Figure5A_horizontal_legend_build/figures').iterdir():shutil.copy2(p,out/'figures'/p.name)
shutil.copytree(current/'captions',out/'captions')
expected=json.loads((current/'Current_v23_3_graphic_hashes.json').read_text());assert len(expected)==29
for name,h in expected.items():assert hashlib.sha256((out/'figures'/name).read_bytes()).hexdigest()==h,name
receipt=dict(status='PASS',reproduced_graphic_assets=29,unchanged_from_v23_2=27,Figure2='Presentation changed; all observed stack segments and panel B values retained',Figure5='A legend entries made horizontal; plotted data and panel B styling retained',numerical_analysis_changes=0,scope='Downstream figure reproduction, not upstream sequencing or assay validation')
(out/'Current_v23_3_reproduction_result.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
