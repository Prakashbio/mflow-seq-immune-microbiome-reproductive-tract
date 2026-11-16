"""Replay existing review checks from an extracted release, with external outputs.

No classifier, simulation generator, new model fit or online action is executed.
Requires the Python/R dependencies recorded in the source environment records.
"""
from pathlib import Path
import argparse,json,shutil,subprocess,sys
R=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True,help='New directory outside the extracted repository')
p.add_argument('--figures',action='store_true',help='Also replay the existing figure scripts; requires Rscript and Poppler')
a=p.parse_args();out=a.output.resolve();root=R.resolve()
assert out!=root and root not in out.parents,'Use an output directory outside the repository'
assert not out.exists(),'Use a new output directory; existing results are never overwritten'
out.mkdir(parents=True);logs=out/'logs';logs.mkdir();base=R/'revisions/revision_v24_8';s=base/'source_materials';steps=[]
def run(name,args,cwd=R):
 log=logs/(name+'.log')
 with log.open('w') as f:r=subprocess.run([sys.executable,*map(str,args)],cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
 steps.append(dict(step=name,exit_code=r.returncode,log=str(log.relative_to(out))))
 assert r.returncode==0, f'{name} failed; see {log}'
run('01_release_integrity',[R/'verify_release_v25_2_2.py'])
run('02_baseline_bundle',[s/'baseline_v20_1_analysis_release/verify_bundle.py'])
run('03_saved_Q15',[s/'baseline_v20_1_analysis_release/analysis_scripts/verify_Q15_joint_factorial_2026_11_05.py',s/'baseline_v20_1_analysis_release/bundle/joint_factorial_q15'])
run('04_reviewer_metadata',[s/'reviewer_metadata_v24_6/verify_reviewer_metadata_v24_6.py',s/'reviewer_metadata_v24_6'])
run('05_score_bridge',[s/'score_provenance_bridge_v24_1/scripts/reproduce_score_bridge_v24_1.py',s,out/'score_bridge'])
own=s/'own_PRE_between_sites_v24_5';work=out/'own_PRE';work.mkdir();shutil.copytree(own/'inputs',work/'inputs')
run('06_own_PRE_replay',[own/'own_PRE_between_sites_v24_5.py','--inputs',work/'inputs','--output',work/'results'])
run('07_own_PRE_independent_check',[own/'verify_own_PRE_between_sites_v24_5.py','--workdir',work])
seq=R/'revisions/revision_v25'
run('08_sequencing_summary',[seq/'scripts/recompute_v25_sequencing_summaries.py','--input',seq/'data','--output',out/'sequencing_summary'])
if a.figures:
 assert shutil.which('Rscript') and shutil.which('pdftoppm'),'Figure replay requires Rscript and pdftoppm on PATH'
 run('09_figure_replay',[base/'reproduce_current_figures.py','--output',out/'figure_replay'])
run('10_source_integrity_after_replay',[R/'verify_release_v25_2_2.py'])
report=dict(status='PASS',steps=steps,figures_requested=a.figures,source_inputs_unchanged=True,
 scope='Replay of existing downstream checks and figures. This does not reproduce the proprietary upstream pipeline or establish experimental validity or exchangeability.')
(out/'Replay_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
