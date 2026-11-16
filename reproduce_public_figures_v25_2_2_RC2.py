"""Rebuild public figure displays; no private raw human FCS read or new model fit."""
from pathlib import Path
import argparse,shutil,subprocess,sys,json,hashlib
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();R=Path(__file__).resolve().parent;out=a.output.resolve();assert not out.exists() and R not in out.parents,'Use a new external output';out.mkdir(parents=True)
source=R/'revisions/revision_v24_8';work=out/'temporary_public_source';shutil.copytree(source,work)
driver=work/'source_materials/current_v23_2_figure_reproduction/scripts/reproduce_v23_2_figures.py';s=driver.read_text();old="evidence/'scripts/plot_v23_2_facs_evidence.py'";assert s.count(old)==1;s=s.replace(old,"evidence/'scripts/plot_v23_2_facs_evidence_from_public_summaries_RC2.py'");driver.write_text(s)
subprocess.run([sys.executable,str(work/'reproduce_current_figures.py'),'--output',str(out/'replay')],check=True)
shutil.copytree(out/'replay/figures',out/'figures');expected=json.loads((source/'Selected_figure_SHA256.json').read_text());assert len(expected)==58
for n,h in expected.items():assert hashlib.sha256((out/'figures'/n).read_bytes()).hexdigest()==h,n
report=dict(status='PASS',graphics_and_captions=58,raw_human_FCS_replayed=False,withheld_human_FCS=34,scope='All selected displays match via saved public bins/native artwork plus existing downstream figures. Private raw-pilot re-binning is not publicly reproducible from this archive.')
(out/'Public_figure_replay_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
