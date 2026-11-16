from pathlib import Path
import argparse,sys,subprocess,shutil,hashlib,json
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=Path(__file__).resolve().parent;s=r/'source_materials';out=a.output.resolve();assert not out.exists(),'Use a new external output directory';out.mkdir(parents=True)
subprocess.run([sys.executable,str(s/'current_v24_figure_reproduction/scripts/reproduce_v24_figures.py'),str(s),str(out/'base')],check=True)
shutil.copytree(out/'base/figures',out/'figures')
for letter,name in zip('ABC',['BSA','PBS1','PBS2']):
 subprocess.run(['pdftoppm','-r','170','-png','-singlefile',str(s/'current_v24_4_control_figures'/(name+'.pdf')),str(out/'figures'/('Supplementary_Figure_S18_'+letter))],check=True)
for f in (r/'selected_figures').glob('*caption.txt'):(out/'figures'/f.name).write_bytes(f.read_bytes())
expected=json.loads((r/'Selected_figure_SHA256.json').read_text())
for name,want in expected.items():assert hashlib.sha256((out/'figures'/name).read_bytes()).hexdigest()==want,name
print('PASS: all',len(expected),'selected current graphics and captions reproduced or source-preserved. No raw-read reprocessing or new model fit.')
