"""Render original operator PDFs at170dpi; no gate/data manipulation."""
from pathlib import Path
import subprocess
r=Path(__file__).resolve().parent
for letter,name in zip('ABC',['BSA','PBS1','PBS2']):
 subprocess.run(['pdftoppm','-r','170','-png','-singlefile',str(r/(name+'.pdf')),str(r/('Supplementary_Figure_S18_'+letter))],check=True)
