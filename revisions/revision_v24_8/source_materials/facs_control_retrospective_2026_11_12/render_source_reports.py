from pathlib import Path
import subprocess
r=Path(__file__).resolve().parent
for letter,name in zip("ABC",["BSA","PBS1","PBS2"]):
 subprocess.run(["pdftoppm","-r","170","-png","-singlefile",str(r/"source_reports"/(name+".pdf")),str(r/"figures"/("Figure_S20"+letter+"_HUT093_"+name+"_retrospective"))],check=True)
