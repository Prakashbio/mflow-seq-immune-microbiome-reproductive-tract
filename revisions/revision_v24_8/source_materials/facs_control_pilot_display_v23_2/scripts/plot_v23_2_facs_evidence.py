"""Display existing control/pilot evidence; no original gate reconstruction.

Usage: python3 plot_v23_2_facs_evidence.py SOURCE_FOLDER OUTPUT_FOLDER
SOURCE_FOLDER holds source_fcs/, source_reports/, source_operator_graphics/.
OUTPUT_FOLDER may equal SOURCE_FOLDER. Only derived figures/captions/audit are
written. Source FCS/EMF/PDF bytes are checked before and after plotting.
"""
from pathlib import Path
import csv,hashlib,io,json,struct,sys,subprocess,tempfile,os,shutil
import numpy as np
from PIL import Image
from pypdf import PdfReader
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

SRC=Path(sys.argv[1]).resolve();OUT=Path(sys.argv[2]).resolve()
fontcache=tempfile.TemporaryDirectory(prefix='mflow_facs_fonts_')
childenv=os.environ.copy();childenv['XDG_CACHE_HOME']=fontcache.name
PDFTOPPM=shutil.which('pdftoppm') or '/opt/homebrew/bin/pdftoppm'
RSCRIPT=shutil.which('Rscript') or '/usr/local/bin/Rscript'
for n in ['figures','captions','audit']:(OUT/n).mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={str(p.relative_to(SRC)):sha(p) for folder in ['source_fcs','source_reports','source_operator_graphics'] for p in (SRC/folder).rglob('*') if p.is_file()}
crop_rows=[];plot_rows=[]
def emf_bitmap(p):
 b=p.read_bytes();off=0
 while off+8<=len(b):
  kind,size=struct.unpack_from('<II',b,off)
  if kind==81:
   rr=b[off:off+size];ob,cb,op,cp=struct.unpack_from('<IIII',rr,48);bmi=rr[ob:ob+cb];bits=rr[op:op+cp]
   im=Image.open(io.BytesIO(b'BM'+struct.pack('<IHHI',14+len(bmi)+len(bits),0,0,14+len(bmi))+bmi+bits)).convert('RGB')
   return im
  if size<8:break
  off+=size
 raise ValueError('No source bitmap in '+str(p))
def crop(im,box,source,part,label):
 z=im.crop(box)
 assert np.array_equal(np.asarray(z),np.asarray(im)[box[1]:box[3],box[0]:box[2]])
 crop_rows.append(dict(part=part,panel=label,source=source,box_pixels=list(box),cropped_pixel_SHA256=hashlib.sha256(np.asarray(z).tobytes()).hexdigest(),pixel_values_unchanged=True))
 return z
def source_plate(stem,title,entries):
 pdf=OUT/'figures'/(stem+'.pdf');ww=5*72;hh=7.5*72
 cc=canvas.Canvas(str(pdf),pagesize=(ww,hh),invariant=1)
 cc.setTitle(title);cc.setFont('Helvetica-Bold',12);cc.drawString(.04*ww,.965*hh,title)
 for rect,im,label in entries:
  x,y,w,h=rect;x*=ww;y*=hh;w*=ww;h*=hh
  cc.setFont('Helvetica',9);cc.drawString(x,y+h+4,label)
  cc.drawImage(ImageReader(im),x,y,width=w,height=h,preserveAspectRatio=True,anchor='c',mask='auto')
 cc.save()
 subprocess.run([PDFTOPPM,'-png','-singlefile','-r','300',str(pdf),str(pdf.with_suffix(''))],check=True,env=childenv)

# Original buffer gates: enlarged native DAPI and IgA/IgM graph panels, no
# recolouring, signal adjustment or retrospective gate changes.
reader=PdfReader(SRC/'source_reports/HUT093_original_report.pdf')
buffer_info=[('A','BSA',5,34858,18086,1170),('B','PBS1',6,35423,10066,1448),('C','PBS2',7,36860,9554,2870)]
for letter,tube,page,p1,parent,igm in buffer_info:
 im=reader.pages[page-1].images[0].image.convert('RGB')
 d=crop(im,(615,288,975,575),'HUT093_original_report.pdf page '+str(page),'S18'+letter,'DAPI/SSC')
 a=crop(im,(0,289,310,575),'HUT093_original_report.pdf page '+str(page),'S18'+letter,'IgA/IgM')
 source_plate('Supplementary_Figure_S18_'+letter,letter+'  '+tube+' buffer only',[
  ([.04,.51,.92,.425],d,'DAPI-labelled gate versus SSC'),
  ([.07,.025,.86,.435],a,'IgA/IgM-labelled gates')])
 caption=(f'Supplementary Figure S18{letter}. Original HUT093 {tube} buffer-control gates. '
          f'No sample, antibodies or DAPI were added. The panels reproduce the original DAPI-labelled/SSC and IgA/IgM-labelled graph regions without re-gating; these labels identify software regions, not staining. '
          f'The acquisition contained 50,000 events; P1 contained {p1:,}, and the DAPI-labelled negative parent contained {parent:,}. '
          f'IgM-region occupancy was {igm:,}/{parent:,} ({100*igm/parent:.2f}%). These are buffer/instrument-background observations, not antibody specificity, viability or clinical false-positive rates. Complete scatter, IgG and population-table panels are retained in the source report (page {page}).')
 (OUT/'captions'/f'Supplementary_Figure_S18_{letter}_caption.txt').write_text(caption+'\n')

def readfcs(p):
 b=p.read_bytes();s=int(b[10:18]);e=int(b[18:26]);v=b[s:e+1].decode('latin1').split(chr(b[s]))[1:-1]
 m={v[i]:v[i+1].strip() for i in range(0,len(v),2)};n=int(m['$PAR']);tot=int(m['$TOT'])
 assert m['$DATATYPE']=='F' and int(m['$ENDDATA'])-int(m['$BEGINDATA'])+1==n*tot*4
 aa=np.frombuffer(b,dtype=('>f4' if m['$BYTEORD']=='4,3,2,1' else '<f4'),offset=int(m['$BEGINDATA']),count=n*tot).reshape(tot,n)
 assert np.isfinite(aa).all()
 return m,aa,[m[f'$P{i}N'] for i in range(1,n+1)]
allfcs=[]
for p in (SRC/'source_fcs').rglob('*.fcs'):
 m,aa,names=readfcs(p);allfcs.append((p,m,aa,names))
def one(date,tube):
 hits=[z for z in allfcs if z[1]['$DATE']==date and z[1].get('TUBE NAME')==tube]
 assert len(hits)==1,(date,tube)
 return hits[0]
def hist_grid(records,shape,stem,titles,channels,supertitle):
 c=150.0;xy=[]
 for p,m,a,nn in records:
  xx=a[:,nn.index(channels[0])].astype(float);yy=a[:,nn.index(channels[1])].astype(float)
  xy.append((np.arcsinh(xx/c),np.arcsinh(yy/c)))
 xmin=min(float(x.min()) for x,y in xy);xmax=max(float(x.max()) for x,y in xy)
 ymin=min(float(y.min()) for x,y in xy);ymax=max(float(y.max()) for x,y in xy)
 xe=np.linspace(xmin-1e-8,xmax+1e-8,151);ye=np.linspace(ymin-1e-8,ymax+1e-8,151)
 hh=[np.histogram2d(x,y,bins=(xe,ye))[0] for x,y in xy];vmax=max(float(h.max()) for h in hh)
 display=[]
 for i,(h,rec,title) in enumerate(zip(hh,records,titles)):
  p,m,a,names=rec
  assert int(h.sum())==len(a)
  for ix,iy in zip(*np.nonzero(h)):
   display.append(dict(panel=i+1,title=title,ix=int(ix)+1,iy=int(iy)+1,events=int(h[ix,iy])))
  plot_rows.append(dict(figure=stem,panel=i+1,source=str(p.relative_to(SRC)),tube=m['TUBE NAME'],date=m['$DATE'],events=len(a),events_in_display_bins=int(h.sum()),gate='None; all recorded events',x=channels[0],y=channels[1],transform='arcsinh(value/150), labelled in exported channel units',bins_per_axis=150,additional_compensation_applied=False))
 csvpath=OUT/'audit'/(stem+'_display_bins.csv')
 with csvpath.open('w',newline='') as ff:
  writer=csv.DictWriter(ff,fieldnames=list(display[0]));writer.writeheader();writer.writerows(display)
 xlabel=('IgA-associated PE detector' if '585/42' in channels[0] else 'SYTO9-associated 530/30 detector')+' (exported units)'
 ylabel=('IgM-associated APC detector' if '660/20' in channels[1] else 'PI-associated 616/23 detector')+' (exported units)'
 args=[RSCRIPT,'--vanilla',str(Path(__file__).with_name('render_v23_2_facs_grids.R')),str(csvpath),str(OUT/'figures'/stem),str(shape[0]),str(shape[1]),str(xe[0]),str(xe[-1]),str(ye[0]),str(ye[-1]),str(vmax),xlabel,ylabel,supertitle]
 rr=subprocess.run(args,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=childenv)
 (OUT/'audit'/(stem+'_render.log')).write_text(rr.stdout)
 if rr.returncode:raise RuntimeError(rr.stdout)
 pdf=OUT/'figures'/(stem+'.pdf')
 subprocess.run([PDFTOPPM,'-png','-singlefile','-r','300',str(pdf),str(pdf.with_suffix(''))],check=True,env=childenv)

# All four recorded indices in each of the three June 13 control groups.
records=[];titles=[]
for i in range(1,5):
 for label,tube in [('Unstained',f'Saliva unstained {i}-4'),('IgA only',f'Saliva IgA+PE {i}-2'),('IgM only',f'Saliva IgM+APC {i}-3')]:
  records.append(one('13-JUN-2023',tube));titles.append(f'{label} · acquisition {i}')
assert all(len(z[2])==50000 for z in records)
hist_grid(records,(4,3),'Supplementary_Figure_S19_A',titles,('2-488-585/42-A','1-633-660/20-A'),'A  Saliva development controls · 13 June 2023')
(OUT/'captions/Supplementary_Figure_S19_A_caption.txt').write_text(
 'Supplementary Figure S19A. Archived unstained and single-colour saliva controls. Columns show unstained, IgA-PE-only and IgM-APC-only preparations; rows retain all four recorded acquisition indices from 13 June 2023. Each acquisition contains 50,000 recorded events. These are 12 acquisition records, not 12 independent biological replicates. Each display bins all recorded exported events, without scatter/DAPI/Ig gating or additional compensation. Axes use arcsinh(value/150), labelled in exported channel units; colour gives events per bin on a logarithmic scale. No DAPI was used in this early series. The displays provide development-control context, not matched clinical FMO coverage or antibody-binding specificity.\n')

# Show actual operator-drawn DAPI gates for the first recorded saliva/vaginal
# preparations, rather than inventing thresholds on the exported event values.
entries=[]
for rect,i,title in [([.04,.51,.92,.415],4,'Saliva 1 · 20 June 2023'),([.04,.025,.92,.455],9,'Vaginal 1 · 7 July 2023')]:
 im=emf_bitmap(SRC/'source_operator_graphics'/f'image{i}.emf')
 z=crop(im,(333,0,min(715,im.width),315),f'operator image{i}.emf','S19B',title)
 entries.append((rect,z,title))
source_plate('Supplementary_Figure_S19_B','B  Original DAPI pilot gates',entries)
(OUT/'captions/Supplementary_Figure_S19_B_caption.txt').write_text(
 'Supplementary Figure S19B. Original operator-drawn DAPI pilot regions. The first numbered saliva preparation from 20 June 2023 and first numbered vaginal preparation from 7 July 2023 are shown using cropped original graph bitmaps, without changed gate positions, event colours or signal adjustment. They contain 50,000 and 20,000 recorded events, respectively. DAPI roles follow the operator’s description. The complete records contain four saliva and two vaginal DAPI/antibody acquisitions; acquisition indices do not establish independent biological n. These regions are not a validated live/dead partition or a quantitative comparison with SYTO9/PI. Full source graphics are retained in the reproduction files.\n')

records=[];titles=[]
for i in range(1,5):records.append(one('20-JUN-2023',f'Saliva Syto9+PI {i}-1'));titles.append(f'Saliva {i} · 20 June')
for i in range(1,3):records.append(one('07-JUL-2023',f'Vaginal swap Syto9+PI {i}'));titles.append(f'Vaginal {i} · 7 July')
assert all(len(z[2])==20000 for z in records)
hist_grid(records,(3,2),'Supplementary_Figure_S19_C',titles,('2-488-530/30-A','2-488-616/23-A'),'C  Candidate SYTO9/PI comparator acquisitions · 2023')
(OUT/'captions/Supplementary_Figure_S19_C_caption.txt').write_text(
 'Supplementary Figure S19C. Available SYTO9/PI-labelled candidate comparator acquisitions. Four saliva records from 20 June 2023 and two vaginal records from 7 July 2023 are displayed from their original exported events; each acquisition contains 20,000 events. Date and acquisition-number correspondence to the DAPI pilots does not verify source-aliquot pairing or same-event staining. All recorded events are shown as two-dimensional bin counts using arcsinh(value/150), with axes labelled in exported channel units and logarithmic colour scaling; no additional compensation or bacterial, live/dead or sorting gate is reconstructed. These displays do not supply a validated comparator denominator, concordance statistic or independent validation sample size.\n')

for name,rows in [('Original_graph_crop_provenance.json',crop_rows),('Retrospective_display_records.json',plot_rows)]:
 (OUT/'audit'/name).write_text(json.dumps(rows,indent=2)+'\n')
with (OUT/'audit/Retrospective_display_records.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(plot_rows[0]));w.writeheader();w.writerows(plot_rows)
assert all(sha(SRC/rel)==h for rel,h in sources.items())
assert len(list((OUT/'figures').glob('*.png')))==6
(OUT/'audit/Immutable_FACS_source_hashes.json').write_text(json.dumps(sources,indent=2)+'\n')
result=dict(status='PASS',source_files_unchanged=len(sources),source_FCS_files=len(allfcs),new_supplementary_figures=2,image_plates=6,original_graph_crops=len(crop_rows),retrospective_event_displays=len(plot_rows),all_recorded_events_retained_in_each_display=True,gate_reconstruction=False,new_inferential_tests=0,scope='Presentation of existing evidence only; no quantitative DAPI validation or clinical binding-specificity result is inferred.')
(OUT/'audit/FACS_display_preservation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
fontcache.cleanup()
