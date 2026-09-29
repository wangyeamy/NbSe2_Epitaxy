from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results';p.mkdir(exist_ok=True)
src=ROOT/'data/inputs/reference_C.in';lines=src.read_text().splitlines();start=lines.index('ATOMIC_POSITIONS angstrom')+1
atoms=[]
for line in lines[start:]:
 t=line.split()
 if len(t)<4:break
 if t[0] in ('Nb','Se'):atoms.append((t[0],np.array(list(map(float,t[1:4])))))
cellstart=next(i for i,x in enumerate(lines) if x.startswith('CELL_PARAMETERS'))
cell=np.array([list(map(float,x.split())) for x in lines[cellstart+1:cellstart+4]])
nb=[r for s,r in atoms if s=='Nb'];se=[r for s,r in atoms if s=='Se']
# Periodic Nb-Se coordination check, independent of the drawing crop.
pse=np.array([r+i*cell[0]+j*cell[1] for i in (-1,0,1) for j in (-1,0,1) for r in se])
counts=[int(np.sum(np.linalg.norm(pse-r,axis=1)<2.95)) for r in nb]
assert counts==[6]*16,counts
centers=[min(nb,key=lambda r:np.linalg.norm(r[:2]-target)) for target in [np.array([2.38,1.38]),np.array([5.93,1.38]),np.array([9.52,1.39])]]
ids=sorted(set(int(k) for r in centers for k in np.where(np.linalg.norm(pse-r,axis=1)<2.95)[0]))
pts=[('Nb',r) for r in centers]+[('Se',pse[k]) for k in ids]
# Oblique view of actual positions, no changes to atom positions or coordination.
def proj(r):return np.array([r[0]+.32*r[1],r[2]+.32*r[1]])
uv=np.array([proj(r) for s,r in pts]);uv-=uv.mean(axis=0)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'axes.titlesize':12,'axes.linewidth':.85,'xtick.labelsize':10,'ytick.labelsize':10,'svg.fonttype':'none','pdf.fonttype':42})
fig,axes=plt.subplots(1,3,figsize=(11.5,3.65),squeeze=False);fig.subplots_adjust(left=.035,right=.99,bottom=.20,top=.83,wspace=.46)
blue='#267F9F';gold='#D5A540';red='#B95848';dark='#253740'
for ax,l,title in zip(axes.flat,'abc',[r'NbSe$_2$ displacement','Projected force','Energy response']):
 ax.set_title(title,loc='left',pad=13);ax.text(-.12,1.075,l,transform=ax.transAxes,fontweight='bold',fontsize=15)
 ax.spines[['top','right']].set_visible(False)
a=axes[0,0];a.axis('off');a.set_aspect('equal',adjustable='datalim')
for i,(s,r) in enumerate(pts):
 if s!='Nb':continue
 for j,(t,v) in enumerate(pts):
  if t=='Se' and np.linalg.norm(r-v)<2.95:a.plot(*np.array([uv[i],uv[j]]).T,color='#A7AFB3',lw=2.2,zorder=1)
for i in sorted(range(len(pts)),key=lambda i:pts[i][1][1],reverse=True):
 s,r=pts[i];a.add_patch(Circle(uv[i],.31 if s=='Nb' else .26,facecolor=blue if s=='Nb' else gold,edgecolor='white',lw=.65,zorder=3+(12-r[1])*.1))
# Arrows denote collective plane displacements for negative q, not forces.
arrow_x=uv[:,0].max()+1.0
for species,upper,label,col in [('Se',True,'Se',gold),('Nb',None,'Nb',blue),('Se',False,'Se',gold)]:
 indices=[i for i,(el,r) in enumerate(pts) if el==species and (upper is None or (r[2]>20)==upper)]
 y=float(uv[indices,1].mean())
 dy=-1.4 if species=='Nb' else .70
 a.annotate('',xy=(arrow_x,y+dy/2),xytext=(arrow_x,y-dy/2),arrowprops=dict(arrowstyle='-|>',color=col,lw=1.8,mutation_scale=8,shrinkA=0,shrinkB=0),zorder=9)
 a.text(arrow_x+.32,y,label,va='center',fontsize=11,color=dark)
a.text(arrow_x,uv[:,1].max()+.65,r'$q<0$',ha='center',fontsize=11)
a.set_xlim(uv[:,0].min()-.65,arrow_x+1.15);a.set_ylim(uv[:,1].min()-.8,uv[:,1].max()+1.2)
v=np.genfromtxt(ROOT/'data/high_moment/displacement.csv',delimiter=',',names=True)
q=v['q_A'];fc=v['F_confined_eV_A'];fu=v['F_unconfined_eV_A'];fi=fc-fu
b=axes[0,1]
for y,c,m,label,lw in [(fc,red,'o','Confined',.85),(fu,blue,'s','Unconfined',.85),(fi,dark,'D',r'$f_{\mathrm{int}}$',1.2)]:b.plot(q,y,color=c,marker=m,ms=5.5,lw=lw,ls='--',label=label)
b.axhline(0,color='#CFD4D7',lw=.65,zorder=0);b.set(xlim=(-.023,.003),ylim=(-.35,.65),xticks=q,xticklabels=['−0.02','0'],xlabel=r'$q$ (Å)',ylabel=r'Force (eV Å$^{-1}$)');b.legend(frameon=False,fontsize=9.5,loc='lower left',handlelength=1.5,labelspacing=.3)
d=axes[0,2]
for y,c,m,label in [(v['deltaE_confined_meV'],red,'o','Confined'),(v['deltaE_unconfined_meV'],blue,'s','Unconfined')]:d.plot(q,y,color=c,marker=m,ms=5.5,lw=.85,ls='--',label=label)
d.axhline(0,color='#CFD4D7',lw=.65,zorder=0);d.set(xlim=(-.023,.003),ylim=(-1.2,7.5),xticks=q,xticklabels=['−0.02','0'],xlabel=r'$q$ (Å)',ylabel=r'$E(q)-E(0)$ (meV cell$^{-1}$)');d.legend(frameon=False,fontsize=9.5,loc='upper right',handlelength=1.5)
for ext in ['png','svg','pdf']:fig.savefig(p/('Figure3.'+ext),dpi=450,facecolor='white')
(p/'provenance.json').write_text(json.dumps({'structure_file':'data/inputs/reference_C.in','sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'Nb_Se_neighbor_cutoff_A':2.95,'coordination_all_Nb':counts,'displayed_atoms':[{'element':s,'xyz_A':r.tolist()} for s,r in pts],'q_endpoints_A':q.tolist(),'confined_projected_force_eV_A':fc.tolist(),'unconfined_projected_force_eV_A':fu.tolist(),'energy_change_meV_cell':6.712504750751146,'negative_work_meV_cell':6.70926934183236},indent=2))
print('Coordination',counts,'selected atoms',len(pts))



