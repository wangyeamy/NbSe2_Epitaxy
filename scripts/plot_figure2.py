from pathlib import Path
import sys,csv,json,hashlib
P=Path(__file__).resolve().parents[1]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
O=P/'results';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':11,'axes.linewidth':.8,'axes.spines.top':False,'axes.spines.right':False,'xtick.labelsize':9,'ytick.labelsize':9,'svg.fonttype':'none','pdf.fonttype':42})
blue='#287DA1';red='#B65C4C';ink='#26353E';gold='#D2A74E'
def save(f,name):
 for ext in ['png','svg','pdf']:f.savefig(O/(name+'.'+ext),dpi=400,facecolor='white')
 plt.close(f)
def panel(ax,s,title):
 ax.text(-.17,1.10,s,transform=ax.transAxes,fontweight='bold',fontsize=14)
 ax.set_title(title,loc='left',pad=15)
src=P/'data/forces/atom_forces.csv'
rows=list(csv.DictReader(src.open(encoding='utf-8-sig')))
keys=['x_angstrom','y_angstrom','z_angstrom','dfx_ev_angstrom','dfy_ev_angstrom','dfz_ev_angstrom']
for r in rows:
 for k in keys:r[k]=float(r[k])
groups=[[r for r in rows if r['layer']==s] for s in ['Se_top','Nb','Se_bottom']]
assert [len(g) for g in groups]==[16]*3
means=[np.mean([r['dfz_ev_angstrom'] for r in g]) for g in groups]
sums=np.array([sum(r['dfz_ev_angstrom'] for r in g) for g in groups]);cancel=1-abs(sums.sum())/abs(sums).sum()
cmap=LinearSegmentedColormap.from_list('force',['#216D9B','#FAFAF7','#BC644E']);norm=TwoSlopeNorm(0,vmin=-.09,vmax=.09)
f=plt.figure(figsize=(8.4,6.6))
a=f.add_axes([.13,.62,.36,.26]);panel(a,'a','Atom-resolved normal force')
a.axvspan(-.1,0,color=blue,alpha=.065);a.axvspan(0,.1,color=red,alpha=.065)
for j,g in enumerate(groups):
 v=np.array([r['dfz_ev_angstrom'] for r in g]);off=np.tile([-.105,-.035,.035,.105],4)
 a.scatter(v,j+off,c=v,cmap=cmap,norm=norm,s=37,edgecolors='white',lw=.5,zorder=3)
 a.plot([v.mean()]*2,[j-.23,j+.23],c=ink,lw=1.3)
 a.text(.098,j+.30,f'{v.mean():+.3f}',ha='right',va='center',fontsize=9)
a.axvline(0,color='#A2A9AE',lw=.7);a.set(xlim=(-.10,.10),ylim=(2.5,-.5),xlabel=r'$\Delta F_z$ (eV Å$^{-1}$)');a.set_yticks(range(3),['Upper Se','Nb','Lower Se']);a.set_xticks([-.1,0,.1]);a.spines['left'].set_visible(False);a.tick_params(axis='y',length=0)
b=f.add_axes([.66,.62,.30,.26]);panel(b,'b','Cancellation between planes')
v=np.r_[sums,sums.sum()];b.barh(range(4),v,color=[blue,red,blue,ink],height=.53)
for j,x in enumerate(v):b.text(x+(.045 if x>0 else -.045),j,f'{x:+.3f}',ha='left' if x>0 else 'right',va='center',fontsize=9)
b.axvline(0,color='#A2A9AE',lw=.7);b.set(xlim=(-1.2,1.58),ylim=(3.65,-.6),xlabel=r'$\sum_i\Delta F_{i,z}$ (eV Å$^{-1}$)');b.set_yticks(range(4),['Upper Se','Nb','Lower Se','Net']);b.spines['left'].set_visible(False);b.tick_params(axis='y',length=0)
for j,g in enumerate(groups):
 a=f.add_axes([.065+j*.32,.18,.285,.30]);panel(a,'cde'[j],['Upper Se','Nb','Lower Se'][j])
 x,y,fx,fy,fz=[np.array([r[k] for r in g]) for k in ['x_angstrom','y_angstrom','dfx_ev_angstrom','dfy_ev_angstrom','dfz_ev_angstrom']]
 sc=a.scatter(x,y,c=fz,cmap=cmap,norm=norm,s=160,edgecolors='#7A858B',lw=.5,zorder=3)
 q=a.quiver(x,y,fx,fy,angles='xy',scale_units='xy',scale=.012,width=.007,color=ink,zorder=4)
 a.set(xlim=(-7,15),ylim=(-2.5,14));a.set_aspect('equal');a.axis('off');a.plot([-5,0],[-1.4,-1.4],c=ink,lw=1.5);a.text(-2.5,-2.1,'5 Å',ha='center',va='top',fontsize=9)
 if j==2:a.quiverkey(q,.67,.073,.02,r'$\Delta F_{xy}$: 0.02 eV Å$^{-1}$',coordinates='figure',labelpos='E',fontproperties={'size':9})
cax=f.add_axes([.15,.066,.30,.017]);cb=f.colorbar(sc,cax=cax,orientation='horizontal',ticks=[-.09,0,.09]);cb.outline.set_visible(False);cb.set_label(r'$\Delta F_z$ (eV Å$^{-1}$)',labelpad=3)
save(f,'Figure2')
