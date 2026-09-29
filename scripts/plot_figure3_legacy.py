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
# Displacement dataset: negative endpoint and reference, separate electronic branch.
q=np.array([-.02,0]);fc=np.array([.5276195975552437,.07662925176554906]);fu=np.array([.1968694592232205,-.2635475440856637]);fi=fc-fu
ec=np.array([6.051404121826115,0]);eu=np.array([-.6611006289250305,0]);di=ec-eu
work=.02*fi.mean()*1000
raw=(P/'data/inputs/reference_C.in').read_text()
atoms=[(r.split()[0],*map(float,r.split()[1:4])) for r in raw.split('ATOMIC_POSITIONS angstrom')[1].split('CELL_PARAMETERS')[0].strip().splitlines()]
selected=[r for r in atoms if r[0] in ['Nb','Se'] and 0<=r[1]<=9 and 2<=r[2]<=10]
f=plt.figure(figsize=(8.4,6.8));gs=f.add_gridspec(2,2,left=.12,right=.97,bottom=.10,top=.90,wspace=.40,hspace=.65)
a=f.add_subplot(gs[0,0]);panel(a,'a','Internal Nb–Se displacement')
pos=np.array([[x+.13*(y-6),z] for s,x,y,z in selected])
for i,u in enumerate(selected):
 for j in range(i):
  v=selected[j]
  if {u[0],v[0]}=={'Nb','Se'} and 2.2<np.linalg.norm(np.array(u[1:])-np.array(v[1:]))<2.85:a.plot(pos[[i,j],0],pos[[i,j],1],c='#A1AAAD',lw=2,zorder=1)
for i,(s,x,y,z) in enumerate(selected):
 dz=-.02/np.sqrt(24) if s=='Nb' else .01/np.sqrt(24)
 a.scatter(*pos[i],s=150 if s=='Nb' else 115,c=blue if s=='Nb' else gold,edgecolors='white',linewidths=.5,zorder=3)
 a.annotate('',xy=pos[i]+[0,dz*250],xytext=pos[i],arrowprops=dict(arrowstyle='-|>',color=ink,lw=1.2),zorder=4)
a.set_aspect('equal');a.set_ylim(pos[:,1].min()-1.9,pos[:,1].max()+1.4);a.axis('off')
a.text(0,-.04,'Nb: −0.00408 Å    Se: +0.00204 Å',transform=a.transAxes,fontsize=9)
a.text(0,-.16,r'$q=-0.02$ Å; displacement arrows ×250',transform=a.transAxes,fontsize=9)
for s,col in [('Nb',blue),('Se',gold)]:a.scatter([],[],s=65,c=col,label=s)
a.legend(loc='upper right',frameon=False,ncol=2,handletextpad=.3,columnspacing=.8)
b=f.add_subplot(gs[0,1]);panel(b,'b','Projection onto the displacement')
for vals,col,marker,label in [(fc,red,'o','Confined'),(fu,blue,'s','Unconfined'),(fi,ink,'D',r'$f_{\mathrm{int}}$')]:b.plot(q,vals,marker=marker,color=col,ls='--',lw=1,ms=6,label=label)
b.axhline(0,color='#C5CBCF',lw=.7);b.set(xlabel=r'$q$ (Å)',ylabel=r'Projected force (eV Å$^{-1}$)',xlim=(-.022,.002),ylim=(-.36,.65));b.set_xticks([-.02,0]);b.legend(frameon=False,fontsize=9,loc='lower left')
c=f.add_subplot(gs[1,1]);panel(c,'d','Confined / unconfined')
for vals,col,m,label in [(ec,red,'o','Confined'),(eu,blue,'s','Unconfined')]:c.plot(q,vals,marker=m,color=col,ls='--',lw=1,ms=6,label=label)
c.axhline(0,color='#C5CBCF',lw=.7);c.set(xlabel=r'$q$ (Å)',ylabel=r'$E(q)-E(0)$ (meV cell$^{-1}$)',xlim=(-.022,.002),ylim=(-1.5,8));c.set_xticks([-.02,0]);c.legend(frameon=False,fontsize=9,loc='upper right')
d=f.add_subplot(gs[1,0]);panel(d,'c','Interfacial contribution to the energy')
qq=np.linspace(-.02,0,100);slope=(fi[1]-fi[0])/.02;ii=-(fi[1]*qq+.5*slope*qq**2)*1000
d.plot(qq,ii,c=ink,lw=1.5,label=r'$-\int_0^q f_{\mathrm{int}}\,dq^\prime$');d.scatter(q,di,s=65,facecolors='white',edgecolors=red,lw=1.5,zorder=4,label='DFT energy difference')
d.set(xlabel=r'$q$ (Å)',ylabel=r'$\Delta I$ (meV cell$^{-1}$)',xlim=(-.022,.002),ylim=(-.6,8.5));d.set_xticks([-.02,0]);d.legend(frameon=False,fontsize=9,loc='upper right');d.text(.98,.51,'At −0.02 Å\nDFT: 6.713 meV\nForce integral: 6.709 meV',transform=d.transAxes,ha='right',fontsize=9)
save(f,'Figure3_legacy_four_panel')
