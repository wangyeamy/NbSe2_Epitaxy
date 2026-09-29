from pathlib import Path
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1];O=P/'results'
means=json.loads((O/'analysis.json').read_text())['high_moment_plane_means']
keys=['Nb','upper_Se','lower_Se']
ys=np.arange(3)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8})
fig,axs=plt.subplots(1,2,figsize=(7.2,3.25),gridspec_kw={'width_ratios':[1.1,1]})
colors=['#287CA2','#BD6551'];ys=np.arange(3)
for j,g in enumerate(['GU','GC']):
 vals=[means[g][k] for k in keys];axs[0].scatter(vals,ys+(-.09 if j==0 else .09),s=45,marker='o' if j==0 else 's',color=colors[j],label=r'$G_'+('U' if j==0 else 'C')+'$',zorder=3)
axs[0].axvline(0,color='.65',lw=.8);axs[0].set(yticks=ys,yticklabels=['Nb','Upper Se','Lower Se'],xlabel=r'Mean $\Delta F_z$ (eV Å$^{-1}$)',xlim=(-.06,.08));axs[0].invert_yaxis();axs[0].legend(frameon=False,loc='lower right');axs[0].set_title('Sublattice response',loc='left',fontsize=11,pad=16)
vals=[.340176795851,-.043773049785,-.007311493738]
axs[1].barh(ys,vals,color=['#BD6551','#287CA2','#287CA2'],height=.45);axs[1].set(yticks=ys,yticklabels=['Nb–Se internal','Rigid z','Se breathing'],xlabel=r'$f_{\mathrm{int}}$ (eV Å$^{-1}$)',xlim=(-.10,.43));axs[1].invert_yaxis();axs[1].axvline(0,color='.65',lw=.8);axs[1].set_title(r'Displacement coupling at $G_C$',loc='left',fontsize=11,pad=16)
for y,v in enumerate(vals):axs[1].text(v-.018 if v>0 else .025,y,f'{v:+.3f}',va='center',ha='right' if v>0 else 'left',fontsize=10,color='white' if v>0 else '#222222')
for a,l in zip(axs,['A','B']):a.text(-.20,1.10,l,transform=a.transAxes,fontweight='bold',fontsize=13)
fig.subplots_adjust(left=.14,right=.98,bottom=.22,top=.82,wspace=.95)

for ext in ['png','svg','pdf']:fig.savefig(O/('FigureS1.'+ext),dpi=400,facecolor='white',bbox_inches='tight')
plt.close(fig)
