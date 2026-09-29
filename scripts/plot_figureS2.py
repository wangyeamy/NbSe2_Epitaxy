from pathlib import Path
import sys,json,csv
P=Path(__file__).resolve().parents[1]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
S=json.loads((P/'data/high_moment/four_state.json').read_text())
D=json.loads((P/'results/endpoint_work_diagnostic.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
f,axs=plt.subplots(2,2,figsize=(8.4,6.7));f.subplots_adjust(left=.18,right=.97,bottom=.12,top=.91,wspace=.60,hspace=.63)
blue='#287DA1';red='#B65C4C';ink='#26353E'
for a,l,t in zip(axs.flat,'abcd',['Four-state energy comparison','Force projections at two geometries','Geometric displacement projections','Endpoint work diagnostic']):
 a.text(-.2,1.12,l,transform=a.transAxes,fontsize=14,fontweight='bold');a.set_title(t,loc='left',pad=17)
a=axs[0,0]
for col,lab,delta in [(blue,'Unconfined',S['energy']['delta_U_GC_minus_GU_ev']),(red,'Confined',S['energy']['delta_C_GC_minus_GU_ev'])]:a.plot([0,1],[0,1000*delta],'o--',c=col,lw=1.2,ms=7,label=lab)
a.set_xticks([0,1],[r'$G_U$',r'$G_C$']);a.set(xlim=(-.13,1.2),ylim=(-4.3,1.2),ylabel=r'$E(G)-E(G_U)$ (meV cell$^{-1}$)');a.legend(frameon=False,fontsize=9,loc='lower left');a.axhline(0,color='#BBBBBB',lw=.6)
keys=['buckling','breathing','rigid_z','shear_x','shear_y'];labels=['Nb–Se internal','Se breathing','Rigid z','Shear x','Shear y'];y=np.arange(5)
a=axs[0,1]
for off,g,col,marker in [(-.12,'GU',blue,'o'),(.12,'GC',red,'s')]:
 vals=[S['forces']['capping_on_substrate_'+g]['modes'][k+'_ev_a_unit_mode'] for k in keys]
 a.scatter(vals,y+off,c=col,marker=marker,s=37,label='$G_'+g[-1]+'$')
a.axvline(0,c='#BBBBBB',lw=.7);a.set_yticks(y,labels);a.set(ylim=(4.6,-.6),xlim=(-.075,.40),xlabel=r'$f_{\mathrm{int},m}$ (eV Å$^{-1}$)');a.legend(frameon=False,fontsize=9,loc='lower right')
a=axs[1,0];vals=[S['GU_to_GC_displacement_modes'][k+'_angstrom_unit_mode']*1000 for k in keys];a.barh(y,vals,color=[red if x>0 else blue for x in vals],height=.55);a.set_yticks(y,labels);a.axvline(0,c='#BBBBBB',lw=.7);a.set(ylim=(4.6,-.6),xlabel=r'$\Delta q_m$ (10$^{-3}$ Å)',xlim=(-12,14))
a=axs[1,1];vals=[D['observed_interaction_energy_change_meV_cell'],D['projected_sum_meV_cell'],D['observed_minus_projected_sum_meV_cell']];a.barh(range(3),vals,color=[ink,blue,red],height=.5)
for j,v in enumerate(vals):a.text(v+(.08 if v>0 else -.08),j,f'{v:+.3f}',ha='left' if v>0 else 'right',va='center',fontsize=9)
a.set_yticks(range(3),['Direct ΔI','Five-mode estimate','Residual']);a.axvline(0,c='#BBBBBB',lw=.7);a.set(ylim=(2.6,-.6),xlim=(-5.4,1.7),xlabel=r'Energy change (meV cell$^{-1}$)')
for ext in ['png','svg','pdf']:f.savefig(P/'results'/('FigureS2.'+ext),dpi=400,facecolor='white')
(P/'results/four_state_summary.json').write_text(json.dumps(S,indent=2));(P/'results/endpoint_work_diagnostic.json').write_text(json.dumps(D,indent=2))
print('Four-state figure complete')
