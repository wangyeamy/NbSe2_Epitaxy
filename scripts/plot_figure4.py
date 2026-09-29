from pathlib import Path
import sys,re,json,csv,hashlib
P=Path(__file__).resolve().parents[1]
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.linewidth':.9,'axes.spines.top':False,'axes.spines.right':False,'xtick.labelsize':9,'ytick.labelsize':9,'svg.fonttype':'none','pdf.fonttype':42})
D=P/'data/dos';info={};data={}
for case in ['stack','uncapped','hbn']:
 path=D/(case+'.dat');header=path.open().readline();ef=float(re.search(r'EFermi\s*=\s*([-+0-9.]+)',header).group(1));v=np.loadtxt(path);assert v.shape[1]==4 and np.isfinite(v).all() and np.all(np.diff(v[:,0])>0)
 x=v[:,0]-ef;y=v[:,1]+v[:,2];assert y.min()>=0;data[case]=(x,y);i=np.argmin(abs(x));info[case]={'fermi_eV':ef,'rows':len(x),'nearest_E_minus_EF_eV':float(x[i]),'DOS_nearest_EF':float(y[i]),'DOS_interpolated_EF':float(np.interp(0,x,y)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
 with (P/'results'/f'{case}_aligned_dos.csv').open('w',newline='') as h:
  w=csv.writer(h);w.writerow(['E_minus_own_EF_eV','DOS_spin_sum_states_per_eV_cell']);w.writerows(zip(x,y))
f,axs=plt.subplots(1,2,figsize=(7.2,3.4),gridspec_kw={'width_ratios':[1.35,1]});f.subplots_adjust(left=.105,right=.98,bottom=.22,top=.85,wspace=.40)
for case,color,label in [('stack','#B65C4C','Confined'),('uncapped','#287DA1','Unconfined')]:
 x,y=data[case];axs[0].plot(x,y,c=color,lw=1.9,label=label)
axs[0].set(xlim=(-1.5,1.5),ylim=(0,max(y[(x>=-1.5)&(x<=1.5)].max() for x,y in [data['stack'],data['uncapped']])*1.16),ylabel=r'DOS (states eV$^{-1}$ cell$^{-1}$)');axs[0].legend(frameon=False,fontsize=9,loc='upper right');axs[0].set_xticks([-1.5,0,1.5]);axs[0].set_facecolor('#F5F4FA')
x,y=data['hbn'];axs[1].plot(x,y,c='#89779F',lw=1.9);mask=abs(x)<=3.5;axs[1].set(xlim=(-3.5,3.5),ylim=(0,y[mask].max()*1.12));axs[1].set_xticks([-3,0,3]);axs[1].set_facecolor('#FCF3ED')
for ax in axs:ax.axvline(0,c='#8B959E',lw=.8,ls='--');ax.set_xlabel(r'$E-E_F$ (eV)')
f.text(.025,.93,'a',fontsize=13,weight='bold');f.text(.105,.93,'Confined / unconfined',fontsize=10)
f.text(.585,.93,'b',fontsize=13,weight='bold');f.text(.665,.93,'Isolated strained hBN',fontsize=10)
f.text(.105,.045,'Spin-summed DOS; each spectrum referenced to its own Fermi energy',fontsize=8.5)
for ext in ['png','svg','pdf']:f.savefig(P/'results'/f'Figure4.{ext}',dpi=400,facecolor='white')
info['relative_change_nearest_percent']=100*(info['stack']['DOS_nearest_EF']/info['uncapped']['DOS_nearest_EF']-1)
info['relative_change_interpolated_percent']=100*(info['stack']['DOS_interpolated_EF']/info['uncapped']['DOS_interpolated_EF']-1)
info['branch']='legacy dense GC; not current high-moment';info['energy_alignment']='each own EF; not absolute band alignment'
(P/'results'/'dos_audit.json').write_text(json.dumps(info,indent=2),encoding='utf-8');print(json.dumps(info,indent=2))
