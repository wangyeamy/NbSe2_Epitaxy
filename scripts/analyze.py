"""Recompute manuscript quantities from distributed data without QE or network access."""
from pathlib import Path
import csv,json,re,math
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
RY_EV=13.605693122994;BOHR_A=0.529177210903

def parse_qe(path,nat):
 text=path.read_text(errors='replace')
 if 'JOB DONE.' not in text or 'convergence has been achieved' not in text:raise ValueError('Incomplete calculation: '+path.name)
 blocks=[]
 for part in text.split('Forces acting on atoms (cartesian axes, Ry/au):')[1:]:
  rows=re.findall(r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)',part.split('Total force')[0])
  if len(rows)>=nat:blocks.append(np.array(rows[:nat],float))
 if not blocks:raise ValueError('Missing complete force block')
 energy=float(re.findall(r'!\s+total energy\s+=\s*([-\d.]+)',text)[-1])
 return energy,blocks[-1]

def metrics(f):return {'rms_eV_A':float(np.sqrt(np.mean(np.sum(f*f,axis=1)))),'max_eV_A':float(np.linalg.norm(f,axis=1).max())}
report={};delta={};mesh={}
for k in ['k2','k3']:
 rec={case:parse_qe(ROOT/f'data/qe_outputs/{k}_{case}.out',nat) for case,nat in [('stack',245),('uncapped',183),('hbn',62)]}
 delta[k]=(rec['stack'][1][:183]-rec['uncapped'][1])*RY_EV/BOHR_A
 mesh[k]={'interaction_eV_cell':(rec['stack'][0]-rec['uncapped'][0]-rec['hbn'][0])*RY_EV,**metrics(delta[k][135:183])}
report['low_moment_mesh']=mesh
report['low_moment_mesh_change']={'interaction_meV_cell':1000*(mesh['k3']['interaction_eV_cell']-mesh['k2']['interaction_eV_cell']),**metrics((delta['k3']-delta['k2'])[135:183])}
rows=list(csv.DictReader((ROOT/'data/forces/atom_forces.csv').open(encoding='utf-8-sig')))
film=[r for r in rows if r['species'] in ['Nb','Se']]
forces=np.array([[float(r['df'+ax+'_ev_angstrom']) for ax in 'xyz'] for r in film])
# Check published exported vector components against independently parsed raw outputs.
np.testing.assert_allclose(forces,delta['k3'][135:183],atol=1e-10,rtol=0)
sums={layer:sum(float(r['dfz_ev_angstrom']) for r in film if r['layer']==layer) for layer in ['Nb','Se_top','Se_bottom']}
report['figure2']={'plane_means_eV_A':{k:v/16 for k,v in sums.items()},'plane_sums_eV_A':sums,'net_eV_A':sum(sums.values()),'cancellation_fraction':1-abs(sum(sums.values()))/sum(map(abs,sums.values())),**metrics(forces)}
S=json.loads((ROOT/'data/high_moment/four_state.json').read_text())
means={}
for geom in ['GU','GC']:
 m=S['forces']['capping_on_substrate_'+geom]['modes']
 B=math.sqrt(24)*m['buckling_ev_a_unit_mode'];R=math.sqrt(48)*m['rigid_z_ev_a_unit_mode'];D=math.sqrt(32)*m['breathing_ev_a_unit_mode']
 nb=(2*B+R)/3;se=(2*R-2*B)/3;up=(se+D)/2;lo=(se-D)/2
 means[geom]={'Nb':nb/16,'upper_Se':up/16,'lower_Se':lo/16}
 np.testing.assert_allclose([(nb-.5*(up+lo))/math.sqrt(24),(nb+up+lo)/math.sqrt(48),(up-lo)/math.sqrt(32)],[m['buckling_ev_a_unit_mode'],m['rigid_z_ev_a_unit_mode'],m['breathing_ev_a_unit_mode']],atol=1e-14)
report['high_moment_plane_means']=means
v=np.genfromtxt(ROOT/'data/high_moment/displacement.csv',delimiter=',',names=True)
q=v['q_A'];fi=v['F_confined_eV_A']-v['F_unconfined_eV_A']
assert np.allclose(q,[-.02,0]);np.testing.assert_allclose(fi,v['f_int_eV_A'])
work=.5*(fi[0]+fi[1])*(q[0]-q[1])*1000
energy=v['deltaE_confined_meV'][0]-v['deltaE_unconfined_meV'][0]
np.testing.assert_allclose(energy,v['deltaI_meV'][0])
np.testing.assert_allclose(fi[1],S['forces']['capping_on_substrate_GC']['modes']['buckling_ev_a_unit_mode'])
report['displacement']={'start_q_A':0,'end_q_A':-.02,'deltaI_meV_cell':float(energy),'work_meV_cell':float(work),'negative_work_meV_cell':float(-work),'energy_minus_negative_work_meV_cell':float(energy+work),'relative_difference_percent':float(abs((energy+work)/energy)*100)}
# Recompute five-mode diagnostic from archived projections, not from a stored total.
terms={}
for mode in ['buckling','breathing','rigid_z','shear_x','shear_y']:
 f=np.mean([S['forces']['capping_on_substrate_'+g]['modes'][mode+'_ev_a_unit_mode'] for g in ['GU','GC']]);dq=S['GU_to_GC_displacement_modes'][mode+'_angstrom_unit_mode']
 terms[mode]=float(-1000*f*dq)
observed=S['energy']['delta_delta_capping_GC_minus_GU_ev']*1000
D={'observed_interaction_energy_change_meV_cell':observed,'endpoint_trapezoid_projected_energy_change_meV_cell':terms,'projected_sum_meV_cell':sum(terms.values()),'observed_minus_projected_sum_meV_cell':observed-sum(terms.values())}
(OUT/'endpoint_work_diagnostic.json').write_text(json.dumps(D,indent=2))
report['four_state_diagnostic']=D
# Primitive-cell benchmark total energies; not an interface-work error estimate.
benchmark=list(csv.DictReader((ROOT/'data/benchmark/results-25009137.csv').open()))
for row in benchmark:row['energy_eV']=float(row['total_energy_ry'])*RY_EV
report['primitive_benchmark']=benchmark
# Independent check that supplied coordinate perturbations implement the stated mode.
def positions(name):
 t=(ROOT/'data/inputs'/name).read_text().split('ATOMIC_POSITIONS angstrom')[1].split('CELL_PARAMETERS')[0]
 r=[line.split() for line in t.strip().splitlines()]
 return [v[0] for v in r],np.array([[float(x) for x in v[1:4]] for v in r])
sc,pc=positions('reference_C.in');su,pu=positions('reference_U.in')
assert sc[:183]==su
np.testing.assert_allclose(pc[:183],pu,atol=1e-9,rtol=0)
for name,ref,symbols in [('C_minus.in',pc,sc),('U_minus.in',pu,su)]:
 sm,pm=positions(name);assert sm==symbols
 expected=np.zeros_like(ref)
 for i,el in enumerate(symbols):
  if el=='Nb':expected[i,2]=-.02/math.sqrt(24)
  elif el=='Se':expected[i,2]=.01/math.sqrt(24)
 np.testing.assert_allclose(pm-ref,expected,atol=1e-9,rtol=0)
report['checks']={'figure2_export_matches_raw_QE':True,'reference_common_atom_coordinates_match':True,'negative_endpoint_mode_matches_inputs':True}
# Table S1: reference choices follow the manuscript, not a fit or extrapolation.
cutref=next(float(r['total_energy_ry']) for r in benchmark if r['case']=='ecut60_k12_n4')
mesh_ref=next(float(r['total_energy_ry']) for r in benchmark if r['case']=='ecut50_k15_n4')
with (OUT/'TableS1.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['series','case','energy_Ry','difference_meV'])
 for r in benchmark:
  if r['case'] in ['ecut40_k12_n4','ecut50_k12_n4','ecut60_k12_n4']:w.writerow(['cutoff',r['case'],r['total_energy_ry'],(float(r['total_energy_ry'])-cutref)*RY_EV*1000])
  if r['case'] in ['ecut50_k6_n4','ecut50_k9_n4','ecut50_k12_n4','ecut50_k15_n4']:w.writerow(['mesh',r['case'],r['total_energy_ry'],(float(r['total_energy_ry'])-mesh_ref)*RY_EV*1000])
with (OUT/'TableS2.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['mesh','interaction_eV_cell','force_RMS_eV_A','force_max_eV_A'])
 for k,r in mesh.items():w.writerow([k,r['interaction_eV_cell'],r['rms_eV_A'],r['max_eV_A']])
with (OUT/'TableS3.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['geometry','Nb_eV_A','upper_Se_eV_A','lower_Se_eV_A'])
 for k,r in means.items():w.writerow([k,r['Nb'],r['upper_Se'],r['lower_Se']])
(OUT/'analysis.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'force_RMS':report['figure2']['rms_eV_A'],'force_cancellation':report['figure2']['cancellation_fraction'],'work_comparison':report['displacement']},indent=2))

