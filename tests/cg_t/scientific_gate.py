"""Projected-cell accuracy/stability gates using actual solver stencils and full fits."""
import argparse,json,sys,subprocess
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from projection_audit import projected
from scientific_gate import stencil,operator,physical_modes,load_material,rk_amplification

def energy_check(coeff,m,direction):
 a,d=projected(coeff,direction,24.1,m)
 points=np.array([[i,j,k] for k in range(2) for j in range(2) for i in range(2)])
 lift=np.exp(-1j*(points@(2*np.pi*d/24.1)));r=lift.conj()[None,:]/8
 def stiffness(bulk,mu):
  c=np.zeros((6,6));c[:3,:3]=bulk-2*mu/3
  c[:3,:3]+=2*mu*np.eye(3);c[3:,3:]=2*mu*np.eye(3)
  return c
 cu=stiffness(m['bulk'],m['shear'])
 dc=[stiffness(m['bulk']*m['wk'][k],m['shear']*m['wm'][k]) for k in range(8)]
 pe=np.kron(np.eye(6),lift[:,None]);re=np.kron(np.eye(6),r)
 cr=np.kron(cu,np.eye(8))-pe@sum(dc)@re
 inv=np.linalg.inv(cr);transfer=np.zeros((6,48));internal=np.zeros((48,48))
 for k in range(8):
  ids=np.arange(6)*8+k;transfer[np.arange(6),ids]=m['tau'][k]
  internal[np.ix_(ids,ids)]=8*m['tau'][k]**2*np.linalg.inv(dc[k])
 gamma=pe@transfer
 g=np.zeros((120,120),complex);g[:24,:24]=np.eye(24)
 g[24:72,24:72]=inv;g[24:72,72:]=-inv@gamma;g[72:,24:72]=-gamma.conj().T@inv
 g[72:,72:]=gamma.conj().T@inv@gamma+internal
 scale=np.ones(120)
 for comp in range(3,6):
  scale[(3+comp)*8:(4+comp)*8]=np.sqrt(2);scale[72+comp*8:72+(comp+1)*8]=np.sqrt(2)
 am=scale[:,None]*a/scale[None,:]
 balance=np.sqrt(np.real(np.diag(g)))
 gb=g/balance[:,None]/balance[None,:]
 diss=(g@am+am.conj().T@g)/balance[:,None]/balance[None,:]
 minimum=float(np.linalg.eigvalsh(gb).min());growth=float(np.linalg.eigvalsh(diss).max())
 assert minimum>0 and growth<1e-8,(minimum,growth)
 return dict(minimum_storage_eigenvalue=minimum,maximum_dissipation_eigenvalue=growth)

def run(folder):
 rows=[];energies=[]
 directions=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1),(1,2,3),(3,1,2),(2,3,1)]
 for name in ['JJU_x6','JJU_x6_upwind','JJU_x6_upwind_drp']:
  coeff=stencil(name)
  for direction in [(1,0,0),(1,2,3)]:
   energy=energy_check(coeff,load_material(folder/'cg8-case-3-2.txt'),direction)
   energies.append(dict(stencil=name,direction=direction,**energy))
  for case in [3,4,5,6]:
   m=load_material(folder/f'cg8-case-{case}-2.txt')
   fmin,fmax=m['config'][6:8]
   frequencies=sorted(set([fmin,fmax,.1,1.,5.,m['config'][3]*m['config'][4],m['config'][3]*m['config'][5]]))
   frequencies=[f for f in frequencies if fmin<=f<=fmax]
   for aspect in [(1,1,1),(.5,1,1),(1,.5,1),(1,1,.5)]:
    for direction in directions:
     for f in frequencies:
      for ppw in [24.1,32,64]:
       a,d=projected(coeff,direction,ppw,m,f,aspect)
       b,_=operator(coeff,direction,ppw,m,f=f,coarse=False,aspect=aspect)
       s,p,g=physical_modes(a,d,f,m['config'][8]);sr,pr,_=physical_modes(b,d,f,m['config'][8])
       assert len(s)==2 and len(p)==1,(name,case,direction,f,ppw,s,p)
       assert all(x['q'] is not None and x['coherence']>.8 for x in s+p)
       qe=max(abs(x['q']/y['q']-1) for x,y in zip(s+p,sr+pr))
       pe=max(abs(x['omega']/y['omega']-1) for x,y in zip(s+p,sr+pr))
       dt=min(.25*min(aspect)/(f*ppw*np.sqrt(m['bulk']+7*m['shear']/3)),2*min(m['tau']))
       amp=float(max(abs(rk_amplification(np.linalg.eigvals(a)*dt))))
       row=dict(stencil=name,case=case,direction=direction,aspect=aspect,f=f,ppw=ppw,q_error=qe,phase_error=pe,
                growth=g,rk=amp,shear=s,pressure=p)
       assert qe<=.03 and pe<=.005,row
       assert g<1e-8 and amp<=1+1e-9,row
       rows.append(row)
  print(name,'completed',flush=True)
 shifts=[]
 for name in ['JJU_x6','JJU_x6_upwind','JJU_x6_upwind_drp']:
  for case in [3,6]:
   m=load_material(folder/f'cg8-case-{case}-2.txt')
   a,d=projected(stencil(name),(1,2,3),24.1,m)
   s,p,_=physical_modes(a,d,cp=m['config'][8]);base=np.array([x['omega'] for x in s+p])
   for origin in [(i,j,k) for i in range(2) for j in range(2) for k in range(2)]:
    a,d=projected(stencil(name),(1,2,3),24.1,m,origin=origin)
    s,p,_=physical_modes(a,d,cp=m['config'][8]);v=np.array([x['omega'] for x in s+p])
    diff=float(max(abs(v/base-1)));assert diff<1e-9,(name,case,origin,diff)
    shifts.append(dict(stencil=name,case=case,origin=origin,frequency_difference=diff))
 research=[]
 for name in ['JJU_x6','JJU_x6_upwind','JJU_x6_upwind_drp']:
  for case in [1,2]:
   m=load_material(folder/f'cg8-case-{case}-2.txt')
   for direction in [(1,0,0),(1,1,1)]:
    for f in [.1,1,5]:
     for ppw in [16,24.1,64]:
      a,d=projected(stencil(name),direction,ppw,m,f)
      s,p,g=physical_modes(a,d,f,m['config'][8])
      assert g<1e-8
      research.append(dict(stencil=name,case=case,direction=direction,f=f,ppw=ppw,growth=g,
        physical_branches=len(s)+len(p)))
 groups={}
 for row in rows:
  for field in ['shear','pressure']:
   key=(row['stencil'],row['case'],row['f'],row['ppw'],tuple(row['aspect']),field)
   groups.setdefault(key,[]).extend(x['q'] for x in row[field])
 spreads=[max(x)/min(x)-1 for x in groups.values()]
 assert max(spreads)<=.02,max(spreads)
 return dict(rows=rows,shifts=shifts,energy_checks=energies,low_q_research=research,summary=dict(cases=len(rows),max_q_error=max(r['q_error'] for r in rows),
   max_phase_error=max(r['phase_error'] for r in rows),max_directional_spread=max(spreads)))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--coefficients-exe',type=Path,required=True)
 p.add_argument('--output',type=Path,default=Path('cgt-bloch-report.json'))
 args=p.parse_args();subprocess.run([str(args.coefficients_exe)],check=True)
 result=run(Path.cwd());args.output.write_text(json.dumps(result,indent=2)+'\n');print(result['summary'])
