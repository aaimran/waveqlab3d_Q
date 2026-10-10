"""Exploratory projected-cell audit for planning; not a production CG-T test."""
import sys,json,argparse
from pathlib import Path
import numpy as np
sys.dont_write_bytecode=True
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'tests/cg8'))
from scientific_gate import stencil,operator,derivatives,physical_modes,rk_amplification,load_material

def projected(coeff,direction,ppw,m,f=1,aspect=(1,1,1),origin=(0,0,0)):
 a,d=operator(coeff,direction,ppw,m,f=f,coarse=True,aspect=aspect)
 h=np.array(aspect)/(f*ppw);theta=2*np.pi*f*d*h
 pts=np.array([[i,j,k] for k in range(2) for j in range(2) for i in range(2)])
 lift=np.exp(-1j*(((pts-np.array(origin))%2)@theta));restrict=lift.conj()/8
 dp=derivatives(coeff,theta,h)
 a[:72,72:]=0;a[72:,:]=0
 for comp in range(6):
  for k in range(8):
   row=72+comp*8+k
   kb=m['bulk']*m['wk'][k];mu=m['shear']*m['wm'][k]
   for axis in range(3):
    g=np.zeros((8,8),complex)
    if comp<3:g=(kb-2*mu/3+(2*mu if axis==comp else 0))*dp[axis]
    else:
     x,y=[(0,1),(0,2),(1,2)][comp-3]
     if axis==x:g=mu*dp[y]
     if axis==y:g=mu*dp[x]
    a[row,axis*8:(axis+1)*8]=restrict@g/m['tau'][k]
   a[row,row]=-1/m['tau'][k]
   a[(3+comp)*8:(4+comp)*8,row]=-lift
 return a,d

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--coefficients-dir',type=Path,default=root/'build/test_runs/cg8_coefficients_unit')
    parser.add_argument('--output',type=Path,default=root/'docs/cg-t-projection-audit.json')
    args=parser.parse_args()
    # A separate weighted inner-product sanity check, including nonuniform weights.
    rng=np.random.default_rng(2409)
    w=rng.uniform(.3,2,8);hf=np.diag(w);hc=np.array([[w.sum()]])
    r=(w/w.sum())[None,:];p_lift=np.linalg.solve(hf,r.T@hc)
    algebra=dict(adjoint_defect=float(np.max(abs(hf@p_lift-r.T@hc))),
                 rp_defect=float(np.max(abs(r@p_lift-1))),
                 minimum_relaxed_energy_eigenvalue=float(np.linalg.eigvalsh(hf-.8*r.T@hc@r).min()))
    assert algebra['adjoint_defect']<1e-14 and algebra['rp_defect']<1e-14
    assert algebra['minimum_relaxed_energy_eigenvalue']>0

    rows=[]
    p=args.coefficients_dir
    for name in ['JJU_x6','JJU_x6_upwind','JJU_x6_upwind_drp']:
     for case in [3,4,5,6]:
      m=load_material(p/f'cg8-case-{case}-2.txt')
      for direction in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]:
       for ppw in [16,24,32,64]:
        for f in [.1,1,5]:
         a,d=projected(stencil(name),direction,ppw,m,f)
         b,_=operator(stencil(name),direction,ppw,m,f=f,coarse=False)
         s,t,g=physical_modes(a,d,f,m['config'][8]);u,v,_=physical_modes(b,d,f,m['config'][8])
         clean=len(s)==2 and len(t)==1 and all(x['q'] is not None for x in s+t)
         qe=max([abs(x['q']/y['q']-1) for x,y in zip(s+t,u+v) if x['q'] is not None],default=None)
         pe=max([abs(x['omega']/y['omega']-1) for x,y in zip(s+t,u+v)],default=None)
         dt=min(.25/(f*ppw*np.sqrt(m['bulk']+7*m['shear']/3)),2*min(m['tau']))
         amp=float(max(abs(rk_amplification(np.linalg.eigvals(a)*dt))))
         rows.append(dict(stencil=name,case=case,direction=direction,ppw=ppw,f=f,clean=clean,q_error=qe,phase_error=pe,growth=g,rk=amp))
    summary={}
    for name in ['JJU_x6','JJU_x6_upwind','JJU_x6_upwind_drp']:
     r=[r for r in rows if r['stencil']==name]
     summary[name]=dict(cases=len(r),unclean=sum(not x['clean'] for x in r),max_q_error=max(x['q_error'] for x in r if x['q_error'] is not None),max_phase_error=max(x['phase_error'] for x in r if x['phase_error'] is not None),max_growth=max(x['growth'] for x in r),max_rk=max(x['rk'] for x in r))
     print(name,summary[name],flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(description='Exploratory projected-cell Bloch prototype; not a production solver validation',algebra=algebra,summary=summary,rows=rows),indent=2)+'\n')

if __name__=='__main__': main()
