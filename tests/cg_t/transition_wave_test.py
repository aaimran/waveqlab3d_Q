"""Independent nine-field full/projected/full packet test, with all eight cell memories."""
import sys,argparse,json,subprocess
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from projection_audit import load_material,stencil
AA=[0,-567301805773/1357537059087,-2404267990393/2016746695238,-3550918686646/2091501179385,-1275806237668/842570457699]
BB=[1432997174477/9575080441755,5161836677717/13612068292357,1720146321549/2090206949498,3134564353537/4481467310338,2277821191437/14882151754819]

def propagate(m,hybrid,name,wave):
 h=.04;nx=600;shape=(nx,2,2);mask=np.zeros(shape,bool)
 if hybrid:mask[200:300,:,:]=True
 cell_mask=mask[::2,0,0]
 bulk=m['bulk'];mu=m['shear'];tau=m['tau'];wk=m['wk'];wm=m['wm']
 dt=min(.25*h/np.sqrt(bulk+7*mu/3),2*min(tau));coeff=stencil(name)
 def derivative(v,axis,forward):
  out=np.zeros_like(v)
  for shift,c in coeff.items():out+=(c if forward else -c)*np.roll(v,-shift if forward else shift,axis=axis)/h
  return out
 def rhs(fields,full,cell):
  grad=[derivative(fields[...,:3],axis,True) for axis in range(3)]
  strain=np.stack([grad[0][...,0],grad[1][...,1],grad[2][...,2],
   (grad[1][...,0]+grad[0][...,1])/2,(grad[2][...,0]+grad[0][...,2])/2,
   (grad[2][...,1]+grad[1][...,2])/2],axis=-1)
  tr=np.sum(strain[...,:3],axis=-1)
  stress_grad=[derivative(fields[...,3:],axis,False) for axis in range(3)]
  df=np.zeros_like(fields)
  df[...,0]=stress_grad[0][...,0]+stress_grad[1][...,3]+stress_grad[2][...,4]
  df[...,1]=stress_grad[0][...,3]+stress_grad[1][...,1]+stress_grad[2][...,5]
  df[...,2]=stress_grad[0][...,4]+stress_grad[1][...,5]+stress_grad[2][...,2]
  df[...,3:]=2*mu*strain-np.sum(full,axis=-2)
  df[...,3:6]+=(bulk-2*mu/3)*tr[...,None]
  feedback=np.repeat(np.sum(cell,axis=1),2,axis=0)[:,None,None,:]
  df[...,3:]-=feedback
  forcing=2*mu*wm[None,None,None,:,None]*strain[...,None,:]
  forcing[...,:3]+=(bulk*wk-2*mu*wm/3)[None,None,None,:,None]*tr[...,None,None]
  dm=(forcing-full)/tau[None,None,None,:,None];dm[mask]=0
  avg=strain.reshape(nx//2,2,2,2,6).mean(axis=(1,2,3));trace=avg[:,:3].sum(axis=1)
  force=2*mu*wm[None,:,None]*avg[:,None,:]
  force[:,:,:3]+=(bulk*wk-2*mu*wm/3)[None,:,None]*trace[:,None,None]
  dc=(force-cell)/tau[None,:,None];dc[~cell_mask]=0
  return df,dm,dc
 fields=np.zeros(shape+(9,));full=np.zeros(shape+(8,6));cell=np.zeros((nx//2,8,6))
 x=np.arange(nx)*h;packet=np.exp(-((x-3)/.6)**2)*np.cos(2*np.pi*(x-3))
 if wave=='S':
  fields[...,1]=packet[:,None,None];fields[...,6]=-np.sqrt(mu)*packet[:,None,None];component=1;end=12
 else:
  cp=np.sqrt(bulk+4*mu/3);lam=bulk-2*mu/3
  fields[...,0]=packet[:,None,None];fields[...,3]=-cp*packet[:,None,None]
  fields[...,4]=-lam/cp*packet[:,None,None];fields[...,5]=-lam/cp*packet[:,None,None]
  component=0;end=6
 rate=np.zeros_like(fields);fr=np.zeros_like(full);cr=np.zeros_like(cell);history=[]
 for step in range(int(end/dt)):
  for aa,bb in zip(AA,BB):
   df,dm,dc=rhs(fields,full,cell)
   rate=aa*rate+dt*df;fr=aa*fr+dt*dm;cr=aa*cr+dt*dc
   fields+=bb*rate;full+=bb*fr;cell+=bb*cr
  history.append([dt*(step+1),float(fields[125,...,component].mean()),float(fields[350,...,component].mean())])
 assert np.isfinite(fields).all() and np.isfinite(full).all() and np.isfinite(cell).all()
 return np.array(history)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--coefficients-exe',type=Path,required=True)
 p.add_argument('--stencil',required=True);p.add_argument('--wave',choices=['P','S'],required=True)
 args=p.parse_args();subprocess.run([str(args.coefficients_exe)],check=True)
 report=[]
 for case in ([3,4] if args.wave=='P' else [5,6]):
  m=load_material(f'cg8-case-{case}-2.txt')
  reference=propagate(m,False,args.stencil,args.wave);hybrid=propagate(m,True,args.stencil,args.wave)
  incident=max(abs(reference[:,1]));assert incident>.1
  window=(reference[:,0]>= (3 if args.wave=='P' else 6))&(reference[:,0]<= (6 if args.wave=='P' else 11))
  reflection=float(max(abs(reference[window,1]-hybrid[window,1]))/incident)
  row=dict(case=case,stencil=args.stencil,wave=args.wave,extra_reflection=reflection,incident_peak=float(incident))
  print(row,flush=True);assert reflection<=1e-3,row;report.append(row)
 Path('cgt-transition-wave-report.json').write_text(json.dumps(report,indent=2)+'\n')
