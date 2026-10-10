"""Independent 3-D periodic-grid wave-packet test of the full/CG switch.

The two transverse directions each contain a period-two cell. All nine fields
and six stress memories evolve; this is not a scalar harmonic-response test.
"""
from pathlib import Path
import argparse,json,subprocess,sys
sys.dont_write_bytecode=True
import numpy as np
from scientific_gate import load_material,stencil

AA=[0,-567301805773/1357537059087,-2404267990393/2016746695238,-3550918686646/2091501179385,-1275806237668/842570457699]
BB=[1432997174477/9575080441755,5161836677717/13612068292357,1720146321549/2090206949498,3134564353537/4481467310338,2277821191437/14882151754819]


def propagate(full,cg,hybrid,name):
    h=.04;nx=600;shape=(nx,2,2)
    points=np.stack(np.meshgrid(np.arange(nx),np.arange(2),np.arange(2),indexing='ij'),axis=-1)
    mask=(points[...,0]>=200)&(points[...,0]<300) if hybrid else np.zeros(shape,bool)
    # Correct global period-two index; the explicit formula avoids x multiples of 2 aliasing y/z.
    k=(points[...,0]%2)+2*(points[...,1]%2)+4*(points[...,2]%2)
    bulk=np.where(mask,cg['bulk'],full['bulk']);mu=np.where(mask,cg['shear'],full['shear'])
    wk=np.broadcast_to(full['wk'],shape+(8,)).copy();wm=np.broadcast_to(full['wm'],shape+(8,)).copy()
    wk[mask]=0;wm[mask]=0
    for mechanism in range(8):
        selected=mask&(k==mechanism);wk[...,mechanism][selected]=cg['wk'][mechanism];wm[...,mechanism][selected]=cg['wm'][mechanism]
    tau=full['tau'];assert np.array_equal(tau,cg['tau'])
    dt=min(.25*h/np.sqrt(np.max(bulk+7*mu/3)),2*np.min(tau))
    # Identical time sampling for reference and hybrid.
    dt=min(dt,.25*h/np.sqrt(max(full['bulk']+7*full['shear']/3,cg['bulk']+7*cg['shear']/3)))
    coefficients=stencil(name)
    def derivative(v,axis,forward):
        out=np.zeros_like(v)
        for shift,c in coefficients.items():
            out+=(c if forward else -c)*np.roll(v,-shift if forward else shift,axis=axis)/h
        return out
    def rhs(fields,memory):
        grad=[derivative(fields[...,:3],axis,True) for axis in range(3)]
        tr=grad[0][...,0]+grad[1][...,1]+grad[2][...,2]
        strain=np.stack([grad[0][...,0]-tr/3,grad[1][...,1]-tr/3,grad[2][...,2]-tr/3,
                        (grad[1][...,0]+grad[0][...,1])/2,(grad[2][...,0]+grad[0][...,2])/2,
                        (grad[2][...,1]+grad[1][...,2])/2],axis=-1)
        stress_grad=[derivative(fields[...,3:],axis,False) for axis in range(3)]
        df=np.empty_like(fields)
        df[...,0]=stress_grad[0][...,0]+stress_grad[1][...,3]+stress_grad[2][...,4]
        df[...,1]=stress_grad[0][...,3]+stress_grad[1][...,1]+stress_grad[2][...,5]
        df[...,2]=stress_grad[0][...,4]+stress_grad[1][...,5]+stress_grad[2][...,2]
        df[...,3:]=2*mu[...,None]*strain-np.sum(memory,axis=-2)
        df[...,3:6]+=bulk[...,None]*tr[...,None]
        forcing=2*mu[...,None,None]*wm[...,None]*strain[...,None,:]
        forcing[...,:3]+=bulk[...,None,None]*wk[...,None]*tr[...,None,None]
        dm=(forcing-memory)/tau[None,None,None,:,None]
        return df,dm
    fields=np.zeros(shape+(9,));memory=np.zeros(shape+(8,6))
    x=np.arange(nx)*h;packet=np.exp(-((x-3)/.6)**2)*np.cos(2*np.pi*(x-3))
    fields[...,1]=packet[:,None,None];fields[...,6]=-np.sqrt(full['shear'])*packet[:,None,None]
    rate=np.zeros_like(fields);memory_rate=np.zeros_like(memory)
    history=[];steps=int(12/dt)
    for step in range(steps):
        for aa,bb in zip(AA,BB):
            df,dm=rhs(fields,memory)
            rate=aa*rate+dt*df;memory_rate=aa*memory_rate+dt*dm
            fields+=bb*rate;memory+=bb*memory_rate
        history.append([dt*(step+1),float(fields[125,...,1].mean()),float(fields[350,...,1].mean())])
    assert np.isfinite(fields).all() and np.isfinite(memory).all()
    return np.array(history)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--coefficients-exe',type=Path,required=True)
    parser.add_argument('--stencil',default='JJU_x4_upwind',choices=['JJU_x4_upwind','JJU_x6_upwind','JJU_x6_upwind_drp'])
    args=parser.parse_args();subprocess.run([str(args.coefficients_exe)],check=True)
    report=[]
    for case in [5,6]:
        cg=load_material(f'cg8-case-{case}-1.txt');full=load_material(f'cg8-case-{case}-2.txt')
        reference=propagate(full,cg,False,args.stencil);hybrid=propagate(full,cg,True,args.stencil)
        assert reference.shape==hybrid.shape
        incident=np.max(np.abs(reference[:,1]));assert incident>.1
        window=(reference[:,0]>=6)&(reference[:,0]<=11)
        extra_reflection=float(np.max(np.abs(reference[window,1]-hybrid[window,1]))/incident)
        assert extra_reflection<=1e-3,(case,extra_reflection)
        report.append(dict(case=case,stencil=args.stencil,extra_reflection=extra_reflection,incident_peak=float(incident)))
        print(report[-1],flush=True)
    Path('cg8-transition-wave-report.json').write_text(json.dumps(report,indent=2)+'\n')
