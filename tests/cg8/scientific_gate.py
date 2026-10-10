"""Actual-stencil periodic-cell audit for eight-mechanism collocated attenuation.

Uses NumPy only. This harness is a scientific gate, not a runtime CG implementation.
"""
from pathlib import Path
import argparse, json, re, subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def stencil(name):
    source = (ROOT / 'src/JU_xJU_yJU_z6.f90').read_text()
    section = source.split('subroutine '+name+'(', 1)[1]
    section = section.split('end subroutine', 1)[0]
    interior = section.split('Jq_xU(1:3) =', 1)[1].split('Jq_xU(4:', 1)[0]
    pattern = r'([+-]?\s*[\d.]+_wp(?:\s*/\s*[\d.]+_wp)?)\s*\*\s*F%F\(x\s*([+-]\s*\d+)?\s*,'
    result = {}
    for coefficient, offset in re.findall(pattern, interior):
        # Restricted numeric expression from the checked-in production stencil.
        terms = coefficient.replace('_wp', '').replace(' ', '').split('/')
        value = float(terms[0]) if len(terms) == 1 else float(terms[0])/float(terms[1])
        result[int((offset or '0').replace(' ', ''))] = value
    if not result or abs(sum(result.values())) > 1e-12 or abs(sum(k*v for k,v in result.items())-1) > 1e-12:
        raise ValueError(f'Cannot recover a consistent production stencil: {name}: {result}')
    if name in ['JJU_x6_upwind','JJU_x6_upwind_drp']:
        optimized_name={'JJU_x6_upwind':'JJU_x6_interior_upwind',
                        'JJU_x6_upwind_drp':'JJU_x6_interior_upwind_drp'}[name]
        optimized=source.split('subroutine '+optimized_name+'(',1)[1].split('end subroutine',1)[0]
        arrays=[]
        for label in ['cf','cb']:
            expression=re.search(label+r'\(1:\d+\)\s*=\s*\(/(.*?)/\)',optimized,re.S).group(1)
            values=[]
            for term in expression.replace('&','').split(','):
                parts=term.replace('_wp','').replace(' ','').replace('\n','').split('/')
                values.append(float(parts[0]) if len(parts)==1 else float(parts[0])/float(parts[1]))
            arrays.append(np.array(values))
        assert np.allclose(arrays[0],[result[k] for k in sorted(result)],rtol=0,atol=3e-15),name
        assert np.allclose(arrays[1],-arrays[0][::-1],rtol=0,atol=3e-15),name
    return result


def derivatives(coefficients,theta,h):
    points=np.array([[i,j,k] for k in range(2) for j in range(2) for i in range(2)])
    matrices=[]
    for axis in range(3):
        d=np.zeros((8,8),complex)
        for row,point in enumerate(points):
            for shift,c in coefficients.items():
                other=point.copy(); other[axis]+=shift
                col=int(other[0]%2+2*(other[1]%2)+4*(other[2]%2))
                # A plane-wave phase at every physical node; the cell amplitudes
                # are periodic. This includes every stencil reach, not just ±1.
                d[row,col]+=c*np.exp(1j*theta[axis]*shift)/np.asarray(h)[axis]
        matrices.append(d)
    return matrices


def load_material(path):
    lines=Path(path).read_text().splitlines()
    config=np.fromstring(lines[0],sep=' ')
    bulk,shear,error=np.fromstring(lines[2],sep=' ')
    return dict(config=config,policy=lines[1].strip(),bulk=bulk,shear=shear,error=error,
                tau=np.fromstring(lines[3],sep=' '),wk=np.fromstring(lines[4],sep=' '),
                wm=np.fromstring(lines[5],sep=' '))


def operator(coefficients,direction,ppw,material,f=1,coarse=True,origin=(0,0,0),aspect=(1,1,1)):
    direction=np.array(direction,dtype=float); direction/=np.linalg.norm(direction)
    h=np.array(aspect,dtype=float)/(f*ppw);theta=direction*(2*np.pi*f)*h
    tau=material['tau'];bulk=material['bulk'];mu=material['shear']
    if coarse:
        points=np.array([[i,j,k] for k in range(2) for j in range(2) for i in range(2)])
        dp=derivatives(coefficients,theta,h); count=8;memories=1
    else:
        points=np.zeros((1,3),int);count=1;memories=8
        dp=[np.array([[sum(c*np.exp(1j*theta[axis]*shift) for shift,c in coefficients.items())/h[axis]]])
            for axis in range(3)]
    dm=[-d.conj().T for d in dp]
    a=np.zeros(((9+6*memories)*count,(9+6*memories)*count),complex)
    def strain_row(comp,node,kb,ms):
        g=np.zeros(3*count,complex)
        if comp<3:
            for axis in range(3):
                g[axis*count:(axis+1)*count]+=(kb-2*ms/3+(2*ms if axis==comp else 0))*dp[axis][node]
        else:
            x,y=[(0,1),(0,2),(1,2)][comp-3]
            g[x*count:(x+1)*count]+=ms*dp[y][node]
            g[y*count:(y+1)*count]+=ms*dp[x][node]
        return g
    for node in range(count):
        for comp in range(3):
            for axis in range(3):
                stress={(0,0):0,(1,1):1,(2,2):2,(0,1):3,(1,0):3,(0,2):4,(2,0):4,(1,2):5,(2,1):5}[comp,axis]
                a[comp*count+node,(3+stress)*count:(4+stress)*count]+=dm[axis][node]
        for comp in range(6):
            stress=(3+comp)*count+node
            a[stress,:3*count]=strain_row(comp,node,bulk,mu)
            ks=[int(((points[node]-origin)%2)@[1,2,4])] if coarse else range(8)
            for slot,k in enumerate(ks):
                memory=(9+6*slot+comp)*count+node
                a[stress,memory]=-1
                a[memory,:3*count]=strain_row(comp,node,bulk*material['wk'][k],mu*material['wm'][k])/tau[k]
                a[memory,memory]=-1/tau[k]
    return a,direction


def physical_modes(a,direction,f=1,cp=2):
    eig,vec=np.linalg.eig(a)
    count=8 if a.shape[0]==120 else 1
    candidates=[]
    for j,z in enumerate(eig):
        if z.imag <= 0: continue
        v=vec[:3*count,j].reshape(3,count); mean=v.mean(axis=1)
        coherence=float(count*np.vdot(mean,mean).real/max(np.vdot(v,v).real,1e-300))
        if coherence < .5: continue
        longitudinal=float(abs(direction@mean)**2/max(np.vdot(mean,mean).real,1e-300))
        q=float(z.imag/(-2*z.real)) if z.real<0 else None
        candidates.append(dict(omega=float(z.imag),q=q,coherence=coherence,longitudinal=longitudinal))
    s=sorted([c for c in candidates if c['longitudinal']<.2],key=lambda c:abs(c['omega']/(2*np.pi*f)-1))[:2]
    p=sorted([c for c in candidates if c['longitudinal']>.8],key=lambda c:abs(c['omega']/(cp*2*np.pi*f)-1))[:1]
    return s,p,float(np.max(eig.real))


def rk_amplification(z):
    # Kennedy-Carpenter coefficients also checked against production source.
    aa=[0,-567301805773/1357537059087,-2404267990393/2016746695238,
        -3550918686646/2091501179385,-1275806237668/842570457699]
    bb=[1432997174477/9575080441755,5161836677717/13612068292357,
        1720146321549/2090206949498,3134564353537/4481467310338,2277821191437/14882151754819]
    u=np.ones_like(z);du=np.zeros_like(z)
    for av,bv in zip(aa,bb):
        du=av*du+z*u;u=u+bv*du
    return u


def run(coefficients_dir):
    rows=[]
    for case in range(1,7):
        cg=load_material(coefficients_dir/f'cg8-case-{case}-1.txt')
        full=load_material(coefficients_dir/f'cg8-case-{case}-2.txt')
        for name in ['JJU_x6','JJU_x4_upwind','JJU_x4_upwind_drp','JJU_x6_upwind','JJU_x6_upwind_drp']:
            coefficients=stencil(name)
            for aspect in [(1,1,1),(.5,1,1),(1,.5,1),(1,1,.5)]:
             for direction in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]:
              for f in [.1,.2,.5,1,2,5]:
               for ppw in [16,32,64]:
                    a,d=operator(coefficients,direction,ppw,cg,f=f,aspect=aspect)
                    reference,_=operator(coefficients,direction,ppw,full,f=f,coarse=False,aspect=aspect)
                    sm,pm,growth=physical_modes(a,d,f,cg['config'][8]);sf,pf,_=physical_modes(reference,d,f,cg['config'][8])
                    # Use full fitted discrete reference to distinguish material Q
                    # from finite-loss modal Q and constitutive dispersion.
                    qerror=max([abs(mode['q']/ref['q']-1) for family,refs in [(sm,sf),(pm,pf)]
                        for mode,ref in zip(family,refs) if mode['q'] is not None],default=None)
                    phase=max([abs(mode['omega']/ref['omega']-1) for family,refs in [(sm,sf),(pm,pf)]
                        for mode,ref in zip(family,refs)],default=None)
                    dt=min(.25*min(aspect)/(f*ppw*np.sqrt(cg['bulk']+7*cg['shear']/3)),2*min(cg['tau']))
                    amp=float(np.max(np.abs(rk_amplification(np.linalg.eigvals(a)*dt))))
                    rows.append(dict(case=case,stencil=name,direction=direction,ppw=ppw,f=f,aspect=aspect,
                        shear=sm,pressure=pm,reference_shear=sf,reference_pressure=pf,
                        q_error=qerror,phase_error=phase,max_growth=growth,rk_max_amplification=amp))
    shifts=[]
    for name in ['JJU_x4_upwind','JJU_x6_upwind','JJU_x6_upwind_drp']:
     for case in [3,6]:
        material=load_material(coefficients_dir/f'cg8-case-{case}-1.txt')
        base,d=operator(stencil(name),(1,1,1),32,material)
        sb,pb,_=physical_modes(base,d,cp=material['config'][8])
        expected=sorted([m['omega'] for m in sb+pb])
        for origin in [(i,j,k) for i in range(2) for j in range(2) for k in range(2)]:
            a,d=operator(stencil(name),(1,1,1),32,material,origin=origin)
            sm,pm,_=physical_modes(a,d,cp=material['config'][8])
            observed=sorted([m['omega'] for m in sm+pm])
            difference=float(np.max(np.abs(np.array(observed)/expected-1)))
            shifts.append(dict(case=case,stencil=name,origin=origin,relative_frequency_difference=difference))
    return dict(description='Actual-stencil Bloch gate using exported production coarse/full fits; reference cs=1',
                rows=rows,pattern_shifts=shifts)


def assert_supported(result):
    supported=[r for r in result['rows'] if r['case']>=3 and r['stencil'] in ['JJU_x4_upwind','JJU_x6_upwind','JJU_x6_upwind_drp']]
    for row in supported:
        assert len(row['shear'])==2 and len(row['pressure'])==1, row
        assert row['q_error']<=.05 and row['phase_error']<=.005, row
        assert row['max_growth']<1e-9 and row['rk_max_amplification']<=1+1e-9, row
    groups={}
    for row in supported:
        for field in ['shear','pressure']:
            key=(row['stencil'],row['case'],row['f'],row['ppw'],tuple(row['aspect']),field)
            groups.setdefault(key,[]).extend(m['q'] for m in row[field])
    spreads=[max(q)/min(q)-1 for q in groups.values()]
    assert max(spreads)<=.02, max(spreads)
    assert max(s['relative_frequency_difference'] for s in result['pattern_shifts'])<1e-9
    # The harness must actually expose the unsupported centered and low-Q cases.
    centered=[r for r in result['rows'] if r['stencil']=='JJU_x6']
    assert any(len(r['shear'])!=2 or len(r['pressure'])!=1 for r in centered)
    assert max(r['q_error'] for r in result['rows'] if r['case']==1 and r['stencil']=='JJU_x4_upwind')>.05
    return dict(cases=len(supported),max_q_error=max(r['q_error'] for r in supported),
                max_phase_error=max(r['phase_error'] for r in supported),max_directional_spread=max(spreads))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--coefficients-dir',type=Path,default=Path.cwd())
    parser.add_argument('--coefficients-exe',type=Path)
    parser.add_argument('--assert-supported',action='store_true')
    args=parser.parse_args()
    if args.coefficients_exe:
        subprocess.run([str(args.coefficients_exe)],cwd=args.coefficients_dir,check=True)
    result=run(args.coefficients_dir)
    if args.assert_supported: result['supported_summary']=assert_supported(result)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'Wrote {len(result["rows"])} actual-stencil Bloch cases to {args.output}')
    if 'supported_summary' in result: print(result['supported_summary'])
