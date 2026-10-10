"""Compare actual receiver samples, not rounded final maxima, across MPI layouts."""
from pathlib import Path
import argparse,json,math,re,subprocess


def launcher(config):
    text=Path(config).read_text()
    def get(name):
        return re.search(r'set\('+name+r' \[==\[(.*?)\]==\]\)',text,re.S)[1]
    return get('MPIEXEC'),get('MPIEXEC_NUMPROC_FLAG'),get('MPIEXEC_PREFLAGS').split(';') if get('MPIEXEC_PREFLAGS') else [],get('MPIEXEC_POSTFLAGS').split(';') if get('MPIEXEC_POSTFLAGS') else []


def read_histories(folder):
    result={}
    for path in sorted((folder/'receivers').rglob('*.dat')):
        result[str(path.relative_to(folder/'receivers'))]=[[float(x.replace('D','E')) for x in line.split()]
            for line in path.read_text().splitlines() if line.strip() and not line.startswith('#')]
    assert result, f'No receiver output in {folder}'
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--exe',required=True)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--mpi-config',required=True)
    parser.add_argument('--fd-type',choices=['upwind','upwind_drp'])
    parser.add_argument('--order',type=int)
    args=parser.parse_args();runner,flag,pre,post=launcher(args.mpi_config)
    text=args.input.read_text().replace('output_seismograms=F','output_seismograms=T')
    if args.fd_type: text=re.sub(r"fd_type='[^']+'",f"fd_type='{args.fd_type}'",text)
    if args.order: text=re.sub(r'order=\d+',f'order={args.order}',text,count=1)
    text=re.sub(r"name='[^']+'","name='cg8-history'",text,count=1)
    text=text.replace('&output_list',"&output_list\n station_file_directory='receivers', station_number_in_list=T, station_number_in_filename=T,")
    # Origin zero cuts period-two cells at the odd MPI split, testing global indexing.
    marker=' boundary_policy=' if ' boundary_policy=' in text else " coefficient_policy="
    text=text.replace("boundary_policy='full-buffer'", "boundary_policy='full-buffer', pattern_origin=0,1,0")
    text+='\n!---begin:station_list---\n1 -0.08d0 0d0 0.04d0\n2 0.08d0 0d0 0.04d0\n3 -0.14d0 0d0 0.04d0\n4 0.14d0 0d0 0.04d0\n!---end:station_list---\n'
    baseline=None;report=[]
    for n in [1,2,3,4]:
        folder=Path.cwd()/f'ranks-{n}';(folder/'data').mkdir(parents=True,exist_ok=True)
        inp=folder/'case.in';inp.write_text(text)
        run=subprocess.run([runner,flag,str(n),*pre,args.exe,*post,str(inp)],cwd=folder,text=True,capture_output=True)
        (folder/'run.log').write_text(run.stdout+'\n'+run.stderr)
        assert run.returncode==0,run.stdout+'\n'+run.stderr
        current=read_histories(folder)
        peak=max(abs(v) for values in current.values() for row in values for v in row[1:])
        assert peak>1e-12 and math.isfinite(peak),'Source did not produce a measurable signal'
        difference=0
        if baseline is not None:
            assert current.keys()==baseline.keys()
            for key,values in current.items():
                assert len(values)==len(baseline[key])
                for row,reference in zip(values,baseline[key]):
                    assert len(row)==len(reference)==4
                    for a,b in zip(row,reference):
                        assert math.isfinite(a)
                        difference=max(difference,abs(a-b))
            assert difference<=1e-10*peak+1e-14,(n,difference,peak)
        else:
            baseline=current
        report.append(dict(ranks=n,receiver_files=len(current),peak=peak,max_absolute_difference=difference))
        print(report[-1],flush=True)
    (Path.cwd()/'mpi-history-report.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__': main()
