"""Exercise split receiver configuration, compatibility, filenames and MPI writers."""
from pathlib import Path
import argparse, os, re, subprocess, json, math, shutil

parser=argparse.ArgumentParser()
parser.add_argument('--exe',type=Path,required=True)
parser.add_argument('--mpiexec',default='mpiexec')
parser.add_argument('--work',type=Path,required=True)
args=parser.parse_args(); exe=args.exe.resolve(); args.work.mkdir(parents=True,exist_ok=True)
root=Path(__file__).resolve().parents[2]
base=(root/'test_problems/test_rup_cart_fric.in').read_text()
base=base.replace("name = 'data/test_rup_cart_fric'","name = 'case'").replace('t_final = 0.5d0','t_final = 0.1d0').replace('w_fault = T','w_fault = F')
base=re.sub(r'&output_list.*?/', '&output_list\n output_exact_moment=F, output_fault_topo=F, output_fields_block1=F,\n output_fields_block2=F, stride_fields=1\n/',base,flags=re.S)
base=base.replace("nqrs = 21, 21, 21", "nqrs = 41, 21, 21")
base=base[:base.index('!---begin:station_listU---')]
settings=""" output_seismograms=T, output_station_info=F, output_station_mapping=F,
 station_xyz_index=F, station_list='infile', station_list_file='',
 station_file_directory='receivers', station_output_order='t vx vy vz',
 station_number_in_list=T, station_number_in_filename=T,
 station_use_block_subdirectories=F, station_add_header=T, station_add_metadata=T,
"""
rows='!---begin:station_list---\n1 -0.6 7.5 0\n2 0.6 7.5 0\n3 0 7.5 0\n!---end:station_list---\n'
def config(side='both',suffix=True,subdirs=False):
 s=settings.replace('station_use_block_subdirectories=F',f'station_use_block_subdirectories={"T" if subdirs else "F"}')
 return base+'\n&output_receiver_stations\n'+s+f" append_block={'T' if suffix else 'F'}, interface_stations='{side}'\n/\n"+rows
report=[]
def run(label,text,n=1,failure=None,preflight=False):
 folder=args.work/label; folder.mkdir(parents=True,exist_ok=True); (folder/'data').mkdir(exist_ok=True)
 if (folder/'receivers').exists(): shutil.rmtree(folder/'receivers')
 (folder/'case.in').write_text(text)
 env=os.environ.copy()
 if preflight: env['WAVEQLAB3D_PREFLIGHT_ONLY']='1'
 result=subprocess.run([args.mpiexec,'-n',str(n),str(exe),'case.in'],cwd=folder,env=env,text=True,capture_output=True,timeout=60)
 log=result.stdout+result.stderr; (folder/'run.log').write_text(log)
 if failure:
  assert result.returncode!=0 and failure in log,(label,result.returncode,log[-2500:])
  assert not list((folder/'receivers').rglob('*.dat')),(label,'opened files before rejection')
 else:
  assert result.returncode==0,(label,log[-3500:])
 report.append({'case':label,'ranks':n,'passed':True})
 print(label,'PASS',flush=True)
 return folder

def histories(folder):
 out={}
 for p in (folder/'receivers').rglob('*.dat'):
  rows=[[float(v.replace('D','E')) for v in l.split()] for l in p.read_text().splitlines() if l and not l.startswith('#')]
  assert rows and all(len(r)==4 and all(math.isfinite(v) for v in r) for r in rows)
  header=next((l for l in p.read_text().splitlines() if l.startswith('# t ') or l.startswith('# vz t ')),None)
  tcol=1 if header and header.startswith('# vz t ') else 0
  assert all(a[tcol]<b[tcol] for a,b in zip(rows,rows[1:]))
  out[str(p.relative_to(folder/'receivers'))]=rows
 return out
both=histories(run('both',config()))
assert len(both)==4,both.keys()
assert any(max(abs(x-y) for x,y in zip(a[1:],b[1:]))>1e-10 for a,b in zip(both['case_station-3_block1.dat'],both['case_station-3_block2.dat'])),'slipping sides must differ'
for side,block in [('block_1',1),('block_2',2)]:
 h=histories(run(side,config(side,False)))
 assert set(h)=={f'case_station-{i}.dat' for i in [1,2,3]},h.keys()
 for i,b in [(1,1),(2,2),(3,block)]: assert h[f'case_station-{i}.dat']==both[f'case_station-{i}_block{b}.dat']
sub=histories(run('subdirs',config(subdirs=True,suffix=False)))
assert len(sub)==4 and len({p.split('/')[0] for p in sub})==2
legacy=base.replace(' stride_fields=1',' stride_fields=1,\n'+settings+" common_stations_blocks='both'")+rows
h=histories(run('legacy',legacy)); assert h==both
for n in [2,3,4]:
 h=histories(run(f'mpi-{n}',config(),n)); assert h.keys()==both.keys()
 old=histories(run(f'legacy-mpi-{n}',legacy,n)); assert h==old
 # The legacy slipping solver has serial/MPI differences; receiver refactoring
 # must match legacy samples at each decomposition, rather than mask them.
 for k in h: assert len(h[k])==len(both[k])
run('collision',config(suffix=False),2,failure='CFG-RECEIVER-COLLISION')
run('bad-side',config('invalid'),2,failure='interface_stations must',preflight=True)
run('conflict',config().replace('stride_fields=1','stride_fields=1, output_seismograms=F'),failure='occurs in &output_list',preflight=True)
run('conflict-after-path',config().replace('stride_fields=1',"stride_fields=1, station_file_directory='a/b', station_add_header=T"),failure='occurs in &output_list',preflight=True)
run('bad-order',config().replace("'t vx vy vz'","'t vx vx vz'"),failure='CFG-RECEIVER-ORDER',preflight=True)
run('duplicate-section',config()+'\n&output_receiver_stations /\n',failure='at most one',preflight=True)
run('malformed',config().replace('append_block=T','unknown_receiver_option=T'),failure='Cannot parse',preflight=True)
run('unmapped',config().replace('2 0.6 7.5 0','2 100 7.5 0'),failure='CFG-RECEIVER-UNMAPPED')
run('duplicate-id',config().replace('2 0.6','1 0.6'),failure='duplicate station number')
external=config().replace("station_list='infile'","station_list='extfile'").replace("station_list_file=''", "station_list_file='stations.txt'")
folder=args.work/'external';folder.mkdir(exist_ok=True);(folder/'stations.txt').write_text(rows)
assert histories(run('external',external))==both
for mode,changes in [('xyz','station_xyz_index=T, station_number_in_filename=F'),('grid','station_number_in_filename=F')]:
 text=config('block_2',False).replace('station_number_in_filename=T',changes)
 h=histories(run(mode,text));assert len(h)==3 and all('_block' not in k for k in h)
folder=args.work/'directory-failure';folder.mkdir(exist_ok=True);(folder/'notadir').write_text('file')
run('directory-failure',config().replace("station_file_directory='receivers'","station_file_directory='notadir/receivers'"),failure='CFG-RECEIVER-DIRECTORY')
spaced=histories(run('spaced-dir',config('block_2',False).replace("station_file_directory='receivers'","station_file_directory='receivers/spaced dir'")))
assert len(spaced)==3
run('empty',config().replace(rows,'!---begin:station_list---\n!---end:station_list---\n'))
run('disabled',config().replace('output_seismograms=T','output_seismograms=F').replace(rows,''))
run('bad-row',config().replace('2 0.6 7.5 0','2 not_a_coordinate 7.5 0'),failure='CFG-STATION-002',preflight=True)
run('missing-ext',config().replace("station_list='infile'","station_list='extfile'"),failure='station_list_file is required',preflight=True)
run('number-without-id',config().replace('station_number_in_list=T','station_number_in_list=F'),failure='requires station_number_in_list',preflight=True)
run('case-normalization',config(' BLOCK_2 ',False))
exact=run('exact',config('block_2',False).replace('output_exact_moment=F','output_exact_moment=T'))
files=list((exact/'receivers').glob('*.dat'))
assert len(files)==6 and sum('_exact_' in p.name for p in files)==3
assert all('# t vx vy vz' in p.read_text() for p in files)
reorder=histories(run('reorder',config('block_2',False).replace("'t vx vy vz'","'vz t vy vx'")))
for i,b in [(1,1),(2,2),(3,2)]:
 for a,ref in zip(reorder[f'case_station-{i}.dat'],both[f'case_station-{i}_block{b}.dat']): assert a==[ref[3],ref[0],ref[2],ref[1]]
curved=(root/'test_problems/test_rup_curv_fric.in').read_text()
curved=curved[:curved.index('&output_list')]
curved=curved.replace("name = 'data/test_rup_curv_fric'","name = 'case'").replace('t_final = 0.5d0','t_final = 0.1d0').replace('w_fault = T','w_fault = F')
# At y=7.5,z=0 the analytical curved interface has x=0.1 km.
curved+='&output_list output_exact_moment=F, output_fault_topo=F /\n'
curved+=config('block_2',False)[len(base):].replace('3 0 7.5 0','3 0.1 7.5 0')
h=histories(run('curved',curved,2));assert len(h)==3
single=(root/'inputfile/test_anelastic_cQ_cg_t_nb1.in').read_text()
single=re.sub(r'&output_list.*?/', '&output_list output_exact_moment=F /',single,flags=re.S)
single+=config('block_2')[len(base):]
run('unavailable-block',single,failure='block_2 requires two blocks',preflight=True)
# Named rows retain the same numerical samples and use the supplied basename.
named=config('block_2',False).replace('station_number_in_list=T','station_number_in_list=F, station_name_in_list=T').replace('station_number_in_filename=T','station_number_in_filename=F')
labels=['body030st-120dp000','body030st120dp000','interface-point']
for i,label in enumerate(labels,1): named=named.replace(f'\n{i} ',f'\n{label} ')
h=histories(run('named',named))
assert set(h)=={label+'.dat' for label in labels}
for i,b in [(1,1),(2,2),(3,2)]: assert h[labels[i-1]+'.dat']==both[f'case_station-{i}_block{b}.dat']
for p in (args.work/'named'/'receivers').glob('*.dat'): assert '# station_name: '+p.stem in p.read_text()
for n in [2,4]:
 h=histories(run(f'named-mpi-{n}',named,n))
 ref=histories(args.work/f'legacy-mpi-{n}')
 for i,b in [(1,1),(2,2),(3,2)]: assert h[labels[i-1]+'.dat']==ref[f'case_station-{i}_block{b}.dat']
nboth=named.replace("interface_stations='block_2'","interface_stations='both'").replace('append_block=F','append_block=T')
h=histories(run('named-both',nboth));assert len(h)==4
assert 'interface-point_block1.dat' in h and 'interface-point_block2.dat' in h
folder=args.work/'named-external';folder.mkdir(exist_ok=True)
(folder/'stations.txt').write_text(named[named.index('!---begin:station_list---'):])
ext=named.replace("station_list='infile'","station_list='extfile'").replace("station_list_file=''","station_list_file='stations.txt'")
assert histories(run('named-external',ext))==histories(args.work/'named')
run('named-duplicate',named.replace(labels[1],labels[0]),failure='Duplicate station name',preflight=True)
run('named-case-duplicate',named.replace(labels[1],labels[0].upper()),failure='Duplicate station name',preflight=True)
for label,bad in [('slash','../outside'),('long','a'*65),('space',"'two words'")]:
 run('named-invalid-'+label,named.replace(labels[0],bad),failure='CFG-STATION-NAME',preflight=True)
run('named-incomplete',named.replace('-0.6 7.5 0','/ 7.5 0'),failure='CFG-STATION-002',preflight=True)
run('named-flag-conflict',named.replace('station_number_in_filename=F','station_number_in_filename=T'),failure='station_name_in_list requires',preflight=True)
run('named-collision',named.replace("interface_stations='block_2'","interface_stations='both'"),2,failure='CFG-RECEIVER-COLLISION')
run('named-exact-collision',named.replace('output_exact_moment=F','output_exact_moment=T').replace(labels[1],labels[0]+'_exact'),failure='CFG-RECEIVER-COLLISION')
tpv5=(root/'test_problems/tpv5.in').read_text().replace('nqrs = 201, 201, 401','nqrs = 41, 41, 41').replace('t_final =  11d0','t_final =  0.1d0')
tpv5=tpv5[:tpv5.index('&output_list')]+"&output_list output_exact_moment=F, output_fault_topo=F /\n"
tpv5+=named[len(base):named.index('!---begin:station_list---')]
tpv5+='!---begin:station_list---\nbody030st-120dp000 3.0 0.0 -12.0\nbody030st120dp000 3.0 0.0 12.0\nbody030st-120dp075 3.0 7.5 -12.0\nbody030st120dp075 3.0 7.5 12.0\n!---end:station_list---\n'
h=histories(run('named-tpv5',tpv5,2))
assert set(h)=={f'body030st{strike}dp{depth}.dat' for strike in ['-120','120'] for depth in ['000','075']}
locked=(root/'inputfile/test_anelastic_cQ_cg_t_traditional6.in').read_text()
locked=re.sub(r"name='[^']+'", "name='case'",locked)
locked=re.sub(r'&output_list.*?/', '&output_list\n output_exact_moment=F, output_fault_topo=F, output_fields_block1=F, output_fields_block2=F, stride_fields=1\n/',locked,flags=re.S)
locked+=config()[len(base):].replace('-0.6 7.5 0','-0.08 0 0.04').replace('0.6 7.5 0','0.08 0 0.04').replace('3 0 7.5 0','3 0 0 0.04')
reference=None
for n in [1,2,3,4]:
 h=histories(run(f'locked-mpi-{n}',locked,n))
 assert len(h)==4
 if reference is None: reference=h
 else:
  assert h.keys()==reference.keys()
  for k in h:
   assert len(h[k])==len(reference[k])
   for a,b in zip(h[k],reference[k]):
    assert all(abs(x-y)<=1e-10*max(1,abs(y)) for x,y in zip(a,b)),(n,k,a,b)
assert max(abs(v) for rows in reference.values() for row in rows for v in row[1:])>1e-12
(args.work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
