program cg8_material_test
 use decomposition_safety,only:get_stencil_requirements,stencil_requirements_t
 use common,only:wp
 use datatypes,only:block_material
 use anelastic_cg8_material
 use anelastic_cg8_layout
 implicit none
 type(block_material)::m
 type(stencil_requirements_t)::req
 real(wp)::dx(3),dy(3),dz(3),rate(9),e(3,3),t(3,3),dev(3,3),forcing(6),before(6),tr
 integer::i,j,k,x,y,z,counts(8),n1,n2,nall
 logical::eligible
 allocate(m%cq8_cg)
 allocate(m%cq8_cg%index(1:2,1:1,1:1),m%cq8_cg%eta_cg(6,1),m%cq8_cg%deta_cg(6,1))
 allocate(m%cq8_cg%eta_full(6,8,1),m%cq8_cg%deta_full(6,8,1))
 m%cq8_cg%index(:,1,1)=[1,-1];m%cq8_cg%origin=1
 m%cq8_cg%coarse%bulk=3;m%cq8_cg%coarse%shear=2
 m%cq8_cg%coarse%bulk_strength=.2_wp;m%cq8_cg%coarse%shear_strength=.1_wp
 m%cq8_cg%coarse%tau=2;m%cq8_cg%eta_cg(:,1)=[(.01_wp*i,i=1,6)];m%cq8_cg%deta_cg=.3_wp
 dx=[.2_wp,.3_wp,.4_wp];dy=[.5_wp,.6_wp,.7_wp];dz=[.8_wp,.9_wp,1.0_wp]
 e(:,1)=dx;e(:,2)=dy;e(:,3)=dz;e=(e+transpose(e))/2
 tr=e(1,1)+e(2,2)+e(3,3);dev=e
 do i=1,3
  dev(i,i)=dev(i,i)-tr/3
 enddo
 t=2*2*.1_wp*dev
 do i=1,3
  t(i,i)=t(i,i)+3*.2_wp*tr
 enddo
 forcing=[t(1,1),t(2,2),t(3,3),t(1,2),t(1,3),t(2,3)]
 before=m%cq8_cg%eta_cg(:,1);rate=1
 call apply_cg8_strain(m,1,1,1,dx,dy,dz,rate)
 if(maxval(abs(rate(4:9)-(1-before)))>1e-13_wp) error stop 'coarse stress subtraction'
 if(maxval(abs(m%cq8_cg%deta_cg(:,1)-(.3_wp+(forcing-before)/2)))>1e-13_wp) &
   error stop 'coarse tensor constitutive update'
 m%cq8_cg%full=m%cq8_cg%coarse;m%cq8_cg%full%bulk_strength=.02_wp
 m%cq8_cg%full%shear_strength=.01_wp;m%cq8_cg%eta_full=.1_wp;m%cq8_cg%deta_full=0
 rate=1
 call apply_cg8_strain(m,2,1,1,dx,dy,dz,rate)
 if(maxval(abs(rate(4:9)-.2_wp))>1e-13_wp) error stop 'full-buffer stress subtraction'
 if(maxval(abs(m%cq8_cg%deta_full(:,1,1)-(.1_wp*forcing-.1_wp)/2))>1e-13_wp) &
   error stop 'full tensor constitutive update'
 before=m%cq8_cg%deta_cg(:,1)
 call scale_cg8_rates(m,.5_wp)
 if(any(m%cq8_cg%deta_cg(:,1)/=before*.5_wp)) error stop 'RK rate scaling'
 before=m%cq8_cg%eta_cg(:,1)
 call update_cg8_memory(m,.25_wp)
 if(maxval(abs(m%cq8_cg%eta_cg(:,1)-before-.25_wp*m%cq8_cg%deta_cg(:,1)))>1e-13_wp) &
   error stop 'RK memory update'
 if(allocated(m%eta4Q8).or.allocated(m%eta4Qf8).or.allocated(m%eta4cQ).or.allocated(m%eta4fQ)) &
   error stop 'legacy state allocated'
 call destroy_cg8_properties(m)
 if(allocated(m%cq8_cg).or.allocated(m%fq8_cg)) error stop 'CG8 destruction'
 call destroy_cg8_properties(m)
 req=get_stencil_requirements('upwind',6)
 if(cg8_minimum_buffer(req%boundary_width,req%halo_width)/=12) error stop 'upwind6 buffer guard'
 req=get_stencil_requirements('upwind_drp',6)
 if(cg8_minimum_buffer(req%boundary_width,req%halo_width)/=16) error stop 'DRP6 buffer guard'
 counts=0
 do z=-1,0;do y=-1,0;do x=-1,0
  k=cg8_mechanism([x,y,z],[-1,-1,-1]);counts(k)=counts(k)+1
 enddo;enddo;enddo
 if(any(counts/=1)) error stop 'parity supercell coverage'
 n1=0;n2=0;nall=0
 do z=1,41;do y=1,41;do x=1,41
  eligible=cg8_cell_eligible([x,y,z],[0,0,0],[41,41,41],[10,10,10],[10,10,10])
  if(eligible) then
    nall=nall+1
    if(x<=21) then
      n1=n1+1
    else
      n2=n2+1
    endif
  endif
 enddo;enddo;enddo
 if(nall/=n1+n2.or.nall/=20**3) error stop 'partition-independent whole cells'
 print *, 'CG8 tensor, RK, compact-state, and layout tests passed'
end program
