program cgt_material_mpi_test
 use common,only:wp
 use mpi
 use mpi3dcomm,only:cartesian3d_t,decompose3d
 use datatypes,only:block_material
 use anelastic_cg_t_material
 use anelastic_cg_t_comm,only:init_cgt_layout
 use anelastic_cg8_layout,only:cg8_mechanism
 implicit none
 type(block_material),target::m
 type(cartesian3d_t)::c
 integer::ierr,x,y,z,j,k,comp,stage,first(3),pt(3),slot,idx,np
 real(wp),allocatable::rate(:,:,:,:)
 real(wp)::dx(3),dy(3),dz(3),e(6),avg(6),force(6),eta(6),de(6),stress(6),err,globalerr
 real(wp),parameter::aa(3)=[.5_wp,-.3_wp,.2_wp],bb(3)=[.1_wp,.03_wp,.08_wp]
 character(64)::method,argument
 call MPI_Init(ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 method='3D'
 call decompose3d(c,19,19,19,1,.false.,.false.,.false.,MPI_COMM_WORLD,method,0,0,0)
 allocate(m%cq_cg_t)
 m%cq_cg_t%origin=1;m%cq_cg_t%exclude_lower=2;m%cq_cg_t%exclude_upper=2
 call get_command_argument(1,argument)
 if(trim(argument)=='--single-cell') then
  m%cq_cg_t%exclude_lower=8;m%cq_cg_t%exclude_upper=8
 endif
 call init_cgt_layout(m%cq_cg_t,c)
 if(np==8) then
  idx=size(m%cq_cg_t%peers)
  call MPI_Allreduce(idx,j,1,MPI_INTEGER,MPI_MAX,c%comm,ierr)
  if(j/=7) error stop 'Eight-rank cell did not exercise corner routing'
 endif
 m%cq_cg_t%coeff%bulk=3;m%cq_cg_t%coeff%shear=2
 do k=1,8
  m%cq_cg_t%coeff%bulk_strength(k)=.02_wp*k/8
  m%cq_cg_t%coeff%shear_strength(k)=.01_wp*k/8
  m%cq_cg_t%coeff%tau(k)=.5_wp+.1_wp*k
 enddo
 do j=1,size(m%cq_cg_t%owner)
  idx=m%cq_cg_t%owned_index(j);if(idx==0) cycle
  do k=1,8
    call initial(m%cq_cg_t%first(:,j),k,eta)
    m%cq_cg_t%eta_cell(:,k,idx)=eta
  enddo
 enddo
 do z=c%ms,c%ps;do y=c%mr,c%pr;do x=c%mq,c%pq
  j=m%cq_cg_t%index(x,y,z);if(j>0) cycle
  do k=1,8
    call initial([x,y,z],k,eta);m%cq_cg_t%eta_full(:,k,-j)=eta
  enddo
 enddo;enddo;enddo
 m%cq_cg_t%deta_cell=.3_wp;m%cq_cg_t%deta_full=.3_wp
 allocate(rate(c%mq:c%pq,c%mr:c%pr,c%ms:c%ps,9));err=0
 do stage=1,3
  rate=1
  call scale_cgt_rates(m,aa(stage));call begin_cgt_stage(m)
  do z=c%ms,c%ps;do y=c%mr,c%pr;do x=c%mq,c%pq
    call gradients([x,y,z],dx,dy,dz,e)
    call apply_cgt_strain(m,x,y,z,dx,dy,dz,rate(x,y,z,:))
  enddo;enddo;enddo
  call finish_cgt_stage(m,rate)
  do z=c%ms,c%ps;do y=c%mr,c%pr;do x=c%mq,c%pq
    j=m%cq_cg_t%index(x,y,z)
    if(j>0) then
      first=m%cq_cg_t%first(:,j);call mean_strain(first,avg)
    else
      first=[x,y,z];call gradients(first,dx,dy,dz,avg)
    endif
    stress=0
    do k=1,8
      call reference(first,k,avg,stage,eta,de)
      stress=stress+eta
      if(j<0) err=max(err,maxval(abs(m%cq_cg_t%deta_full(:,k,-j)-de)))
    enddo
    err=max(err,maxval(abs(rate(x,y,z,4:9)-(1-stress))))
  enddo;enddo;enddo
  do j=1,size(m%cq_cg_t%owner)
    idx=m%cq_cg_t%owned_index(j);if(idx==0) cycle
    call mean_strain(m%cq_cg_t%first(:,j),avg)
    do k=1,8
      call reference(m%cq_cg_t%first(:,j),k,avg,stage,eta,de)
      err=max(err,maxval(abs(m%cq_cg_t%deta_cell(:,k,idx)-de)))
    enddo
  enddo
  call update_cgt_memory(m,bb(stage))
 enddo
 call MPI_Allreduce(err,globalerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,c%comm,ierr)
 if(globalerr>1e-12_wp) error stop 'Projected/full tensor and multistage RK reference mismatch'
 if(allocated(m%eta4Q8).or.allocated(m%eta4cQ).or.allocated(m%cq8_cg)) error stop 'Legacy state allocated'
 call destroy_cgt_properties(m);call destroy_cgt_properties(m)
 if(allocated(m%cq_cg_t)) error stop 'Projected state cleanup'
 if(c%rank==0) print *, 'CG-T full/cell tensor, nonzero memory, RK, split-cell MPI passed:',np,globalerr
 call MPI_Finalize(ierr)
contains
 subroutine gradients(point,dx,dy,dz,e)
  integer,intent(in)::point(3)
  real(wp),intent(out)::dx(3),dy(3),dz(3),e(6)
  dx=[.02_wp*point(1),.03_wp*point(2),.04_wp*point(3)]
  dy=[.05_wp*point(1),.06_wp*point(2),.07_wp*point(3)]
  dz=[.08_wp*point(1),.09_wp*point(2),.10_wp*point(3)]
  e=[dx(1),dy(2),dz(3),(dy(1)+dx(2))/2,(dz(1)+dx(3))/2,(dz(2)+dy(3))/2]
 end subroutine
 subroutine mean_strain(first,e)
  integer,intent(in)::first(3)
  real(wp),intent(out)::e(6)
  real(wp)::v(6),a(3),b(3),d(3)
  integer::i,p(3)
  e=0
  do i=0,7
    p=first+[modulo(i,2),modulo(i/2,2),i/4]
    call gradients(p,a,b,d,v);e=e+v
  enddo
  e=e/8
 end subroutine
 subroutine initial(first,k,eta)
  integer,intent(in)::first(3),k
  real(wp),intent(out)::eta(6)
  integer::i
  eta=[(.0001_wp*(sum(first)+k+i),i=1,6)]
 end subroutine
 subroutine reference(first,k,e,stage,eta,de)
  integer,intent(in)::first(3),k,stage
  real(wp),intent(in)::e(6)
  real(wp),intent(out)::eta(6),de(6)
  real(wp)::tensor(3,3),dev(3,3),force(6),tr,tau
  integer::i,t
  tensor=reshape([e(1),e(4),e(5),e(4),e(2),e(6),e(5),e(6),e(3)],[3,3])
  tr=tensor(1,1)+tensor(2,2)+tensor(3,3);dev=tensor
  do i=1,3
    dev(i,i)=dev(i,i)-tr/3
  enddo
  tensor=4*.01_wp*k/8*dev
  do i=1,3
    tensor(i,i)=tensor(i,i)+3*.02_wp*k/8*tr
  enddo
  force=[tensor(1,1),tensor(2,2),tensor(3,3),tensor(1,2),tensor(1,3),tensor(2,3)]
  tau=.5_wp+.1_wp*k;call initial(first,k,eta);de=.3_wp
  do t=1,stage
    de=aa(t)*de+(force-eta)/tau
    if(t<stage) eta=eta+bb(t)*de
  enddo
 end subroutine
end program
