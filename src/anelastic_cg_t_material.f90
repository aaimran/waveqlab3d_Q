module anelastic_cg_t_material
 use common,only:wp
 use datatypes,only:block_material,block_grid_t
 use anelastic_cg8_model
 use anelastic_cg_t_model,only:build_cgt_coefficients
 use anelastic_cg_t_types
 use anelastic_cg8_layout,only:cg8_minimum_buffer,cg8_mechanism
 use anelastic_cg_t_comm,only:init_cgt_layout,exchange_cgt_stage,cgt_fail
 use decomposition_safety,only:get_stencil_requirements,stencil_requirements_t
 use mpi
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::init_cgt_properties,apply_cgt_strain,begin_cgt_stage,finish_cgt_stage
 public::scale_cgt_rates,update_cgt_memory,destroy_cgt_properties,cgt_stats
contains
 subroutine error(message,routine)
  character(*),intent(in)::message,routine
  call cgt_fail(trim(routine)//': '//trim(message))
 end subroutine
 subroutine init_cgt_properties(m,g,p,id,frequency,fd_type,order,pml_lower,pml_upper,npml)
  type(block_material),intent(inout),target::m
  type(block_grid_t),intent(in)::g
  type(cg8_parameters),intent(in)::p
  integer,intent(in)::id,order,npml
  logical,intent(in)::frequency,pml_lower(3),pml_upper(3)
  character(*),intent(in)::fd_type
  type(cgt_state),pointer::s
  type(stencil_requirements_t)::req
  real(wp)::v(3),lo(3),hi(3),global_lo(3),global_hi(3),rho,cs,cp,bulk,shear,spacing(3),min_cs
  real(wp)::geom_lo(10),geom_hi(10),geom(10),global_geom_lo(10),global_geom_hi(10)
  integer::x,y,z,i,status,guard,ierr
  integer(kind=8)::bytes,global_bytes,map_bytes,global_map_bytes
  character(len=256)::message
  if(allocated(m%cq_cg_t).or.allocated(m%fq_cg_t)) call error('CG-T already initialized','init_cgt_properties')
  if(frequency) then
    allocate(m%fq_cg_t);s=>m%fq_cg_t
  else
    allocate(m%cq_cg_t);s=>m%cq_cg_t
  endif
  req=get_stencil_requirements(fd_type,order);guard=cg8_minimum_buffer(req%boundary_width,req%halo_width)
  if(p%buffer_layers>=0) then
    if(p%buffer_layers<guard) call error('CG-T buffer_layers is below the stencil minimum','init_cgt_properties')
    guard=2*((p%buffer_layers+1)/2)
  endif
  s%origin=p%pattern_origin;s%exclude_lower=guard;s%exclude_upper=guard
  where(pml_lower) s%exclude_lower=guard+npml
  where(pml_upper) s%exclude_upper=guard+npml
  lo=huge(1.0_wp);hi=-huge(1.0_wp)
  geom_lo=huge(1.0_wp);geom_hi=-huge(1.0_wp)
  do z=g%C%ms,g%C%ps;do y=g%C%mr,g%C%pr;do x=g%C%mq,g%C%pq
    v=m%m(x,y,z,1:3)
    if(.not.all(ieee_is_finite(v))) call error('CG-T material must be finite','init_cgt_properties')
    lo=min(lo,v);hi=max(hi,v)
    geom=[g%metricx(x,y,z,:),g%metricy(x,y,z,:),g%metricz(x,y,z,:),g%J(x,y,z)]
    if(.not.all(ieee_is_finite(geom))) call error('CG-T grid metrics must be finite','init_cgt_properties')
    geom_lo=min(geom_lo,geom);geom_hi=max(geom_hi,geom)
  enddo;enddo;enddo
  call MPI_Allreduce(lo,global_lo,3,MPI_DOUBLE_PRECISION,MPI_MIN,g%C%comm,ierr)
  call MPI_Allreduce(hi,global_hi,3,MPI_DOUBLE_PRECISION,MPI_MAX,g%C%comm,ierr)
  if(any(abs(global_hi-global_lo)>1.0e-12_wp*max(1.0_wp,abs(global_hi)))) &
    call error('CG-T currently requires uniform material within each Cartesian block','init_cgt_properties')
  call MPI_Allreduce(geom_lo,global_geom_lo,10,MPI_DOUBLE_PRECISION,MPI_MIN,g%C%comm,ierr)
  call MPI_Allreduce(geom_hi,global_geom_hi,10,MPI_DOUBLE_PRECISION,MPI_MAX,g%C%comm,ierr)
  if(any(abs(global_geom_hi-global_geom_lo)>1.0e-10_wp*max(1.0_wp,abs(global_geom_hi))).or. &
     any(abs(global_geom_hi([2,3,4,6,7,8]))>1.0e-10_wp).or.any(global_geom_lo([1,5,9,10])<=0)) &
    call error('CG-T requires an axis-aligned uniform Cartesian metric','init_cgt_properties')
  rho=global_lo(3);shear=global_lo(2);bulk=global_lo(1)+2*shear/3
  if(rho<=0.or.shear<=0.or.bulk<=0) call error('CG-T material requires positive moduli and density','init_cgt_properties')
  cs=sqrt(shear/rho);cp=sqrt((bulk+4*shear/3)/rho)
  if(cp/cs<sqrt(3.0_wp)-1.0e-10_wp.or.cp/cs>2.0_wp+1.0e-10_wp) &
    call error('CG-T currently validates sqrt(3) <= Vp/Vs <= 2','init_cgt_properties')
  call build_cgt_coefficients(p,id,rho,cs,cp,s%coeff,status,message)
  if(status/=0) call error(trim(message),'init_cgt_properties')
  ! Cartesian metric spacing is checked below using inverse mapping metrics.
  spacing=[g%hq/abs(g%metricx(g%C%mq,g%C%mr,g%C%ms,1)), &
           g%hr/abs(g%metricy(g%C%mq,g%C%mr,g%C%ms,2)), &
           g%hs/abs(g%metricz(g%C%mq,g%C%mr,g%C%ms,3))]
  if(minval(spacing)/maxval(spacing)<0.5_wp-1.0e-10_wp) &
    call error('CG-T currently validates grid-spacing aspect ratios between 0.5 and 1','init_cgt_properties')
  ! A passive additive shear modulus has phase speed >= its relaxed speed.
  ! This bound avoids missing an extremum between sampled frequencies.
  min_cs=sqrt(s%coeff%shear*(1-sum(s%coeff%shear_strength))/rho)
  if(min_cs/(p%fmax*maxval(spacing))<24.1_wp) &
    call error('CG-T fit band requires at least 24.1 points per shortest S wavelength','init_cgt_properties')
  if(any(s%coeff%bulk_strength<0).or.any(s%coeff%shear_strength<0).or. &
     sum(s%coeff%bulk_strength)>=1.or.sum(s%coeff%shear_strength)>=1) &
    call error('Nonpassive additive coefficient spectrum','init_cgt_properties')
  call init_cgt_layout(s,g%C)
  m%m(:,:,:,1)=s%coeff%bulk-2*s%coeff%shear/3
  m%m(:,:,:,2)=s%coeff%shear
  ! Uniform loaded density also defines physical/ghost values for this block.
  m%m(:,:,:,3)=rho
  bytes=size(s%strain,kind=8)*8_8+size(s%feedback,kind=8)*8_8
  do i=1,size(s%peers)
    bytes=bytes+8_8*(int(size(s%peers(i)%strain_send),8)+int(size(s%peers(i)%strain_recv),8)+ &
      int(size(s%peers(i)%stress_send),8)+int(size(s%peers(i)%stress_recv),8))
  enddo
  call MPI_Allreduce(bytes,global_bytes,1,MPI_INTEGER8,MPI_SUM,g%C%comm,ierr)
  if(ierr/=MPI_SUCCESS) call cgt_fail('MPI workspace accounting failed')
  if(g%C%rank==0) then
    write(*,'(A,I0)') 'CG-T strain/feedback/communication workspace bytes=',global_bytes
    write(*,'(A,ES12.4)') 'CG-T additive material fit error=',s%coeff%max_error
  endif
  map_bytes=int(storage_size(0)/8,8)*(size(s%index,kind=8)+size(s%first,kind=8)+ &
    size(s%owner,kind=8)+size(s%owned_index,kind=8))+ &
    int(storage_size(.false.)/8,8)*(size(s%seen,kind=8)+size(s%full_seen,kind=8))
  do i=1,size(s%peers)
    map_bytes=map_bytes+int(storage_size(0)/8,8)*(1_8+size(s%peers(i)%send_cell,kind=8)+ &
      size(s%peers(i)%send_slot,kind=8)+size(s%peers(i)%recv_cell,kind=8)+ &
      size(s%peers(i)%recv_slot,kind=8)+size(s%peers(i)%stress_send_cell,kind=8)+ &
      size(s%peers(i)%stress_recv_cell,kind=8))
  enddo
  call MPI_Allreduce(map_bytes,global_map_bytes,1,MPI_INTEGER8,MPI_SUM,g%C%comm,ierr)
  if(ierr/=MPI_SUCCESS) call cgt_fail('MPI map accounting failed')
  if(g%C%rank==0) write(*,'(A,I0)') 'CG-T index/route/coverage array bytes=',global_map_bytes
  nullify(s)
 end subroutine
 subroutine begin_cgt_stage(m)
  type(block_material),intent(inout)::m
  if(allocated(m%cq_cg_t)) call begin_state(m%cq_cg_t)
  if(allocated(m%fq_cg_t)) call begin_state(m%fq_cg_t)
 end subroutine
 subroutine begin_state(s)
  type(cgt_state),intent(inout)::s
  s%strain=0;s%seen=.false.;s%full_seen=.false.
 end subroutine
 subroutine apply_cgt_strain(m,x,y,z,dx,dy,dz,rate)
  type(block_material),intent(inout)::m
  integer,intent(in)::x,y,z
  real(wp),intent(in)::dx(:),dy(:),dz(:)
  real(wp),intent(inout)::rate(:)
  if(allocated(m%cq_cg_t)) call apply_state(m%cq_cg_t,x,y,z,dx,dy,dz,rate)
  if(allocated(m%fq_cg_t)) call apply_state(m%fq_cg_t,x,y,z,dx,dy,dz,rate)
 end subroutine
 subroutine tensor_forcing(c,k,strain,forcing)
  type(cg8_coefficients),intent(in)::c
  integer,intent(in)::k
  real(wp),intent(in)::strain(6)
  real(wp),intent(out)::forcing(6)
  real(wp)::tr,bulk,shear
  tr=sum(strain(1:3));bulk=c%bulk*c%bulk_strength(k);shear=c%shear*c%shear_strength(k)
  forcing=2*shear*strain;forcing(1:3)=forcing(1:3)+(bulk-2*shear/3)*tr
 end subroutine
 subroutine apply_state(s,x,y,z,dx,dy,dz,rate)
  type(cgt_state),intent(inout)::s
  integer,intent(in)::x,y,z
  real(wp),intent(in)::dx(:),dy(:),dz(:)
  real(wp),intent(inout)::rate(:)
  real(wp)::strain(6),forcing(6)
  integer::j,k,slot
  strain=[dx(1),dy(2),dz(3),(dy(1)+dx(2))/2,(dz(1)+dx(3))/2,(dz(2)+dy(3))/2]
  j=s%index(x,y,z)
  if(j>0) then
    slot=cg8_mechanism([x,y,z],s%origin)
    if(s%seen(slot,j)) call cgt_fail('Duplicate projected strain write in derivative regions')
    s%seen(slot,j)=.true.;s%strain(:,slot,j)=strain
  else
    j=-j
    if(s%full_seen(j)) call cgt_fail('Duplicate full-buffer memory update')
    s%full_seen(j)=.true.
    do k=1,8
      rate(4:9)=rate(4:9)-s%eta_full(:,k,j)
      call tensor_forcing(s%coeff,k,strain,forcing)
      s%deta_full(:,k,j)=s%deta_full(:,k,j)+(forcing-s%eta_full(:,k,j))/s%coeff%tau(k)
    enddo
  endif
 end subroutine
 subroutine finish_cgt_stage(m,rate)
  type(block_material),intent(inout)::m
  real(wp),intent(inout)::rate(:,:,:,:)
  if(allocated(m%cq_cg_t)) call finish_state(m%cq_cg_t,rate)
  if(allocated(m%fq_cg_t)) call finish_state(m%fq_cg_t,rate)
 end subroutine
 subroutine finish_state(s,rate)
  type(cgt_state),intent(inout)::s
  ! The caller supplies the owned-node subarray, so its indices start at one.
  real(wp),intent(inout)::rate(:,:,:,:)
  integer::j,i,k,x,y,z,mx,my,mz
  real(wp)::strain(6),forcing(6)
  if(.not.all(s%full_seen)) call cgt_fail('Missing full-buffer derivative update')
  call exchange_cgt_stage(s)
  do j=1,size(s%owner)
    i=s%owned_index(j);if(i==0) cycle
    if(.not.all(s%seen(:,j))) call cgt_fail('Incomplete eight-node cell strain')
    strain=0
    do k=1,8
      strain=strain+s%strain(:,k,j)
    enddo
    strain=strain/8
    do k=1,8
      call tensor_forcing(s%coeff,k,strain,forcing)
      s%deta_cell(:,k,i)=s%deta_cell(:,k,i)+(forcing-s%eta_cell(:,k,i))/s%coeff%tau(k)
    enddo
  enddo
  mx=lbound(s%index,1);my=lbound(s%index,2);mz=lbound(s%index,3)
  do z=lbound(s%index,3),ubound(s%index,3)
   do y=lbound(s%index,2),ubound(s%index,2)
    do x=lbound(s%index,1),ubound(s%index,1)
      j=s%index(x,y,z);if(j<=0) cycle
      rate(x-mx+1,y-my+1,z-mz+1,4:9)=rate(x-mx+1,y-my+1,z-mz+1,4:9)-s%feedback(:,j)
    enddo
   enddo
  enddo
 end subroutine
 subroutine scale_cgt_rates(m,a)
  type(block_material),intent(inout)::m
  real(wp),intent(in)::a
  if(allocated(m%cq_cg_t)) then
    m%cq_cg_t%deta_cell=a*m%cq_cg_t%deta_cell;m%cq_cg_t%deta_full=a*m%cq_cg_t%deta_full
  endif
  if(allocated(m%fq_cg_t)) then
    m%fq_cg_t%deta_cell=a*m%fq_cg_t%deta_cell;m%fq_cg_t%deta_full=a*m%fq_cg_t%deta_full
  endif
 end subroutine
 subroutine update_cgt_memory(m,dt)
  type(block_material),intent(inout)::m
  real(wp),intent(in)::dt
  if(allocated(m%cq_cg_t)) then
    m%cq_cg_t%eta_cell=m%cq_cg_t%eta_cell+dt*m%cq_cg_t%deta_cell
    m%cq_cg_t%eta_full=m%cq_cg_t%eta_full+dt*m%cq_cg_t%deta_full
  endif
  if(allocated(m%fq_cg_t)) then
    m%fq_cg_t%eta_cell=m%fq_cg_t%eta_cell+dt*m%fq_cg_t%deta_cell
    m%fq_cg_t%eta_full=m%fq_cg_t%eta_full+dt*m%fq_cg_t%deta_full
  endif
 end subroutine
 subroutine cgt_stats(m,maximum,finite)
  type(block_material),intent(in)::m
  real(wp),intent(out)::maximum
  logical,intent(out)::finite
  maximum=0;finite=.true.
  if(allocated(m%cq_cg_t)) call state_stats(m%cq_cg_t,maximum,finite)
  if(allocated(m%fq_cg_t)) call state_stats(m%fq_cg_t,maximum,finite)
 end subroutine
 subroutine state_stats(s,maximum,finite)
  type(cgt_state),intent(in)::s
  real(wp),intent(out)::maximum
  logical,intent(out)::finite
  maximum=0
  if(size(s%eta_cell)>0) maximum=max(maximum,maxval(abs(s%eta_cell)))
  if(size(s%eta_full)>0) maximum=max(maximum,maxval(abs(s%eta_full)))
  finite=all(ieee_is_finite(s%eta_cell)).and.all(ieee_is_finite(s%eta_full)).and. &
    all(ieee_is_finite(s%deta_cell)).and.all(ieee_is_finite(s%deta_full))
 end subroutine
 subroutine destroy_cgt_properties(m)
  type(block_material),intent(inout)::m
  if(allocated(m%cq_cg_t)) deallocate(m%cq_cg_t)
  if(allocated(m%fq_cg_t)) deallocate(m%fq_cg_t)
 end subroutine
end module
