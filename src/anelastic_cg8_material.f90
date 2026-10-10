module anelastic_cg8_material
 use common,only:wp
 use datatypes,only:block_material,block_grid_t
 use anelastic_cg8_model
 use anelastic_cg8_types
 use anelastic_cg8_layout
 use decomposition_safety,only:get_stencil_requirements,stencil_requirements_t
 use mpi
 use, intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::init_cg8_properties,apply_cg8_strain,scale_cg8_rates,update_cg8_memory,destroy_cg8_properties,cg8_stats
contains
 subroutine error(message,routine)
  use diagnostics,only:fatal_local
  character(*),intent(in)::message,routine
  call fatal_local('RUN-CG8-INIT',message,routine)
 end subroutine

 subroutine init_cg8_properties(m,g,p,id,frequency,fd_type,order,pml_lower,pml_upper,npml)
  type(block_material),intent(inout),target::m
  type(block_grid_t),intent(in)::g
  type(cg8_parameters),intent(in)::p
  integer,intent(in)::id,order,npml
  logical,intent(in)::frequency,pml_lower(3),pml_upper(3)
  character(*),intent(in)::fd_type
  type(cg8_state),pointer::s
  type(stencil_requirements_t)::req
  real(wp)::v(3),lo(3),hi(3),global_lo(3),global_hi(3),rho,cs,cp,bulk,shear,spacing(3),min_cs,freq
  real(wp)::geom_lo(10),geom_hi(10),geom(10),global_geom_lo(10),global_geom_hi(10)
  integer::x,y,z,i,j,k,ncg,nfull,status,guard,ierr,shape(3)
  integer(kind=8)::bytes,global_bytes,counts(2),global_counts(2)
  complex(wp)::rs
  character(len=256)::message
  if(allocated(m%cq8_cg).or.allocated(m%fq8_cg)) call error('CG8 already initialized','init_cg8_properties')
  if(frequency) then
    allocate(m%fq8_cg);s=>m%fq8_cg
  else
    allocate(m%cq8_cg);s=>m%cq8_cg
  endif
  req=get_stencil_requirements(fd_type,order);guard=cg8_minimum_buffer(req%boundary_width,req%halo_width)
  if(p%buffer_layers>=0) then
    if(p%buffer_layers<guard) call error('CG8 buffer_layers is below the stencil minimum','init_cg8_properties')
    guard=2*((p%buffer_layers+1)/2)
  endif
  s%origin=p%pattern_origin;s%exclude_lower=guard;s%exclude_upper=guard
  where(pml_lower) s%exclude_lower=guard+npml
  where(pml_upper) s%exclude_upper=guard+npml
  shape=[g%C%nq,g%C%nr,g%C%ns]
  lo=huge(1.0_wp);hi=-huge(1.0_wp)
  geom_lo=huge(1.0_wp);geom_hi=-huge(1.0_wp)
  do z=g%C%ms,g%C%ps;do y=g%C%mr,g%C%pr;do x=g%C%mq,g%C%pq
    v=m%m(x,y,z,1:3)
    if(.not.all(ieee_is_finite(v))) call error('CG8 material must be finite','init_cg8_properties')
    lo=min(lo,v);hi=max(hi,v)
    geom=[g%metricx(x,y,z,:),g%metricy(x,y,z,:),g%metricz(x,y,z,:),g%J(x,y,z)]
    if(.not.all(ieee_is_finite(geom))) call error('CG8 grid metrics must be finite','init_cg8_properties')
    geom_lo=min(geom_lo,geom);geom_hi=max(geom_hi,geom)
  enddo;enddo;enddo
  call MPI_Allreduce(lo,global_lo,3,MPI_DOUBLE_PRECISION,MPI_MIN,g%C%comm,ierr)
  call MPI_Allreduce(hi,global_hi,3,MPI_DOUBLE_PRECISION,MPI_MAX,g%C%comm,ierr)
  if(any(abs(global_hi-global_lo)>1.0e-12_wp*max(1.0_wp,abs(global_hi)))) &
    call error('CG8 currently requires uniform material within each Cartesian block','init_cg8_properties')
  call MPI_Allreduce(geom_lo,global_geom_lo,10,MPI_DOUBLE_PRECISION,MPI_MIN,g%C%comm,ierr)
  call MPI_Allreduce(geom_hi,global_geom_hi,10,MPI_DOUBLE_PRECISION,MPI_MAX,g%C%comm,ierr)
  if(any(abs(global_geom_hi-global_geom_lo)>1.0e-10_wp*max(1.0_wp,abs(global_geom_hi))).or. &
     any(abs(global_geom_hi([2,3,4,6,7,8]))>1.0e-10_wp).or.any(global_geom_lo([1,5,9,10])<=0)) &
    call error('CG8 requires an axis-aligned uniform Cartesian metric','init_cg8_properties')
  rho=global_lo(3);shear=global_lo(2);bulk=global_lo(1)+2*shear/3
  if(rho<=0.or.shear<=0.or.bulk<=0) call error('CG8 material requires positive moduli and density','init_cg8_properties')
  cs=sqrt(shear/rho);cp=sqrt((bulk+4*shear/3)/rho)
  if(cp/cs<sqrt(3.0_wp)-1.0e-10_wp.or.cp/cs>2.0_wp+1.0e-10_wp) &
    call error('CG8 currently validates sqrt(3) <= Vp/Vs <= 2','init_cg8_properties')
  call build_cg8_coefficients(p,id,rho,cs,cp,.true.,s%coarse,status,message)
  if(status/=0) call error(trim(message),'init_cg8_properties')
  call build_cg8_coefficients(p,id,rho,cs,cp,.false.,s%full,status,message)
  if(status/=0) call error(trim(message),'init_cg8_properties')
  ! Cartesian metric spacing is checked below using inverse mapping metrics.
  spacing=[g%hq/abs(g%metricx(g%C%mq,g%C%mr,g%C%ms,1)), &
           g%hr/abs(g%metricy(g%C%mq,g%C%mr,g%C%ms,2)), &
           g%hs/abs(g%metricz(g%C%mq,g%C%mr,g%C%ms,3))]
  if(minval(spacing)/maxval(spacing)<0.5_wp-1.0e-10_wp) &
    call error('CG8 currently validates grid-spacing aspect ratios between 0.5 and 1','init_cg8_properties')
  min_cs=huge(1.0_wp)
  do i=1,257
    freq=p%fmin*(p%fmax/p%fmin)**(real(i-1,wp)/256)
    rs=s%coarse%shear*cg8_response(freq,s%coarse%tau,s%coarse%shear_strength,.true.)
    min_cs=min(min_cs,1/real(sqrt(rho/rs),wp))
    rs=s%full%shear*cg8_response(freq,s%full%tau,s%full%shear_strength,.false.)
    min_cs=min(min_cs,1/real(sqrt(rho/rs),wp))
  enddo
  if(min_cs/(p%fmax*maxval(spacing))<16.1_wp) &
    call error('CG8 fit band requires at least 16.1 points per shortest S wavelength','init_cg8_properties')
  allocate(s%index(g%C%mq:g%C%pq,g%C%mr:g%C%pr,g%C%ms:g%C%ps));ncg=0;nfull=0
  do z=g%C%ms,g%C%ps;do y=g%C%mr,g%C%pr;do x=g%C%mq,g%C%pq
    if(cg8_cell_eligible([x,y,z],s%origin,shape,s%exclude_lower,s%exclude_upper)) then
      ncg=ncg+1;s%index(x,y,z)=ncg
    else
      nfull=nfull+1;s%index(x,y,z)=-nfull
    endif
  enddo;enddo;enddo
  counts=[int(ncg,8),int(nfull,8)]
  call MPI_Allreduce(counts,global_counts,2,MPI_INTEGER8,MPI_SUM,g%C%comm,ierr)
  if(global_counts(1)==0) call error('CG8 domain/buffers leave no complete coarse interior cell','init_cg8_properties')
  allocate(s%eta_cg(6,ncg),s%deta_cg(6,ncg),s%eta_full(6,8,nfull),s%deta_full(6,8,nfull))
  s%eta_cg=0;s%deta_cg=0;s%eta_full=0;s%deta_full=0
  ! Material ghosts use the same global mask; memory values remain owned-node only.
  do z=lbound(m%m,3),ubound(m%m,3);do y=lbound(m%m,2),ubound(m%m,2);do x=lbound(m%m,1),ubound(m%m,1)
    if(cg8_cell_eligible([x,y,z],s%origin,shape,s%exclude_lower,s%exclude_upper)) then
      bulk=s%coarse%bulk;shear=s%coarse%shear
    else
      bulk=s%full%bulk;shear=s%full%shear
    endif
    m%m(x,y,z,1)=bulk-2*shear/3;m%m(x,y,z,2)=shear
  enddo;enddo;enddo
  bytes=96_8*(int(ncg,8)+8_8*int(nfull,8))
  call MPI_Allreduce(bytes,global_bytes,1,MPI_INTEGER8,MPI_SUM,g%C%comm,ierr)
  if(g%C%rank==0) then
    write(*,'(A,I0,A,2I12)') 'CG8 block ',id,': coarse/full nodes=',global_counts
    write(*,'(A,I0)') 'CG8 memory-variable bytes=',global_bytes
    write(*,'(A,2ES12.4)') 'CG8 coarse/full material fit errors=',s%coarse%max_error,s%full%max_error
  endif
  nullify(s)
 end subroutine

 subroutine apply_state(s,x,y,z,dx,dy,dz,rate)
  type(cg8_state),intent(inout)::s
  integer,intent(in)::x,y,z
  real(wp),intent(in)::dx(:),dy(:),dz(:)
  real(wp),intent(inout)::rate(:)
  integer::idx,k
  real(wp)::strain(6),tr,forcing(6),bulk,shear
  tr=dx(1)+dy(2)+dz(3)
  strain=[dx(1)-tr/3,dy(2)-tr/3,dz(3)-tr/3,(dy(1)+dx(2))/2,(dz(1)+dx(3))/2,(dz(2)+dy(3))/2]
  idx=s%index(x,y,z)
  if(idx>0) then
    k=cg8_mechanism([x,y,z],s%origin)
    bulk=s%coarse%bulk*s%coarse%bulk_strength(k);shear=s%coarse%shear*s%coarse%shear_strength(k)
    forcing=2*shear*strain;forcing(1:3)=forcing(1:3)+bulk*tr
    rate(4:9)=rate(4:9)-s%eta_cg(:,idx)
    s%deta_cg(:,idx)=s%deta_cg(:,idx)+(forcing-s%eta_cg(:,idx))/s%coarse%tau(k)
  else
    idx=-idx;rate(4:9)=rate(4:9)-sum(s%eta_full(:,:,idx),dim=2)
    do k=1,8
      bulk=s%full%bulk*s%full%bulk_strength(k);shear=s%full%shear*s%full%shear_strength(k)
      forcing=2*shear*strain;forcing(1:3)=forcing(1:3)+bulk*tr
      s%deta_full(:,k,idx)=s%deta_full(:,k,idx)+(forcing-s%eta_full(:,k,idx))/s%full%tau(k)
    enddo
  endif
 end subroutine

 subroutine apply_cg8_strain(m,x,y,z,dx,dy,dz,rate)
  type(block_material),intent(inout)::m
  integer,intent(in)::x,y,z
  real(wp),intent(in)::dx(:),dy(:),dz(:)
  real(wp),intent(inout)::rate(:)
  if(allocated(m%cq8_cg)) call apply_state(m%cq8_cg,x,y,z,dx,dy,dz,rate)
  if(allocated(m%fq8_cg)) call apply_state(m%fq8_cg,x,y,z,dx,dy,dz,rate)
 end subroutine

 subroutine scale_cg8_rates(m,a)
  type(block_material),intent(inout)::m
  real(wp),intent(in)::a
  if(allocated(m%cq8_cg)) then
    m%cq8_cg%deta_cg=a*m%cq8_cg%deta_cg;m%cq8_cg%deta_full=a*m%cq8_cg%deta_full
  endif
  if(allocated(m%fq8_cg)) then
    m%fq8_cg%deta_cg=a*m%fq8_cg%deta_cg;m%fq8_cg%deta_full=a*m%fq8_cg%deta_full
  endif
 end subroutine

 subroutine update_cg8_memory(m,dt)
  type(block_material),intent(inout)::m
  real(wp),intent(in)::dt
  if(allocated(m%cq8_cg)) then
    m%cq8_cg%eta_cg=m%cq8_cg%eta_cg+dt*m%cq8_cg%deta_cg
    m%cq8_cg%eta_full=m%cq8_cg%eta_full+dt*m%cq8_cg%deta_full
  endif
  if(allocated(m%fq8_cg)) then
    m%fq8_cg%eta_cg=m%fq8_cg%eta_cg+dt*m%fq8_cg%deta_cg
    m%fq8_cg%eta_full=m%fq8_cg%eta_full+dt*m%fq8_cg%deta_full
  endif
 end subroutine

 subroutine cg8_stats(m,maximum,finite)
  type(block_material),intent(in)::m
  real(wp),intent(out)::maximum
  logical,intent(out)::finite
  maximum=0;finite=.true.
  if(allocated(m%cq8_cg)) call state_stats(m%cq8_cg,maximum,finite)
  if(allocated(m%fq8_cg)) call state_stats(m%fq8_cg,maximum,finite)
 end subroutine
 subroutine state_stats(s,maximum,finite)
  type(cg8_state),intent(in)::s
  real(wp),intent(out)::maximum
  logical,intent(out)::finite
  maximum=0
  if(size(s%eta_cg)>0) maximum=max(maximum,maxval(abs(s%eta_cg)))
  if(size(s%eta_full)>0) maximum=max(maximum,maxval(abs(s%eta_full)))
  finite=all(ieee_is_finite(s%eta_cg)).and.all(ieee_is_finite(s%eta_full)).and. &
    all(ieee_is_finite(s%deta_cg)).and.all(ieee_is_finite(s%deta_full))
 end subroutine
 subroutine destroy_cg8_properties(m)
  type(block_material),intent(inout)::m
  if(allocated(m%cq8_cg)) deallocate(m%cq8_cg)
  if(allocated(m%fq8_cg)) deallocate(m%fq8_cg)
 end subroutine
end module
