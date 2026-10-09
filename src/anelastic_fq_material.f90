module anelastic_fq_material

  use common, only : wp
  use datatypes, only : block_material, block_grid_t
  use anelastic_fq_model, only : fq_parameters, build_fq_coefficients, fq_max_relative_error
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: init_anelastic_fq_properties, destroy_anelastic_fq_properties
  public :: apply_anelastic_fq_strain

contains

  subroutine apply_anelastic_fq_strain(M,x,y,z,Dx,Dy,Dz,rate)
    type(block_material), intent(inout) :: M
    integer, intent(in) :: x,y,z
    real(wp), intent(in) :: Dx(:),Dy(:),Dz(:)
    real(wp), intent(inout) :: rate(:)
    integer :: i,n
    real(wp) :: tr,mu2,pm,sm,bulk
    n=M%n_mechanism_fQ
    rate(4)=rate(4)-sum(M%eta4fQ(x,y,z,1:n)); rate(5)=rate(5)-sum(M%eta5fQ(x,y,z,1:n))
    rate(6)=rate(6)-sum(M%eta6fQ(x,y,z,1:n)); rate(7)=rate(7)-sum(M%eta7fQ(x,y,z,1:n))
    rate(8)=rate(8)-sum(M%eta8fQ(x,y,z,1:n)); rate(9)=rate(9)-sum(M%eta9fQ(x,y,z,1:n))
    tr=Dx(1)+Dy(2)+Dz(3); mu2=2.0_wp*M%M(x,y,z,2)
    do i=1,n
       sm=M%strength_s_fQ(i)
       pm=M%strength_p_fQ(i)
       bulk=(M%M(x,y,z,1)+mu2)*pm-mu2*sm
       M%Deta4fQ(x,y,z,i)=M%Deta4fQ(x,y,z,i)+(mu2*sm*Dx(1)+bulk*tr-M%eta4fQ(x,y,z,i))/M%tau_fQ(i)
       M%Deta5fQ(x,y,z,i)=M%Deta5fQ(x,y,z,i)+(mu2*sm*Dy(2)+bulk*tr-M%eta5fQ(x,y,z,i))/M%tau_fQ(i)
       M%Deta6fQ(x,y,z,i)=M%Deta6fQ(x,y,z,i)+(mu2*sm*Dz(3)+bulk*tr-M%eta6fQ(x,y,z,i))/M%tau_fQ(i)
       M%Deta7fQ(x,y,z,i)=M%Deta7fQ(x,y,z,i)+(M%M(x,y,z,2)*sm*(Dy(1)+Dx(2))-M%eta7fQ(x,y,z,i))/M%tau_fQ(i)
       M%Deta8fQ(x,y,z,i)=M%Deta8fQ(x,y,z,i)+(M%M(x,y,z,2)*sm*(Dz(1)+Dx(3))-M%eta8fQ(x,y,z,i))/M%tau_fQ(i)
       M%Deta9fQ(x,y,z,i)=M%Deta9fQ(x,y,z,i)+(M%M(x,y,z,2)*sm*(Dz(2)+Dy(3))-M%eta9fQ(x,y,z,i))/M%tau_fQ(i)
    end do
  end subroutine apply_anelastic_fq_strain

  subroutine init_anelastic_fq_properties(M,G,parameters,block_id)
    use mpi3dcomm, only : allocate_array_body
    use mpi3dbasic, only : error
    type(block_material), intent(inout) :: M
    type(block_grid_t), intent(in) :: G
    type(fq_parameters), intent(in) :: parameters
    integer, intent(in) :: block_id
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp) :: wref,val_s,val_p,vs,vp,mu_s,mu_p,max_s,max_p
    integer :: i,j,k,l,n,status
    character(len=256) :: fit_message

    if (M%anelastic_fQ .or. allocated(M%eta4fQ)) &
         call error('anelastic-fQ material is already initialized','init_anelastic_fq_properties')
    n=parameters%n_mechanisms
    M%anelastic_fQ=.true.; M%n_mechanism_fQ=n
    M%fref_fQ=parameters%fref; M%fmin_fQ=parameters%fmin; M%fmax_fQ=parameters%fmax
    M%Qs0_fQ=parameters%Qs0(block_id); M%Qp0_fQ=parameters%Qp0(block_id)
    M%coefficient_policy_fQ=parameters%coefficient_policy
    M%gamma_fQ=parameters%gamma; M%f_transition_fQ=parameters%f_transition
    allocate(M%tau_fQ(n),M%strength_s_fQ(n),M%strength_p_fQ(n))
    call build_fq_coefficients(parameters,block_id,M%tau_fQ,M%strength_s_fQ,M%strength_p_fQ,status,fit_message)
    if(status /= 0) call error(trim(fit_message),'init_anelastic_fq_properties')
    if (any(M%tau_fQ <= 0.0_wp) .or. any(M%strength_s_fQ < 0.0_wp) .or. &
        any(M%strength_p_fQ < 0.0_wp)) &
         call error('anelastic-fQ produced invalid coefficients','init_anelastic_fq_properties')

    call fq_max_relative_error(M%Qs0_fQ,parameters,M%tau_fQ,M%strength_s_fQ,max_s)
    call fq_max_relative_error(M%Qp0_fQ,parameters,M%tau_fQ,M%strength_p_fQ,max_p)
    if (max(max_s,max_p) > parameters%max_fit_error) then
       write(fit_message,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'anelastic-fQ fitted response exceeds max_fit_error: S=',max_s, &
            ', P=',max_p,', limit=',parameters%max_fit_error
       call error(trim(fit_message),'init_anelastic_fq_properties')
    end if

    call allocate_array_body(M%eta4fQ,G%C,n,ghost_nodes=.true.); M%eta4fQ=0.0_wp
    call allocate_array_body(M%Deta4fQ,G%C,n,ghost_nodes=.true.); M%Deta4fQ=0.0_wp
    call allocate_array_body(M%eta5fQ,G%C,n,ghost_nodes=.true.); M%eta5fQ=0.0_wp
    call allocate_array_body(M%Deta5fQ,G%C,n,ghost_nodes=.true.); M%Deta5fQ=0.0_wp
    call allocate_array_body(M%eta6fQ,G%C,n,ghost_nodes=.true.); M%eta6fQ=0.0_wp
    call allocate_array_body(M%Deta6fQ,G%C,n,ghost_nodes=.true.); M%Deta6fQ=0.0_wp
    call allocate_array_body(M%eta7fQ,G%C,n,ghost_nodes=.true.); M%eta7fQ=0.0_wp
    call allocate_array_body(M%Deta7fQ,G%C,n,ghost_nodes=.true.); M%Deta7fQ=0.0_wp
    call allocate_array_body(M%eta8fQ,G%C,n,ghost_nodes=.true.); M%eta8fQ=0.0_wp
    call allocate_array_body(M%Deta8fQ,G%C,n,ghost_nodes=.true.); M%Deta8fQ=0.0_wp
    call allocate_array_body(M%eta9fQ,G%C,n,ghost_nodes=.true.); M%eta9fQ=0.0_wp
    call allocate_array_body(M%Deta9fQ,G%C,n,ghost_nodes=.true.); M%Deta9fQ=0.0_wp

    wref=2.0_wp*pi*parameters%fref
    if (any(.not.ieee_is_finite(M%M)) .or. any(M%M(:,:,:,2) <= 0.0_wp) .or. &
        any(M%M(:,:,:,3) <= 0.0_wp) .or. any(M%M(:,:,:,1)+2.0_wp*M%M(:,:,:,2) <= 0.0_wp)) &
         call error('anelastic-fQ requires finite positive density and P/S moduli','init_anelastic_fq_properties')
    do i=G%C%mq,G%C%pq; do j=G%C%mr,G%C%pr; do k=G%C%ms,G%C%ps
       if (M%M(i,j,k,2) <= 0.0_wp .or. M%M(i,j,k,3) <= 0.0_wp) &
            call error('mu and density must be positive for anelastic-fQ','init_anelastic_fq_properties')
       val_s=0.0_wp; val_p=0.0_wp
       do l=1,n
          val_s=val_s+M%strength_s_fQ(l)/((wref*wref*M%tau_fQ(l)**2+1.0_wp))
          val_p=val_p+M%strength_p_fQ(l)/((wref*wref*M%tau_fQ(l)**2+1.0_wp))
       end do
       if (val_s >= 1.0_wp .or. val_p >= 1.0_wp) &
            call error('invalid anelastic-fQ modulus correction','init_anelastic_fq_properties')
       vs=sqrt(M%M(i,j,k,2)/M%M(i,j,k,3))
       vp=sqrt((M%M(i,j,k,1)+2.0_wp*M%M(i,j,k,2))/M%M(i,j,k,3))
       mu_s=M%M(i,j,k,3)*vs*vs/(1.0_wp-val_s)
       mu_p=M%M(i,j,k,3)*vp*vp/(1.0_wp-val_p)
       M%M(i,j,k,2)=mu_s
       M%M(i,j,k,1)=mu_p-2.0_wp*mu_s
    end do; end do; end do

    if (G%C%rank == 0) then
       write(*,'(A,I0,A,I0)') 'anelastic-fQ block ',block_id,': mechanisms=',n
       write(*,'(A,ES12.4,A,ES12.4)') '  Qs0=',M%Qs0_fQ,', Qp0=',M%Qp0_fQ
       write(*,'(A,A)') '  coefficient policy=',trim(parameters%coefficient_policy)
       write(*,'(A,ES12.4,A,ES12.4)') '  gamma=',parameters%gamma,', transition Hz=',parameters%f_transition
       write(*,'(A,A)') '  relaxation policy=',trim(parameters%relaxation_policy)
       write(*,'(A,A)') '  transition policy=',trim(parameters%transition_policy)
       if(parameters%transition_policy == 'smooth') write(*,'(A,2ES12.4)') &
            '  transition lower/upper Hz=',parameters%f_transition*parameters%transition_lower_ratio, &
            parameters%f_transition*parameters%transition_upper_ratio
       write(*,'(A,F8.3,A,F8.3,A)') '  max relative Q error: S=',100.0_wp*max_s, &
            ' %, P=',100.0_wp*max_p,' %'
    end if
  end subroutine init_anelastic_fq_properties


  subroutine destroy_anelastic_fq_properties(M)
    type(block_material), intent(inout) :: M
    if (allocated(M%tau_fQ)) deallocate(M%tau_fQ,M%strength_s_fQ,M%strength_p_fQ)
    if (allocated(M%eta4fQ)) then
       deallocate(M%eta4fQ,M%eta5fQ,M%eta6fQ,M%eta7fQ,M%eta8fQ,M%eta9fQ)
       deallocate(M%Deta4fQ,M%Deta5fQ,M%Deta6fQ,M%Deta7fQ,M%Deta8fQ,M%Deta9fQ)
    end if
    M%anelastic_fQ=.false.; M%n_mechanism_fQ=0
  end subroutine destroy_anelastic_fq_properties

end module anelastic_fq_material
