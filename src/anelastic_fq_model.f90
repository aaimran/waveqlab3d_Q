module anelastic_fq_model
  use common, only : wp
  use withers_tables, only : get_relaxation_times
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  type, public :: fq_parameters
     real(wp) :: Qs0(2)=-1.0_wp, Qp0(2)=-1.0_wp
     real(wp) :: gamma=0.0_wp, f_transition=1.0_wp
     real(wp) :: fref=1.0_wp, fmin=0.05_wp, fmax=20.0_wp
     integer :: n_mechanisms=8, nnls_samples=256, nnls_max_iterations=200000
     character(len=32) :: coefficient_policy='nnls-block-ps'
     character(len=32) :: relaxation_policy='band', nnls_objective='relative-q'
     real(wp) :: nnls_tolerance=1.0e-10_wp, max_fit_error=0.10_wp
  end type
  public :: read_fq_parameters, build_fq_coefficients, fq_relaxation_dt_limit
  public :: fq_target_q, fq_response, fq_max_relative_error
contains
  subroutine read_fq_parameters(infile,p,status,message)
    integer, intent(in) :: infile
    type(fq_parameters), intent(out) :: p
    integer, intent(out) :: status
    character(*), intent(out) :: message
    real(wp) :: Qs0(2),Qp0(2),gamma,f_transition,fref,fmin,fmax,nnls_tolerance,max_fit_error
    integer :: n_mechanisms,nnls_samples,nnls_max_iterations,stat
    character(len=32) :: coefficient_policy,relaxation_policy,nnls_objective
    namelist /anelastic_fQ_list/ Qs0,Qp0,gamma,f_transition,fref,fmin,fmax,n_mechanisms, &
         coefficient_policy,relaxation_policy,nnls_objective,nnls_samples,nnls_tolerance, &
         nnls_max_iterations,max_fit_error
    Qs0=p%Qs0; Qp0=p%Qp0; gamma=p%gamma; f_transition=p%f_transition
    fref=p%fref; fmin=p%fmin; fmax=p%fmax; n_mechanisms=p%n_mechanisms
    coefficient_policy=p%coefficient_policy; relaxation_policy=p%relaxation_policy
    nnls_objective=p%nnls_objective; nnls_samples=p%nnls_samples
    nnls_tolerance=p%nnls_tolerance; nnls_max_iterations=p%nnls_max_iterations; max_fit_error=p%max_fit_error
    rewind(infile)
    read(infile,nml=anelastic_fQ_list,iostat=stat)
    status=1; message='anelastic-fQ requires a valid &anelastic_fQ_list'
    if(stat /= 0) return
    p%Qs0=Qs0; p%Qp0=Qp0; p%gamma=gamma; p%f_transition=f_transition
    p%fref=fref; p%fmin=fmin; p%fmax=fmax; p%n_mechanisms=n_mechanisms
    p%coefficient_policy=trim(adjustl(coefficient_policy)); p%relaxation_policy=trim(adjustl(relaxation_policy))
    p%nnls_objective=trim(adjustl(nnls_objective)); p%nnls_samples=nnls_samples
    p%nnls_tolerance=nnls_tolerance; p%nnls_max_iterations=nnls_max_iterations; p%max_fit_error=max_fit_error
    call validate_fq(p,status,message)
  end subroutine

  subroutine validate_fq(p,status,message)
    type(fq_parameters), intent(in) :: p
    integer, intent(out) :: status
    character(*), intent(out) :: message
    real(wp) :: tau(8)
    status=1
    message='anelastic-fQ requires two finite Qs0 and Qp0 values >= 15'
    if(any(.not.ieee_is_finite(p%Qs0)) .or. any(.not.ieee_is_finite(p%Qp0))) return
    if(any(p%Qs0 < 15.0_wp) .or. any(p%Qp0 < 15.0_wp)) return
    message='anelastic-fQ requires finite gamma in [0,0.9] and positive f_transition'
    if(.not.ieee_is_finite(p%gamma) .or. .not.ieee_is_finite(p%f_transition)) return
    if(p%gamma < 0.0_wp .or. p%gamma > 0.9_wp .or. p%f_transition <= 0.0_wp) return
    message='anelastic-fQ requires 0 < fmin < fmax and fmin <= fref <= fmax'
    if(.not.all(ieee_is_finite([p%fmin,p%fmax,p%fref]))) return
    if(p%fmin <= 0.0_wp .or. p%fmax <= p%fmin .or. p%fref < p%fmin .or. p%fref > p%fmax) return
    message='anelastic-fQ n_mechanisms must be 4..8 and nnls_samples >= n_mechanisms'
    if(p%n_mechanisms < 4 .or. p%n_mechanisms > 8 .or. p%nnls_samples < p%n_mechanisms) return
    message='anelastic-fQ supports coefficient_policy=nnls-block-ps and nnls_objective=relative-q'
    if(p%coefficient_policy /= 'nnls-block-ps' .or. p%nnls_objective /= 'relative-q') return
    message='anelastic-fQ relaxation_policy must be band, or fq8-table with n_mechanisms=8'
    if(p%relaxation_policy /= 'band' .and. p%relaxation_policy /= 'fq8-table') return
    if(p%relaxation_policy == 'fq8-table' .and. p%n_mechanisms /= 8) return
    message='anelastic-fQ requires positive finite tolerances and nnls_max_iterations >= 1'
    if(.not.all(ieee_is_finite([p%nnls_tolerance,p%max_fit_error]))) return
    if(p%nnls_tolerance <= 0.0_wp .or. p%max_fit_error <= 0.0_wp .or. p%nnls_max_iterations < 1) return
    call relaxation_times(p,tau(1:p%n_mechanisms))
    message='anelastic-fQ frequencies produce invalid relaxation times'
    if(.not.all(ieee_is_finite(tau(1:p%n_mechanisms)))) return
    if(any(tau(1:p%n_mechanisms) <= 0.0_wp)) return
    status=0; message=''
  end subroutine

  pure real(wp) function fq_target_q(f,q0,gamma,transition) result(q)
    real(wp), intent(in) :: f,q0,gamma,transition
    q=q0
    if(f > transition) q=q0*(f/transition)**gamma
  end function

  subroutine relaxation_times(p,tau)
    type(fq_parameters), intent(in) :: p
    real(wp), intent(out) :: tau(:)
    real(wp), parameter :: pi=3.141592653589793_wp
    integer :: k,n
    n=size(tau)
    if(p%relaxation_policy == 'fq8-table') then
       call get_relaxation_times(p%gamma,tau)
       tau=tau/p%f_transition
    else
       do k=1,n
          tau(k)=exp(log(1.0_wp/(2.0_wp*pi*p%fmax))+ &
               real(k-1,wp)/real(n-1,wp)*log(p%fmax/p%fmin))
       end do
    endif
  end subroutine

  subroutine build_fq_coefficients(p,block_id,tau,ss,sp,status,message)
    type(fq_parameters), intent(in) :: p
    integer, intent(in) :: block_id
    real(wp), intent(out) :: tau(:),ss(:),sp(:)
    integer, intent(out) :: status
    character(*), intent(out) :: message
    real(wp) :: es,ep
    call validate_fq(p,status,message)
    if(status /= 0) return
    status=1; message='anelastic-fQ invalid block id or coefficient array extent'
    if(block_id < 1 .or. block_id > 2) return
    if(size(tau) /= p%n_mechanisms .or. size(ss) /= size(tau) .or. size(sp) /= size(tau)) return
    call relaxation_times(p,tau)
    call fit_strengths(p,p%Qs0(block_id),tau,ss,status,message)
    if(status /= 0) return
    call fit_strengths(p,p%Qp0(block_id),tau,sp,status,message)
    if(status /= 0) return
    status=1; message='anelastic-fQ invalid coefficients or non-positive relaxed modulus'
    if(.not.all(ieee_is_finite(tau)) .or. any(tau <= 0.0_wp)) return
    if(.not.all(ieee_is_finite(ss)) .or. .not.all(ieee_is_finite(sp))) return
    if(any(ss < 0.0_wp) .or. any(sp < 0.0_wp) .or. sum(ss) >= 1.0_wp .or. sum(sp) >= 1.0_wp) return
    call fq_max_relative_error(p%Qs0(block_id),p%gamma,p%f_transition,tau,ss,p%fmin,p%fmax,es)
    call fq_max_relative_error(p%Qp0(block_id),p%gamma,p%f_transition,tau,sp,p%fmin,p%fmax,ep)
    if(max(es,ep) > p%max_fit_error) then
       write(message,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'anelastic-fQ fitted response exceeds max_fit_error: S=',es,', P=',ep,', limit=',p%max_fit_error
       return
    endif
    status=0; message=''
  end subroutine

  subroutine fit_strengths(p,q0,tau,strength,status,message)
    type(fq_parameters), intent(in) :: p
    real(wp), intent(in) :: q0,tau(:)
    real(wp), intent(out) :: strength(:)
    integer, intent(out) :: status
    character(*), intent(out) :: message
    real(wp) :: a(p%nnls_samples,size(tau)),w(size(tau)),residual(p%nnls_samples)
    real(wp) :: f,q,x,change,old,new,delta,denom
    real(wp), parameter :: pi=3.141592653589793_wp
    integer :: i,k,iter
    do i=1,p%nnls_samples
       f=p%fmin*(p%fmax/p%fmin)**(real(i-1,wp)/real(p%nnls_samples-1,wp))
       q=fq_target_q(f,q0,p%gamma,p%f_transition)
       do k=1,size(tau)
          x=2.0_wp*pi*f*tau(k)
          ! w=q0*strength. At gamma=0 this is the cQ relative-Q system.
          a(i,k)=(q*x+1.0_wp)/(q0*(1.0_wp+x*x))
       enddo
    enddo
    status=1; message='anelastic-fQ target law produces non-finite fitting values'
    if(.not.all(ieee_is_finite(a))) return
    w=0.0_wp; residual=1.0_wp
    status=1; message='anelastic-fQ NNLS did not converge within nnls_max_iterations'
    do iter=1,p%nnls_max_iterations
       change=0.0_wp
       do k=1,size(tau)
          old=w(k); denom=max(dot_product(a(:,k),a(:,k)),tiny(1.0_wp))
          new=max(0.0_wp,old+dot_product(a(:,k),residual)/denom)
          delta=new-old
          residual=residual-delta*a(:,k)
          w(k)=new; change=max(change,abs(delta))
       enddo
       if(change <= p%nnls_tolerance) then
          status=0; message=''; exit
       endif
    enddo
    strength=w/q0
  end subroutine

  pure complex(wp) function fq_response(f,tau,strength) result(response)
    real(wp), intent(in) :: f,tau(:),strength(:)
    real(wp), parameter :: pi=3.141592653589793_wp
    integer :: k
    response=cmplx(1.0_wp,0.0_wp,wp)
    do k=1,size(tau)
       response=response-strength(k)/cmplx(1.0_wp,2.0_wp*pi*f*tau(k),wp)
    enddo
  end function

  pure subroutine fq_max_relative_error(q0,gamma,transition,tau,strength,fmin,fmax,max_error)
    real(wp), intent(in) :: q0,gamma,transition,tau(:),strength(:),fmin,fmax
    real(wp), intent(out) :: max_error
    real(wp) :: f,q,error
    complex(wp) :: response
    integer :: i
    max_error=0.0_wp
    ! Dense validation independent of the number of fitting samples.
    do i=1,1025
       f=fmin*(fmax/fmin)**(real(i-1,wp)/1024.0_wp)
       response=fq_response(f,tau,strength)
       if(aimag(response) <= tiny(1.0_wp)) then
          max_error=huge(1.0_wp); return
       endif
       q=real(response,wp)/aimag(response)
       error=abs(q/fq_target_q(f,q0,gamma,transition)-1.0_wp)
       max_error=max(max_error,error)
    enddo
    if(transition >= fmin .and. transition <= fmax) then
       response=fq_response(transition,tau,strength)
       if(aimag(response) <= tiny(1.0_wp)) then
          max_error=huge(1.0_wp); return
       endif
       max_error=max(max_error,abs(real(response,wp)/aimag(response)/q0-1.0_wp))
    endif
  end subroutine

  real(wp) function fq_relaxation_dt_limit(p) result(limit)
    type(fq_parameters), intent(in) :: p
    real(wp) :: tau(p%n_mechanisms)
    call relaxation_times(p,tau)
    limit=2.0_wp*minval(tau)
  end function
end module anelastic_fq_model
