program fq_transition_test
  use common, only : wp
  use anelastic_fq_model
  use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
  implicit none
  type(fq_parameters) :: p,sharp
  real(wp), allocatable :: tau(:),ss(:),sp(:),ts(:),ssharp(:),psharp(:)
  real(wp) :: f,lower,upper,q,previous,h,slope_left,slope_right,es,ep
  integer :: i,j,n,b,status
  character(len=256) :: message
  p%Qs0=[40.0_wp,120.0_wp]; p%Qp0=[69.3_wp,155.9_wp]
  p%gamma=0.6_wp; p%f_transition=2.5_wp
  if(p%transition_policy /= 'sharp') error stop 'sharp must remain the default'
  ! The unchanged sharp formula is the backward-compatibility oracle.
  do i=1,101
     f=0.01_wp*10000.0_wp**(real(i-1,wp)/100.0_wp)
     q=40.0_wp
     if(f > p%f_transition) q=40.0_wp*(f/p%f_transition)**p%gamma
     if(fq_target_q(f,40.0_wp,p) /= q) error stop 'sharp target changed'
  enddo
  p%transition_policy='smooth'
  lower=p%f_transition*p%transition_lower_ratio; upper=p%f_transition*p%transition_upper_ratio
  if(fq_target_q(lower,40.0_wp,p) /= 40.0_wp) error stop 'lower join value'
  q=40.0_wp*p%transition_upper_ratio**p%gamma
  if(abs(fq_target_q(upper,40.0_wp,p)/q-1.0_wp) > 1.0e-14_wp) error stop 'upper join value'
  if(fq_target_q(p%f_transition,40.0_wp,p) <= 40.0_wp) error stop 'smooth center should exceed plateau'
  ! Verify the log-log slope on BOTH sides of each join without duplicating the cubic formula.
  h=1.0e-5_wp
  slope_left=log(fq_target_q(lower,40.0_wp,p)/fq_target_q(lower*exp(-h),40.0_wp,p))/h
  slope_right=log(fq_target_q(lower*exp(h),40.0_wp,p)/fq_target_q(lower,40.0_wp,p))/h
  if(abs(slope_left) > 1.0e-4_wp .or. abs(slope_right) > 1.0e-4_wp) error stop 'lower slope continuity'
  slope_left=log(fq_target_q(upper,40.0_wp,p)/fq_target_q(upper*exp(-h),40.0_wp,p))/h
  slope_right=log(fq_target_q(upper*exp(h),40.0_wp,p)/fq_target_q(upper,40.0_wp,p))/h
  if(abs(slope_left-p%gamma) > 1.0e-4_wp .or. abs(slope_right-p%gamma) > 1.0e-4_wp) &
       error stop 'upper slope continuity'
  do j=0,3
     p%gamma=0.3_wp*real(j,wp)
     previous=40.0_wp
     do i=1,1001
        f=p%f_transition*0.01_wp*10000.0_wp**(real(i-1,wp)/1000.0_wp)
        q=fq_target_q(f,40.0_wp,p)
        if(q < previous-1.0e-12_wp) error stop 'smooth target must be monotone'
        if(j == 0 .and. q /= 40.0_wp) error stop 'gamma zero must remain constant'
        previous=q
     enddo
  enddo
  ! Check a non-default interval to ensure ratios affect the target.
  p%gamma=0.6_wp; p%transition_lower_ratio=0.7_wp; p%transition_upper_ratio=1.5_wp
  lower=p%f_transition*p%transition_lower_ratio; upper=p%f_transition*p%transition_upper_ratio
  if(fq_target_q(lower,40.0_wp,p) /= 40.0_wp) error stop 'custom lower join'
  if(abs(fq_target_q(upper,40.0_wp,p)/(40.0_wp*1.5_wp**p%gamma)-1.0_wp) > 1.0e-14_wp) &
       error stop 'custom upper join'
  p%transition_lower_ratio=0.8_wp; p%transition_upper_ratio=1.2_wp; p%f_transition=1.0_wp
  p%max_fit_error=0.30_wp
  do n=4,8
     p%n_mechanisms=n; p%relaxation_policy='band'
     allocate(tau(n),ss(n),sp(n),ts(n),ssharp(n),psharp(n))
     sharp=p; sharp%transition_policy='sharp'
     do b=1,2
        call build_fq_coefficients(p,b,tau,ss,sp,status,message)
        if(status /= 0) then
           print *,trim(message); error stop 'smooth band fit'
        endif
        call fq_max_relative_error(p%Qs0(b),p,tau,ss,es)
        call fq_max_relative_error(p%Qp0(b),p,tau,sp,ep)
        if(max(es,ep) > p%max_fit_error) error stop 'smooth error bound'
        call build_fq_coefficients(sharp,b,ts,ssharp,psharp,status,message)
        if(status /= 0) error stop 'sharp comparison fit'
        if(any(tau /= ts)) error stop 'smoothing must not change relaxation times'
        if(maxval(abs(ss-ssharp)) < 1.0e-8_wp) error stop 'fit did not use smooth target'
     enddo
     if(n == 8) then
        p%relaxation_policy='fq8-table'
        call build_fq_coefficients(p,2,tau,ss,sp,status,message)
        if(status /= 0) then
           print *,trim(message); error stop 'smooth table fit'
        endif
        ! Transition boundaries outside the fitting band are allowed.
        p%relaxation_policy='band'; p%f_transition=0.01_wp
        call build_fq_coefficients(p,1,tau,ss,sp,status,message)
        if(status /= 0) then
           print *,trim(message); error stop 'low out-of-band transition fit'
        endif
        p%f_transition=100.0_wp
        call build_fq_coefficients(p,1,tau,ss,sp,status,message)
        if(status /= 0) then
           print *,trim(message); error stop 'high out-of-band transition fit'
        endif
        p%f_transition=1.0_wp; p%relaxation_policy='band'; p%gamma=0.0_wp
        sharp=p; sharp%transition_policy='sharp'
        call build_fq_coefficients(p,1,tau,ss,sp,status,message)
        if(status /= 0) error stop 'smooth constant-Q fit'
        call build_fq_coefficients(sharp,1,ts,ssharp,psharp,status,message)
        if(status /= 0 .or. any(ss /= ssharp) .or. any(sp /= psharp)) error stop 'constant-Q fit changed'
        p%gamma=0.6_wp
        p%transition_policy='unknown'
        call expect_invalid('transition_policy')
        p%transition_policy='smooth'; p%transition_lower_ratio=0.0_wp
        call expect_invalid('transition ratios')
        p%transition_lower_ratio=0.8_wp; p%transition_upper_ratio=1.0_wp
        call expect_invalid('transition ratios')
        p%transition_upper_ratio=1.2_wp; p%transition_lower_ratio=0.2_wp
        call expect_invalid('monotone')
        p%transition_lower_ratio=ieee_value(1.0_wp,ieee_quiet_nan)
        call expect_invalid('transition ratios')
     endif
     deallocate(tau,ss,sp,ts,ssharp,psharp)
  enddo
  print *, 'Frequency-Q smooth transition tests passed'
contains
  subroutine expect_invalid(expected)
    character(*), intent(in) :: expected
    call build_fq_coefficients(p,1,tau,ss,sp,status,message)
    if(status == 0 .or. index(message,expected) == 0) then
       print *,trim(message); error stop 'expected transition validation failure'
    endif
  end subroutine
end program
