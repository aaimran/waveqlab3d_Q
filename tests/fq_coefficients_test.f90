program fq_coefficients_test
  use common, only : wp
  use anelastic_fq_model
  use anelastic_cq_model, only : cq_parameters,build_cq_coefficients
  use anelastic_fq8_model, only : fq8_parameters,build_fq8_coefficients,fq8_realized_q
  implicit none
  type(fq_parameters) :: p
  type(cq_parameters) :: c
  type(fq8_parameters) :: legacy
  real(wp), allocatable :: tau(:),ss(:),sp(:),ct(:),cs(:),cp(:)
  real(wp) :: es,ep,f,q,qlegacy,lt(8),ls(8),lp(8)
  complex(wp) :: r
  integer :: n,b,status,i
  character(len=256) :: message
  p%Qs0=[40.0_wp,120.0_wp]; p%Qp0=[69.3_wp,155.9_wp]; p%gamma=0.6_wp
  p%max_fit_error=0.30_wp
  do n=4,8
     p%n_mechanisms=n
     allocate(tau(n),ss(n),sp(n))
     do b=1,2
        call build_fq_coefficients(p,b,tau,ss,sp,status,message)
        if(status /= 0) then
           print *,n,b,trim(message); error stop 'frequency-Q fit failed'
        endif
        if(any(tau <= 0.0_wp) .or. any(ss < 0.0_wp) .or. any(sp < 0.0_wp)) error stop 'invalid coefficients'
        if(sum(ss) >= 1.0_wp .or. sum(sp) >= 1.0_wp) error stop 'invalid relaxed modulus'
        if(maxval(abs(ss-sp)) < 1.0e-8_wp) error stop 'P/S fits must differ'
        call fq_max_relative_error(p%Qs0(b),p,tau,ss,es)
        call fq_max_relative_error(p%Qp0(b),p,tau,sp,ep)
        print *, 'mechanisms/block/error:',n,b,es,ep
        if(max(es,ep) > p%max_fit_error) error stop 'fit accuracy'
     enddo
     deallocate(tau,ss,sp)
  enddo
  n=8; p%n_mechanisms=n
  allocate(tau(n),ss(n),sp(n),ct(n),cs(n),cp(n))
  p%gamma=0.0_wp; p%max_fit_error=0.10_wp
  c%Qs0=p%Qs0; c%Qp0=p%Qp0
  do b=1,2
     call build_fq_coefficients(p,b,tau,ss,sp,status,message)
     if(status /= 0) error stop 'gamma zero fit failed'
     call build_cq_coefficients(c,b,ct,cs,cp)
     if(maxval(abs(tau-ct)) > 1.0e-12_wp) error stop 'constant-Q relaxation mismatch'
     do i=1,101
        f=p%fmin*(p%fmax/p%fmin)**(real(i-1,wp)/100.0_wp)
        r=fq_response(f,tau,ss); q=real(r,wp)/aimag(r)
        r=fq_response(f,ct,cs/c%Qs0(b)); qlegacy=real(r,wp)/aimag(r)
        if(abs(q/qlegacy-1.0_wp) > 1.0e-3_wp) error stop 'constant-Q S response mismatch'
        r=fq_response(f,tau,sp); q=real(r,wp)/aimag(r)
        r=fq_response(f,ct,cp/c%Qp0(b)); qlegacy=real(r,wp)/aimag(r)
        if(abs(q/qlegacy-1.0_wp) > 1.0e-3_wp) error stop 'constant-Q P response mismatch'
     enddo
  enddo
  ! Match the existing conventional full-layout fQ8 times and fitting band.
  legacy%coefficient_method='conventional-nnls'; legacy%coarse_grain=0
  legacy%Qs0=50.0_wp; legacy%Qp0=80.0_wp; legacy%gamma=0.6_wp
  legacy%f_transition=1.0_wp; legacy%fref=1.0_wp
  call build_fq8_coefficients(legacy,lt,ls,lp,status,message)
  if(status /= 0) error stop 'legacy coefficient failure'
  p%Qs0=legacy%Qs0; p%Qp0=legacy%Qp0; p%gamma=legacy%gamma
  p%relaxation_policy='fq8-table'; p%fmin=0.1_wp; p%fmax=10.0_wp
  p%nnls_tolerance=1.0e-12_wp; p%nnls_max_iterations=1000000
  call build_fq_coefficients(p,1,tau,ss,sp,status,message)
  if(status /= 0) then
     print *,trim(message); error stop 'fQ8 comparison fit failed'
  endif
  if(maxval(abs(tau-lt)) > 1.0e-12_wp) error stop 'fQ8 relaxation mismatch'
  do i=1,101
     f=0.1_wp*100.0_wp**(real(i-1,wp)/100.0_wp)
     r=fq_response(f,tau,ss); q=real(r,wp)/aimag(r)
     qlegacy=fq8_realized_q(f,lt,ls)
     if(abs(q/qlegacy-1.0_wp) > 1.0e-4_wp) error stop 'fQ8 S response mismatch'
     r=fq_response(f,tau,sp); q=real(r,wp)/aimag(r)
     qlegacy=fq8_realized_q(f,lt,lp)
     if(abs(q/qlegacy-1.0_wp) > 1.0e-4_wp) error stop 'fQ8 P response mismatch'
  enddo
  p%nnls_max_iterations=1
  call build_fq_coefficients(p,1,tau,ss,sp,status,message)
  if(status == 0 .or. index(message,'did not converge') == 0) error stop 'missing convergence failure'
  p%nnls_max_iterations=1000000; p%max_fit_error=1.0e-8_wp
  call build_fq_coefficients(p,1,tau,ss,sp,status,message)
  if(status == 0 .or. index(message,'max_fit_error') == 0) error stop 'missing fit failure'
  p%relaxation_policy='band'; p%max_fit_error=0.1_wp
  if(abs(fq_relaxation_dt_limit(p)-2.0_wp/(2.0_wp*acos(-1.0_wp)*p%fmax)) > 1.0e-12_wp) &
       error stop 'relaxation dt limit'
  print *, 'Independent frequency-Q coefficient tests passed'
end program
