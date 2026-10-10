module anelastic_cg_t_model
 use common,only:wp
 use anelastic_cg8_model,only:cg8_parameters,cg8_coefficients,build_cg8_coefficients,cg8_response,cg8_target
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::build_cgt_coefficients
contains
 subroutine build_cgt_coefficients(p,id,rho,cs,cp,c,status,message)
  type(cg8_parameters),intent(in)::p
  integer,intent(in)::id
  real(wp),intent(in)::rho,cs,cp
  type(cg8_coefficients),intent(out)::c
  integer,intent(out)::status
  character(*),intent(out)::message
  real(wp)::f,extra(4),err,qs,qp,reference_s,reference_p
  complex(wp)::rs,rp
  integer::i,n
  status=1;message='CG-T max_fit_error must be finite, positive, and <= 0.02'
  if(.not.ieee_is_finite(p%max_fit_error)) return
  if(p%max_fit_error<=0.or.p%max_fit_error>0.02_wp) return
  ! The additive branch is the only permitted spectrum. No nodal harmonic fit.
  call build_cg8_coefficients(p,id,rho,cs,cp,.false.,c,status,message)
  if(status/=0) return
  status=1;message='CG-T additive spectrum is not passive'
  if(any(c%bulk_strength<0).or.any(c%shear_strength<0).or. &
     sum(c%bulk_strength)>=1.or.sum(c%shear_strength)>=1.or.any(c%tau<=0)) return
  n=max(1025,4*p%fit_samples+1);err=0
  extra=[p%fref,p%f_transition,p%f_transition*p%transition_lower_ratio,p%f_transition*p%transition_upper_ratio]
  do i=1,n+4
    if(i<=n) then
      f=p%fmin*(p%fmax/p%fmin)**(real(i-1,wp)/real(n-1,wp))
    else
      f=extra(i-n);if(f<p%fmin.or.f>p%fmax) cycle
    endif
    rs=c%shear*cg8_response(f,c%tau,c%shear_strength,.false.)
    rp=c%bulk*cg8_response(f,c%tau,c%bulk_strength,.false.)+4*rs/3
    if(aimag(rs)<=0.or.aimag(rp)<=0) return
    qs=real(rs,wp)/aimag(rs);qp=real(rp,wp)/aimag(rp)
    err=max(err,abs(qs/cg8_target(f,p%Qs0(id),p)-1),abs(qp/cg8_target(f,p%Qp0(id),p)-1))
  enddo
  c%max_error=err
  message='CG-T additive fit failed independent dense-frequency error check'
  if(.not.ieee_is_finite(err).or.err>p%max_fit_error) return
  rs=c%shear*cg8_response(p%fref,c%tau,c%shear_strength,.false.)
  rp=c%bulk*cg8_response(p%fref,c%tau,c%bulk_strength,.false.)+4*rs/3
  reference_s=1/real(sqrt(rho/rs),wp);reference_p=1/real(sqrt(rho/rp),wp)
  message='CG-T reference phase velocity normalization failed'
  if(abs(reference_s/cs-1)>1e-9_wp.or.abs(reference_p/cp-1)>1e-9_wp) return
  status=0;message=''
 end subroutine
end module
