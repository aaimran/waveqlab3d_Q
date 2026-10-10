module anelastic_cg8_model
  use common, only: wp
  use anelastic_fq_model, only: fq_parameters,fq_target_q
  use withers_tables, only: get_relaxation_times
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type, public :: cg8_parameters
    real(wp) :: Qs0(2)=-1.0_wp,Qp0(2)=-1.0_wp
    real(wp) :: fref=1.0_wp,fmin=0.05_wp,fmax=20.0_wp,gamma=0.0_wp,f_transition=1.0_wp
    character(len=32) :: transition_policy='sharp',relaxation_policy='band'
    character(len=32) :: coefficient_policy='constrained-effective-fit',boundary_policy='full-buffer'
    real(wp) :: transition_lower_ratio=0.8_wp,transition_upper_ratio=1.2_wp
    real(wp) :: fit_tolerance=1.0e-9_wp,max_fit_error=0.05_wp
    integer :: fit_samples=256,fit_max_iterations=500,buffer_layers=-1,pattern_origin(3)=1
  end type
  type, public :: cg8_coefficients
    real(wp) :: tau(8)=0.0_wp,bulk_strength(8)=0.0_wp,shear_strength(8)=0.0_wp
    real(wp) :: bulk=0.0_wp,shear=0.0_wp,max_error=0.0_wp
    logical :: coarse=.true.
  end type
  public :: validate_cg8,build_cg8_coefficients,cg8_response,cg8_target,cg8_times
contains
  pure real(wp) function cg8_target(f,q0,p) result(q)
    real(wp),intent(in):: f,q0
    type(cg8_parameters),intent(in):: p
    type(fq_parameters):: target
    target%gamma=p%gamma;target%f_transition=p%f_transition;target%transition_policy=p%transition_policy
    target%transition_lower_ratio=p%transition_lower_ratio;target%transition_upper_ratio=p%transition_upper_ratio
    q=fq_target_q(f,q0,target)
  end function

  subroutine validate_cg8(p,nblocks,status,message)
    type(cg8_parameters),intent(in):: p
    integer,intent(in):: nblocks
    integer,intent(out):: status
    character(*),intent(out):: message
    status=1; message='CG8 requires one or two blocks with finite Qs0/Qp0 >= 15'
    if(nblocks<1.or.nblocks>2) return
    if(.not.all(ieee_is_finite(p%Qs0(1:nblocks))).or..not.all(ieee_is_finite(p%Qp0(1:nblocks)))) return
    if(any(p%Qs0(1:nblocks)<15.0_wp).or.any(p%Qp0(1:nblocks)<15.0_wp)) return
    message='CG8 requires 0 < fmin <= fref <= fmax with fmin < fmax'
    if(.not.all(ieee_is_finite([p%fmin,p%fmax,p%fref]))) return
    if(p%fmin<=0.or.p%fmax<=p%fmin.or.p%fref<p%fmin.or.p%fref>p%fmax) return
    message='CG8 requires gamma in [0,0.9] and positive f_transition'
    if(.not.all(ieee_is_finite([p%gamma,p%f_transition]))) return
    if(p%gamma<0.or.p%gamma>0.9_wp.or.p%f_transition<=0) return
    message='CG8 transition policy must be sharp or smooth with valid monotone ratios'
    if(p%transition_policy/='sharp'.and.p%transition_policy/='smooth') return
    if(.not.all(ieee_is_finite([p%transition_lower_ratio,p%transition_upper_ratio]))) return
    if(p%transition_lower_ratio<=0.or.p%transition_lower_ratio>=1.or.p%transition_upper_ratio<=1) return
    if(p%transition_policy=='smooth') then
      if(log(p%transition_lower_ratio)+2*log(p%transition_upper_ratio)<0) return
    endif
    message='CG8 supports band or withers-times relaxation policies and constrained-effective-fit coefficients'
    if(p%relaxation_policy/='band'.and.p%relaxation_policy/='withers-times') return
    if(p%coefficient_policy/='constrained-effective-fit') return
    message='CG8 requires full-buffer boundaries, buffer_layers=-1 or nonnegative, and valid fit controls'
    if(p%boundary_policy/='full-buffer'.or.p%buffer_layers< -1.or.p%fit_samples<32.or.p%fit_max_iterations<1) return
    if(.not.all(ieee_is_finite([p%fit_tolerance,p%max_fit_error]))) return
    if(p%fit_tolerance<=0.or.p%max_fit_error<=0.or.p%max_fit_error>0.05_wp) return
    status=0;message=''
  end subroutine

  subroutine cg8_times(p,tau)
    type(cg8_parameters),intent(in)::p
    real(wp),intent(out)::tau(8)
    integer::k
    if(p%relaxation_policy=='withers-times') then
      call get_relaxation_times(p%gamma,tau);tau=tau/p%f_transition
    else
      ! Absorption cutoffs extend beyond the fit band, as in the paper.
      do k=1,8
        tau(k)=exp(log(1.0_wp/(8.0_wp*acos(-1.0_wp)*p%fmax))+ &
          real(k-1,wp)/7.0_wp*log(16.0_wp*p%fmax/p%fmin))
      enddo
    endif
  end subroutine

  pure complex(wp) function cg8_response(f,tau,w,coarse) result(r)
    real(wp),intent(in)::f,tau(8),w(8)
    logical,intent(in)::coarse
    complex(wp)::local(8)
    local=1.0_wp-w/cmplx(1.0_wp,2*acos(-1.0_wp)*f*tau,wp)
    if(coarse) then
      r=8.0_wp/sum(1.0_wp/local)
    else
      r=1.0_wp-sum(1.0_wp-local)
    endif
  end function

  subroutine normalize(tau,w,coarse,fref,rho,cs,cp,bulk,shear,ok)
    real(wp),intent(in)::tau(8),w(16),fref,rho,cs,cp
    logical,intent(in)::coarse
    real(wp),intent(out)::bulk,shear
    logical,intent(out)::ok
    complex(wp)::rs,rk,rp
    real(wp)::lo,hi,mid,velocity
    integer::i
    rs=cg8_response(fref,tau,w(9:16),coarse);rk=cg8_response(fref,tau,w(1:8),coarse)
    shear=rho*cs*cs*real(1/sqrt(rs),wp)**2
    lo=0;hi=rho*cp*cp
    rp=(4.0_wp/3.0_wp)*shear*rs
    velocity=1/real(sqrt(rho/rp),wp)
    ok=velocity<cp.and.ieee_is_finite(shear).and.shear>0
    if(.not.ok) return
    do i=1,60
      rp=hi*rk+(4.0_wp/3.0_wp)*shear*rs
      velocity=1/real(sqrt(rho/rp),wp)
      if(velocity>=cp) exit
      hi=2*hi
    enddo
    if(velocity<cp.or..not.ieee_is_finite(hi)) then
      ok=.false.;return
    endif
    do i=1,60
      mid=(lo+hi)/2
      rp=mid*rk+(4.0_wp/3.0_wp)*shear*rs
      velocity=1/real(sqrt(rho/rp),wp)
      if(velocity<cp) then
        lo=mid
      else
        hi=mid
      endif
    enddo
    bulk=(lo+hi)/2
    ok=bulk>0.and.ieee_is_finite(bulk).and.abs(velocity/cp-1)<1.0e-8_wp
  end subroutine

  subroutine residual(p,id,tau,w,coarse,rho,cs,cp,f,r,bulk,shear,ok)
    type(cg8_parameters),intent(in)::p
    integer,intent(in)::id
    real(wp),intent(in)::tau(8),w(16),rho,cs,cp,f(:)
    logical,intent(in)::coarse
    real(wp),intent(out)::r(2*size(f)),bulk,shear
    logical,intent(out)::ok
    complex(wp)::rs,rk,rp
    integer::i,n
    n=size(f)
    ok=all(w>=0).and.all(w<1.0_wp-1.0e-10_wp)
    if(.not.coarse) ok=ok.and.sum(w(1:8))<1.and.sum(w(9:16))<1
    if(.not.ok) return
    call normalize(tau,w,coarse,p%fref,rho,cs,cp,bulk,shear,ok)
    if(.not.ok) return
    do i=1,n
      rs=shear*cg8_response(f(i),tau,w(9:16),coarse)
      rk=bulk*cg8_response(f(i),tau,w(1:8),coarse);rp=rk+(4.0_wp/3.0_wp)*rs
      if(aimag(rs)<=tiny(1.0_wp).or.aimag(rp)<=tiny(1.0_wp)) then
        ok=.false.;return
      endif
      r(i)=real(rs,wp)/aimag(rs)/cg8_target(f(i),p%Qs0(id),p)-1
      r(n+i)=real(rp,wp)/aimag(rp)/cg8_target(f(i),p%Qp0(id),p)-1
    enddo
    ok=all(ieee_is_finite(r))
  end subroutine

  subroutine solve(a,b,x,ok)
    real(wp),intent(in)::a(16,16),b(16)
    real(wp),intent(out)::x(16)
    logical,intent(out)::ok
    real(wp)::m(16,17),row(17),factor
    integer::i,j,pivot
    m(:,1:16)=a;m(:,17)=b;ok=.false.
    do i=1,16
      pivot=i-1+maxloc(abs(m(i:16,i)),dim=1)
      if(abs(m(pivot,i))<tiny(1.0_wp)) return
      row=m(i,:);m(i,:)=m(pivot,:);m(pivot,:)=row
      do j=i+1,16
        factor=m(j,i)/m(i,i);m(j,i:17)=m(j,i:17)-factor*m(i,i:17)
      enddo
    enddo
    do i=16,1,-1
      x(i)=(m(i,17)-sum(m(i,i+1:16)*x(i+1:16)))/m(i,i)
    enddo
    ok=all(ieee_is_finite(x))
  end subroutine

  subroutine seed_strengths(p,q0,tau,coarse,w)
    type(cg8_parameters),intent(in)::p
    real(wp),intent(in)::q0,tau(8)
    logical,intent(in)::coarse
    real(wp),intent(out)::w(8)
    real(wp)::a(p%fit_samples,8),r(p%fit_samples),f,q,x,old,new,delta,change
    integer::i,k,it
    do i=1,p%fit_samples
      f=p%fmin*(p%fmax/p%fmin)**(real(i-1,wp)/real(p%fit_samples-1,wp))
      q=cg8_target(f,q0,p)
      do k=1,8
        x=2*acos(-1.0_wp)*f*tau(k);a(i,k)=(q*x+1)/(1+x*x)
      enddo
    enddo
    w=0;r=1
    do it=1,20000
      change=0
      do k=1,8
        old=w(k);new=max(0.0_wp,old+dot_product(a(:,k),r)/dot_product(a(:,k),a(:,k)))
        delta=new-old;r=r-delta*a(:,k);w(k)=new;change=max(change,abs(delta))
      enddo
      if(change<1.0e-12_wp) exit
    enddo
    if(coarse) w=8*w
    ! Only an initial feasible guess; the effective nonlinear fit follows.
    w=max(1.0e-8_wp,min(.95_wp,w))
  end subroutine

  subroutine build_cg8_coefficients(p,id,rho,cs,cp,coarse,c,status,message)
    type(cg8_parameters),intent(in)::p
    integer,intent(in)::id
    real(wp),intent(in)::rho,cs,cp
    logical,intent(in)::coarse
    type(cg8_coefficients),intent(out)::c
    integer,intent(out)::status
    character(*),intent(out)::message
    real(wp),allocatable::f(:),r(:),rt(:),rm(:),jac(:,:),fv(:),rv(:)
    real(wp)::w(16),trial(16),a(16,16),b(16),delta(16),eps,alpha,cost,bulk,shear,bk,sm,damping,qk,loss,loss_scale
    integer::i,j,it,line,n
    logical::ok,moved,converged
    status=1;message='CG8 requires finite positive density, shear speed, and bulk modulus'
    if(.not.all(ieee_is_finite([rho,cs,cp]))) return
    if(rho<=0.or.cs<=0.or.cp*cp<=4*cs*cs/3) return
    if(id<1.or.id>2) return
    call validate_cg8(p,id,status,message)
    if(status/=0) return
    status=1;message='CG8 incompatible P/S targets: inferred bulk dissipation is negative'
    loss=cp*cp/p%Qp0(id)-(4.0_wp/3)*cs*cs/p%Qs0(id)
    loss_scale=max(cp*cp/p%Qp0(id),(4.0_wp/3)*cs*cs/p%Qs0(id))
    if(loss < -32*epsilon(1.0_wp)*loss_scale) return
    c%coarse=coarse;call cg8_times(p,c%tau)
    if(.not.all(ieee_is_finite(c%tau)).or.any(c%tau<=0)) then
      message='CG8 invalid relaxation times';return
    endif
    n=p%fit_samples;allocate(f(n),r(2*n),rt(2*n),jac(2*n,16),rm(2*n))
    do i=1,n
      f(i)=p%fmin*(p%fmax/p%fmin)**(real(i-1,wp)/real(n-1,wp))
    enddo
    if(loss>32*epsilon(1.0_wp)*loss_scale) then
      qk=(cp*cp-4*cs*cs/3)/loss
      call seed_strengths(p,qk,c%tau,coarse,w(1:8))
    else
      w(1:8)=1.0e-12_wp
    endif
    call seed_strengths(p,p%Qs0(id),c%tau,coarse,w(9:16))
    damping=1.0e-6_wp;converged=.false.
    do it=1,p%fit_max_iterations
      call residual(p,id,c%tau,w,coarse,rho,cs,cp,f,r,bulk,shear,ok)
      if(.not.ok) exit
      cost=sum(r*r)
      do j=1,16
        trial=w;eps=1.0e-6_wp;trial(j)=trial(j)+eps
        call residual(p,id,c%tau,trial,coarse,rho,cs,cp,f,rt,bk,sm,ok)
        if(.not.ok) exit
        jac(:,j)=(rt-r)/eps
        if(w(j)>eps) then
          trial=w;trial(j)=trial(j)-eps
          call residual(p,id,c%tau,trial,coarse,rho,cs,cp,f,rm,bk,sm,ok)
          if(ok) jac(:,j)=(rt-rm)/(2*eps)
        endif
        ok=.true.
      enddo
      if(.not.ok) exit
      a=matmul(transpose(jac),jac);b=-matmul(transpose(jac),r)
      ! A scaled projected gradient gives a coordinatewise KKT convergence check.
      do j=1,16
        trial(j)=max(1.0e-12_wp,min(1.0_wp-1.0e-9_wp,w(j)+b(j)/max(a(j,j),tiny(1.0_wp))))
      enddo
      if(maxval(abs(trial-w))<p%fit_tolerance) then
        converged=.true.;exit
      endif
      do j=1,16
        a(j,j)=a(j,j)+damping
      enddo
      do j=1,16
        if((w(j)<=1.0e-9_wp.and.b(j)<0).or.(w(j)>=1-1.0e-8_wp.and.b(j)>0)) then
          a(:,j)=0;a(j,:)=0;a(j,j)=1;b(j)=0
        endif
      enddo
      call solve(a,b,delta,ok)
      if(.not.ok) exit
      moved=.false.;alpha=1
      do line=1,24
        trial=max(1.0e-12_wp,min(1.0_wp-1.0e-9_wp,w+alpha*delta))
        call residual(p,id,c%tau,trial,coarse,rho,cs,cp,f,rt,bk,sm,ok)
        if(ok) then
          if(sum(rt*rt)<cost) then
            moved=.true.;exit
          endif
        endif
        alpha=alpha/2
      enddo
      if(moved) then
        if(maxval(abs(trial-w))<p%fit_tolerance) converged=.true.
        w=trial;damping=max(1.0e-10_wp,damping/2)
      else
        ! Projected stationarity: no feasible improving step at the search resolution.
        if(maxval(abs(alpha*delta))<p%fit_tolerance) converged=.true.
        damping=damping*10
      endif
      if(converged) exit
    enddo
    message='CG8 constrained effective fit did not converge'
    if(.not.converged) return
    allocate(fv(1028),rv(2056))
    do i=1,1025
      fv(i)=p%fmin*(p%fmax/p%fmin)**(real(i-1,wp)/1024.0_wp)
    enddo
    fv(1026:1028)=max(p%fmin,min(p%fmax,p%f_transition* &
      [p%transition_lower_ratio,1.0_wp,p%transition_upper_ratio]))
    call residual(p,id,c%tau,w,coarse,rho,cs,cp,fv,rv,c%bulk,c%shear,ok)
    message='CG8 invalid constrained material response'
    if(.not.ok) return
    c%bulk_strength=w(1:8);c%shear_strength=w(9:16);c%max_error=maxval(abs(rv))
    if(c%max_error>p%max_fit_error) then
      write(message,'(A,ES12.4,A,ES12.4)') 'CG8 effective fit error=',c%max_error,', limit=',p%max_fit_error
      return
    endif
    status=0;message=''
  end subroutine
end module
