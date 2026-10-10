program test
 use common,only:wp
 use anelastic_cg8_model
 use anelastic_cq8_cg_model
 use anelastic_fq8_cg_model
 implicit none
 type(cg8_parameters)::p
 type(cg8_coefficients)::c
 type(cq8_cg_parameters)::pc
 type(fq8_cg_parameters)::pf
 integer::status,i,j,unit
 character(len=256)::message,path
 complex(wp)::rs,rp
 real(wp)::cs,cp,reference_cp
 do j=1,6
  p=cg8_parameters();reference_cp=2
  if(j>=5) reference_cp=sqrt(3.0_wp)
  if(j==1) then
   p%Qs0=20.0_wp;p%Qp0=20.0_wp
  elseif(j==2) then
   p%Qs0=50.0_wp;p%Qp0=50.0_wp
  elseif(j==3) then
   p%Qs0=400.0_wp;p%Qp0=400.0_wp
  else
   p%Qs0=400.0_wp;p%Qp0=600.0_wp;p%gamma=.3_wp;p%transition_policy='smooth'
   p%fmin=.1_wp;p%fmax=10.0_wp;p%relaxation_policy='withers-times'
   if(j==5) then
     p%gamma=0;p%relaxation_policy='band'
   endif
  endif
  do i=1,2
   call build_cg8_coefficients(p,1,1.0_wp,1.0_wp,reference_cp,i==1,c,status,message)
   print *, 'case=',j,'coarse=',i==1,'status=',status,'error=',c%max_error,trim(message)
   if(status/=0) error stop 'CG8 fit failed'
   if(c%bulk<=0.or.c%shear<=0) error stop 'invalid moduli'
   if(any(c%bulk_strength<0).or.any(c%shear_strength<0)) error stop 'negative relaxation'
   rs=c%shear*cg8_response(p%fref,c%tau,c%shear_strength,c%coarse)
   rp=c%bulk*cg8_response(p%fref,c%tau,c%bulk_strength,c%coarse)+(4.0_wp/3)*rs
   cs=1/real(sqrt(1/rs),wp);cp=1/real(sqrt(1/rp),wp)
   if(abs(cs-1)>1e-8_wp.or.abs(cp/reference_cp-1)>1e-8_wp) error stop 'reference phase speed mismatch'
   write(path,'(A,I0,A,I0,A)') 'cg8-case-',j,'-',i,'.txt'
   open(newunit=unit,file=trim(path),status='replace')
   write(unit,*) p%Qs0(1),p%Qp0(1),p%gamma,p%f_transition,p%transition_lower_ratio,p%transition_upper_ratio, &
     p%fmin,p%fmax,reference_cp
   write(unit,*) trim(p%transition_policy)
   write(unit,*) c%bulk,c%shear,c%max_error
   write(unit,*) c%tau
   write(unit,*) c%bulk_strength
   write(unit,*) c%shear_strength
   close(unit)
  enddo
 enddo
 ! Independent readers: cQ and fQ with gamma=0 must agree.
 open(newunit=unit,status='scratch',form='formatted')
 write(unit,'(A)') '&anelastic_cQ8_cg_list Qs0=50, Qp0=50 /'
 call read_cq8_cg_parameters(unit,1,pc,status,message)
 if(status/=0) error stop 'constant reader'
 close(unit)
 open(newunit=unit,status='scratch',form='formatted')
 write(unit,'(A)') '&anelastic_fQ8_cg_list Qs0=50, Qp0=50, gamma=0 /'
 call read_fq8_cg_parameters(unit,1,pf,status,message)
 if(status/=0) error stop 'frequency reader'
 close(unit)
 if(any(pc%settings%Qs0/=pf%settings%Qs0).or.pc%settings%gamma/=pf%settings%gamma) &
   error stop 'gamma zero configuration mismatch'
 p%Qs0=15;p%Qp0=1000;p%gamma=0
 call build_cg8_coefficients(p,1,1.0_wp,1.0_wp,2.0_wp,.true.,c,status,message)
 if(status==0.or.index(message,'incompatible')==0) error stop 'incompatible targets accepted'
 p%Qs0=50;p%Qp0=50;p%fit_max_iterations=1
 call build_cg8_coefficients(p,1,1.0_wp,1.0_wp,2.0_wp,.true.,c,status,message)
 if(status==0.or.index(message,'converge')==0) error stop 'missing convergence failure'
end program
