program cgt_coefficients_test
 use common,only:wp
 use anelastic_cg8_model,only:cg8_parameters,cg8_coefficients
 use anelastic_cg_t_model,only:build_cgt_coefficients
 use anelastic_cq_cg_t_model
 use anelastic_fq_cg_t_model
 implicit none
 type(cq_cg_t_parameters)::cq
 type(fq_cg_t_parameters)::fq
 type(cg8_parameters)::p
 type(cg8_coefficients)::a,b
 integer::u,status,kind
 character(512)::message
 open(newunit=u,status='scratch')
 write(u,'(A)') '&anelastic_cQ_cg_t_list Qs0=400,Qp0=600,fmin=.1,fmax=10,n_mechanisms=8/'
 call read_cq_cg_t_parameters(u,1,cq,status,message)
 if(status/=0) error stop 'cQ projected reader'
 close(u)
 open(newunit=u,status='scratch')
 write(u,'(A)') '&anelastic_fQ_cg_t_list Qs0=400,Qp0=600,fmin=.1,fmax=10,gamma=0/'
 call read_fq_cg_t_parameters(u,1,fq,status,message)
 if(status/=0) error stop 'fQ projected reader'
 close(u)
 call build_cgt_coefficients(cq%settings,1,1.0_wp,1.0_wp,2.0_wp,a,status,message)
 if(status/=0) then
  print *,trim(message);error stop 'cQ additive fit'
 endif
 call build_cgt_coefficients(fq%settings,1,1.0_wp,1.0_wp,2.0_wp,b,status,message)
 if(status/=0) error stop 'fQ gamma zero fit'
 if(any(a%tau/=b%tau).or.any(a%bulk_strength/=b%bulk_strength).or. &
    any(a%shear_strength/=b%shear_strength).or.a%bulk/=b%bulk.or.a%shear/=b%shear) error stop 'gamma zero identity'
 if(a%coarse.or.b%coarse) error stop 'Harmonic nodal coefficient branch used'
 do kind=1,3
  open(newunit=u,status='scratch')
  select case(kind)
  case(1)
   write(u,'(A)') '&anelastic_cQ_cg_t_list Qs0=400,Qp0=600,n_mechanisms=4/'
  case(2)
   write(u,'(A)') "&anelastic_cQ_cg_t_list Qs0=400,Qp0=600,coefficient_policy='constrained-effective-fit'/"
  case(3)
   write(u,'(A)') '&anelastic_cQ_cg_t_list Qs0=400,Qp0=600,max_fit_error=.03/'
  end select
  call read_cq_cg_t_parameters(u,1,cq,status,message)
  if(status==0) error stop 'Unsupported projected coefficient setting accepted'
  close(u)
 enddo
 p=fq%settings;p%gamma=.3_wp;p%transition_policy='smooth';p%relaxation_policy='withers-times';p%max_fit_error=.02_wp
 call build_cgt_coefficients(p,1,1.0_wp,1.0_wp,sqrt(3.0_wp),b,status,message)
 if(status/=0) then
  print *,trim(message);error stop 'Frequency-dependent projected dense fit'
 endif
 p%Qp0(1)=4000
 call build_cgt_coefficients(p,1,1.0_wp,1.0_wp,2.0_wp,b,status,message)
 if(status==0) error stop 'Incompatible P/S projected fit accepted'
 p=fq%settings;p%fit_max_iterations=1
 call build_cgt_coefficients(p,1,1.0_wp,1.0_wp,2.0_wp,b,status,message)
 if(status==0) error stop 'Failed projected fit convergence accepted'
 print *, 'CG-T readers, additive fit, dense verification, gamma-zero, and rejection tests passed'
end program
