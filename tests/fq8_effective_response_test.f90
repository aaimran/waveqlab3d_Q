program fq8_effective_response_test
 use common,only:wp
 use anelastic_fq8_model
 use withers_tables,only:get_relaxation_times,get_withers_weights
 implicit none
 type(fq8_parameters)::p
 real(wp)::tau(8),ws(8),wpv(8),q
 integer::unit,status
 character(len=256)::message
 ! Published-table harmonic response remains mathematical test data, not an fQ8 runtime mode.
 call get_relaxation_times(.6_wp,tau)
 call get_withers_weights(.6_wp,50.0_wp,ws)
 q=fq8_effective_q(1.0_wp,tau,ws)
 if(abs(q-54.8324966843126_wp)>1e-10_wp) error stop 'table harmonic reference changed'
 if(abs(ws(1)+.000688_wp)>1e-14_wp) error stop 'negative table reference changed'
 if(abs(fq8_phase_velocity_ratio(1.0_wp,1.0_wp,tau,ws)-1)>1e-13_wp) error stop 'harmonic phase reference'
 p%Qs0=50;p%Qp0=80;p%gamma=.6_wp;p%f_transition=1;p%fref=1
 call build_fq8_coefficients(p,tau,ws,wpv,status,message)
 if(status/=0.or.any(ws<0).or.any(wpv<0).or.sum(ws)>=1.or.sum(wpv)>=1) error stop 'full-layout fit'
 open(newunit=unit,status='scratch',action='readwrite')
 write(unit,'(A)') '&anelastic_fQ8_list Qs0=50, Qp0=80, gamma=.6, f_transition=1, fref=1, coarse_grain=0 /'
 call read_fq8_parameters(unit,p,status,message);close(unit)
 if(status/=0) error stop 'legacy explicit full mode compatibility'
 open(newunit=unit,status='scratch',action='readwrite')
 write(unit,'(A)') '&anelastic_fQ8_list Qs0=50, Qp0=80, gamma=.6, f_transition=1, fref=1, coarse_grain=2 /'
 call read_fq8_parameters(unit,p,status,message);close(unit)
 if(status==0.or.index(message,'anelastic-fQ8-cg')==0) error stop 'legacy coarse migration missing'
 p%coefficient_method='withers-2015'
 call build_fq8_coefficients(p,tau,ws,wpv,status,message)
 if(status==0.or.index(message,'anelastic-fQ8-cg')==0) error stop 'old Withers runtime accepted'
 print *, 'Full fQ8 and published-reference tests passed'
end program
