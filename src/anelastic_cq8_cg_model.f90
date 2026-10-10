module anelastic_cq8_cg_model
  use common, only: wp
  use anelastic_cg8_model, only: cg8_parameters,validate_cg8
  implicit none
  private
  type, public :: cq8_cg_parameters
    type(cg8_parameters) :: settings
  end type
  public :: read_cq8_cg_parameters
contains
  subroutine read_cq8_cg_parameters(infile,nblocks,p,status,message)
    integer,intent(in)::infile,nblocks
    type(cq8_cg_parameters),intent(out)::p
    integer,intent(out)::status
    character(*),intent(out)::message
    real(wp),allocatable::Qs0(:),Qp0(:)
    real(wp)::fref,fmin,fmax,fit_tolerance,max_fit_error
    integer::fit_samples,fit_max_iterations,buffer_layers,pattern_origin(3),stat
    character(len=32)::relaxation_policy,coefficient_policy,boundary_policy
    namelist /anelastic_cQ8_cg_list/ Qs0,Qp0,fref,fmin, &
      fmax,relaxation_policy,coefficient_policy,fit_samples, &
      fit_tolerance,fit_max_iterations,max_fit_error,boundary_policy, &
      buffer_layers,pattern_origin
    status=1;message='CG8 requires one or two blocks'
    if(nblocks<1.or.nblocks>2) return
    allocate(Qs0(nblocks),Qp0(nblocks));Qs0=-1;Qp0=-1
    fref=p%settings%fref
    fmin=p%settings%fmin
    fmax=p%settings%fmax
    relaxation_policy=p%settings%relaxation_policy
    coefficient_policy=p%settings%coefficient_policy
    fit_samples=p%settings%fit_samples
    fit_tolerance=p%settings%fit_tolerance
    fit_max_iterations=p%settings%fit_max_iterations
    max_fit_error=p%settings%max_fit_error
    boundary_policy=p%settings%boundary_policy
    buffer_layers=p%settings%buffer_layers
    pattern_origin=p%settings%pattern_origin
    rewind(infile)
    read(infile,nml=anelastic_cQ8_cg_list,iostat=stat)
    message='A valid &anelastic_cQ8_cg_list is required'
    if(stat/=0) return
    p%settings%Qs0(1:nblocks)=Qs0;p%settings%Qp0(1:nblocks)=Qp0
    p%settings%fref=fref
    p%settings%fmin=fmin
    p%settings%fmax=fmax
    p%settings%relaxation_policy=relaxation_policy
    p%settings%coefficient_policy=coefficient_policy
    p%settings%fit_samples=fit_samples
    p%settings%fit_tolerance=fit_tolerance
    p%settings%fit_max_iterations=fit_max_iterations
    p%settings%max_fit_error=max_fit_error
    p%settings%boundary_policy=boundary_policy
    p%settings%buffer_layers=buffer_layers
    p%settings%pattern_origin=pattern_origin
    call validate_cg8(p%settings,nblocks,status,message)
  end subroutine
end module
