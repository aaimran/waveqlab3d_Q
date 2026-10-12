module interrupt_control
#ifdef WQL3D_SIGNAL_INTEL
  use ifport, only: signal
#elif defined(WQL3D_SIGNAL_GNU)
  use, intrinsic :: iso_c_binding, only: c_int
#endif
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi
  implicit none
  private
  public :: install_interrupt_handler, interruption_requested, finish_interrupted
  integer, save :: stop_reason = 0
  ! Independent flags avoid losing a signal if another handler interrupts it.
  integer, volatile, save :: int_received = 0, term_received = 0, usr1_received = 0
  integer, parameter :: sigint_number = 2, sigterm_number = 15
  integer, parameter :: sigusr1_number = WQL3D_SIGUSR1
#ifdef WQL3D_SIGNAL_INTEL
  public :: int_received, term_received, usr1_received
  public :: sigint_number, sigterm_number, sigusr1_number
  integer(4), external :: wql3d_signal_handler
#endif
contains
  subroutine install_interrupt_handler()
    integer :: ierr, rank
    logical :: failed, any_failed
#ifdef WQL3D_SIGNAL_INTEL
    failed = signal(sigint_number, wql3d_signal_handler, -1) < 0
    if (signal(sigterm_number, wql3d_signal_handler, -1) < 0) failed = .true.
    if (sigusr1_number > 0) then
      if (signal(sigusr1_number, wql3d_signal_handler, -1) < 0) failed = .true.
    end if
#elif defined(WQL3D_SIGNAL_GNU)
    integer :: status
    intrinsic :: signal
    call signal(sigint_number, record_interrupt, status)
    failed = status == -1
    call signal(sigterm_number, record_interrupt, status)
    if (status == -1) failed = .true.
    if (sigusr1_number > 0) then
      call signal(sigusr1_number, record_interrupt, status)
      if (status == -1) failed = .true.
    end if
#else
    failed = .false.
#endif
    call MPI_Allreduce(failed, any_failed, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (any_failed) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      if (rank == 0) write(error_unit, '(A)') 'Unable to install stop signal handlers.'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
  end subroutine

  ! Collective: call on every rank at the same point, outside signal handlers.
  logical function interruption_requested() result(requested)
    integer :: local_reason, global_reason, ierr
    local_reason = 0
    if (int_received /= 0) local_reason = 1
    if (usr1_received /= 0) local_reason = 2
    if (term_received /= 0) local_reason = 3
    call MPI_Allreduce(local_reason, global_reason, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    stop_reason = max(stop_reason, global_reason)
    requested = stop_reason /= 0
  end function

  ! Call collectively after closing any initialized simulation output files.
  subroutine finish_interrupted()
    integer :: rank, ierr, exit_status
#if defined(WQL3D_SIGNAL_GNU) || defined(WQL3D_SIGNAL_INTEL)
    intrinsic :: exit
#endif
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (rank == 0) then
      select case (stop_reason)
      case (3)
        write(error_unit, '(A)') 'Simulation stopped after SIGTERM (scheduler termination or cancellation).'
      case (2)
        write(error_unit, '(A)') 'Simulation stopped after SIGUSR1 (scheduler warning).'
      case default
        write(error_unit, '(A)') 'Simulation interrupted by Ctrl+C.'
      end select
      flush(error_unit)
    end if
    call MPI_Finalize(ierr)
    select case (stop_reason)
    case (3)
      exit_status = 128 + sigterm_number
    case (2)
      exit_status = 128 + sigusr1_number
    case default
      exit_status = 128 + sigint_number
    end select
#if defined(WQL3D_SIGNAL_INTEL) || defined(WQL3D_SIGNAL_GNU)
    ! Compiler runtime EXIT avoids per-rank Fortran STOP messages.
    call exit(exit_status)
#else
    stop 1
#endif
  end subroutine

#ifdef WQL3D_SIGNAL_GNU
  ! GNU SIGNAL invokes a C-convention callback implemented entirely in Fortran.
  subroutine record_interrupt(signum) bind(C)
    integer(c_int), value :: signum
    if (signum == sigint_number) int_received = 1
    if (signum == sigterm_number) term_received = 1
    if (sigusr1_number > 0 .and. signum == sigusr1_number) usr1_received = 1
  end subroutine
#endif
end module

#ifdef WQL3D_SIGNAL_INTEL
! IFPORT requires an external INTEGER(4) function with a reference argument.
integer(4) function wql3d_signal_handler(signum) result(status)
  use interrupt_control, only: int_received, term_received, usr1_received, &
       sigint_number, sigterm_number, sigusr1_number
  implicit none
  integer(4), intent(in) :: signum
  !DIR$ ATTRIBUTES DEFAULT :: wql3d_signal_handler
  if (signum == sigint_number) int_received = 1
  if (signum == sigterm_number) term_received = 1
  if (sigusr1_number > 0 .and. signum == sigusr1_number) usr1_received = 1
  status = 0
end function
#endif
