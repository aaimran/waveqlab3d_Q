program interrupt_mpi_test
  use interrupt_control
  use mpi
#ifdef WQL3D_SIGNAL_INTEL
  use ifport, only: getpid
#endif
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  integer :: rank, ierr
  double precision :: started
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call install_interrupt_handler()
  write(output_unit, '(A,1X,I0,1X,I0)') 'READY', rank, getpid()
  flush(output_unit)
  started = MPI_Wtime()
  do
    if (interruption_requested()) call finish_interrupted()
    if (MPI_Wtime() - started > 10.0d0) call MPI_Abort(MPI_COMM_WORLD, 99, ierr)
  end do
end program
