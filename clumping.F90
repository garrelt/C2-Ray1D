module clumping_factor

  use precision
  use file_admin, only: stdinput
  use my_mpi

  implicit none

  real(kind=dp),public :: clumping
  
contains

  subroutine clumping_init()

#ifdef MPI       
    integer :: ierror
#endif
    
    if (rank == 0) then
       write(*,'(A,$)') 'Enter clumping factor: '
       read(stdinput,*) clumping
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(clumping,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
#endif

  end subroutine clumping_init

end module clumping_factor
