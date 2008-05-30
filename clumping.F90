!>
!! \brief This module contains clumping data and routines
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:


module clumping_factor

  use precision
  use file_admin, only: stdinput
  use my_mpi

  implicit none

  !> Global clumping factor
  real(kind=dp),public :: clumping
  
contains

  !> Initializes the global clumping factor
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
