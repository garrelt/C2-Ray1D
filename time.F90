!>
!! \brief This module handles the time variables
!!
!! Module for C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2008-05-29 (no date before)
!<
module times

  ! This module handles the time variables

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput
  use astroconstants, only: YEAR

  implicit none

  private

  real(kind=dp),public :: end_time !< End time of simulation
  real(kind=dp),public :: dt !< Time step
  real(kind=dp),public :: output_time !< Time interval between outputs
  
contains

  ! =======================================================================

  !>
  !!  Initializes the module variables
  !<
  subroutine time_ini( )
    
    ! Initializes number of time steps per frame (integration and output)

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 10-Mar-2004)

    real(kind=dp) :: end_time_yrs,output_time_yrs
    integer :: number_timesteps

#ifdef MPI
    integer :: ierror
#endif

    if (rank == 0) then
       ! Ask for end time
       write(*,'(A,$)') 'Enter end time of calculation (years): '
       read(stdinput,*) end_time_yrs

       !  Convert to seconds
       end_time=end_time_yrs*YEAR

       ! Ask for number of time steps
       write(*,'(A,$)') 'Enter number of time steps: '
       read(stdinput,*) number_timesteps

       ! Ask for interval between outputs
       write(*,'(A,$)') 'Enter time interval between outputs (years): '
       read(stdinput,*) output_time_yrs
       ! Convert to seconds
       output_time=output_time_yrs*YEAR

    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(end_time,1,MPI_DOUBLEPRECISION,0,&
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(output_time,1,MPI_DOUBLEPRECISION,0,MPI_COMM_NEW,ierror)
#endif
    
    ! Set value of time step
    dt=end_time/real(number_timesteps)

  end subroutine time_ini

end module times
