!>
!! \brief This module contains parameters having to do
!!  with file in/ouput
!!
!! Module for C2-Ray/Capreole (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2008-05-27
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)
!<
module file_admin

  implicit none

  private
  
  !> set to not 5 if input via file (not implemented in all codes, check main program!)  
  integer,public,parameter :: stdinput=5 
  integer,public,parameter :: logf=30 !< unit number of log file(s)
  integer,public,parameter :: ah3=40  !< unit number of ah3 files (if applicable)

end module file_admin
