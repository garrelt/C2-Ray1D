!>
!! \brief This module contains data and routines for cosmological problems
!!
!! Module for Capreole / C2-Ray (F90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2010-03-08 (but older than that)
!!
!! This module initializes the cosmology and contains functions for changing
!! redshift to time and vice versa. It also contains cosmological cooling
!! routines. Note that it does not evolve the proper lengths, volumes and
!! densities. This is done in the cosmological_evolution module.
!!
!! List of routines:
!! - cosmology_init: initializes cosmological time and sets lengths
!!            and volumes from comoving to proper scaling.
!! - time2zred: convert time to redshift
!! - zred2time: convert redshift to time
!! - cosmo_cool: cosmological adiabatic cooling rate
!! - compton_cool: Compton cooling wrt the CMB.
!!
!! This module keeps track of the current redshift (zred) and the flag
!! for cosmological calculations (cosmological).
!!

module cosmology


  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, file_input

  ! use cosmology_parameters

  implicit none

  real(kind=dp),parameter :: cmbtemp=2.726 !< CMB temperature

  logical :: cosmological !< true if cosmological simulation
  real(kind=dp) :: h !< Hubble constant (in 100 km/s/Mpc)
  real(kind=dp) :: Omega0 !< Total matter density (in critical density)
  ! Derived parameters
  real(kind=dp) :: H0 !< Hubble constant (cgs)
 
  real(kind=dp) :: zred_t0 !< initial redshift
  real(kind=dp) :: t0      !< time of initial redshift
  real(kind=dp) :: zred    !< current redshift
  real(kind=dp) :: zfactor !< scaling factor between two redshifts (used in cosmological_evolution)

contains
  ! =======================================================================

  !> Initialize cosmology and cosmological constants
  subroutine cosmology_init (cosmo_switch)
    
    use astroconstants, only: Mpc

    logical, intent(in) :: cosmo_switch !< cosmological problem or not?

    real(kind=dp) :: zred0

    ! Activate cosmological stuff
    cosmological = cosmo_switch

    if (cosmological) then
       ! Ask for cosmological parameters
       if (rank == 0) then
          if (.not.file_input) write(*,'(A,$)') 'Initial redshift?'
          read(stdinput,*) zred0
          if (.not.file_input) write(*,'(A,$)') 'Hubble constant? (0-1)'
          read(stdinput,*) h         
          if (.not.file_input) write(*,'(A,$)') 'Density parameter Omega0?'
          read(stdinput,*) Omega0
       endif
#ifdef MPI
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(zred0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
       call MPI_BCAST(h,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
       call MPI_BCAST(Omega0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif
       
       ! Derived parameters
       H0=h*100.0*1e5/Mpc ! Hubble constant (cgs)
       
       ! Cosmological time corresponding to (initial) redshift zred0
       t0 = 2.*(1.+zred0)**(-1.5)/(3.*H0*sqrt(Omega0))
       
       ! Initialize redshift
       zred_t0=zred0 ! keep initial redshift
       ! zred is the master redshift variable, but 
       ! needs to be zero so that grid and material variables will be
       ! initialized as comoving (z=0) and then changed to proper 
       ! (z=zred_t0).
       zred=0.0 
       
    endif

  end subroutine cosmology_init

  ! =======================================================================

  !> Calculates the cosmological redshift for a given time
  function time2zred (time)

    ! Calculates the cosmological redshift for a given time

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 21-May-2005)
    
    ! Version: f90

    ! History: - 20-Aug-2006, conversion to f90

    real(kind=dp) :: time2zred
    real(kind=dp),intent(in) :: time

    ! Calculate the redshift
    time2zred = -1+(1.+zred_t0)*(t0/(t0+time))**(2./3.)

  end function time2zred

  ! =======================================================================

  !> Calculates the time for a given cosmological redshift
  function zred2time (zred1)

    ! Calculates the time for a given cosmological redshift

    ! Author: Garrelt Mellema

    ! Date: 30-Sep-2006
    
    ! Version: f90

    ! History: 

    real(kind=dp) :: zred2time
    real(kind=dp),intent(in) :: zred1

    ! Calculate the redshift
    zred2time = t0*( ((1.0+zred_t0)/(1.0+zred1))**1.5 - 1.0 )

  end function zred2time

  ! =======================================================================

  !> Calculates the cosmological redshift from time
  !! and the scale factor zfactor for use in cosmo_evol
  subroutine redshift_evol (time)

    ! Calculates the cosmological redshift from time
    ! and the scale factor zfactor for use in cosmo_evol
  
    ! Author: Garrelt Mellema
    ! Date: 10-Mar-2006
    ! Version: F90

    ! History:
    ! - 19-Nov-2004: first version f77

    real(kind=dp),intent(in) :: time

    real(kind=dp) :: zred_prev

    ! Calculate the change since the previous redshift.
    ! Note: the initial redshift should be ZERO since
    ! the variables are initialized as comoving!
    zred_prev=zred
    zred = -1+(1.+zred_t0)*((t0+time)/t0)**(-2./3.)

    ! Take the average zfactor between zred_prev and zred
    zfactor=(1.0+zred_prev)/(1.+zred)

  end subroutine redshift_evol

  ! =======================================================================

  !> Calculates the cosmological adiabatic cooling
  real(kind=dp) function cosmo_cool (e_int)

    ! Calculates the cosmological adiabatic cooling

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: f90

    ! History:
    ! 19-Nov-2004: first version f77

    real(kind=dp),intent(in) :: e_int

    real(kind=dp) :: dzdt

    ! Cosmological cooling rate:
    ! 2*(da/dt)/a*e
    ! or
    ! 2*(dz/dt)/(1+z)*e
    ! with a the cosmological scale factor

    ! dz/dt (for flat LambdaCDM)
    dzdt=H0*(1.+zred)*sqrt(Omega0*(1.+zred)**3+1.-Omega0)

    !Cooling rate
    cosmo_cool=e_int*2.0/(1.0+zred)*dzdt

  end function cosmo_cool

  ! =======================================================================

  !> Calculates the (cosmological) Compton cooling rate (against the CMB)
  real(kind=dp) function compton_cool (temper,eldens)
    
    ! Calculates the (cosmological) Compton cooling rate
    
    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: first version
    
    ! History:
    ! 16-May-2005: first version f77
    

    ! parameter reference?

    real(kind=dp),intent(in) :: temper ! temperature
    real(kind=dp),intent(in) :: eldens ! electron density
    
    !Cooling rate
    compton_cool=5.65e-36*eldens*(1.0+zred)**4* &
         (temper-cmbtemp*(1.0+zred))
    
  end function compton_cool
  
end module cosmology
