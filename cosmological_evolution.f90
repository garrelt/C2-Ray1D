!>
!! \brief This module contains routines for cosmological evolution
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 04-Mar-2006
!!
!! The density and lengths and volumes evolve due to the expansion
!! of the Universe. The cosmo_evol routines evolves the geometrical
!! variables (r, dr, vol) and the material density variables (ndens)
!! according to the current redshift (zfactor from the cosmology module).
!!
!! \b Programming \b note: this is probably a separate module from the cosmology
!! module because of
!! dependency problems. For the cosmo_evol routines, the grid and material
!! variables have to be available, but for those the cosmology module has
!! to be available. As cosmo_evol really belongs within the cosmology module
!! a different solution would be desirable (GM, 20100308).

module cosmological_evolution

  ! This module contains routines having to do with the cosmological 
  ! evolution of state and grid variables
  
  ! - cosmo_evol: cosmological evolution of space, density

  use precision, only: dp
  use cosmology, only: zfactor
  use grid, only: r,dr,vol
  use material, only: ndens
    
  implicit none

contains

  ! =======================================================================

  !> Calculates the cosmological evolution of space and densities\n
  !! Changes variables: ndens, r, dr, vol
  subroutine cosmo_evol ()

    ! Calculates the cosmological evolution of space and densities

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: F90 first version

    ! History:
    ! - 19-Nov-2004: first version f77

    real(kind=dp) :: zfactor3

    zfactor3=zfactor*zfactor*zfactor

    ! Change the grid coordinates
    r(:)=r(:)*zfactor

    dr=dr*zfactor
    
    vol(:)=vol(:)*zfactor3

    ! Change the densities
    ndens(:)=ndens(:)/zfactor3

  end subroutine cosmo_evol

end module cosmological_evolution
