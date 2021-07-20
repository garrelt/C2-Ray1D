!>
!! \brief This module contains astronomical constants and units
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2003-06-01
!!
!! \b Version: cgs units
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)

module astroconstants

  use precision, only: dp

  implicit none

  ! A collection of astronomical units and conversion factors
  ! Units: cgs
  
  real(kind=dp),parameter :: R_SOLAR=6.9599e10_dp !< Solar radius
  real(kind=dp),parameter :: L_SOLAR=3.826e33_dp !< Solar luminosity
  real(kind=dp),parameter :: M_SOLAR=1.98892d33 !< Solar mass
  
  real(kind=dp),parameter :: YEAR=3.15576E+07_dp !< Julian year

  real(kind=dp),parameter :: pc=3.086e18_dp !< parsec
  real(kind=dp),parameter :: kpc=1e3_dp*pc !< kiloparsec
  real(kind=dp),parameter :: Mpc=1e6_dp*pc !< megaparsec
  real(kind=dp),parameter :: lightyear=9.463e17_dp !< lightyear
  real(kind=dp),parameter :: AU=1.49597870E+13_dp !< Astronomical Unit

end module astroconstants
