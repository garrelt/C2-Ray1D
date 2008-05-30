!>
!! \brief This module contains parameters specific for C2-Ray
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! Note: this file contains parameters for different versions
!!       of C2-Ray (dimensions, single/multiple sources), so 
!!       not all parameters will be used.

module c2ray_parameters

  ! This module collects parameters needed by C2-Ray

  use precision
  use astroconstants

  implicit none

  !> Which fraction of the cells can be left unconverged in order
  !! to improve performance (used in rad_evolve3d)
  real(kind=dp),parameter :: convergence_fraction=1.5e-5

  ! Set to true to let C2-Ray not change the temperature
  ! set interactively in mat_ini
  !logical,parameter :: isothermal=.true.

  !> A really small number
  real(kind=dp),parameter :: epsilon=1e-40_dp

  !> Convergence criterion for per source calculation (evolve0d)
  real(kind=dp),parameter :: convergence1=1.0e-3

  !> Convergence criterion for global calculation (evolve0d)
  real(kind=dp),parameter :: convergence2=1.0e-2

  !> Parameters for nominal SED
  real(kind=dp),parameter :: teff_nominal=0.0
  real(kind=dp),parameter :: s_star_nominal=1e48_dp
  !real(kind=dp),parameter :: s_star_nominal=1e50_dp
  
  !> Subgrid clumping\n
  !! 1: constant clumping (with clumping_factor)\n
  !! 2: 3.5Mpc PM, WMAP1 clumping\n
  !! 3: 3.5Mpc PM, WMAP3 clumping\n
  !! 4: 1 Mpc P3M
  integer,parameter :: type_of_clumping=3
  !> Clumping factor if constant
  real,parameter :: clumping_factor=1.0  

  ! Cosmological cooling
  ! Set in cosmology module!!
  !logical,parameter :: cosmological=.true.

  !> Thermal: minimum temperature
  real(kind=dp),parameter :: minitemp=1.0 ! minimum temperature
  !> Thermal: fraction of the cooling time step below which no iteration is done
  real(kind=dp),parameter :: relative_denergy=0.1

  ! 
  !> Source properties: Photon per atom for high mass sources
  real,parameter :: phot_per_atom1=250.0
  !> Source properties: Photon per atom for low mass sources
  real,parameter :: phot_per_atom2=250.0
  !> Source properties: Life time of sources (if set at compile time)
  real,parameter :: lifetime=20e6*YEAR
  !> Source properties: Lower limit neutral fraction for suppression criterion
  real,parameter :: StillNeutral=0.9 ! lower limit of neutral criterium
  !> Source properties: Upper limit for low mass sources (not used)
  real,parameter :: LowMassLimit=1e9 ! in solar masses

end module c2ray_parameters
