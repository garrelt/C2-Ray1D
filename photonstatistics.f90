!>
!! \brief This module handles the calculation of the photon statistics
!!
!! Module for C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 26-Sep-2006
!!
!!

module photonstatistics
  
  ! This module handles the calculation of the photon statistics
  ! For: C2-Ray

  ! Author: Garrelt Mellema

  ! Date: 26-Sep-2006

  ! Photon statistics
  ! photon_loss is a sum over all sources, summing is done in evolve0d.
  ! For parallelization consider a reduction, or may also be dropped
  ! entirely.
  ! Photon_losst contains the threaded values

  use precision, only: dp
  use cgsconstants, only: albpow,bh00,colh0,temph0
  use sizes, only: mesh
  use grid, only: vol
  use material, only: ndens, xh, temper, clumping
  use tped, only: electrondens
  !use subgrid_clumping, only: clumping

  logical,parameter :: do_photonstatistics=.true.
  real(kind=dp) :: totrec !< total number of recombinations
  real(kind=dp) :: totcollisions !< total number of collisional ionizations
  real(kind=dp) :: dh0 !< change in the number of neutral Hs
  real(kind=dp) :: total_ion !< total number of ionizations 
  real(kind=dp) :: grtotal_ion !< total number of ionization since start
  real(kind=dp) :: photon_loss !< photons lost from the grid

  real(kind=dp),private :: h0_before !< number of neutral H at begin of time step
  real(kind=dp),private :: h0_after !< number of neutral H at end of time step
  real(kind=dp),private :: h1_before !< number of ionized H at begin of time step
  real(kind=dp),private :: h1_after !< number of ionized H at end of time step
  integer,private :: i !< mesh loop variable
  integer,private :: j !< mesh loop variable
  integer,private :: k !< mesh loop variable

contains

  !> initializes the grand total of ionization to zero
  subroutine initialize_photonstatistics ()

    ! set total number of ionizing photons used to zero
    grtotal_ion=0.0

  end subroutine initialize_photonstatistics

  !> Calculates the photon statistics
  subroutine calculate_photon_statistics (dt)

    real(kind=dp),intent(in) :: dt !< time step

    ! Call the individual routines needed for this calculation

    call state_after () ! number of neutrals after integration
    call total_rates (dt) ! total photons used in balancing recombinations etc.
    call total_ionizations () ! final statistics
    
  end subroutine calculate_photon_statistics

  !> Calculates the number of neutrals and ions at the start of the time step
  subroutine state_before ()

    ! Photon statistics: calculate the number of neutrals before integration
    h0_before=0.0
    h1_before=0.0
    do i=1,mesh
       h0_before=h0_before+vol(i)*ndens(i)*xh(i,0)
       h1_before=h1_before+vol(i)*ndens(i)*xh(i,1)
    enddo
    
  end subroutine state_before

  !> Calculates total number of recombinations and collisions
  subroutine total_rates(dt)

    real(kind=dp),intent(in) :: dt !< time step

    real(kind=dp),dimension(0:1) :: yh
 
    ! Photon statistics: Determine total number of recombinations/collisions
    ! Should match the code in doric_module

    totrec=0.0
    totcollisions=0.0
    do i=1,mesh
       yh(0)=xh(i,0)
       yh(1)=xh(i,1)
       totrec=totrec+vol(i)*ndens(i)*xh(i,1)*    &
            electrondens(ndens(i),yh)*  &
            clumping*bh00*(temper(i)/1e4)**albpow
       totcollisions=totcollisions+vol(i)*ndens(i)*   &
            xh(i,0)*electrondens(ndens(i),yh)* &
            colh0*sqrt(temper(i))*exp(-temph0/temper(i))
    enddo

    totrec=totrec*dt
    totcollisions=totcollisions*dt

  end subroutine total_rates
  
  !> Calculates the number of neutrals and ions at the end of the time step
  subroutine state_after()
    
    ! Photon statistics: Calculate the number of neutrals after the integration
    h0_after=0.0
    h1_after=0.0
    do i=1,mesh
       h0_after=h0_after+vol(i)*ndens(i)*xh(i,0)
       h1_after=h1_after+vol(i)*ndens(i)*xh(i,1)
    enddo
    
  end subroutine state_after
  
  !> Calculates the total number of ionizations during a time step
  subroutine total_ionizations ()
    
    ! Photon statistics: Total number of new ionizations
    dh0=(h0_before-h0_after)
    total_ion=totrec+dh0
    
  end subroutine total_ionizations

end module photonstatistics
