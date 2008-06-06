!>
!! \brief This module contains routines for calculating the ionization and temperature evolution of the entire grid (1D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 1D version similar to the 3D version.

module evolve

  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of the entire grid (1D).
     
  ! Version:
  ! 1D version similar to the 3D version.

  use precision, only: dp
  use file_admin, only: logf
  use my_mpi ! supplies all the MPI definitions
  use sizes, only: Ndim, mesh
  use grid, only: r,vol,dr
  use material, only: ndens, xh, temper
  use photonstatistics, only: state_before, calculate_photon_statistics, &
       photon_loss

  implicit none

  private

  public :: evolve1D !< evolve 1D grid

  !> H column density at the back end of the cell
  real(kind=dp),dimension(mesh) :: coldensh_out 

contains

  ! =======================================================================
  
  !> Calculates the evolution of the hydrogen ionization state
  !! and temperature for the whole 1D grid\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 21-Aug-2006 (f77/OMP: 13-Jun-2005)\n
  !! \b Version: Simple 1D
    
  subroutine evolve1D (dt)

    ! Calculates the evolution of the hydrogen ionization state
     
    ! Author: Garrelt Mellema
     
    ! Date: 21-Aug-2006 (f77/OMP: 13-Jun-2005)

    ! Version: Simple 1D
    
    ! History:
    ! 11-Jun-2004 (GM) : grid arrays now passed via common (in grid.h)
    !    and material arrays also (in material.h).
    ! 11-Jun-2004 (GM) : adapted for multiple sources.
    !  3-Jan-2005 (GM) : reintegrated with updated Ifront3D
    ! 20-May-2005 (GM) : split original eveolve0D into two routines
    ! 13-Jun-2005 (HM) : OpenMP version : Hugh Merz
    ! 31-Mar-2008 (GM) : clean up.

    !> The time step
    real(kind=dp),intent(in) :: dt

    !> Will contains the integer position of the cell being treated
    integer,dimension(Ndim) :: pos
      
    ! Loop variables
    integer :: i,j,k,l,nx,niter

    ! Flag variable (passed back from evolve0D_global)
    integer :: conv_flag

    ! End of declarations

    ! Initial state (for photon statistics)
    call state_before ()

    ! reset photon loss counter
    photon_loss=0.0
       
    ! Put column density to zero
    coldensh_out(:)=0.0

    ! Loop through grid
    do i=1,mesh
       pos(1)=i
       call evolve0D(dt,pos)
    end do

    ! Calculate photon statistics
    call calculate_photon_statistics (dt)

    return
  end subroutine evolve1D


  !=======================================================================

  !> Calculates the evolution of the hydrogen ionization state in a
  !! single cell.\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 21-Aug-2006 (20-May-2005, 5-Jan-2005, 02 Jun 2004)\n
  !! \b Version: single source

  subroutine evolve0D(dt,rtpos)
    
    ! Calculates the evolution of the hydrogen ionization state in a
    ! single cell.
    
    ! Author: Garrelt Mellema
    
    ! Date: 21-Aug-2006 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: single sources, fixed temperature
    
    ! Multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to phih_grid, but the ionization fractions are not updated. 

    use mathconstants, only: pi
    use tped, only: electrondens
    use doric_module, only: doric, coldens, coldensh_bndry
    use radiation, only: photoion, photrates
    use c2ray_parameters, only: epsilon,convergence1,convergence2
    use material, only: isothermal
    use thermalevolution, only: thermal

    implicit none

    real(kind=dp),parameter :: max_coldensh=2e19 !< column density for stopping chemisty
    
    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    !> The time step
    real(kind=dp),intent(in) :: dt 
    !> mesh position of cell being done
    integer,dimension(Ndim),intent(in) :: rtpos 
    
    logical :: finalpass
    integer :: nx,nd,nit,idim ! loop counters
    integer,dimension(Ndim) :: pos
    integer,dimension(Ndim) :: srcpos1
    real(kind=dp) :: coldensh_in
    real(kind=dp) :: coldensh_cell
    real(kind=dp) :: path
    real(kind=dp) :: de
    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0
    real(kind=dp) :: ndens_p
    real(kind=dp) :: avg_temper
    
    real(kind=dp) :: dist,vol_ph
    real(kind=dp) :: temper0,temper1,temper2
    real(kind=dp) :: yh_av0
    real(kind=dp) :: convergence

    type(photrates) :: phi

    convergence=convergence1

    pos(1)=rtpos(1)
    ! Initialize local ionization states to the global ones
    do nx=0,1
       yh(nx)=xh(pos(1),nx)
       yh0(nx)=xh(pos(1),nx)
       yh_av(nx)=xh(pos(1),nx) ! use calculated xh_av
    enddo
    ! Initialize local temperature and density
    avg_temper=temper(pos(1))
    ! Initialize scalars temper0 and temper1
    temper0=temper(pos(1))    ! original temperature as scalar
    temper1=temper0           ! will contain the new temperature

    ndens_p=ndens(pos(1))

    ! Find the column density at the entrance point of the cell (short
    ! characteristics)
    if (pos(1).eq.1) then
       coldensh_in = coldensh_bndry()
    else
       coldensh_in = coldensh_out(pos(1)-1)
    endif
    path=dr

    ! Find the distance to the source
    dist=r(pos(1))

    ! Find the volume of the shell this cell is part of (dilution factor)
    vol_ph=vol(pos(1))

    ! Iterate to get mean ionization state (column density / optical depth) 
    ! in cell
    nit=0
    do 
       nit=nit+1
       
       ! Save the value of yh_av found in the previous iteration
       yh_av0=yh_av(0)
       temper2=temper1

       ! Calculate (time averaged) column density of cell
       coldensh_cell=coldens(path,yh_av(0),ndens_p)
       
       ! Calculate (photon-conserving) photo-ionization rate
       call photoion(phi,coldensh_in,coldensh_in+coldensh_cell, &
            vol_ph)
       phi%h=phi%h/(yh_av(0)*ndens_p)

       ! Restore yh to initial values (for doric)
       yh(:)=yh0(:)
       
       ! Calculate (mean) electron density
       de=electrondens(ndens_p,yh_av)
       
       ! Calculate the new and mean ionization states (yh and yh_av)
       call doric(dt,avg_temper,de,ndens_p,yh,yh_av,phi)
       
       temper1=temper0       ! set temper1 to the original temperature

       if (nit.eq.1.and.abs(yh0(0)-yh(0)).gt.0.1) then
          continue
       else
          if (.not.isothermal) &
               call thermal(dt,temper1,avg_temper,de,ndens_p, &
               yh,yh_av,yh0,phi)
       endif
          
       ! Test for convergence on ionization fraction and temperature
       if ((abs((yh_av(0)-yh_av0)/yh_av(0)).lt.convergence &
            .or. &
            (yh_av(0).lt.1e-12)) &
            .and. &
            abs(temper1-temper2)/temper1.lt.convergence) exit
       
       ! Warn about non-convergence
       if (nit.gt.5000) then
          write(logf,*) 'Convergence failing'
          write(logf,*) 'xh: ',yh_av(0),yh_av0
          write(logf,*) 'temper: ',temper1,temper2
          exit
       endif
          
    enddo ! end of iteration
    
    ! Copy ionic abundances back
    xh(pos(1),:)=yh(:)
      
    ! Copy temperature back
    temper(pos(1))=temper1
      
    ! Add the (time averaged) column density of this cell
    ! to the total column density (for this source)
    coldensh_out(pos(1))=coldensh_in+ &
         coldens(path,yh_av(0),ndens_p)

    ! Photon statistics: register number of photons leaving the grid
    if (pos(1).eq.mesh) &
         photon_loss=photon_loss+ &
         phi%h_out*vol(pos(1))/vol_ph

  end subroutine evolve0D

end module evolve
