!>
!! \brief This module contains data and routines for handling the material properties on the grid (1D)
!!
!! These properties are; density, temperature, clumping, ionization fractions
!! 
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!
!! \b Version:  1D test problems\n
!! Problem 1: constant density (Strömgren problem)\n
!! Problem 2: 1/r density\n
!! Problem 3: 1/r^2 density (with core of radius r_core)\n
!! Problem 4: cosmological constant density (Shapiro & Giroux problem)\n


module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start

  ! Version: 1D test problems

  ! Problem 1: constant density (Strömgren problem)
  ! Problem 2: 1/r density
  ! Problem 3: 1/r^2 density (with core of radius r_core)
  ! Problem 4: cosmological constant density (Shapiro & Giroux problem)

  use precision, only: dp,si
  use cgsconstants, only: bh00
  use astroconstants, only: YEAR
  use sizes, only: mesh
  use file_admin, only: stdinput
  use my_mpi
  use grid, only: r,vol
  use cgsconstants, only: m_p
  use abundances, only: mu
  use cosmology, only: cosmology_init,H0,t0,zred_t0

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! temper - temperature (K) of a cell
  ! xh - ionization fractions for one cell
  real(kind=dp) :: ndens(mesh) !< number density (cm^-3) of a cell
  real(kind=dp) :: temper(mesh) !< temperature (K) of a cell
  real(kind=dp) :: xh(mesh,0:1) !< ionization fractions for one cell
  real(kind=dp) :: clumping !< global clumping factor
  real(kind=dp) :: r_core !< core radius (for problems 2 and 3) 
  real(kind=dp) :: dens_core !< core density (for problems 2 and 3)
  integer :: testnum !< number of test problem (1 to 4)
  logical :: isothermal !< is the run isothermal?
  ! needed for analytical solution of cosmological Ifront
  real(kind=dp) :: t1 !< parameter for analytical solution of test 4 
  real(kind=dp) :: eta !< parameter for analytical solution of test 4  

#ifdef MPI
  integer,private :: ierror !< MPI error flag
#endif

contains

  ! ============================================================================

  !> Initializes material properties on grid\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))\n
  !! \b Version:
  !! - 1D\n
  !! - Four different test problems\n
  !! - Initially completely neutral\n

  subroutine mat_ini (restart)

    ! Initializes material properties on grid

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))

    ! Version: 
    ! - 1D
    ! - Four different test problems
    ! - Initially completely neutral

    integer,intent(out) :: restart !< will be /= 0 if a restart is intended

    integer :: i,n ! loop counters
    real(kind=dp) :: dens_val
    real(kind=dp) :: temper_val
    real(kind=dp) :: alpha
    character(len=1) :: answer

    ! restart
    restart=0 ! no restart by default

    ! Ask for input
    if (rank == 0) then
       write(*,'(A,$)') 'Which test? (1-4): '
       read(stdinput,*) testnum
#ifdef MPI
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(testnum,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
    endif

    ! Set alpha according to test problem
    select case (testnum)
    case(1,4) 
       alpha=0.0
    case(2) 
       alpha=-1.0
    case(3) 
       alpha=-2.0
    end select

    if (rank == 0) then
       if (testnum.eq.1.or.testnum.eq.4) then
          write(*,'(A,$)') 'Enter density (cm^-3): '
          read(stdinput,*) dens_val
       elseif (testnum.eq.2.or.testnum.eq.3) then
          write(*,'(A,$)') 'Enter reference (core) radius (cm): '
          read(stdinput,*) r_core
          write(*,'(A,$)') 'Enter density at reference (core)', &
               ' radius(cm^-3): '
          read(stdinput,*) dens_val
       endif
       
       write(*,'(A,$)') 'Enter clumping factor: '
       read(stdinput,*) clumping
       write(*,'(A,$)') 'Enter initial temperature (K): '
       read(stdinput,*) temper_val
       write(*,'(A,$)') 'Isothermal? (y/n): '
       read(stdinput,*) answer
       ! Isothermal?
       if (answer.eq.'y'.or.answer.eq.'Y') then
          isothermal=.true.
       else
          isothermal=.false.
       endif
    endif
#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(dens_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    if (testnum.eq.2.or.testnum.eq.3) &
         call MPI_BCAST(r_core,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(clumping,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(isothermal,1,MPI_LOGICAL,0,MPI_COMM_NEW,ierror)
#endif

       
    ! For test problem 4: cosmological parameters
    if (testnum.eq.4) then
       call cosmology_init (.true.)
    else
       call cosmology_init (.false.)
    endif

    ! Assign density and temperature to grid
    
    select case (testnum)
    case(1)
       do i=1,mesh
          ndens(i)=dens_val
          temper(i)=temper_val
       enddo
       
    case(2,3)
       dens_core=dens_val
       do i=1,mesh
          !     This is an attempt to make the initial conditions more
          !     correct: the density of a cell is the integrated density
          !     (mass) divided by the volume. Seems to give a worse fit
          !     to the analytical solution (slightly)
          !              rl=r(i)-0.5*dr
          !              rr=r(i)+0.5*dr
          !              ndens(i)=4.0*pi*dens_val*r_core**(-alpha)/vol(i)*
          !     $            (rr**(3.0+alpha)-rl**(3.0+alpha))/(3.0+alpha)
          
          !     This is just a straight sampling of the density distribution
          !     using the value at the cell centre.
          if (testnum.eq.3.and.r(i).le.r_core) then
             ! Flat core for test 3
             ndens(i)=dens_val
          else
             ndens(i)=dens_val*(r(i)/r_core)**alpha
          endif
          temper(i)=temper_val
       enddo
       
    case(4)
       ! For cosmological simulations, mean IGM
       dens_core=dens_val    ! save initial density value
       ! Parameters needed for the analytical solution
       ! Recombination time at z0
       t1 = 1./(bh00*clumping*dens_core) 
       eta = t0/t1*(1.+zred_t0)**3
       
       ! Set the density to the comoving value 
       ! (as is the spatial coordinate)
       ! evol_cosmo will set it to proper values.
       do i=1,mesh
          ndens(i)=dens_val
          temper(i)=temper_val
       enddo
       dens_val=dens_val*(1.+zred_t0)**3 !otherwise recombination time 
       !scale below would be wrong
    end select
    
    
    ! Assign ionization fractions (completely neutral)
    do i=1,mesh
       xh(i,0)=1.0-1.2e-3
       xh(i,1)=1.2e-3
    enddo
    
    write(*,'(A,1pe10.3,A)') 'Recombination time scale: ', &
         1.0/(dens_val*clumping*bh00*YEAR),' years'
    
  end subroutine mat_ini

end module material
