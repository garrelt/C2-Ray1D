!>
!! \brief Main program for C2Ray-1D
!!
!! \b Author: Garrelt Mellema \n
!! \b Date: 23-Sep-2006
!<
Program C2Ray

  ! Author: Garrelt Mellema

  ! Date: 23-Sep-2006

  ! Goal:
  ! One dimensional photo-ionization calculation for a series of
  ! test problems.
  
  ! Does not include hydrodynamics
  ! Assumes constant time step

  ! Needs following `modules'
  ! c2ray_parameters : all the tunable parameters for the code
  ! my_mpi : sets up the MPI parallel environment
  ! output_module : output routines
  ! grid : sets up the grid
  ! radiation : radiation tools
  ! pmfast : interface to PMFAST output
  ! cosmology : cosmological utilities
  ! material : material properties
  ! times : time and time step utilities
  ! sourceprops : source properties
  ! evolve : evolve grid in time

  use precision, only: dp
  use astroconstants, only: YEAR
  use my_mpi, only: mpi_setup, mpi_end, rank
  use output_module, only: setup_output,output,close_down
  use grid, only: grid_ini
  use radiation, only: rad_ini
  use cosmology, only: cosmology_init, redshift_evol, &
       time2zred, zred2time, zred, cosmological
  use cosmological_evolution, only: cosmo_evol
  use material, only: mat_ini, testnum
  use times, only: time_ini
  use evolve, only: evolve1D
  use file_admin, only:stdinput, flag_for_file_input
  use times, only: end_time,dt,output_time

#ifdef XLF
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

  ! CPU time variables
  real :: tstart !< Start time for CPU report
  real :: tend !< End time for CPU report

  ! Wall clock time variables
  integer :: cntr1 !< Start time wall clock
  integer :: cntr2 !< End time wall clock
  integer :: countspersec !< counts per second (for wall clock time)

  ! 
  integer :: nstep
  integer :: restart,nz,flag

  ! Time variables
  real(kind=dp) :: time !< actual time (s)
  real(kind=dp) :: next_output_time !< time of next output (s)
  real(kind=dp) :: actual_dt !< actual time step (s)

  !> Input file
  character(len=512) :: inputfile

  ! Initialize cpu timer
  call cpu_time(tstart)

  ! Initialize wall cock times
  call system_clock(cntr1)

  ! Set up MPI structure
  call mpi_setup()

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (iargc() > 0) then
     call getarg(1,inputfile)
     if (rank == 0) then
        write(*,*) 'reading input from ',trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
     endif
  endif

  ! Initialize output
  call setup_output()

  ! Initialize grid
  call grid_ini()

  ! Initialize the material properties
  call mat_ini (restart)

  ! Initialize photo-ionization calculation
  call rad_ini( )

  ! Initialize time step parameters
  call time_ini ()

  ! Set time to zero
  time=0.0
  next_output_time=0.0

  ! Initialize cosmology
  if (cosmological) then
     call redshift_evol(time)
     call cosmo_evol( )
     write(*,*) zred
  endif
  !call cosmology_init(zred_array(1),time)

  ! Loop until end time is reached
  nstep=0
  do
  
     ! Write output
     if (abs(time-next_output_time).le.1e-6*time) then
        call output(time,dt,end_time)
        next_output_time=next_output_time+output_time
     endif
     
     ! Make sure you produce output at the correct time
     ! dt=YEAR*10.0**(min(5.0,(-2.0+real(nstep)/1e5*10.0)))
     actual_dt=min(next_output_time-time,dt)
     nstep=nstep+1

     ! Report time and time step
     write(30,'(A,2(1pe10.3,1x),A)') 'Time, dt:', &
          time/YEAR,actual_dt/YEAR,' (years)'
     
     ! For cosmological simulations evolve proper quantities
     if (cosmological) then
        call redshift_evol(time+0.5*actual_dt)
        call cosmo_evol()
     endif

     ! Take one time step
     call evolve1D(actual_dt)

     ! Update time
     time=time+actual_dt
            
     if (abs(time-end_time).lt.1e-6*end_time) exit
  enddo

  ! stop
  ! Scale to the current redshift
  if (cosmological) then
     call redshift_evol(time)
     call cosmo_evol()
  endif

  ! Write final output
  call output(time,dt,end_time)

  ! Clean up some stuff
  call close_down ()
  
  ! Find out CPU time
  call cpu_time(tend)
  call system_clock(cntr2,countspersec)

  write(30,*) 'CPU time: ',tend-tstart,' s'
  write(30,*) 'Wall clock time: ',(cntr2-cntr1)/countspersec,' s'

  ! End the run
  call mpi_end ()

end Program C2Ray
