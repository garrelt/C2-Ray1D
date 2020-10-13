!>
!! \brief This module data and routines which deal with radiative
!!     effects. 
!!
!!  Its main part deal with photo-ionizing radiation, but it
!!     also initializes other radiative properties, such as cooling (which
!!     are contained in different modules).
!!     It can be used in hydrodynamic or stand-alone radiative transfer 
!!     calculations.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2010-Mar-08 (but older)
!!
!! \b Version: 1D version similar to the 3D version.
!!
!! <b>Programming note:</b> This version for the 1D code also sets and
!! contains the source properties. In the 3Dm version there is a separate
!! sourceprops module for that. It would be nice to be consistent and also
!! have a source properties module for the 1D version. However, for the 1D
!! version the effective temperature is a real input variable and the
!! photo-ionization rates cannot be calculated before this variable is
!! set. This means that sourceprops would have to be called before rad_ini
!! (unlike in the 3Dm version where it has to be called after). Also, both
!! sourceprops and radiation would use the romberg integrator, so it should
!! be clear who initializes it (unless the romberg module checks for this
!! itself).
!!

module radiation
  
  !     This module contains data and routines which deal with radiative
  !     effects. Its main part deal with photo-ionizing radiation, but it
  !     also initializes other radiative properties, such as cooling (which
  !     are contained in different modules).
  !     It can be used in hydrodynamic or stand-alone radiative transfer 
  !     calculations.
  !
  !     Author: Garrelt Mellema
  ! 
  !     Date: 02-Nov-2012 (31-Jan-2008 (02-Jun-2004 (04-Mar-2004)
  
  ! Version
  ! Simplified version
  ! - Only hydrogen
  ! - Option for grey photo-ionization cross section
  ! - MPI enabled (broadcasts of radiative parameters to all nodes).

  ! Notes:
  ! - the initialization of the radiative cooling does not really belong
  !   here.
  ! - isothermal is sometimes an input parameter, and sometimes a compile
  !   time parameter. This needs to be streamlined. Probably along similar
  !   lines as the stellar parameters are dealt with.

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf, file_input
  use mathconstants, only: pi
  use cgsconstants, only: sigmasb, hplanck, kb, tpic2
  use cgsphotoconstants, only: frth0, frtop1, frtop2, sh0, betah0, sigh
  use astroconstants, only: R_SOLAR, L_SOLAR
  use romberg, only: scalar_romberg,vector_romberg,romberg_initialisation
  use c2ray_parameters, only: teff_nominal, S_star_nominal!, isothermal
  use material, only: isothermal

  implicit none

  !-----------------------------------------------------------------------
  !     NumFreq - Number of integration points in one of the three 
  !               frequency interval.
  !     NumTau -  Number of table points for the optical depth.
  !     NumFreqBnd - Number of frequency bands (1 for hydrogen only)
  !-----------------------------------------------------------------------

  !> Number of integration points in one of the three frequency interval.
  integer,parameter :: NumFreq=128
  !> Number of table points for the optical depth.
  integer,parameter :: NumTau=2000
  !> Number of frequency bands (1 for hydrogen only)
  integer,parameter :: NumFreqBnd=1
  
  !> This parameter sets the optical depth at the entrance of the grid.
  !> It can be used if radiation enters the simulation volume from the
  !> outside.
  real(kind=dp) :: tauHI=0.0

  ! Parameters defining the optical depth entries in the table.
  ! minlogtau is log10(lowest optical depth) (table position 1)
  ! maxlogtau is log10(highest optical depth) (table position NumTau)
  ! dlogtau is the step size in log10(tau) between table entries
  real(kind=dp),parameter :: minlogtau=-20.0 !< log10(lowest optical depth) 
  real(kind=dp),parameter :: maxlogtau=4.0 !< log10(highest optical depth)
  !> step size in log10(tau) between table entries
  real(kind=dp),parameter :: dlogtau=(maxlogtau-minlogtau)/real(NumTau)
  !> Optical depth array
  real(kind=dp) :: tau(0:NumTau)

  !> Logical that determines the use of grey opacities
  logical,parameter :: grey=.false. ! use grey opacities?

  ! stellar properties
  real(kind=dp) :: teff !< Black body effective temperature
  real(kind=dp) :: rstar !< Black body radius
  real(kind=dp) :: lstar !< Black body luminosity
  real(kind=dp) :: S_star !< Black body ionizing photons rate

  !> Frequency steps for integration
  real(kind=dp),dimension(NumFreqBnd) :: steph0 

  ! Photo-ionization integrals (rates)
  !> photo-ionization integral for H0 (optically thick case)
  real(kind=dp),dimension(:,:),allocatable  :: hphot
  !> photo-ionization heating integral for H0 (optically thick case)
  real(kind=dp),dimension(:,:),allocatable  :: hheat
  !> photo-ionization integral for H0 (optically thin case)
  real(kind=dp),dimension(:,:),allocatable  :: hphot1
  !> photo-ionization heating integral for H0 (optically thin case)
  real(kind=dp),dimension(:,:),allocatable  :: hheat1

  !> This type contains all the photo-ionization rates
  !> The in and out rates are used to ensure photon-conservation.
  !> See the C2-Ray paper.
  type photrates
     real(kind=dp) :: h        !< total H ionizing rate
     real(kind=dp) :: hv_h     !< total H heating rate
     real(kind=dp) :: h_in     !< in-rate
     real(kind=dp) :: hv_h_in  !< in-heating rate
     real(kind=dp) :: h_out    !< out-rate
     real(kind=dp) :: hv_h_out !< out-heating rate
  end type photrates

  ! photo-ionization rates (disabled as they are passed as arguments)
  !real(kind=dp),public :: phih,hvphih
  !real(kind=dp),public :: phih_in,phih_out
  !real(kind=dp),public :: hvphih_in,hvphih_out

#ifdef MPI       
    integer,private :: ierror
#endif

contains

!=======================================================================

  !> initializes constants and tables for radiation processes (heating, cooling and ionization)
  subroutine rad_ini ()

    ! initializes constants and tables for radiation processes
    ! (heating, cooling and ionization)

    use radiative_cooling, only: setup_cool

    ! Initialize integration routines
    call romberg_initialisation(NumFreq)

    ! Ask for the parameters of the spectrum
    call spectrum_parameters ()

    ! Determine spectrum diagnostics
    call spectrum_diagnostics ()

    ! Find the photo-ionization integrals for this spectrum
    call integrate_spectrum ()

    ! Set the radiative boundary conditions
    !call rad_boundary() ! NO LONGER NEEDED

    ! Set source position
    ! call source_position() CALLED ELSEWHERE

    ! Setup cooling
    if (.not.isothermal) call setup_cool () ! SHOULD BE CALLED ELSEWHERE

  end subroutine rad_ini

  !=======================================================================

  !> Input routine: establish the ionizing spectrum
  subroutine spectrum_parameters

    ! Input routine: establish the ionizing spectrum
     
    ! Author: Garrelt Mellema
    ! Update: 18-Feb-2004

    use file_admin, only: stdinput, file_input
    
    integer :: nchoice
    real(kind=dp) :: totflux

    ! Ask for input

    ! a) Effective temperature

    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
    ! Note that it is assumed that if teff_nominal is set, 
    ! S_star_nominal is ALSO set.
    if (rank == 0 .and. teff_nominal == 0.0) then
       if (.not.file_input) write(*,'(A)') ' '
       teff=0.0
       do while (teff < 2000.0 .or. teff > 200000.) 
          if (.not.file_input) write(*,'(A,$)') 'Give black body effective temperature: '
          read(stdinput,*) teff
          if (.not.file_input) write(*,*)
          if (teff < 2000.0 .or. teff > 200000.) then
             write(*,*) 'Error: Effective temperature out of range. Try again'
             write(*,*) 'Valid range: 2000 to 200,000'
          endif
       enddo

       ! Find total flux (Stefan-Boltzmann law)
       totflux=sigmasb*teff**4
       
       ! b) Luminosity, radius, or ionizing photon rate?
       if (.not.file_input) then
          write(*,'(A)') ' '
          write(*,'(A)') 'You can specify' 
          write(*,'(A)') ' 1) a stellar radius'
          write(*,'(A)') ' 2) a luminosity'
          write(*,'(A)') ' 3) Total number of ionizing photons'
       endif
       nchoice=0
       do while (nchoice <= 0 .or. nchoice > 3)
          if (.not.file_input) write(*,'(A,$)') 'Preferred option (1, 2 or 3): '
          read(stdinput,*) nchoice
          if (nchoice <= 0 .or. nchoice > 3) then
             write(*,*) 'Error: Choose between 1 2 or 3'
          endif
       enddo
       if (nchoice.eq.1) then
          if (.not.file_input) write(*,'(A,$)') 'Give radius in solar radii: '
          read(stdinput,*) rstar
          rstar=rstar*r_solar
          lstar=rstar*rstar*(4.0d0*pi*totflux)
          ! Number of photo-ionizing photons set to zero
          ! determined in spec_diag routine
          S_star=0.0
       elseif (nchoice .eq. 2) then
          if (.not.file_input) write(*,'(A,$)') 'Give luminosity in solar luminosities: '
          read(stdinput,*) lstar
          lstar=lstar*l_solar
          rstar=dsqrt(lstar/(4.0d0*pi*totflux))
          ! Number of photo-ionizing photons set to zero
          ! determined in spec_diag routine
          S_star=0.0
       else
          if (.not.file_input) write(*,'(A,$)') 'Give S_* (ionizing photons s^-1): '
          read(stdinput,*) S_star
          ! Assign some fiducial values, these are scaled to correspond 
          ! to S_star in routine spec_diag
          rstar=r_solar
          lstar=rstar*rstar*(4.0d0*pi*totflux)
       endif
    else
       ! teff and S_star are assumed to have been set in the c2ray_parameter 
       ! module
       teff=teff_nominal
       S_star=S_star_nominal
       totflux=sigmasb*teff**4
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       rstar=r_solar
       lstar=rstar*rstar*(4.0d0*pi*totflux)
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(teff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
    call MPI_BCAST(rstar,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
    call MPI_BCAST(lstar,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
#endif
    
  end subroutine spectrum_parameters

  !=======================================================================
  
  !> Calculates properties of the black body spectrum
  subroutine spectrum_diagnostics ()

    ! Calculates properties of spectrum
    ! This version: number of ionizing photons, S*, which can be
    ! used to calculate the Stromgren radius and other photon-statistics

    ! Author: Garrelt Mellema
    ! Update: 18-Feb-2004

    ! Tested against numbers listed on 
    ! http://nimbus.pa.uky.edu/plasma2000/input_for_nebular_models.htm
    ! (19 Feb 2004)

    integer :: i
    real(kind=dp) :: rfr,frmax,stepfl,flux
    real(kind=dp) :: fr(0:NumFreq),weight(0:NumFreq),bb(0:NumFreq)
    real(kind=dp) :: S_star_unscaled,scaling
    
    ! This is h/kT (unit 1/Hz, or sec)
    rfr=hplanck/(kb*teff)

    ! Upper limit of frequency integration
    frmax=min(frtop1,10.0*frtop2)

    ! Frequency step
    stepfl=(frmax-frth0)/real(NumFreq)

    ! Fill the arrays (frequency, weight, spectrum)
    do i=0,NumFreq
       fr(i)=frth0+stepfl*real(i)
       weight(i)=stepfl
       bb(i)=tpic2*fr(i)*fr(i)/(exp(fr(i)*rfr)-1.0)
    enddo

    ! Find ionizing flux by integrating over spectrum
    flux=scalar_romberg(bb,weight,NumFreq,NumFreq,0)

    ! Find out what is the S_star for the radius supplied.
    S_star_unscaled=4.0*pi*rstar*rstar*flux

    ! If S_star is zero, it is set here.
    if (S_star == 0.0) then
       S_star=S_star_unscaled
    else
       ! Find out the factor by which to change the radius
       ! and luminosity to get the required S_star.
       scaling=S_star/S_star_unscaled
       rstar=sqrt(scaling)*rstar
       lstar=scaling*lstar
    endif
    
    ! Report back
    if (rank == 0) then
       write(logf,'(/a)')           'Using a black body with'
       write(logf,'(a,1pe10.3,a)')   ' Teff=       ',teff,' K'
       write(logf,'(a,1pe10.3,a)')   ' Radius=     ',rstar/r_solar, &
            ' R_solar'
       write(logf,'(a,1pe10.3,a)')   ' Luminosity= ',lstar/l_solar, &
            ' L_solar'
       write(logf,'(A,1PE10.3,A//)') ' Number of H ionizing photons: ', &
            S_star,' s^-1'
    endif

  end subroutine spectrum_diagnostics
  
  !=======================================================================

  !> Calculates photo-ionization tables by integrating over frequency
  subroutine integrate_spectrum ()

    ! Calculates photo-ionization tables by integrating over frequency

    ! Author: Garrelt Mellema
    ! Date: 01-Nov-2012 (19-Feb-2004)
    ! Version: Single integration step

    ! Note 1: we calculate two integrals over each rate: one for optically
    ! thick cells (ensuring photon-conservation for those cells), and one 
    ! for optically thin cells. The latter are marked with 1.

    integer :: i,n
    real(kind=dp) :: frmax
    real(kind=dp) :: rfr
    real(kind=dp) :: fr(0:NumFreq)
    real(kind=dp) :: h0cross_freqdep(0:NumFreq)
    real(kind=dp),dimension(:),allocatable :: integrand
    
    ! Allocate the integrand needed for making the tables
    allocate(integrand(0:NumFreq))

    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do n=1,NumTau
       tau(n)=10.0**(minlogtau+dlogtau*real(n-1))
    enddo
    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! Warn about grey opacities:
    if (grey .and. rank == 0) write(logf,*) 'WARNING: Using grey opacities'

    ! Allocate photo-ionization tables
    allocate(hphot(0:NumTau,NumFreqBnd))
    allocate(hphot1(0:NumTau,NumFreqBnd))
    if (.not.isothermal) then
       allocate(hheat(0:NumTau,NumFreqBnd))
       allocate(hheat1(0:NumTau,NumFreqBnd))
    endif

    ! This is h/kT
    rfr=hplanck/(kb*teff)

    ! frequency band 1
    ! (there is space for NumFreqBnd frequency bands, only
    ! one is used here).
    if (frth0 < frtop1) then
       
       ! Upper limit of frequency integration
       frmax=min(frtop1,10.0*frtop2)

       ! Step size in frequency 
       steph0(1)=(frmax-frth0)/real(NumFreq)

       do i=0,NumFreq
          fr(i)=frth0+steph0(1)*real(i)
          
          ! Frequency dependence of the absorption
          ! cross section:
          if (grey) then
             h0cross_freqdep(i)=1.0
          else
             h0cross_freqdep(i)=(betah0*(fr(i)/frth0)**(-sh0)+ &
                  (1.0-betah0)*(fr(i)/frth0)**(-sh0-1.0))
          endif
       enddo

       ! Photo-ionization rate table (optically thick)
       do i=0,NumFreq
          integrand(i)=tpic2*fr(i)*fr(i)/(exp(fr(i)*rfr)-1.0)
       enddo
       hphot=make_table(fr,steph0(1),integrand,h0cross_freqdep,1)

       ! Photo-ionization rate table (optically thin)
       do i=0,NumFreq
          integrand(i)=tpic2*fr(i)*fr(i)*h0cross_freqdep(i)/(exp(fr(i)*rfr)-1.0)
       enddo
       hphot1=make_table(fr,steph0(1),integrand,h0cross_freqdep,1)

       ! Photo-ionization heating table (optically thick)
       if (.not.isothermal) then
          do i=0,NumFreq
             integrand(i)=hplanck*(fr(i)-frth0)*tpic2*fr(i)*fr(i)/ &
                  (exp(fr(i)*rfr)-1.0)
          enddo
          hheat=make_table(fr,steph0(1),integrand,h0cross_freqdep,1)

          ! Photo-ionization heating table (optically thin)
          do i=0,NumFreq
             integrand(i)=hplanck*(fr(i)-frth0)*tpic2*fr(i)*fr(i)*h0cross_freqdep(i)/ &
                  (exp(fr(i)*rfr)-1.0)
          enddo
          hheat1=make_table(fr,steph0(1),integrand,h0cross_freqdep,1)
       endif
    endif

  end subroutine integrate_spectrum
  
  ! =======================================================================

  !> Function to integrate the different rates.
  function make_table(fr,deltafr,integrand,absfr,IFreqBand)

    ! Version: Black body spectrum

    real(kind=dp),dimension(0:NumFreq),intent(in) :: fr !< frequency array
    real(kind=dp),intent(in) :: deltafr !< step size of frequency array
    !> core function to be integrated
    real(kind=dp),dimension(0:NumFreq),intent(in) :: integrand 
    !> frequency dependence of absorption coefficient
    real(kind=dp),dimension(0:NumFreq),intent(in) :: absfr 
    integer,intent(in) :: IFreqBand !< which frequency band
    !> result: photoionization table for different optical depths
    real(kind=dp),dimension(0:NumTau,NumFreqBnd) :: make_table
    
    real(kind=dp),dimension(0:NumFreq,0:NumTau) :: weight(0:NumFreq,0:NumTau)
    real(kind=dp),dimension(0:NumFreq,0:NumTau) :: func
    real(kind=dp),dimension(0:NumTau) :: integral_result
    integer :: i, n

    do i=0,NumFreq
       do n=0,NumTau
          weight(i,n)=deltafr
          ! Protect against floating point errors
          ! This needs to be checked. I remember that
          ! -700 is the minimum exponent allowed for
          ! doubleprecision...
          if (tau(n)*absfr(i) < 700.0) then
             func(i,n)=integrand(i)*exp(-tau(n)*absfr(i))
          else
             func(i,n)=0.0
          endif
       enddo
    enddo

    call vector_romberg (func,weight,NumFreq,NumFreq,NumTau,integral_result)
    do n=0,NumTau
       make_table(n,IFreqBand)=4.0*pi*rstar*rstar*integral_result(n)
    enddo

  end function make_table

  ! =======================================================================
  
  !> Calculates photo-ionization rates by looking them up in the tables
  subroutine photoion (phi,hcolum_in,hcolum_out,vol)
    
    ! Calculates photo-ionization rates
    
    ! Author: Garrelt Mellema
    ! Date: 28-Sep-2008 (11-May-2005 (f90) (18 feb 2004)
    
    ! Version:
    ! Simplified version derived from Coral version.
    ! Only hydrogen is dealt with, and one frequency band is used.

    !use sourceprops, only: NormFlux

    type(photrates),intent(out) :: phi !< result of the routine
    real(kind=dp),intent(in) :: hcolum_in !< H0 column density at front side
    real(kind=dp),intent(in) :: hcolum_out !< H0 column density at back side
    real(kind=dp),intent(in) :: vol !< volume of shell cell is part of
    !integer,intent(in) :: nsrc !< number of the source

    real(kind=dp) :: tauh_in,tauh_out
    real(kind=dp) ::  tau1,odpos1,dodpo1
    integer :: iodpo1,iodp11
    
    ! find the optical depths (in and outgoing)
    tauh_in=sigh*hcolum_in
    tauh_out=sigh*hcolum_out
    
    ! find the table positions for the optical depth (ingoing)
    tau1=log10(max(1.0e-20_dp,tauh_in))
    ! odpos1=min(1.0d0*NumTau,max(0.0d0,1.0d0+(tau1-minlogtau)/
    odpos1=min(real(NumTau,dp),max(0.0_dp,1.0+(tau1-minlogtau)/dlogtau))
    iodpo1=int(odpos1)
    dodpo1=odpos1-real(iodpo1,dp)
    iodp11=min(NumTau,iodpo1+1)
    
    ! Find the hydrogen photo-ionization rate (ingoing)
    ! Since all optical depths are hydrogen, we can use
    ! tau1 for all.
    phi%h_in=(hphot(iodpo1,1)+ &
         (hphot(iodp11,1)-hphot(iodpo1,1))*dodpo1)
    if (.not.isothermal) phi%hv_h_in= &
         (hheat(iodpo1,1)+(hheat(iodp11,1)-hheat(iodpo1,1))*dodpo1)

    ! Test for optically thick/thin case
    if (abs(tauh_out-tauh_in) > 1e-2) then 
       
       ! find the table positions for the optical depth (outgoing)
       tau1=log10(max(1.0e-20_dp,tauh_out))
       ! odpos1=min(1.0d0*NumTau,max(0.0d0,1.0d0+(tau1-minlogtau)/
       odpos1=min(real(NumTau,dp),max(0.0_dp,1.0+(tau1-minlogtau)/dlogtau))
       iodpo1=int(odpos1)
       dodpo1=odpos1-real(iodpo1)
       iodp11=min(NumTau,iodpo1+1)
        
       ! find the hydrogen photo-ionization rate (outgoing)
       phi%h_out=(hphot(iodpo1,1)+ &
            (hphot(iodp11,1)-hphot(iodpo1,1))*dodpo1)
       if (.not.isothermal) phi%hv_h_out= &
            (hheat(iodpo1,1)+(hheat(iodp11,1)-hheat(iodpo1,1))*dodpo1)
       
       ! The photon conserving photo-ionization rate is the difference between
       ! the one coming in, and the one going out.
       phi%h=(phi%h_in-phi%h_out)/vol
       if (.not.isothermal) phi%hv_h=(phi%hv_h_in-phi%hv_h_out)/vol
       
    else ! optically thin case
       
       ! Find the hydrogen photo-ionization rate for the optically thin
       ! case, and from this derive the outgoing rate.
       ! Since all optical depths are hydrogen, we can use
       ! tau1 for all.
       phi%h=(tauh_out-tauh_in)*( &
            hphot1(iodpo1,1)+(hphot1(iodp11,1)-hphot1(iodpo1,1))*dodpo1)/vol
       phi%h_out=phi%h_in-phi%h*vol

       if (.not.isothermal) then
          phi%hv_h=(tauh_out-tauh_in)*( &
               hheat1(iodpo1,1)+(hheat1(iodp11,1)-hheat1(iodpo1,1))*dodpo1)/vol
          phi%hv_h_out=phi%hv_h_in-phi%hv_h*vol
       endif

    endif

  end subroutine photoion

end module radiation

