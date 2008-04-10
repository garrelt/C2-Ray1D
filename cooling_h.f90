module radiative_cooling

  ! This module file contains subroutines for radiative cooling
  ! - coolin     - calculate cooling rate
  ! - setup_cool - setup cooling table

  ! Version: H-only cooling

  use precision, only: dp
  
  implicit none

  integer,parameter,private :: temppoints=81
  real(kind=dp),dimension(temppoints),private :: h0_cool,h1_cool
  real(kind=dp),private :: mintemp, dtemp

contains
  
  !===========================================================================

  function coolin(nucldens,eldens,xh,temp0)
    
    real(kind=dp) :: coolin

    real(kind=dp),intent(in) :: nucldens,eldens,xh(0:1),temp0

    real(kind=dp) :: tpos, dtpos
    integer :: itpos,itpos1

    tpos=(log10(temp0)-mintemp)/dtemp+1.0d0
    itpos=min(temppoints-1,max(1,int(tpos)))
    dtpos=tpos-real(itpos)
    itpos1=min(temppoints,itpos+1)

    ! Cooling curve
    coolin=nucldens*eldens*( &
         xh(0)*(h0_cool(itpos)+ &
         (h0_cool(itpos1)-h0_cool(itpos))*dtpos)+ &
         xh(1)*(h1_cool(itpos)+ &
         (h1_cool(itpos1)-h1_cool(itpos))*dtpos))
    
    return
  end function coolin

  !===========================================================================

  subroutine setup_cool ()

    real(kind=dp),dimension(temppoints) :: temp
    integer :: itemp
    integer :: element,ion,nchck

    ! Open cooling table (H0)
    open(unit=22,file='tables/H0-cool-B.tab',status='old')
    ! Read the cooling data
    read(22,*) element,ion,nchck
    if (nchck.eq.0) then
       do itemp=1,temppoints
          read(22,*) temp(itemp),h0_cool(itemp)
       enddo
    else
       write(*,*) 'Error reading cooling tables'
    endif
    close(22)
    
    mintemp=temp(1)
    dtemp=temp(2)-temp(1)
    ! not needed: maxtemp=temp(temppoints)
    
    ! Open cooling table (H1)
    open(unit=22,file='tables/H1-cool-B.tab',status='old')
    ! Read the cooling data
    read(22,*) element,ion,nchck
    if (nchck.eq.0) then
       do itemp=1,temppoints
          read(22,*) temp(itemp),h1_cool(itemp)
       enddo
    else
       write(*,*) 'Error reading cooling tables'
    endif
    close(22)
    
    ! Convert cooling to from log to linear 
    do itemp=1,temppoints
       h0_cool(itemp)=10.0d0**h0_cool(itemp)
       h1_cool(itemp)=10.0d0**h1_cool(itemp)
    enddo

    return
  end subroutine setup_cool
  
end module radiative_cooling
