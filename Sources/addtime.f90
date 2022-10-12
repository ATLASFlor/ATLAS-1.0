subroutine addtime(iyr0,imo0,idy0,ihr0,iyr,imo,idy,ihr,imi,ise,time)
  !*********************************************************************
  !*
  !*    Adds time seconds to the initial date YYYY-MM-DD-HH
  !*
  !*    INPUTS:
  !*      IYR0    - integer - Current year
  !*      IMO0    - integer - Current month (1-12)
  !*      IDY0    - integer - Current day
  !*      IHR0    - integer - Current hour (0-23)
  !*      TIME    - real*8  - time increment in seconds
  !*
  !**********************************************************************
  use KindType
  implicit none
  integer(ip) :: iyr0,imo0,idy0,ihr0,iyr,imo,idy,ihr,imi,ise
  integer(ip) :: nhincr,ijl
  real   (rp) :: time,work
  !
  !***  Initialization
  !
  iyr = iyr0
  imo = imo0
  idy = idy0
  ihr = ihr0
  !
  !***  Time < 1 min
  !
  if(time < 60.0_rp) then
     imi = 0
     ise = int(time)
     return
     !
     !***  Time < 1 h
     !
  else if(time < 3600.0_rp) then
     imi = int(time/60.0_rp)
     ise = int(time)-60*imi
     return
  else
     !
     !***  Time > 1 h
     !
     nhincr = int(time/3600.0_rp)
     work   = time - nhincr*3600.0_rp
     imi    = int(work/60.0_rp)
     ise    = int(work)-60*imi

     call julday  (iyr,imo,idy,ijl)      ! Computes the Julian day ijl
     call timeincr(iyr,ijl,ihr,nhincr)   ! Updates the  Julian day
     call grday   (iyr,ijl,imo,idy)      ! Converts back to Gregorian
     return
  end if
end subroutine addtime


subroutine julday(iyr,imo,iday,ijuldy)
  !*********************************************************************
  !*
  !*    Computes the Julian day number from the Gregorian date
  !*    (month, day)
  !*
  !*********************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: iyr,imo,iday,ijuldy,ierr
  integer(ip), dimension(12) :: kday =(/0,31,59,90,120,151,181,212,243,273,304,334/)
  !
  ierr=0
  !
  !***  Check for valid month
  !
  if(imo < 1.or.imo > 12) ierr=1
  !
  !***  Check for valid day in 30-day months
  !
  if(imo == 4.or.imo == 6.or.imo == 9.or.imo == 11)then
     if(iday > 30)ierr=1
  else if(imo == 2)then
     if(mod(iyr,4) == 0)then  ! February in a leap year
        if(iday > 29)ierr=1
     else                     ! February in a non-leap year
        if(iday > 28)ierr=1
     endif
  else
     !
     !***  Check for valid day in 31-day months
     !
     if(iday > 31) ierr=1
  endif
  !
  if(ierr == 1)call runend('Wrong computation of the julian day')
  !
  !***  Compute the Julian day
  !
  ijuldy=kday(imo)+iday
  if(imo <= 2) return
  if(mod(iyr,4) == 0) ijuldy=ijuldy+1
  return
end subroutine julday


subroutine grday(iyr,ijul,imo,iday)
  !**********************************************************************
  !*
  !*    Compute the Gregorian date (month, day) from the Julian day
  !*
  !**********************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: iyr,ijul,imo,iday,ileap,i
  integer(ip), dimension(12,2) :: kday =reshape((/31,59,90,120,151,181,212,243,273,304,334,365,  &
       31,60,91,121,152,182,213,244,274,305,335,366/),shape(kday))
  !
  ileap=1
  if(mod(iyr,4) == 0) ileap=2
  if(ijul < 1.or.ijul > kday(12,ileap))go to 11
  !
  do i=1,12
     if(ijul > kday(i,ileap)) cycle
     imo=i
     iday=ijul
     if(imo /= 1)iday=ijul-kday(imo-1,ileap)
     return
  end do

11 call runend('Wrong computation of the gegorian day')
end subroutine grday


subroutine timeincr(iyr,ijul,ihr,nhrinc)
  !**********************************************************************
  !*
  !*    Increments the time and date by "NHRINC" hours
  !*
  !*
  !*    INPUTS:
  !*       IYR    - integer - Current year
  !*       IJUL   - integer - Current Julian day
  !*       IHR    - integer - Current hour (00-23)
  !*       NHRINC - integer - Time increment (hours)
  !*
  !*       NOTE: "NHRINC" must >= -24
  !*             Hour is between 00-23
  !*
  !*    OUTPUT:
  !*       IYR    - integer - Updated year
  !*       IJUL   - integer - Updated Julian day
  !*       IHR    - integer - Updated hour (00-23)
  !*
  !**********************************************************************
  use KindType
  implicit none
  !
  integer(ip) :: iyr,ijul,ihr,nhrinc
  integer(ip) :: nleft,ninc,ileap
  !
  !***  Check nhrinc
  !
  if(nhrinc < 0) call runend('TIMEINCR invalid time increment')
  !
  !***  Save increment remaining (needed if nhrinc > 8760)
  !
  nleft=nhrinc
  !
  !***  Process change in hour
  !
  if(nhrinc > 0)then
     !
10   ninc=MIN0(nleft,8760)
     nleft=nleft-ninc
     !
     !***  Increment time
     !
     ihr=ihr+ninc
     if(ihr <= 23) return
     !
     !***  Increment day
     !
     ijul=ijul+ihr/24
     ihr=mod(ihr,24)
     !
     !***  ILEAP = 0 (non-leap year) or 1 (leap year)
     !
     if(mod(iyr,4) == 0)then
        ileap=1
     else
        ileap=0
     endif
     !
     !***  Update year
     !
     if(ijul > 365+ileap) then
        iyr=iyr+1
        ijul=ijul-(365+ileap)
     endif
     !
     !***  Repeat if more hours need to be added
     !
     if(nleft > 0) goto 10
     !
  elseif(nhrinc < 0)then
     !***     Decrement time
     ihr=ihr+nhrinc
     if(ihr < 0)then
        ihr=ihr+24
        ijul=ijul-1
        if(ijul < 1)then
           iyr=iyr-1
           if(mod(iyr,4) == 0) then
              ijul=366
           else
              ijul=365
           endif
        endif
     endif
  endif
  !
  return
end subroutine timeincr
