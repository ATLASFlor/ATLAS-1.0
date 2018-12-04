subroutine outrest
  !******************************************************************************************
  !*
  !* write in a text output file the final particles location with their characteritics
  !* 
  !*
  !******************************************************************************************

  use InpOut
  use Master
  use KindType                                     
  implicit none
  integer (ip) :: iphase, ipart
  integer (ip) :: partcount ! To enumerate particles in the file
  real    (rp) :: syear,smonth,sday,shour,sminute
  !
  !*** Opens the file
  !
  open(56,FILE=TRIM(frest),status='REPLACE')
  !
  !*** Calculates the time end: year, month, day, hour, minutes and write it in the file
  !
  syear   = INT(time_end/100000000)                                       !Restart dates
  smonth  = INT((time_end-syear*100000000)/1000000)
!  sday    = INT((time_end-syear*100000000-smonth*1000000)/10000)
!  shour   = INT((time_end-syear*100000000-smonth*1000000-sday*10000)/100)
!  sminute = INT((time_end-syear*100000000-smonth*1000000-sday*10000-shour*100))
  sday    = 21.
  shour   = 00.
  sminute = 00.

  write(56,564)'END_YEAR', syear, 'END_MONTH',smonth, 'END_DAY',sday, &
       'END_HOUR',shour,'END_MINUTE',sminute
  write(56,565)'particle','state','bin','rho','diam','mass','sphe','psi','lon','lat','z'
564 FORMAT(a,1X,f5.0,/,a,1X,f3.0,/,a,1X,f3.0,/,a,1X,f3.0,/,a,1X,f3.0)
565 FORMAT(a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a)
  !
  !*** Loop over particles and writes information on the file
  !
  partcount=0
  do iphase=1,nphases
     do ipart=phase(iphase)%firstpart,phase(iphase)%lastpart
        if ((part(ipart)%state==1).or.(part(ipart)%state==-1)) then
           partcount=partcount+1
           write(56,570)partcount,part(ipart)%state,part(ipart)%ibin,part(ipart)%rho,&
                part(ipart)%diam,part(ipart)%mass,part(ipart)%sphe,part(ipart)%psi,&
                part(ipart)%lon,part(ipart)%lat,part(ipart)%z
570        FORMAT(i9,1X,i2,1X,i3,1X,f11.6,1X,f11.6,1X,f24.6,1X,f11.6,1X,f11.6,1X,f11.6,1X,f11.6,1X,f15.6)
        end if
     end do
  end do
  write(56,571)'TOTAL_PARTICLES',partcount
571 FORMAT(a,1X,i9)
  !
  !*** Close the file
  !
  close(56)
  return
end subroutine outrest
