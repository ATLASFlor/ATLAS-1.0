subroutine runend(message)
  !*************************************************************************
  !*
  !*    This routine stops the run and writes termination status
  !*
  !*************************************************************************
  use kindType 
  use InpOut
  implicit none
  !
  integer(ip)      :: iwarn
  character(len=*) :: message
  !
  character(len=10) :: rdate
  character(len=8 ) :: rtime
  !
  !***  Warning list
  !
  write(lulog,1) nwarn
  if(out_screen) write(*,1) nwarn
1 format(/,'---> Number of warnings :',i2)
  do iwarn = 1,nwarn
     write(lulog,2) iwarn,warning(iwarn)
     if(out_screen) write(*,2) iwarn,warning(iwarn)
  end do
2 format('     ',i2,' : ',a)
  !
  !*** Write message and stop the run.
  !
  if(TRIM(message).ne.'OK') then             ! Abnormal termination
     !  
     write(lulog,10) TRIM(message)
     if(out_screen) write(*,10) TRIM(message)
10   format('---> Number of errors   : 1',/, &
          '---> ',a)
     !
     write(lulog,11)
     write(lulog,12)
     if(out_screen) write(*,11) 
     if(out_screen) write(*,12) 
11   format(/,'---> Program ATLAS ends ABNORMALLY') 
12   format(  '---> Please check the log file for details',/)
     stop 1                         ! environment 1
     !
  else                                     ! Normal termination
     !   
     call datem(rdate,rtime)
     write(lulog,20)
     write(lulog,12)
     write(lulog,21)rdate, rtime
     if(out_screen) write(*,20) 
     if(out_screen) write(*,12) 
     if(out_screen) write(*,21) rdate, rtime
20   format('---> Number of errors   : 0',/, &
          '---> Program ATLAS ends NORMALLY')
21   format('Ending Date        :',a10,/, &
          'Ending Time        :' ,a8,/)
     stop 0                         ! environment 0
     !
  end if
  !
end subroutine runend
!
!
!
subroutine wriwar(message)
  !*************************************************************************
  !*
  !*    This routine writes a warning message to the warnings list
  !*
  !*************************************************************************
  use InpOut
  implicit none
  character(len=*) :: message
  !
  nwarn = nwarn + 1
  nwarn = MIN(nwarn,maxwarn)
  warning(nwarn) = message(1:LEN_TRIM(message))
  !
  return      
end subroutine wriwar
