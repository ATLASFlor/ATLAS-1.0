subroutine runend(message)
  !*************************************************************************
  !*
  !*    This routine stops the run and writes termination status
  !*
  !*************************************************************************
  implicit none
  !
  integer(4)      :: iwarn, nwarn
  character(len=*) :: message
  character(len=50):: warning(20) 
  !
  character(len=10) :: rdate
  character(len=8 ) :: rtime
  !
  !***  Warning list
  !
if(nwarn.gt.0)then
 write(*,1) nwarn
1 format(/,'---> Number of warnings :',i2)
  do iwarn = 1,nwarn
 write(*,2) iwarn,warning(iwarn)
  end do
2 format('     ',i2,' : ',a)
end if
  !
  !*** Write message and stop the run.
  !
  if(TRIM(message).ne.'OK') then             ! Abnormal termination
     !  
write(*,10) TRIM(message)
10   format('---> Number of errors   : 1',/, &
          '---> ',a)
     !
 write(*,11) 
! write(*,12) 
11   format(/,'---> Program erainNC ends ABNORMALLY') 
!12   format(  '---> Please check the log file for details',/)
     stop 1                         ! environment 1
     !
  else                                     ! Normal termination
     !   
write(*,20) 
!     if(out_screen) write(*,12) 
20   format('---> Number of errors   : 0',/, &
          '---> Program erainNC ends NORMALLY')
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
  implicit none
  character(len=*) :: message
  integer(4)       :: nwarn =0
  integer(4)       :: maxwarn=20
  character(len=50):: warning(20) 
  !
  nwarn = nwarn + 1
  nwarn = MIN(nwarn,maxwarn)
  warning(nwarn) = message(1:LEN_TRIM(message))
  !
  return      
end subroutine wriwar
