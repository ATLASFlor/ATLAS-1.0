subroutine openlog
  !**************************************************************
  !*
  !*    Opens the control file
  !*
  !**************************************************************
  use InpOut
  use Master
  implicit none
  !
  character(len=10) :: rdate
  character(len=8 ) :: rtime
  !
  !***  Clock
  !
  call datem(rdate,rtime)
  !
  !***  Opens and writes the log file
  !
  open(lulog,file=TRIM(flog),status='unknown')
  write(lulog,1)            version,rdate,rtime
  if(out_screen) write(*,1) version,rdate,rtime
  !
1 format(/, &
       '--------------------------------------------------',/,  &
       '                  Program ATLAS                   ',/,  &
       '--------------------------------------------------',/,  &
       '  Version               : ',f5.2,/, &
       '  Starting date         : ',a10,/, &
       '  Starting time         : ','  ',a8,/)
  !             
  return
  !
end subroutine openlog
!

!
subroutine datem(rdate,rtime)
  !***********************************************************************
  !*
  !*    Get system date and time from system clock
  !*
  !*    rtime*8    HH:MM:SS
  !*    rdate*10   MM-DD-YYYY
  !*
  !***********************************************************************
  implicit none
  !
  character(len=8)  :: rtime,sdate
  character(len=10) :: rdate,stime
  !
  call DATE_AND_TIME(sdate,stime)
  ! 
  rdate='  -  -    '
  rdate(1:2 )=sdate(5:6)
  rdate(4:5 )=sdate(7:8)
  rdate(7:10)=sdate(1:4)
  rtime='  :  :  '
  rtime(1:2 )=stime(1:2)
  rtime(4:5 )=stime(3:4)
  rtime(7:8 )=stime(5:6)
  !
  return
end subroutine datem
