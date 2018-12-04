subroutine readmeteo(nm,itime)
  !******************************************************************************************
  !*    READMETEO
  !*    Open the corresponding Meteorological file and get variables values.
  !*    nm: correspond to the number nm meteorological file introduced in the .inp
  !*    itime: Correspond to the time step to be read.
  !*
  !******************************************************************************************
  use InpOut
  use Master
  use KindType
  use Wrf 
  use GFS 
  use Debug                                        
  implicit none
  !
  integer (ip)          :: nm, itime
  character(len=s_mess) :: str
  character(len=1     ) :: ext   
  character(len=19    ) :: string_time   
  !
  write(ext,'(i1)') nm
  str = 'METEO_DATA_'//ext
  !
  !*** Open and read meteo model and time-dependent variables
  !
  if(TRIM(meteo(nm)%modeltype)=='WRF') then
     ! 
     call readwrf(nm,itime)
     !
  else if (TRIM(meteo(nm)%modeltype)=='GFS') then
     !   
     call readgfs(nm,itime)
     !
  else if (TRIM(meteo(nm)%modeltype)=='DEBUG')then
     !
     call readDEBUG(nm,itime)
     ! 
     !
  else
     call runend('The model Type is not recognize')
  end if
  !
  !*** Write in log file 
  !
  write(string_time,'(F13.0)')meteo(nm)%time(itime)
  string_time=string_time(1:4)//'/'//string_time(5:6)//'/'//string_time(7:8)//'-'//string_time(9:10)&
       //':'//string_time(11:12)//':00'
  write(lulog,5) TRIM(meteo(nm)%modeltype), & 
       TRIM(str), &
       string_time!meteo(nm)%time(itime) f16.0
5 format(/, &
       'Read ',a,' data file : ',a,' at time : ',a,/)
  if(out_screen) write(*,1) TRIM(meteo(nm)%modeltype), & 
       TRIM(str), &
       string_time!meteo(nm)%time(itime) f16.0
1 format(/, &
       'Read ',a,' data file : ',a,' at time : ',a,/)
  !      
  return
end subroutine readmeteo
