!***************************************************************
!*
!*		Module for input/output
!*
!***************************************************************
MODULE InpOut
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Type of wgrib files
  !
  logical :: wgrib_bin = .true.  ! set to .true./.false. for wgrib -bin or -text
  !
  !***  File logical units
  !
  integer(ip), parameter  :: lulog = 11    ! log file
  integer(ip), parameter  :: lucnf = 12    ! netCDF file
  integer(ip), parameter  :: lugrb = 13    ! data file (grib decoded dataset)
  integer(ip), parameter  :: lunc  = 14    ! netCDF file
  !
  !***  File names
  !
  character(len=s_file) :: lulogname,lucnfname,lugrbname,luncname,lugrbbase
  integer(ip)           :: iarg
  !
  !***  List of Warnings
  !
  integer(ip), parameter :: maxwarn = 100
  integer(ip)            :: nwarn = 0
  character(len=s_mess)  :: warning(maxwarn)
  !
END MODULE InpOut
