!***************************************************************
!*
!*		Module for master operations
!*
!***************************************************************
MODULE Master
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Version of input file
  !
  character(len=3) :: version
  !
  !*** Grib mesoscal/global model
  !
  character(len=25) :: date
  character(len=25) :: model = ' '
  integer(ip) :: ibyr,ibmo,ibdy,ibhr,dt
  real   (rp) :: cen_lon,cen_lat,missing_value
  !
  !*** Mesh
  !
  logical     :: invert_x,invert_y
  integer(ip) :: nx,ny,np,nt,n2d,n3d
  integer(ip) ::          it,i2d,i3d,itt,ipos
  real   (rp) :: lonmin,lonmax,latmin,latmax
  !
  !*** netcdf
  !
  integer(ip)  :: ncID
  integer(ip)  :: nx_nc_ID,ny_nc_ID,np_nc_ID,nt_nc_ID
  integer(ip)  :: var2d_ID(100),var3d_ID(100)
  !
  character(len=20) :: var2d(100),var3d(100)
  !
  !*** Attributes
  !
  character(len=50 ) :: att2d_u(100),att3d_u(100)
  character(len=100) :: att2d_d(100),att3d_d(100)
  !
  !*** nc dimensions (=coordinate variables) names and ID
  !
  integer(ip)  :: lon_nc_ID
  integer(ip)  :: lat_nc_ID
  integer(ip)  :: pre_nc_ID
  integer(ip)  :: time_nc_ID
  !
  character(len=25) :: lon_nc_name   = 'lon'
  character(len=25) :: lat_nc_name   = 'lat'
  character(len=25) :: pre_nc_name   = 'pres'
  character(len=25) :: time_nc_name  = 'time'
  !
  !*** Variables
  !
  real(rp), allocatable :: pres(:)
  real(rp), allocatable :: timesec(:)
  real(rp), allocatable :: lon(:)
  real(rp), allocatable :: lat(:)
  real(rp), allocatable :: work(:,:,:)
  real(rp), allocatable :: work1d(:)
  real(rp), allocatable :: work2d(:,:)
  real(4 ), allocatable :: work4(:,:)
  !
END MODULE Master
