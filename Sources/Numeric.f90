!***************************************************************
!*
!*		Module for numeric operations
!*
!***************************************************************
MODULE Numeric
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !***  Constants
  !
  real(rp) :: pi     = 3.14159265358979323846_rp   ! pi
  real(rp) :: Rearth = 6356d3                      ! Earth's radius
  ! 
  !
  !***  System of coordinates and transformations
  !
  character(len=10)  :: coord_sys
  real(rp)           :: dX1,dX2
  !
  !
  !***  Mesh related variables
  !
  integer(ip) :: nx,ny,nz
  real   (rp) :: lonmin,lonmax,latmin,latmax
  real   (rp) :: xmin  ,xmax  ,ymin  ,ymax
  real   (rp) :: xorigr,yorigr,dx,dy,tpgmax,tpgmin,ztop
  !
END MODULE Numeric
