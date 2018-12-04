   subroutine inival
  !***********************************************************************
  !*
  !*     Initializes the mesoscale grid values. Values read from the
  !*     model-dependent configuration file
  !*
  !***********************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  integer(ip) :: stoi1,info,i,ix,iy
  real(rp)    :: dlon,dlat
  character(len=s_mess) :: cvoid
  logical     :: found
  !
  !***   Read dimensions
  !
  open(lucnf,FILE=TRIM(lucnfname),status='old')
  read(lucnf,*,iostat=info) cvoid
  read(lucnf,*,iostat=info) cvoid
  read(lucnf,*,iostat=info) cvoid
  if( info /= 0 ) call runend('Config file not found or incorrect header')
  read(lucnf,*) cvoid,lonmin
  read(lucnf,*) cvoid,lonmax
  read(lucnf,*) cvoid,latmin
  read(lucnf,*) cvoid,latmax
  read(lucnf,*) cvoid,nx
  read(lucnf,*) cvoid,ny
  read(lucnf,*) cvoid,np
  read(lucnf,*) cvoid,nt
  read(lucnf,*) cvoid,dt
  read(lucnf,*) cvoid,n2d
  read(lucnf,*) cvoid,n3d
  read(lucnf,*) cvoid,cen_lon
  read(lucnf,*) cvoid,cen_lat
  read(lucnf,*) cvoid,missing_value
  !
  !***  Allocates memory
  !
  allocate(pres(np))
  allocate(timesec(0:nt))
  allocate(lon(nx))
  allocate(lat(ny))
  allocate(work  (nx,ny,np))
  allocate(work2d(nx,ny   ))
  allocate(work4 (nx,ny   ))
  !
  !***  Reads the mesh and sets pressure levels
  !
  read(lucnf,*) cvoid
  do i=1,np
     read(lucnf,*) pres(i)
  end do
  !
  !***  Reads variables
  !
  read(lucnf,*) cvoid
  do i=1,n2d
     read(lucnf,*) cvoid,var2d(i)
     read(lucnf,*) cvoid,att2d_d(i)
     read(lucnf,*) cvoid,att2d_u(i)
  end do
  !
  read(lucnf,*) cvoid
  do i=1,n3d
     read(lucnf,*) cvoid,var3d(i)
     read(lucnf,*) cvoid,att3d_d(i)
     read(lucnf,*) cvoid,att3d_u(i)
  end do
  close(lucnf)
  !
  !***  Computes ibyr,ibmo,ibdy,ibhr.
  !***  date is in the form YYMMDDHH
  !
  ibyr = stoi1(date(1:1))*10 + stoi1(date(2:2))
  if(ibyr.lt.40) then
     ibyr = 2000 + ibyr
  else
     ibyr = 1900 + ibyr
  end if
  ibmo = stoi1(date(3:3))*10 + stoi1(date(4:4))
  ibdy = stoi1(date(5:5))*10 + stoi1(date(6:6))
  ibhr = stoi1(date(7:7))*10 + stoi1(date(8:8))
  !
  !***  Computes timesec.  Seconds after 0000UTC for the day
  !***  ibyr,ibmo,ibdy
  !
  do it=0,nt
     timesec(it) = (ibhr + it*dt)*3600
  end do
  !
  !***  Computes lon/lat
  !
  dlon = (lonmax-lonmin)/(nx-1)
  do ix = 1,nx
     lon(ix) = lonmin + (ix-1)*dlon
  end do
  dlat = (latmax-latmin)/(ny-1)
  do iy = 1,ny
     lat(iy) = latmin + (iy-1)*dlat
  end do
  !
  !***  Model specific conversion
  !
  invert_x = .false.
  invert_y = .false.
  if(TRIM(model).eq.'ecmwf'.or.TRIM(model).eq.'ECMWF') then
     invert_y = .true.
  end if
  if(TRIM(model).eq.'eraIn'.or.TRIM(model).eq.'ERAIN') then
     invert_y = .true.
  end if
  if(TRIM(model).eq.'gfs05deg'.or.TRIM(model).eq.'GFS05DEG') then
     invert_y = .true.
  end if
 if(TRIM(model).eq.'gfs1deg'.or.TRIM(model).eq.'GFS1DEG') then
     invert_y = .true.
  end if
  if(TRIM(model).eq.'ncepFNL'.or.TRIM(model).eq.'NCEPFNL') then
     invert_y = .true.
  end if
  if(TRIM(model).eq.'ncep1'.or.TRIM(model).eq.'NCEP1') then
     invert_y = .true.
  end if
  !
  !*** Checks if longitude has to change from (0,360) to (-180,180)
  !
  found = .false.
  do ix = 1,nx
     if( (lon(ix).ge.180.0).and.(.not.found) ) then
        ipos = ix
        found = .true.
     end if
  end do
  !
  if(found) invert_x = .true. 
  !
  if(invert_x) then
     allocate(work1d(nx))
     work1d(1:nx) = lon(1:nx)       
     !
     i = 0
     do ix = ipos,nx
        i = i + 1
        lon(i) = work1d(ix) - 360.0
     end do
     do ix = 1,ipos-1
        i = i + 1
        lon(i) = work1d(ix)
     end do
     !
     lonmin = lon(1)
     lonmax = lon(nx)
     !
   end if
   !
   return
   end subroutine inival
!
!
!
   integer function stoi1(string1)
  !**************************************************************
  !*
  !*    Decodes a character*1 string
  !*
  !**************************************************************
  implicit none
  character(len=1) :: string1
  !
  if(string1.eq.'0') then
     stoi1 = 0
  else if(string1.eq.'1') then
     stoi1 = 1
  else if(string1.eq.'2') then
     stoi1 = 2
  else if(string1.eq.'3') then
     stoi1 = 3
  else if(string1.eq.'4') then
     stoi1 = 4
  else if(string1.eq.'5') then
     stoi1 = 5
  else if(string1.eq.'6') then
     stoi1 = 6
  else if(string1.eq.'7') then
     stoi1 = 7
  else if(string1.eq.'8') then
     stoi1 = 8
  else if(string1.eq.'9') then
     stoi1 = 9
  end if
  return
end function stoi1
