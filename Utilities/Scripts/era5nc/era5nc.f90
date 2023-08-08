program erainNC
!###############################################################################################
!# This Program converts two netcdf files downloaded by ECMWF Era5 meteorological file
!# to one netcdf file capable to be readed by ATLAS atmospheric dispersion model.
!#
!# The netcdf files are referred to Invariant variables (Geopotential, Land-sea mask)
!# Surface variables (Boundary layer height, 10 metre U wind component, 10 metre V wind component, 2 metre 
!# temperature (analysis: 00:00:00, 06:00:00, 12:00:00, 18:00:00). step=all)
!# Pressure variables (Geopotential, Relative humidity, Temperature, U component of wind, V component of wind, Vertical velocity)
!###
!# Author: Florencia Reckziegel
!# 
!###############################################################################################
!
!***Defines variables
!
use netcdf
implicit none
!
real(4), parameter :: g      = 9.81          ! gravity
character(len=80):: file_invariant, file_surface, file_surface2, file_pressure
character(len=50) :: attr_desc, attr_units
character(len=10) :: attr_title
! 
real(4)    :: lonmax, lonmin, latmax, latmin, res
integer(4) :: year, month, day, endday,hour
integer(4) :: i,ip
integer(4) :: cen_lon,cen_lat,time_incr
real(8)    :: missing_value
real   (4) :: add_offset, scale_factor  ! attributes for variables to read
! NETCDF parameters
integer(4) :: ncID                                      ! netcdf file id
!ID
integer(4) :: nxdimid,nydimid,npdimid,ntdimid           ! dimensions pressure
integer(4) :: lonvarid,latvarid,presvarid,tid           ! variables pressure
integer(4) :: tvarid,geopid,relhumid,tempid,uid,vid,wid ! variables pressure
integer(4) :: lsmid, sgeopid                            ! variables invariant
integer(4) :: tpid,blhid                                ! variables surface
integer(4) :: u10id,v10id,t2id                          ! variables surface2
! Pressure variables and dimensions
integer(4) :: nx,ny,np                              !Dimensions of the new netcdf file
! Time
integer(4) :: ntm,nt!,ntw                            ! time steps of the total file(ntm), and time of user's interest (nt). 
!integer(4) :: nt1,nt2                               ! initial and final time indices of user's interest inside the total time steps of ecmwf files.
real   (4),allocatable :: lon (:)                   ! Lon vector corresponding to the entire ecmwf downloaded files.
real   (4),allocatable :: lat (:)                   ! Lat vector corresponding to the entire ecmwf downloaded files.
real   (4),allocatable :: pres(:)                   ! Pres vector corresponding to the entire ecmwf downloaded files.
!real   (4),allocatable :: time(:)                   ! tiempo TEMPORAL time(ntm)
integer(4),allocatable :: timeecmwf(:)              ! Time vector corresponding to the entire ecmwf downloaded files.
integer(4),allocatable :: timeu(:)                  ! time vector with the user's interest times
real   (4),allocatable :: geopotential(:,:,:,:)     ! geopotential(lon,lat,pressure,time)
real   (4),allocatable :: relhum(:,:,:,:)           ! relative humidity(lon,lat,pressure,time)
real   (4),allocatable :: temperature(:,:,:,:)      ! temperature(lon,lat,pressure,time)
real   (4),allocatable :: uwind(:,:,:,:)            ! horizontal u wind(lon,lat,pressure,time)
real   (4),allocatable :: vwind(:,:,:,:)            ! horizonat crosswind v (lon,lat,pressure,time)
real   (4),allocatable :: wwind(:,:,:,:)            ! vertical wind w (lon,lat,pressure,time)
! Invariant variables
real   (4),allocatable :: surfgeop(:,:,:)           ! surface geopotential (lon,lat,time)
real   (4),allocatable :: lsm(:,:,:)                ! land sea mask (lon,lat,time)
! Surface variables
real   (4),allocatable :: work(:,:,:)               ! Auxiliar variable(lon,lat,time)
real   (4),allocatable :: work1(:,:,:)              ! Auxiliar variable(lon,lat,1)
real   (4),allocatable :: work2(:,:,:,:)            ! Auxiliar variable(lon,lat,pres,time)
real   (4),allocatable :: p    (:,:,:,:)            ! Auxiliar variable pressure(lon,lat,pres,time)
real   (4),allocatable :: tp(:,:,:)                 ! Total precipitation (lon,lat,time)
real   (4),allocatable :: blh(:,:,:)                ! boundary layer height (lon,lat,time)
! Surface2 variables
real   (4),allocatable :: u10(:,:,:)                ! 10 metres u wind (lon,lat, time)
real   (4),allocatable :: v10(:,:,:)                ! 10 metres v wind (lon,lat,time)
real   (4),allocatable :: t2(:,:,:)                 ! 2 metres temperature (lon,lat,time)
! Auxiliar
real   (4),allocatable :: tv(:,:,:,:)                 ! virtual temperature (lon,lat,time)
real   (4),allocatable :: tpo(:,:,:,:)                 ! potential temperature (lon,lat,time)
real   (4),allocatable :: ro(:,:,:,:)                 ! density (lon,lat,time)
! Output netcdf 
integer(4) :: nlon_id,nlat_id,npres_id,nt_id ! dimensions
integer(4) :: lon_id,lat_id,pres_id,time_id ! variables
integer(4) :: zsfc_id,u10sfc_id,v10sfc_id,T2sfc_id, lsmsfc_id, blhsfc_id ! variables
integer(4) :: T_id, U_id, V_id, W_id, R_id, Z_id ! variables
!
!***Dictionary
! 
!**Pressure
character(len=10):: nx_lon     = 'longitude'
character(len=10):: ny_lat     = 'latitude'
character(len=10):: np_pres    = 'level'
character(len=10):: nt_time    = 'time'
character(len=10):: lonname    = 'longitude'
character(len=10):: latname    = 'latitude'
character(len=10):: presname   = 'level'
character(len=10):: timename   = 'time'
character(len=10):: geopname   = 'z'
character(len=10):: relhumname = 'r'
character(len=10):: tempname   = 't'
character(len=10):: uname      = 'u'
character(len=10):: vname      = 'v'
character(len=10):: wname      = 'w'
!**Invariant
character(len=10):: sgeopname = 'z'
character(len=10):: lsmname   = 'lsm'
!**Surface
character(len=10):: tpname  = 'tp'
character(len=10):: blhname = 'blh'
!**Surface 2
character(len=10):: u10name = 'u10'
character(len=10):: v10name = 'v10'
character(len=10):: t2name  = 't2m'
!**Output netcdf
character(len=25):: outputfile
! 
! 
! 
character(len=10)::nxname, nyname,npname,ntname !dimensions

!##############################################################################################
!***User parameters
! 
! Set files path
file_surface   = "//ATLAS/ATLAS_MPI/Data/output.nc"	!  Write between quotation marks the netcdf surface file name
file_pressure  = "//ATLAS/ATLAS_MPI/Data/outputPress.nc"	!  Write between quotation marks the netcdf pressure file name
! 
! Set area NOT USED FOR NOW
! 
lonmax = 180.
lonmin = -180.
latmax = 90.
latmin = -90.
res    = 0.5 !#Resolution recomended (default in ECMWF EraInterim files)
! 
! Set time
! 
year  = 2022!!
month = 8
day   = 7!19! !# the meteo starts at 00:00 UTC of this day. Must be the first day of your meteo downloaded.
endday=  9!22! !# the meteo ends at 00:00 UTC of this day.
!
! Set name of output file
!
outputfile = "prueba.era5.nc"
!
!### The user must not change anythings belong this line
!###############################################################################################
!
!*** Open files and extract the necessary information to save later in the new netcdf file
!*** PRESSURE file
!
if(NF90_OPEN(TRIM(file_pressure),NF90_NOWRITE,ncID)/=0) &
   call runend('Error in  nf90_open pressure file')
!
!*** Extract the necessary information, and check consistency
!Por ahora no chequeo nada. asumo que lon lat y nx ny son consistentes entre todos los archivos. Despues hay que chequear
    !
    !*** Gets the Id number for the Nx Dimension. and read it
    !
    if(nf90_inq_dimid(ncID, nx_lon, nxdimid)/=NF90_NOERR) &
         call runend('Error inquiring nx dimension id in pressure file') 
    if(nf90_inquire_dimension(ncID, nxdimid, NXname, nx)/=NF90_NOERR) &
         call runend('Error inquiring NX dimension in pressure file')
    !
    !*** Gets the Id number for the Ny Dimension. and read it
    !
    if(nf90_inq_dimid(ncID, ny_lat, nydimid)/=NF90_NOERR) &
         call runend('Error inquiring ny dimension id in pressure file') 
    if(nf90_inquire_dimension(ncID, nydimid, NYname, ny)/=NF90_NOERR) &
         call runend('Error inquiring NY dimension in pressure file')
    !
    !*** Gets the Id number for the Np Dimension. and read it
    !
    if(nf90_inq_dimid(ncID, np_pres, npdimid)/=NF90_NOERR) &
         call runend('Error inquiring np dimension id in pressure file') 
    if(nf90_inquire_dimension(ncID, npdimid, npname, np)/=NF90_NOERR) &
         call runend('Error inquiring NP dimension in pressure file')
    !
    !*** Gets the Id number for the Time. Read the Time dimension. 
    !
    if(nf90_inq_dimid(ncID, nt_time, ntdimid)/=NF90_NOERR) &
         call runend('Error inquiring nt dimension id in pressure file') 
    if(nf90_inquire_dimension(ncID, ntdimid, ntname, ntm)/=NF90_NOERR) &
         call runend('Error inquiring nt dimension in pressure file')
allocate(lon(nx))
allocate(lat(ny))
allocate(pres(np))

! 
! 
! 
! 
!allocate(time(ntm))
    !*** time
!    if(nf90_inq_varid(ncID, 'time', tid)/=NF90_NOERR) &
!         call runend('Error getting time variable id in pressure file')
!    if(nf90_get_var(ncID,tid,time)/=NF90_NOERR) &
!         call runend('Error getting variable time in pressure file')
!write(*,*)ntm,time
    !
    !*** Get time-independent variables
    !*** lon
    if(nf90_inq_varid(ncID, lonname, lonvarID)/=NF90_NOERR) &
         call runend('Error getting Longitude variable id in pressure file')
    if(nf90_get_var(ncID,lonvarID,lon)/=NF90_NOERR) &
         call runend('Error getting variable Longitude in pressure file')
    if (lon(1).gt.180 .and. lon(nx).gt.180)then
       lon(:)=lon(:)-360
    end if
    !*** lat
    if(nf90_inq_varid(ncID, latname, latvarID)/=NF90_NOERR) &
         call runend('Error getting Latitude variable id in pressure file')
    if(nf90_get_var(ncID,latvarID,lat)/=NF90_NOERR) &
         call runend('Error getting variable Latitudein pressure file')
    !*** pres
    if(nf90_inq_varid(ncID, presname, presvarID)/=NF90_NOERR) &
         call runend('Error getting Times pres level id in pressure file')
    if(nf90_get_var(ncID,presvarID,pres)/=NF90_NOERR) &
         call runend('Error getting variable pres in pressure file')
    pres=pres!*101.3 ! convert mb to Pa (in atlas program)
! 

nt=(endday-day)*24 !nt2-nt1+1
allocate(timeecmwf(ntm)) !readed time
allocate(timeu(nt))      ! Time to put in the netcdf file
    !
    !*** Get times (Variable)
    !
    if(nf90_inq_varid(ncID, timename, tvarID)/=NF90_NOERR) &
         call runend('Error getting Times variable id in pressure file')
    if(nf90_get_var(ncID,tvarID,timeecmwf)/=NF90_NOERR) &
         call runend('Error getting variable Time in pressure file')
    ! 
    !*** Time vector in seconds since 00:00 of the start day. To put in the new netcdf file.
    ! 
    do i=1,nt
       timeu(i)=(i-1)*3600 
    end do
    ! 
    !*** Allocate time dependant variables
    ! 
    allocate(geopotential(nx,ny,np,ntm))
    allocate(relhum(nx,ny,np,ntm))
    allocate(temperature(nx,ny,np,ntm))
    allocate(uwind(nx,ny,np,ntm))
    allocate(vwind(nx,ny,np,ntm))
    allocate(wwind(nx,ny,np,ntm))
    ! Invariant
    allocate(work1(nx,ny,1))
    allocate(surfgeop(nx,ny,ntm))
    allocate(lsm(nx,ny,ntm))
    ! Surface
    allocate(tp(nx,ny,ntm))
    allocate(blh(nx,ny,ntm))
    allocate(work(nx,ny,ntm))
    allocate(u10(nx,ny,ntm))
    allocate(v10(nx,ny,ntm))
    allocate(t2(nx,ny,ntm))
    ! 
    allocate(tv(nx,ny,np,ntm))
    allocate(tpo(nx,ny,np,ntm))
    allocate(ro(nx,ny,np,ntm))
    !
    !*** Gets Variables ID
    !
    if(nf90_inq_varid(ncID, geopname, geopID)/=NF90_NOERR) &
         call runend('Error getting geopotential variable id in pressure file')
    !
    if(nf90_inq_varid(ncID, relhumname, relhumID)/=NF90_NOERR) &
         call runend('Error getting relative humidity variable id in pressure file')
    !
    if(nf90_inq_varid(ncID, tempname, tempID)/=NF90_NOERR) &
         call runend('Error getting temperature variable id in pressure file')
    !
    if(nf90_inq_varid(ncID, uname, uID)/=NF90_NOERR) &
         call runend('Error getting u wind variable id in pressure file')
    !
    if(nf90_inq_varid(ncID, vname, vID)/=NF90_NOERR) &
         call runend('Error getting v wind variable id in pressure file')
    !
    if(nf90_inq_varid(ncID, wname, wID)/=NF90_NOERR) &
         call runend('Error getting w wind variable id in pressure file')
    ! 
    !*** Read Geopotential variable
    !
    if( nf90_get_var(ncID,geopid,geopotential,start=(/1,1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading geopotential variable in pressure file')
    if(nf90_get_att(ncID, geopid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for geopotential variable')
    if(nf90_get_att(ncID, geopid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for geopotential variable')
    geopotential=geopotential*scale_factor+add_offset
    !geopotential=geopotential/g (in atlas program)
    ! 
    !*** Read Relative humidity variable
    !
    if( nf90_get_var(ncID,relhumid,relhum,start=(/1,1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading relative humidity variable in pressure file')
    if(nf90_get_att(ncID, relhumid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for relative humidity variable')
    if(nf90_get_att(ncID, relhumid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for relative humidity variable')
    relhum=relhum*scale_factor+add_offset 
    ! 
    !*** Read Temperature variable
    !
    if( nf90_get_var(ncID,tempid,temperature,start=(/1,1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading temperature variable in pressure file') 
    if(nf90_get_att(ncID, tempid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for temperature variable')
    if(nf90_get_att(ncID, tempid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for temperature variable')
    temperature=temperature*scale_factor+add_offset
    ! 
    !*** Read U variable
    !
    if( nf90_get_var(ncID,uid,uwind,start=(/1,1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading u wind variable in pressure file') 
    if(nf90_get_att(ncID, uid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for u wind variable')
    if(nf90_get_att(ncID, uid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for u wind variable')
    uwind=uwind*scale_factor+add_offset
    ! 
    !*** Read V variable
    !
    if( nf90_get_var(ncID,vid,vwind,start=(/1,1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading v wind variable in pressure file') 
    if(nf90_get_att(ncID, vid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for v wind variable')
    if(nf90_get_att(ncID, vid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for v wind variable')
    vwind=vwind*scale_factor+add_offset
    ! 
    !*** Read W variable
    !
    if( nf90_get_var(ncID,wid,wwind,start=(/1,1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading w wind variable in pressure file') 
    if(nf90_get_att(ncID, wid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for vertical wind variable')
    if(nf90_get_att(ncID, wid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for vertical wind variable')
    wwind=wwind*scale_factor+add_offset 
    do ip=1,np
       wwind(:,:,ip,:)=wwind(:,:,ip,:)!/pres(ip) !given in pressure levels. division in atlas program
    end do
    ! 
    !*** Close pressure file
    !  
if( nf90_close(ncID) /= NF90_NOERR ) &
    call runend('reading pressure : error closing netcdf file')
!
!*** INVARIANT file
!
if(NF90_OPEN(TRIM(file_surface),NF90_NOWRITE,ncID)/=0) &
   call runend('Error in  nf90_open surface file') 
    !
    !*** Gets Variables ID
    !
    if(nf90_inq_varid(ncID, sgeopname, sgeopID)/=NF90_NOERR) &
         call runend('Error getting geopotential variable id in invariant file')
    !
    if(nf90_inq_varid(ncID, lsmname, lsmID)/=NF90_NOERR) &
         call runend('Error getting relative humidity variable id in invariant file')
         ! 
    if(nf90_inq_varid(ncID, blhname, blhid)/=NF90_NOERR) &
         call runend('Error getting boundary layer height variable id in surface file')

    ! 
    !*** Read surface geopotential variable
    !
    if( nf90_get_var(ncID,sgeopid,work1,start=(/1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading surface geopotential variable in invariant file') 
    if(nf90_get_att(ncID, sgeopid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for surface geopotential variable')
    if(nf90_get_att(ncID, sgeopid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for surface geopotential variable')
    work1=work1*scale_factor+add_offset
    !
    do i=1,ntm!2
      surfgeop(:,:,i)=work1(:,:,1) !1
    end do
    ! 
    !*** Read land sea mask variable
    !
    if( nf90_get_var(ncID,lsmid,work1,start=(/1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading land sea mask variable in invaariant file') 
    if(nf90_get_att(ncID, lsmid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for land sea mask variable')
    if(nf90_get_att(ncID, lsmid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for land sea mask variable')
    work1=work1*scale_factor+add_offset
    !
    do i=1,ntm !2
      lsm(:,:,i)=work1(:,:,1)!1
    end do 
        ! 
    !*** Read Boundary layer height variable
    !
    if( nf90_get_var(ncID,blhid,work,start=(/1,1,1/),count=(/nx,ny,ntm/)) /=NF90_NOERR )&
         call runend('Error reading boundary layer height variable in surface file') 
    if(nf90_get_att(ncID, blhid, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for boundary layer height variable')
    if(nf90_get_att(ncID, blhid, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for boundary layer height variable')
    work=work*scale_factor+add_offset
    !work(:,:,1)=work(:,:,2)
!
	do i=1,ntm
 	  blh(:,:,i)=work(:,:,i)
	end do
! 
deallocate(work)
deallocate(work1) 
    ! 
    !*** Gets Variables ID
    !
    if(nf90_inq_varid(ncID, u10name, u10ID)/=NF90_NOERR) &
         call runend('Error getting 10 metres u wind variable id in surface file')
    !
    if(nf90_inq_varid(ncID, v10name, v10ID)/=NF90_NOERR) &
         call runend('Error getting total 10 metre v wind variable id in surface file')
    !
    if(nf90_inq_varid(ncID, t2name, t2ID)/=NF90_NOERR) &
         call runend('Error getting 2 metre temperature variable id in surface file')
    ! 
    !*** Read 10 metres u wind variable
    !
    if( nf90_get_var(ncID,u10id,u10,start=(/1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading 10 metres u wind variable in surface file') 
    if(nf90_get_att(ncID, u10id, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for 10 metres u wind variable')
    if(nf90_get_att(ncID, u10id, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for 10 metres u wind variable')
    u10=u10*scale_factor+add_offset
    ! 
    !*** Read 10 metres v wind variable
    !
    if( nf90_get_var(ncID,v10id,v10,start=(/1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading 10 metres v wind variable in surface file') 
    if(nf90_get_att(ncID, v10id, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for 10 metres v wind variable')
    if(nf90_get_att(ncID, v10id, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for 10 metres v wind variable')
    v10=v10*scale_factor+add_offset
    ! 
    !*** Read 2 metres temperature variable
    !
    if( nf90_get_var(ncID,t2id,t2,start=(/1,1,1/)) /=NF90_NOERR )&
         call runend('Error reading 2 metres temperature variable in surface file')  
    if(nf90_get_att(ncID, t2id, 'scale_factor', scale_factor) /=NF90_NOERR )&
         call runend('Error reading scale_factor attribute for 2 metres temperature variable')
    if(nf90_get_att(ncID, t2id, 'add_offset', add_offset) /=NF90_NOERR )&
         call runend('Error reading add_offset attribute for 2 metres temperature variable')
    t2=t2*scale_factor+add_offset
! 
! Close the file
if( nf90_close(ncID) /= NF90_NOERR ) &
    call runend('reading surface : error closing netcdf file')
! 
!*** Writes in a new netcdf file
     !*** Create a netcdf file 
     !
     if( nf90_create(TRIM(outputfile), NF90_CLOBBER, ncid)/=NF90_NOERR) &
          call runend('Error creting netcdf file')
     !
     !*** Define Dimensions
     !
     if( nf90_def_dim(ncid, 'lon', nx, nlon_id)/=NF90_NOERR)&
          call runend('Error creting netcdf nx dimension')
     !   
     if( nf90_def_dim(ncid, 'lat', ny, nlat_id)/=NF90_NOERR)&
          call runend('Error creting netcdf ny dimension')
     !
     if( nf90_def_dim(ncid, 'pres', np, npres_id)/=NF90_NOERR)&
          call runend('Error creting netcdf nz dimension')
     !
     if( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, nt_id)/=NF90_NOERR)&
          call runend('Error creting netcdf nt dimension')
     ! 
     !*** Define time-independent variables
     !
     if( nf90_def_var(ncID, 'lon' ,NF90_FLOAT, (/nlon_id/), lon_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable lon')
     !
     if( nf90_def_var(ncID, 'lat' ,NF90_FLOAT, (/nlat_id/), lat_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable lat')
     !
     if( nf90_def_var(ncID, 'pres' ,NF90_FLOAT, (/npres_id/), pres_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable pres')
     !
     if( nf90_def_var(ncID, 'time' ,NF90_FLOAT, (/nt_id/), time_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable time')
     !
     !*** Define Time Dependent Variables
     !
     if( nf90_def_var(ncID, 'Z:sfc' ,NF90_FLOAT, (/nlon_id,nlat_id,nt_id/), zsfc_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable geopotential at surface')
     !
     if( nf90_def_var(ncID, 'U10:sfc' ,NF90_FLOAT, (/nlon_id,nlat_id,nt_id/), u10sfc_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable u-component of wind at 10m')
     !
     if( nf90_def_var(ncID, 'V10:sfc' ,NF90_FLOAT, (/nlon_id,nlat_id,nt_id/), v10sfc_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable v-component (meridional) of wind at 10m')
     ! 
     if( nf90_def_var(ncID, 'T2:sfc' ,NF90_FLOAT, (/nlon_id,nlat_id,nt_id/), t2sfc_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable temperature at 2m')
     !
     if( nf90_def_var(ncID, 'LSM:sfc' ,NF90_FLOAT, (/nlon_id,nlat_id,nt_id/), lsmsfc_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable Land mask')
     !
     if( nf90_def_var(ncID, 'BLH:sfc' ,NF90_FLOAT, (/nlon_id,nlat_id,nt_id/), blhsfc_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable Planetary Boundary Layer Height')
     ! 
     if( nf90_def_var(ncID, 'T' ,NF90_FLOAT, (/nlon_id,nlat_id,npres_id,nt_id/), T_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable temperature"')
     !
     if( nf90_def_var(ncID, 'U' ,NF90_FLOAT, (/nlon_id,nlat_id,npres_id,nt_id/), U_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable u velocity')
     !
     if( nf90_def_var(ncID, 'V' ,NF90_FLOAT, (/nlon_id,nlat_id,npres_id,nt_id/), V_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable v velocity')
     !
     if( nf90_def_var(ncID, 'W' ,NF90_FLOAT, (/nlon_id,nlat_id,npres_id,nt_id/), W_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable Vertical omega velocity')
     !
     if( nf90_def_var(ncID, 'R' ,NF90_FLOAT, (/nlon_id,nlat_id,npres_id,nt_id/), R_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable Relative humidity')
     !
     if( nf90_def_var(ncID, 'Z' ,NF90_FLOAT, (/nlon_id,nlat_id,npres_id,nt_id/), Z_ID) /= NF90_NOERR ) &
          call runend('Error in nf90_def_variable for variable Geopotential height')
     !
     !*** Assign attribute values
     !
     attr_desc  = 'Geopotential height at surface'
     attr_units = 'gpm'
     if( nf90_put_att(ncID,zsfc_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for geopotential at surface Variable')
     if( nf90_put_att(ncID, zsfc_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for geopotential at surface Variable')
     !
     attr_desc  = 'u-component of wind at 10m'
     attr_units = 'm/s'
     if( nf90_put_att(ncID,u10sfc_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for u-component of wind at 10m Variable')
     if( nf90_put_att(ncID, u10sfc_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for u-component of wind at 10m Variable')
     !
     attr_desc  = 'v-component (meridional) of wind at 10m'
     attr_units = 'm/s'
     if( nf90_put_att(ncID, v10sfc_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for v-component (meridional) of wind at 10m Variable')
     if( nf90_put_att(ncID, v10sfc_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for v-component (meridional) of wind at 10m Variable')
     ! 
     attr_desc  = 'temperature at 2m'
     attr_units = 'K'
     if( nf90_put_att(ncID,t2sfc_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for temperature at 2m Variable')
     if( nf90_put_att(ncID, t2sfc_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for temperature at 2m Variable')
     !
     attr_desc  = 'Land mask'
     attr_units = '-'
     if( nf90_put_att(ncID,lsmsfc_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Land mask Variable')
     if( nf90_put_att(ncID, lsmsfc_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Land mask Variable')
     !
     attr_desc  = 'Planetary Boundary Layer Height'
     attr_units = 'm'
     if( nf90_put_att(ncID, blhsfc_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Planetary Boundary Layer Height Variable')
     if( nf90_put_att(ncID, blhsfc_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Planetary Boundary Layer Height Variable')
     ! 
     attr_desc  = 'temperature'
     attr_units = 'K'
     if( nf90_put_att(ncID, T_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for temperature Variable')
     if( nf90_put_att(ncID, T_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for temperature Variable')
     ! 
     attr_desc  = 'u velocity'
     attr_units = 'm/s'
     if( nf90_put_att(ncID, U_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for u velocity Variable')
     if( nf90_put_att(ncID, U_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for u velocity Variable')
     ! 
     attr_desc  = 'v velocity'
     attr_units = 'm/s'
     if( nf90_put_att(ncID, V_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for v velocity Variable')
     if( nf90_put_att(ncID, V_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for v velocity Variable')
     ! 
     attr_desc  = 'Vertical omega velocity'
     attr_units = 'Pa/s'
     if( nf90_put_att(ncID, W_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Vertical omega velocity Variable')
     if( nf90_put_att(ncID, W_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Vertical omega velocity Variable')
     ! 
     attr_desc  = 'Relative humidity'
     attr_units = '%'
     if( nf90_put_att(ncID, R_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Relative humidity Variable')
     if( nf90_put_att(ncID, R_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Relative humidity Variable')
     ! 
     attr_desc  = 'Geopotential height'
     attr_units = 'gpm'
     if( nf90_put_att(ncID, Z_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Geopotential height Variable')
     if( nf90_put_att(ncID, Z_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('Error in nf90_put_att for Geopotential height Variable')
     ! 
     !*** Put global attributes
     ! 
     attr_title = 'erainNC-1'
     if(nf90_put_att(ncID, NF90_GLOBAL, 'TITLE', attr_title)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for tittle')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'LONMIN', lonmin)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'LONMAX', lonmax)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'LATMIN', latmin)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'LATMAX', latmax)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'NX', nx)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'NY', ny)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'NP', np)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'NT', nt)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'YEAR', year)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'MONTH', month )/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     if(nf90_put_att(ncID, NF90_GLOBAL, 'DAY', day)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     hour  = 00
     if(nf90_put_att(ncID, NF90_GLOBAL, 'HOUR', hour)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     time_incr = 3600
     if(nf90_put_att(ncID, NF90_GLOBAL, 'TIME_INCR', time_incr)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     cen_lon = 0.
     if(nf90_put_att(ncID, NF90_GLOBAL, 'CEN_LON', cen_lon)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')
     ! 
     cen_lat = 0.
     if(nf90_put_att(ncID, NF90_GLOBAL, 'CEN_LAT', cen_lat)/= NF90_NOERR ) &
          call runend('Error in nf90_put_att global for year')

     !
     !*** Leave define mode 
     !
     if( nf90_enddef(ncid) /= NF90_NOERR ) call runend('Error in nf90_enddef')
     ! 
     !*** Fill values for time-independent variables
     ! 
     if( nf90_put_var(ncID, lon_ID,lon , start=(/1/), count=(/nx/) ) /= NF90_NOERR ) &
          call wriwar('Error in nf90_put_var for variable lon')
     !
     if( nf90_put_var(ncID, lat_ID,lat , start=(/1/), count=(/ny/) ) /= NF90_NOERR ) &
          call wriwar('Error in nf90_put_var for variable lat')
     !
     if( nf90_put_var(ncID, pres_ID,pres , start=(/1/), count=(/np/) ) /= NF90_NOERR ) &
          call wriwar('Error in nf90_put_var for variable pres')
     !
     if( nf90_put_var(ncID, time_ID,timeu , start=(/1/), count=(/nt/) ) /= NF90_NOERR ) &
          call wriwar('Error in nf90_put_var for variable time')
     ! 
     !*** Fill values for time-dependent variables
     ! 
  if(nf90_put_var(ncID, zsfc_ID,surfgeop(:,:,:),start=(/1,1,1/),count=(/(nx),(ny),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable geopotential surface')
! 
  if(nf90_put_var(ncID, u10sfc_ID,u10(:,:,:),start=(/1,1,1/),count=(/(nx),(ny),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable u velocity at 10m ')
! 
  if(nf90_put_var(ncID, v10sfc_ID,v10(:,:,:),start=(/1,1,1/),count=(/(nx),(ny),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable uv velocity at 10m ')
! 
  if(nf90_put_var(ncID, T2sfc_ID,t2(:,:,:),start=(/1,1,1/),count=(/(nx),(ny),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable temperature at 2m ')
! 
  if(nf90_put_var(ncID, lsmsfc_ID,lsm(:,:,:),start=(/1,1,1/),count=(/(nx),(ny),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable land sea mask')
! 
  if(nf90_put_var(ncID, blhsfc_ID,blh(:,:,:),start=(/1,1,1/),count=(/(nx),(ny),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable boundary layer height ')
! 
  if(nf90_put_var(ncID, T_ID,temperature(:,:,:,:),start=(/1,1,1,1/),count=(/(nx),(ny),(np),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable Temperature')
! 
  if(nf90_put_var(ncID, U_ID,uwind(:,:,:,:),start=(/1,1,1,1/),count=(/(nx),(ny),(np),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable U wind')
! 
  if(nf90_put_var(ncID, V_ID,vwind(:,:,:,:),start=(/1,1,1,1/),count=(/(nx),(ny),(np),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable V wind')
! 
  if(nf90_put_var(ncID, W_ID,wwind(:,:,:,:),start=(/1,1,1,1/),count=(/(nx),(ny),(np),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable W wind')
! 
  if(nf90_put_var(ncID, R_ID,relhum(:,:,:,:),start=(/1,1,1,1/),count=(/(nx),(ny),(np),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable Relative humidity')
! 
  if(nf90_put_var(ncID, Z_ID,geopotential(:,:,:,:),start=(/1,1,1,1/),count=(/(nx),(ny),(np),(nt)/))/= NF90_NOERR) &
       call wriwar('Error in nf90_put_var for variable geopotential')
  ! 
  !*** CLOSE
  ! 
  if( nf90_close(ncid)/=NF90_NOERR) call runend('Error closing netcdf file')
call runend("OK")
end program erainNC
