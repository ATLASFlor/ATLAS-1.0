!***************************************************************
!*
!*              Module for WRF Meteo
!* 
!***************************************************************
MODULE Wrf
  use KindType
  use Master
  use InpOut
  IMPLICIT NONE
  SAVE

  !
  !*** WRF DICTIONARY
  !
  ! WRF DIMENSIONS:
  character(len=25) :: nx_wrf_name      = 'west_east'
  character(len=25) :: nx_stag_wrf_name = 'west_east_stag'
  character(len=25) :: ny_wrf_name      = 'south_north'
  character(len=25) :: ny_stag_wrf_name = 'south_north_stag'
  character(len=25) :: np_wrf_name      = 'bottom_top'
  character(len=25) :: np_stag_wrf_name = 'bottom_top_stag'
  character(len=25) :: nt_wrf_name      = 'Time'
  !
  !WRF VARIABLES:
  character(len=25) :: time_wrf_name    = 'Times'
  character(len=25) :: lon_wrf_name     = 'XLONG'
  character(len=25) :: lat_wrf_name     = 'XLAT'
  character(len=25) :: hgt_wrf_name     = 'HGT'
  character(len=25) :: landmask_wrf_name= 'LANDMASK'
  character(len=25) :: ph_wrf_name      = 'PH'
  character(len=25) :: phb_wrf_name     = 'PHB'
  character(len=25) :: u_wrf_name       = 'U'
  character(len=25) :: v_wrf_name       = 'V'
  character(len=25) :: w_wrf_name       = 'W'
  character(len=25) :: T_wrf_name       = 'T'
  character(len=25) :: p_wrf_name       = 'P'
  character(len=25) :: pb_wrf_name      = 'PB'
  character(len=25) :: qv_wrf_name      = 'QVAPOR'
  character(len=25) :: pblh_wrf_name    = 'PBLH'
  character(len=25) :: ust_wrf_name     = 'UST'
  character(len=25) :: hfx_wrf_name     = 'HFX'
  !
  !WRF ATTRIBUTES:
  character(len=25) :: DX_wrf_name      = 'DX'
  character(len=25) :: DY_wrf_name      = 'DY'
  !
CONTAINS
  !
  !
  integer(ip) function stoi1(string1)
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
  !
  subroutine openwrf(nm)
    !******************************************************************************************
    !*    OPENWRF
    !*    Opens the Meteorological file to create the meteo mesh
    !*
    !*    nm       Correspond to the meteo file readed by input file  (Input value)
    !*
    !******************************************************************************************
    use Master
    use KindType
    use netcdf
    use TimeFun
    use Elsest_mod
    implicit none
    !
    integer  (ip)         :: nm
    integer  (ip)         :: ielem,ix,iy,ipoin,it,it_WRF
    integer  (ip)         :: ibyr_WRF,ibmo_WRF,ibdy_WRF,ibhr_WRF,ibmi_WRF, ibse_WRF
    integer  (ip)         :: iyr_WRF,imo_WRF,idy_WRF,ihr_WRF,imi_WRF
    integer  (ip)         :: iyr,imo,idy,ihr,imi, ise
    integer  (ip)         :: iyr2,imo2,idy2,ihr2
    !
    integer  (ip)         :: ncid                   ! Id Number of the meteo mesh
    integer  (ip)         :: NXdimid                ! NX Dimension Id. 
    integer  (ip)         :: NXstag_dimid           ! NX_stag Dimension Id. 
    integer  (ip)         :: NYdimid                ! NY Dimension Id. 
    integer  (ip)         :: NYstag_dimid           ! NY_stag Dimension Id. 
    integer  (ip)         :: NPdimid                ! NP Dimension Id (pressure levels)
    integer  (ip)         :: NPstag_dimid           ! NP_stag Dimension Id (pressure levels). 
    integer  (ip)         :: Ntdimid                ! Dt Dimension Id. 
    !
    integer  (ip)         :: tvarid                 ! Times    Variable Id.
    integer  (ip)         :: longid                 ! Longitud Variable Id.
    integer  (ip)         :: latid                  ! Latitud  Variable Id.
    integer  (ip)         :: hgtid                  ! Height   Variable Id.
    integer  (ip)         :: landmaskid             ! Height   Variable Id.
    !
    real     (rp)         :: DXattva                ! DX attribute. 
    real     (rp)         :: DYattva                ! DY attribute. 
    !
    !
    character(len=s_mess) :: NXname, NYname,NPname    ! Names of dimension NX, NY and NP
    character(len=s_mess) :: NXstag_name, NYstag_name, NPstag_name ! Names of dimension NX_stag, NY_stag
    character(len=s_mess) :: Ntname                   ! Name of time step dimension.
    character(len=s_mess) :: str
    character(len=1     ) :: ext                      ! 
    !
    character(len=19), allocatable :: time_WRF_string(:)
    real(rp),          allocatable :: WRF_coord(:,:)
    real(rp), allocatable          :: work2d (:,:)  ! Local array used for topography variable
    real(rp), allocatable          :: lwork2d(:,:)  ! Local array used for Landmask variable
    !
    !    Time Variables
    !  
    logical     :: found
    real   (rp) :: dt_WRF
    !
    !    Elsest variables
    !
    real   (rp) :: xcoor(3), coloc(3)
    real   (rp) :: shapf(4), deriv(2,4), relse(20)
    integer(ip) :: ielse(10), ltopo, ifoun, i
    !
    !*** Elsest initializations
    !    
    ielse(1)      = 10                                     ! number nx boxes (search)
    ielse(2)      = 10                                     ! ny
    ielse(3)      = 10                                     ! nz
    ielse(4)      = 1                                      ! data format (0=type,1=linked list)
    ielse(5)      = 10                                     ! Maximum number of possible meshes
    ielse(6)      = 1                                      ! Not used
    ielse(7)      = 0                                      ! Output unit
    ielse(8)      = 0                                      ! Search strategy (0=bin,1=Quad)
    !
    relse(1)      = 0.01_rp                                ! Tolerance for iteration
    relse(3)      =-1e8_rp                                 ! Patch bounding box
    relse(4)      = 1e8_rp                                   
    relse(5)      =-1e8_rp
    relse(6)      = 1e8_rp
    relse(7)      =-1e8_rp
    relse(8)      = 1e8_rp
    !
    ltopo = 0 ! (0 quads hexahedra , 1 for triangle tetra, 2 Pentahedra, 3 Pyramid)
    !
    !*** Open the WRF model. Gets the WRF file ID 
    !
    if(NF90_OPEN(meteo(nm)%file_path,NF90_NOWRITE,ncID)/=0) &
         call runend('Error in  nf90_open wrf file') 
    !
    write(ext,'(i1)') nm
    str = 'METEO_DATA_'//ext
    !
    !*** Gets the Id number for the Nx Dimension. and read it
    !
    if(nf90_inq_dimid(ncID, nx_wrf_name, NXdimid)/=NF90_NOERR) &
         call runend('Error inquiring nx dimension id in wrf file') 
    if(nf90_inquire_dimension(ncID, NXdimid, NXname, meteo(nm)%nx)/=NF90_NOERR) &
         call runend('Error inquiring NX dimension in wrf file')
    !
    !*** Gets the Id number for the Nx_stag Dimension. and read it
    !
    if(nf90_inq_dimid(ncID, nx_stag_wrf_name, NXstag_dimid)/=NF90_NOERR) &
         call runend('Error inquiring nx_stag dimension id in wrf file') 
    if(nf90_inquire_dimension(ncID, NXstag_dimid, NXstag_name, meteo(nm)%nx_stag)/=NF90_NOERR) &
         call runend('Error inquiring NX_stag dimension in wrf file')
    !
    !*** Gets the Id number for the Ny Dimension.
    !
    if(nf90_inq_dimid(ncID, ny_wrf_name, NYdimid)/=NF90_NOERR) &
         call runend('Error inquiring ny dimension id in wrf file') 
    if(nf90_inquire_dimension(ncID, NYdimid, NYname, meteo(nm)%ny)/=NF90_NOERR) &
         call runend('Error inquiring NY dimension in wrf file')
    !
    !*** Gets the Id number for the Ny_stag Dimension.
    !
    if(nf90_inq_dimid(ncID, ny_stag_wrf_name, NYstag_dimid)/=NF90_NOERR) &
         call runend('Error inquiring ny_stag dimension id in wrf file') 
    if(nf90_inquire_dimension(ncID, NYstag_dimid, NYstag_name, meteo(nm)%ny_stag)/=NF90_NOERR) &
         call runend('Error inquiring NY_stag dimension in wrf file')
    !     
    !*** Gets the Id number for the number of pressure levels Dimension and read it. (np)
    !
    if(nf90_inq_dimid(ncID, np_wrf_name, NPdimid)/=NF90_NOERR) &
         call runend('Error inquiring np dimension id in wrf file') 
    if(nf90_inquire_dimension(ncID, NPdimid, NPname, meteo(nm)%np)/=NF90_NOERR) &
         call runend('Error inquiring NP dimension in wrf file')
    !     
    !*** Gets the Id number for the number of staggered pressure levels Dimension and read it. (npstag)
    !
    if(nf90_inq_dimid(ncID, np_stag_wrf_name, NPstag_dimid)/=NF90_NOERR) &
         call runend('Error inquiring np dimension id in wrf file') 
    if(nf90_inquire_dimension(ncID, NPstag_dimid, NPstag_name, meteo(nm)%np_stag)/=NF90_NOERR) &
         call runend('Error inquiring NP dimension in wrf file')
    !
    !*** Gets the ID number for Time dimension. Read nt
    ! 
    if(nf90_inq_dimid(ncID, nt_wrf_name, Ntdimid)/=NF90_NOERR) &
         call runend('Error inquiring nt dimension id in wrf file') 
    if(nf90_inquire_dimension(ncID, Ntdimid, Ntname, meteo(nm)%nt)/=NF90_NOERR) &
         call runend('Error inquiring Nt dimension in wrf file')
    !
    !*** Allocate memory for module arrays
    !
    meteo(nm)%nelem2d = (meteo(nm)%nx - 1)*(meteo(nm)%ny - 1)
    meteo(nm)%npoin2d = meteo(nm)%nx*meteo(nm)%ny
    !
    allocate (meteo(nm)%lon    (meteo(nm)%nx,meteo(nm)%ny))
    allocate (meteo(nm)%lat    (meteo(nm)%nx,meteo(nm)%ny))
    allocate (meteo(nm)%time   (meteo(nm)%nt             ))
    allocate (meteo(nm)%timesec(meteo(nm)%nt            ))
    !
    allocate (meteo(nm)%hgt     (meteo(nm)%npoin2d))
    allocate (meteo(nm)%landmask(meteo(nm)%npoin2d))
    !
    allocate (meteo(nm)%lnods(meteo(nm)%nnode,meteo(nm)%nelem2d))
    !
    allocate (meteo(nm)%z (meteo(nm)%npoin2d,meteo(nm)%np))   ! stored in 2D planes to interpolate with Elsest
    allocate (meteo(nm)%u (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%v (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%w (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%T (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%p (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%Tv(meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%ro(meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%qv(meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%tp(meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%pblh(meteo(nm)%npoin2d))
    allocate (meteo(nm)%hfx(meteo(nm)%npoin2d))
    allocate (meteo(nm)%ust(meteo(nm)%npoin2d))
    !
    !***  Allocate memory for local arrays (need to be deallocated later)
    !
    allocate(time_WRF_string(meteo(nm)%nt))
    allocate(WRF_coord(2,meteo(nm)%npoin2d))
    allocate(work2d (meteo(nm)%nx,meteo(nm)%ny))
    allocate(lwork2d(meteo(nm)%nx,meteo(nm)%ny))
    !
    !*** Get times
    !
    if(nf90_inq_varid(ncID, time_wrf_name, tvarID)/=NF90_NOERR) &
         call runend('Error getting Times variable id in openwrf')
    if(nf90_get_var(ncID,tvarID,time_WRF_string)/=NF90_NOERR) &
         call runend('Error getting variable Time in openwrf')
    !  
    !*** Convert the string yyyy-mm-dd_hh:mm:ss to the integer YYYYMMDDHHMM for the first time
    !
    ibyr_WRF = stoi1(time_WRF_string(1)(1 :1 ))*1000 + &
         stoi1(time_WRF_string(1)(2 :2 ))*100  + &
         stoi1(time_WRF_string(1)(3 :3 ))*10   + &
         stoi1(time_WRF_string(1)(4 :4 ))
    ibmo_WRF = stoi1(time_WRF_string(1)(6 :6 ))*10   + &
         stoi1(time_WRF_string(1)(7 :7 ))
    ibdy_WRF = stoi1(time_WRF_string(1)(9 :9 ))*10   + &
         stoi1(time_WRF_string(1)(10:10))
    ibhr_WRF = stoi1(time_WRF_string(1)(12:12))*10   + &
         stoi1(time_WRF_string(1)(13:13))
    ibmi_WRF = stoi1(time_WRF_string(1)(15:15))*10   + &
         stoi1(time_WRF_string(1)(16:16))
    !
    meteo(nm)%time(1) = 1d8*ibyr_WRF + 1d6*ibmo_WRF + 1d4*ibdy_WRF + 1d2*ibhr_WRF + ibmi_WRF
    !
    !*** Averiguates the WRF time step by guess (1h increment iteration) and
    !*** comparing with the second value
    !
    if(meteo(nm)%nt.gt.1) then
       iyr2 = stoi1(time_WRF_string(2)(1 :1 ))*1000 + &
            stoi1(time_WRF_string(2)(2 :2 ))*100  + &
            stoi1(time_WRF_string(2)(3 :3 ))*10   + &
            stoi1(time_WRF_string(2)(4 :4 ))
       imo2 = stoi1(time_WRF_string(2)(6 :6 ))*10   + &
            stoi1(time_WRF_string(2)(7 :7 ))
       idy2 = stoi1(time_WRF_string(2)(9 :9 ))*10   + &
            stoi1(time_WRF_string(2)(10:10))
       ihr2 = stoi1(time_WRF_string(2)(12:12))*10   + &
            stoi1(time_WRF_string(2)(13:13))
       !
       dt_WRF = 0.
       found = .false.
       do while(.not.found)
          dt_WRF = dt_WRF + 3600.
          call addtime(ibyr_WRF,ibmo_WRF,ibdy_WRF,ibhr_WRF, &
               iyr,imo,idy,ihr,imi,ise,dt_WRF)
          if( (iyr.eq.iyr2).and.(imo.eq.imo2).and.(idy.eq.idy2).and.(ihr.eq.ihr2) ) found = .true.
          if( dt_WRF.gt.(24.*3600.) ) call runend('read_WRF_grid : unable to find dt_WRF')
       end do
    else
       dt_WRF = 0.
    end if
    !
    !*** Calculates timesec_WRF(:) in sec after 0000UTC
    !
    meteo(nm)%timesec(1) = 3600.*ibhr_WRF + 60.*ibmi_WRF + ibse_WRF
    do it_WRF = 2,meteo(nm)%nt
       meteo(nm)%timesec(it_WRF) = meteo(nm)%timesec(1) + (it_WRF-1)*dt_WRF
    end do
    !
    !*** Convert the string yyyy-mm-dd_hh:mm:ss to the integer YYYYMMDDHHMM
    !
    do it = 2,meteo(nm)%nt 
       iyr_WRF = stoi1(time_WRF_string(it)(1 :1 ))*1000 + &
            stoi1(time_WRF_string(it)(2 :2 ))*100  + &
            stoi1(time_WRF_string(it)(3 :3 ))*10   + &
            stoi1(time_WRF_string(it)(4 :4 ))
       imo_WRF = stoi1(time_WRF_string(it)(6 :6 ))*10   + &
            stoi1(time_WRF_string(it)(7 :7 ))
       idy_WRF = stoi1(time_WRF_string(it)(9 :9 ))*10   + &
            stoi1(time_WRF_string(it)(10:10))
       ihr_WRF = stoi1(time_WRF_string(it)(12:12))*10   + &
            stoi1(time_WRF_string(it)(13:13))
       imi_WRF = stoi1(time_WRF_string(it)(15:15))*10   + &
            stoi1(time_WRF_string(it)(16:16))
       !
       meteo(nm)%time(it) = 1d8*iyr_WRF + 1d6*imo_WRF + 1d4*idy_WRF + 1d2*ihr_WRF + imi_WRF
       !
    end do
    !
    !***  Calculates time lag by iteration (that is, the time in seconds between the WRF origin
    !***  and the Simulation origin). Both origins are referred to 0000UTC, but may belong to different
    !***  days
    !
    found = .false.
    meteo(nm)%time_lag = 0.0_rp
    do while(.not.found)
       call addtime(ibyr_WRF,ibmo_WRF,ibdy_WRF,0, &
            iyr,imo,idy,ihr,imi,ise,meteo(nm)%time_lag)
       !
       if( (iyr.eq.ibyr).and.(imo.eq.ibmo).and.(idy.eq.ibdy) ) then
          found = .true.
       else
          meteo(nm)%time_lag = meteo(nm)%time_lag + 86400.0_rp
       end if
       if(meteo(nm)%time_lag.gt.meteo(nm)%timesec(meteo(nm)%nt))&
            call runend('Time lag not found. Check interval consistency')
    end do
    !
    !*** Checks that the required time is within the bounds
    !
    if(sidt.gt.0)then
       if((meteo(nm)%time_lag+(sist )).lt.meteo(nm)%timesec(1           )) &
            call runend('Minimum datatime is after required time:'//TRIM(str))
       if((meteo(nm)%time_lag+(siend)).gt.meteo(nm)%timesec(meteo(nm)%nt)) &
            call runend('Maximum datatime is before required time:'//TRIM(str))
    else if(sidt.lt.0)then
       if((meteo(nm)%time_lag+(siend)).lt.meteo(nm)%timesec(1           )) &
            call runend('Minimum datatime is after required time:'//TRIM(str))
       if((meteo(nm)%time_lag+(sist )).gt.meteo(nm)%timesec(meteo(nm)%nt)) &
            call runend('Maximum datatime is before required time:'//TRIM(str))
    end if
    !
    !*** Redifine the timesec adding the time_lag. 
    !*** So, the meteo file and the simulation have the same origin.
    !
    meteo(nm)%timesec(:) = meteo(nm)%timesec(:) - meteo(nm)%time_lag
    !
    !*** Get wrf lon at box center
    !
    if(nf90_inq_varid(ncID, lon_wrf_name, longID)/=NF90_NOERR) &
         call runend('Error getting XLONG variable id in openwrf')
    if(nf90_get_var(ncID,longID,meteo(nm)%lon)/=NF90_NOERR) &
         call runend('Error getting variable XLONG in openwrf')
    !
    !*** Get wrf lat at box center
    !
    if(nf90_inq_varid(ncID, lat_wrf_name, latID)/=NF90_NOERR) &
         call runend('Error getting XLAT variable id in openwrf')
    if(nf90_get_var(ncID,latID,meteo(nm)%lat)/=NF90_NOERR) &
         call runend('Error getting variable XLAT in openwrf')
    !
    !*** Get Topography for wrf file at box center
    !
    if(nf90_inq_varid(ncID, hgt_wrf_name, hgtID)/=NF90_NOERR) &
         call runend('Error getting HGT variable id in openwrf')
    if(nf90_get_var(ncID,hgtID,work2d)/=NF90_NOERR) &
         call runend('Error getting variable HGT in openwrf')
    !
    !*** Get LANDMASK for wrf file at box center. 1 for LAND, 0 for WATER.
    !
    if(nf90_inq_varid(ncID, landmask_wrf_name, landmaskID)/=NF90_NOERR) &
         call runend('Error getting LANDMASK variable id in openwrf')
    if(nf90_get_var(ncID,landmaskID,lwork2d)/=NF90_NOERR) &
         call runend('Error getting variable LANDMASK in openwrf')
    work2d=work2d*lwork2d

    ipoin = 0
    do iy = 1,meteo(nm)%ny
       do ix = 1,meteo(nm)%nx
          ipoin = ipoin + 1
          meteo(nm)%hgt     (ipoin) = work2d (ix,iy)
          meteo(nm)%landmask(ipoin) = lwork2d(ix,iy)
       end do
    end do

    !
    !*** Get the Dx and Dy Attributes to calculate the model resolution
    !
    if(nf90_get_att(ncID,NF90_GLOBAL,DX_wrf_name, DXattva)/=NF90_NOERR) &
         call runend('Error getting attribute DX in openwrf')
    !
    if(nf90_get_att(ncID,NF90_GLOBAL,DY_wrf_name, DYattva)/=NF90_NOERR) &
         call runend('Error getting attribute DY in openwrf')
    !
    meteo(nm)%res  = DYattva*360/(2*pi*Rearth)
    !
    !*** Builds mesh structure. 2d
    !
    ! Associate each element with 4 nods
    ielem = 0    
    do iy = 1,meteo(nm)%ny-1
       do ix = 1,meteo(nm)%nx-1
          ielem = ielem + 1
          meteo(nm)%lnods(1,ielem) = (iy-1)*meteo(nm)%nx + ix
          meteo(nm)%lnods(2,ielem) = (iy-1)*meteo(nm)%nx + ix + 1
          meteo(nm)%lnods(3,ielem) = (iy  )*meteo(nm)%nx + ix + 1
          meteo(nm)%lnods(4,ielem) = (iy  )*meteo(nm)%nx + ix
       end do
    end do
    !
    ! Associate each point with their coordinates. 1: lon. 2: lat.
    ipoin = 0
    do iy = 1,meteo(nm)%ny
       do ix = 1,meteo(nm)%nx
          ipoin = ipoin + 1
          WRF_coord(1,ipoin) = meteo(nm)%lon(ix,iy)
          WRF_coord(2,ipoin) = meteo(nm)%lat(ix,iy)
       end do
    end do
    !
    !*** Computes the interpolation functions in 2D for all points of the background mesh ON A PLANE 
    !*** to interpolate WRF variables (later on)
    !
    do ipoin = 1,grid%npoin2d            ! loop over background mesh points in 2D 
       !
       xcoor(1) = grid%coord(1,ipoin)
       xcoor(2) = grid%coord(2,ipoin)
       !
       ! elses 2_ search element.  
       !  returns:
       !    ielem > 0  element number
       !    ielem < 0  is not in the domain  
       ! 
       call elsest(&
            2_ip,ielse,meteo(nm)%nnode,2_ip,meteo(nm)%npoin2d,meteo(nm)%nelem2d,&
            meteo(nm)%lnods,ltopo,WRF_coord,xcoor,relse,&
            ielem,shapf,deriv,coloc)
       !
       ! Restrict the shape factor to [0,1]
       !
       do i=1,4
          shapf(i) = min(shapf(i),1.0)
          shapf(i) = max(shapf(i),0.0)
       end do
       !
       !*** Put the flag 'nm' in the background grid elements which are in this meteo domain and have better resolution than others. 
       !
       if(ielem.gt.0) then
          if(grid%model(ipoin).eq.-1) then
             grid%model  (ipoin)   = nm
             grid%element(ipoin)   = ielem
             grid%shape(1:4,ipoin) = shapf(1:4)
          else if(meteo(grid%model(ipoin))%res.gt.meteo(nm)%res) then
             grid%model  (ipoin)   = nm
             grid%element(ipoin)   = ielem
             grid%shape(1:4,ipoin) = shapf(1:4)
          end if
       end if
       !
    end do
    !
    !*** Deallocates 
    !
    deallocate(time_WRF_string)
    deallocate(WRF_coord)
    deallocate(work2d)
    deallocate(lwork2d)
    !
    !***    Release elsest memory
    !
    call elsest_deallo()
    !
    !*** Close the netcdf file
    !
    if( nf90_close(ncID) /=NF90_NOERR) call runend('read_WRF : Error in nf90_close')
    !
    !
    !*** Print Information in log file
    !
    !***  Writes the log file
    !
    write(lulog,2)  'Opened WRF file and readed necesary Time Independent variables. Created mesh structure.'
    if(out_screen) write(*,2) 'Opened WRF file and readed necesary Time Independent variables. Created mesh structure.'
    !
2   format(/,  a /)
    !
    return
  end subroutine openwrf
  !
  !
  !
  subroutine readwrf(nm,itime)
    !******************************************************************************************
    !*    READWRF
    !*    Opens the Meteorological file and read time dependent variables for time itime
    !*    nm: correspond to the number nm meteorological file introduced in the .inp
    !*    itime: Correspond to the time step to be read.
    !*
    !******************************************************************************************
    use Master
    use KindType
    use netcdf
    use TimeFun
    implicit none
    !
    integer  (ip)         :: nm, itime
    integer  (ip)         :: ipoin, iz, ix, iy, i, ielem
    integer  (ip)         :: jpoin, jz, kz
    real     (rp)         :: zshape
    logical               :: found
    character(len=s_mess) :: str
    character(len=1     ) :: ext  
    !
    integer(ip) :: ncID         ! Netcdf file ID
    integer(ip) :: phbid        ! base-state Geopotential id
    integer(ip) :: phid         ! Perturbation Geopotential id
    integer(ip) :: uid          ! Horizontal winds id
    integer(ip) :: vid          ! Vertical winds id
    integer(ip) :: wid          ! Omega
    integer(ip) :: tid          ! Potencial Temperature
    integer(ip) :: pid          ! perturbation pressure id
    integer(ip) :: pbid         ! base state pressure id
    integer(ip) :: qvid         ! base state pressure id
    integer(ip) :: pblhid       ! planetary boundary layer heights id
    integer(ip) :: hfxid        ! Surface heat flux id
    integer(ip) :: ustid        ! U* id
    !
    real    (rp), allocatable :: p(:,:)
    real    (rp), allocatable :: work_WRF (:,:,:)
    real    (rp), allocatable :: work_WRF1(:,:,:)
    real    (rp), allocatable :: work_WRF0(:,:,:)
    real    (rp), allocatable :: work_WRF2(:,:)
    real    (rp), allocatable :: work_WRF3(:,:)
    !
    !*** Initializations
    !
    write(ext,'(i1)') nm
    str = 'METEO_DATA_'//ext
    !
    !*** Open the WRF model. Gets the WRF file ID 
    !
    if(NF90_OPEN(meteo(nm)%file_path,NF90_NOWRITE,ncID)/=NF90_NOERR) &
         call runend('Error in  nf90_open wrf file') 
    !
    !*** Gets Variables ID
    !
    if( nf90_inq_varid(ncID,phb_WRF_name,phbID) /= NF90_NOERR)&
         call runend('Error getting phbID in '//TRIM(str))
    !
    if( nf90_inq_varid(ncID,ph_WRF_name,phID) /= NF90_NOERR) &
         call runend('Error getting phID in '//TRIM(str) )
    !
    if(nf90_inq_varid(ncID, u_wrf_name, uID)/=NF90_NOERR) &
         call runend('Error getting U variable id in readwrf')
    !
    if(nf90_inq_varid(ncID, v_wrf_name, vID)/=NF90_NOERR) &
         call runend('Error getting V variable id in readwrf')
    !
    if(nf90_inq_varid(ncID, w_wrf_name, wID)/=NF90_NOERR) &
         call runend('Error getting W variable id in readwrf')
    !
    !*** Allocate memory
    !
    allocate(work_WRF(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np_stag))
    !
    !*** Read Heights. This is necessary because pressure levels change with time.
    !*** Reads geopotentials PHB, PH and computes z and zstag at time
    !
    allocate (work_WRF1(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np_stag))   ! PHB(itime, bottom_top_stag, south_north, west_east)
    !
    if( nf90_get_var(ncID,phbID,work_WRF,start=(/1,1,1,itime/)) /=NF90_NOERR ) &
         call runend('Error reading phb in '//TRIM(str)) 
    !
    if( nf90_get_var(ncID,phID,work_WRF1,start=(/1,1,1,itime/)) /=NF90_NOERR ) &
         call runend('Error reading ph in '//TRIM(str))
    work_WRF1 = ( work_WRF1 + work_WRF ) / 9.81_rp           ! zstag = (ph+phb)/g at time   height in the staggered grid
    !
    work_WRF1 = 0.5_rp*( work_WRF1(:,:,1:meteo(nm)%np_stag-1) &
         +     work_WRF1(:,:,2:meteo(nm)%np_stag) )   ! z at time
    !
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%z(ipoin,iz)= work_WRF1(ix,iy,iz) - meteo(nm)%hgt(ipoin)   ! Height above the terrain (in m) 
          end do
       end do
    end do
    if(MAXVAL(meteo(nm)%z(:,:)).lt.(grid%ztop)) &
         call wriwar('Input max vertical layer is higher than meteo max layer. Values extrapolated in '//TRIM(str))
    !
    deallocate(work_WRF1)
    deallocate(work_WRF )
    !
    !*** Read U Variable: Horizontal wind component 
    !
    allocate(work_WRF1(meteo(nm)%nx_stag,meteo(nm)%ny,meteo(nm)%np))
    !
    if( nf90_get_var(ncID,uID,work_WRF1,start=(/1,1,1,itime/)) /=NF90_NOERR )&
         call runend('Error reading u in '//TRIM(str))
    !
    work_WRF1  = 0.5_rp*(work_WRF1(1:meteo(nm)%nx_stag-1,1:meteo(nm)%ny,1:meteo(nm)%np)+ &
         work_WRF1(2:meteo(nm)%nx_stag  ,1:meteo(nm)%ny,1:meteo(nm)%np))      ! u box at time 
    !
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%u(ipoin,iz)= work_WRF1(ix,iy,iz)
          end do
       end do
    end do
    deallocate(work_WRF1)
    !
    !*** Read V Variable: Vertical wind component
    !
    allocate(work_WRF1(meteo(nm)%nx,meteo(nm)%ny_stag,meteo(nm)%np))
    !
    if( nf90_get_var(ncID,vID,work_WRF1,start=(/1,1,1,itime/)) /=NF90_NOERR )&
         call runend('Error reading v in '//TRIM(str))
    !
    work_WRF1 = 0.5_rp*(work_WRF1(1:meteo(nm)%nx,1:meteo(nm)%ny_stag-1,1:meteo(nm)%np)+ &
         work_WRF1(1:meteo(nm)%nx,2:meteo(nm)%ny_stag  ,1:meteo(nm)%np))      ! v box at time 
    !
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%v(ipoin,iz)= work_WRF1(ix,iy,iz)
          end do
       end do
    end do
    deallocate(work_WRF1)
    !
    !*** Read W Variable: 
    !
    allocate(work_WRF1(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np_stag))
    !
    if( nf90_get_var(ncID,wID,work_WRF1,start=(/1,1,1,itime/)) /=NF90_NOERR )&
         call runend('Error reading w in '//TRIM(str))
    !
    work_WRF1 = 0.5_rp*(work_WRF1(1:meteo(nm)%nx,1:meteo(nm)%ny,1:meteo(nm)%np_stag-1)+ &
         work_WRF1(1:meteo(nm)%nx,1:meteo(nm)%ny,2:meteo(nm)%np_stag))      ! v box at time 
    !
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%w(ipoin,iz)= work_WRF1(ix,iy,iz)
          end do
       end do
    end do
    deallocate(work_WRF1)
    !
    allocate(work_WRF (meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np))
    allocate(work_WRF1(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np))
    !
    !*** p(nx_WRF,ny_WRF,nz_WRF) --> p(nx,ny,nz) Pressure (reference plus perturbed)
    !
    if( nf90_inq_varid(ncID,p_WRF_name,pID) /= 0) call runend('readWRF : Error getting p ID')
    if( nf90_get_var(ncID,pID,work_WRF,start=(/1,1,1,itime/)) /=0 ) call runend('read_WRF: Error reading p')
    !
    if( nf90_inq_varid(ncID,pb_WRF_name,pbID) /= 0) call runend('read_WRF : Error getting pbID')
    if( nf90_get_var(ncID,pbID,work_WRF1,start=(/1,1,1,itime/)) /=0 ) call runend('read_WRF: Error reading pb')
    work_WRF = work_WRF + work_WRF1  ! p+pb 
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%p(ipoin,iz)= work_WRF(ix,iy,iz)
          end do
       end do
    end do
    !
    !***  meteo(nm)%T(npoin2d,meteo(nm)%np) --> Tp(npoin2d,meteo(nm)%np) potential temperature and T(nx,ny,nz) temperature
    !***  It is also necessary to use pressure (stored in work_WRF) because T is given as theta temperature
    !***  (minus a reference value = 300)
    !***
    !***  T(:,:)=(Tp(:,:)+300)*((p(:,:)+pb(:,:))/pref_theta)**(r/cp)
    !
    if( nf90_inq_varid(ncID,t_WRF_name,tID) /= 0) call runend('readWRF : Error getting t ID')
    if( nf90_get_var(ncID,tID,work_WRF1,start=(/1,1,1,itime/)) /=0 ) call runend('readWRF: Error reading T')
    !
    work_WRF =  (work_WRF1+300.0_rp)*((work_WRF/1e5_rp)**(0.285_rp))    ! T box
    !
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%T(ipoin,iz) = work_WRF(ix,iy,iz)
             meteo(nm)%tp(ipoin,iz)= work_WRF1(ix,iy,iz)
          end do
       end do
    end do

    !
    !*** qv(nx_WRF,ny_WRF,nz_WRF) --> qv(nx,ny,nz) specific humidity. WRF gives the Water vapor mixing ratio (e),
    !*** so the conversion q = e/1+e is needed
    !
    if( nf90_inq_varid(ncID,qv_WRF_name,qvID) /= 0) call runend('read_WRF : Error getting qvID')
    if( nf90_get_var(ncID,qvID,work_WRF1,start=(/1,1,1,itime/)) /=0 ) call runend('read_WRF: Error reading qv')
    work_WRF1 = work_WRF1/(1.0_rp + work_WRF1)                 ! specific humidity
    !
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%qv(ipoin,iz)= work_WRF1(ix,iy,iz)
          end do
       end do
    end do
    !
    !*** 2D scalar variables, given at th box center
    !
    allocate(work_WRF3(meteo(nm)%nx,meteo(nm)%ny))
    !
    !*** pblh(nx_WRF,ny_WRF) --> pblh(nx,ny) Atmospheric Boundary Layer
    !
    if( nf90_inq_varid(ncID,pblh_WRF_name,pblhID) /= 0) call runend('read_WRF : Error getting pblhID')
    if( nf90_get_var(ncID,pblhID,work_WRF3,start=(/1,1,itime/)) /=0 ) call runend('read_WRF: Error reading pblh')
    !
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%pblh(ipoin)= work_WRF3(ix,iy)
          end do
       end do
    !
    !***hfx(nx_WRF,ny_WRF) --> hfx(nx,ny) Surface heat flux
    !
    if( nf90_inq_varid(ncID,hfx_WRF_name,hfxID) /= 0) call runend('read_WRF : Error getting hfxID')
    if( nf90_get_var(ncID,hfxID,work_WRF3,start=(/1,1,itime/)) /=0 ) call runend('read_WRF: Error reading hfx')
    !
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%hfx(ipoin)= work_WRF3(ix,iy)
          end do
       end do
    !
    !***ust(nx_WRF,ny_WRF) --> ust(nx,ny) U* Friction velocity
    !
    if( nf90_inq_varid(ncID,ust_WRF_name,ustID) /= 0) call runend('read_WRF : Error getting ustID')
    if( nf90_get_var(ncID,ustID,work_WRF3,start=(/1,1,itime/)) /=0 ) call runend('read_WRF: Error reading ust')
    !
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             meteo(nm)%ust(ipoin)= work_WRF3(ix,iy)
          end do
       end do
    !
    deallocate(work_WRF1)
    deallocate(work_WRF)
    deallocate(work_WRF3)
    !
    allocate(work_WRF2(meteo(nm)%npoin2d,meteo(nm)%np))
    !
    !*** Tv(nx,ny,nz) virtual temperature
    !
    work_WRF2 = meteo(nm)%qv(:,:)/(1.0_rp-meteo(nm)%qv(:,:))      ! Mixing ratio (=humidity ratio)
    !
    meteo(nm)%Tv(:,:) = meteo(nm)%T(:,:)*(1.0_rp + 1.6077_rp*work_WRF2(:,:))/(1.0_rp+work_WRF2(:,:))
    !
    !*** ro(nx,ny,nz) density. Gas law using virtual temperature (account for water in air)
    !
    meteo(nm)%ro(:,:) = meteo(nm)%p(:,:)/(287.06*meteo(nm)%Tv(:,:))
    !
    deallocate(work_WRF2)
    !
    !*** close file
    !
    if( nf90_close(ncID) /= NF90_NOERR ) &
         call runend('readwrf : error closing netcdf file')
    !
    !*** Interpolate to the background mesh
    !
    do iz = 1, grid%nz
       ipoin = 0
       do iy = 1, grid%ny
          do ix = 1, grid%nx
             ipoin = ipoin+1
             !
             if(grid%model(ipoin).eq.nm)then
                !
                !*** Store previous values
                ! 
                grid%u(ix,iy,iz,1) = grid%u(ix,iy,iz,2)
                grid%u(ix,iy,iz,2) = 0.0
                !
                grid%v(ix,iy,iz,1) = grid%v(ix,iy,iz,2)
                grid%v(ix,iy,iz,2) = 0.0
                !
                grid%w(ix,iy,iz,1) = grid%w(ix,iy,iz,2)
                grid%w(ix,iy,iz,2) = 0.0
                !
                grid%T(ix,iy,iz,1) = grid%T(ix,iy,iz,2)
                grid%T(ix,iy,iz,2) = 0.0
                !
                grid%ro(ix,iy,iz,1) = grid%ro(ix,iy,iz,2)
                grid%ro(ix,iy,iz,2) = 0.0
                !
                grid%drhodz(ix,iy,iz,1) = grid%drhodz(ix,iy,iz,2)
                grid%drhodz(ix,iy,iz,2) = 0.0
                !
                grid%qv(ix,iy,iz,1) = grid%qv(ix,iy,iz,2)
                grid%qv(ix,iy,iz,2) = 0.0
                !
                grid%tp(ix,iy,iz,1) = grid%tp(ix,iy,iz,2)
                grid%tp(ix,iy,iz,2) = 0.0
                ! 
                !*** Interpolate IN HEIGHT and respect the box(ielem) position of the particle.
                ! 
                ielem = grid%element(ipoin)   ! meteo model element (2d) 
                do i = 1,4
                   jpoin = meteo(nm)%lnods(i,ielem)
                   jz = 0
                   found = .false.
                   do while(.not.found)
                      jz = jz + 1
                      ! 
                      if (grid%z(iz).lt.meteo(nm)%z(jpoin,1  ) ) then
                         kz = 1
                         zshape = 1.0
                         found  = .true.
                      else if((grid%z(iz).ge.meteo(nm)%z(jpoin,jz  )).and. & 
                           (grid%z(iz).lt.meteo(nm)%z(jpoin,jz+1)) ) then
                         kz = jz
                         zshape = (meteo(nm)%z(jpoin,kz+1)-grid%z(iz))/(meteo(nm)%z(jpoin,kz+1)-meteo(nm)%z(jpoin,kz))
                         found  = .true.
                      else if(jz.eq.(meteo(nm)%np-1)) then
                         kz = jz
                         zshape = 0.0
                         found = .true.     ! top not found. Use the highest value
                      end if
                      ! 
                   end do
                   !       
                   ! Interpolate U
                   ! 
                   grid%u(ix,iy,iz,2) = grid%u(ix,iy,iz,2) + & 
                        zshape *(grid%shape(i,ipoin)*meteo(nm)%u(jpoin,kz))+&
                        (1-zshape)*(grid%shape(i,ipoin)*meteo(nm)%u(jpoin,kz+1))
                   !       
                   ! Interpolate V
                   ! 
                   grid%v(ix,iy,iz,2) = grid%v(ix,iy,iz,2) + & 
                        zshape *(grid%shape(i,ipoin)*meteo(nm)%v(jpoin,kz))+&
                        (1-zshape)*(grid%shape(i,ipoin)*meteo(nm)%v(jpoin,kz+1))
                   !       
                   ! Interpolate W
                   ! 
                   grid%w(ix,iy,iz,2) = grid%w(ix,iy,iz,2) + & 
                        zshape *(grid%shape(i,ipoin)*meteo(nm)%w(jpoin,kz))+&
                        (1-zshape)*(grid%shape(i,ipoin)*meteo(nm)%w(jpoin,kz+1))
                   !       
                   ! Interpolate T
                   ! 
                   grid%T(ix,iy,iz,2) = grid%T(ix,iy,iz,2) + & 
                        zshape *(grid%shape(i,ipoin)*meteo(nm)%T(jpoin,kz))+&
                        (1-zshape)*(grid%shape(i,ipoin)*meteo(nm)%T(jpoin,kz+1))
                   !       
                   ! Interpolate RO
                   ! 
                   grid%ro(ix,iy,iz,2) = grid%ro(ix,iy,iz,2) + & 
                        zshape *(grid%shape(i,ipoin)*meteo(nm)%ro(jpoin,kz))+&
                        (1-zshape)*(grid%shape(i,ipoin)*meteo(nm)%ro(jpoin,kz+1))
                   !       
                   ! Interpolate QV
                   ! 
                   grid%qv(ix,iy,iz,2) = grid%qv(ix,iy,iz,2) + & 
                        zshape *(grid%shape(i,ipoin)*meteo(nm)%qv(jpoin,kz))+&
                        (1-zshape)*(grid%shape(i,ipoin)*meteo(nm)%qv(jpoin,kz+1))
                   !       
                   ! Interpolate QV
                   ! 
                   grid%tp(ix,iy,iz,2) = grid%tp(ix,iy,iz,2) + & 
                        zshape *(grid%shape(i,ipoin)*meteo(nm)%tp(jpoin,kz))+&
                        (1-zshape)*(grid%shape(i,ipoin)*meteo(nm)%tp(jpoin,kz+1))
                   ! 
                end do  ! i = 1,4
                ! 
                !Change the sign in backward case
                if(sidt.lt.0)then
                   grid%u (ix,iy,iz,2)=grid%u (ix,iy,iz,2)*(-1)
                   grid%v (ix,iy,iz,2)=grid%v (ix,iy,iz,2)*(-1)
                end if
             end if
          end do
       end do
    end do
    ! 
    !***interpolate in the grid for 2D variables.
    ! 
    ipoin = 0
    do iy = 1, grid%ny
       do ix = 1, grid%nx
         ipoin = ipoin+1
             if(grid%model(ipoin).eq.nm)then
                !
                !*** Store previous values
                ! 
                grid%pblh(ix,iy,1) = grid%pblh(ix,iy,2)
                grid%pblh(ix,iy,2) = 0.0
                !
                grid%hfx(ix,iy,1) = grid%hfx(ix,iy,2)
                grid%hfx(ix,iy,2) = 0.0
                !
                grid%ust(ix,iy,1) = grid%ust(ix,iy,2)
                grid%ust(ix,iy,2) = 0.0
                !
                grid%rmonin(ix,iy,1) = grid%rmonin(ix,iy,2)
                grid%rmonin(ix,iy,2) = 0.0
                !
                grid%wst(ix,iy,1) = grid%wst(ix,iy,2)
                grid%wst(ix,iy,2) = 0.0
                ! 
                !*** Interpolate respect the box(ielem) position of the particle.
                ! 
                ielem = grid%element(ipoin)   ! meteo model element (2d) 
                do i = 1,4
                   jpoin = meteo(nm)%lnods(i,ielem)
                   !       
                   ! Interpolate PBLH
                   ! 
                    grid%pblh(ix,iy,2) = grid%pblh(ix,iy,2) + & 
                                         grid%shape(i,ipoin)*meteo(nm)%pblh(jpoin)
                   !       
                   ! Interpolate HFX
                   ! 
                    grid%hfx(ix,iy,2) = grid%hfx(ix,iy,2) + & 
                                        grid%shape(i,ipoin)*meteo(nm)%hfx(jpoin)
                   !       
                   ! Interpolate U*
                   ! 
                    grid%ust(ix,iy,2) = grid%ust(ix,iy,2) + & 
                                        grid%shape(i,ipoin)*meteo(nm)%ust(jpoin)
                end do
              end if
              ! 
       end do
    end do
    !
    !***  Computes other variables not given by WRF
    !

    do iy = 1,grid%ny
       do ix = 1,grid%nx
          !
          !***  rmonin(nx,ny)
          !
          grid%rmonin(ix,iy,2) =  - (grid%ust(ix,iy,2)*grid%ust(ix,iy,2)*grid%ust(ix,iy,2)*grid%tp(ix,iy,1,2)&
               *grid%ro(ix,iy,1,2)*1006.0)/(0.4*9.81*grid%hfx(ix,iy,2))
          grid%rmonin(ix,iy,2) = min(grid%rmonin(ix,iy,2), 1e3_rp)
          grid%rmonin(ix,iy,2) = max(grid%rmonin(ix,iy,2),-1e3_rp)
          !
          !*** Calculate the convective velocity scale wst
          !
          if((grid%hfx(ix,iy,2)).lt.0.) then       !if surface heat flux is less than 0.
                   grid%wst(ix,iy,2)=-grid%pblh(ix,iy,2)*g/grid%tp(ix,iy,1,2)*grid%hfx(ix,iy,2)/cpa
                   if(grid%wst(ix,iy,2).ge.0)then
                      grid%wst(ix,iy,2)=(grid%wst(ix,iy,2))**0.333 
                   else
                      grid%wst(ix,iy,2)=-(-grid%wst(ix,iy,2))**0.333 
                   end if
              ! grid%wst(ix,iy,2)=(-grid%pblh(ix,iy,2)*g/grid%tp(ix,iy,1,2)*grid%hfx(ix,iy,2)/cpa)**0.333   
              ! g:gravity, cpa:specific heat for dry air
              ! excess=-bs*grid%hfx/cpa/wst
              ! if (iter.lt.itmax) goto 30
          else
             grid%wst(ix,iy,2)=0.
          endif
          !
          !*** Calculates drhodz: air density vertical gradient
          !
          grid%drhodz(ix,iy,1,2)=(grid%ro(ix,iy,2,2)-grid%ro(ix,iy,1,2))/(grid%z(2)-grid%z(1))
          do iz=2,grid%nz-1
            grid%drhodz(ix,iy,iz,2)=(grid%ro(ix,iy,iz+1,2)-grid%ro(ix,iy,iz-1,2))/ &
            (grid%z(iz+1)-grid%z(iz-1))
          end do
          grid%drhodz(ix,iy,grid%nz,2)=grid%drhodz(ix,iy,grid%nz-1,2)
       end do
    end do

    return
  end subroutine readwrf
  !
end MODULE Wrf
