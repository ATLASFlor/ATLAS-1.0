!***************************************************************
!*
!*              Module for GFS Meteo
!*
!***************************************************************
MODULE GFS
   use KindType
   use Master
   use InpOut
   IMPLICIT NONE
   SAVE

   !
   !*** GFS DICTIONARY
   !
   ! GFS DIMENSIONS:
   character(len=25) :: nx_gfs_name = 'lon'
   character(len=25) :: ny_gfs_name = 'lat'
   character(len=25) :: np_gfs_name = 'pres'
   character(len=25) :: nt_gfs_name = 'time'
   !
   !GFS VARIABLES:
   character(len=25) :: time_gfs_name = 'time'
   character(len=25) :: lon_gfs_name = 'lon'
   character(len=25) :: lat_gfs_name = 'lat'
   character(len=25) :: hgt_gfs_name = 'HGT:surface'
   character(len=25) :: landmask_gfs_name = 'LAND'
   character(len=25) :: fi_GFS_name = 'HGT'
   character(len=25) :: u_gfs_name = 'UGRD'
   character(len=25) :: v_gfs_name = 'VGRD'
   character(len=25) :: w_gfs_name = 'VVEL'
   character(len=25) :: t_gfs_name = 'TMP'
   character(len=25) :: q_gfs_name = 'RH'
   character(len=25) :: p_gfs_name = 'pres'
   character(len=25) :: pblh_gfs_name = 'HPBL'
   !
   !GFS ATTRIBUTES:
   character(len=25) :: year_GFS_name = 'YEAR'
   character(len=25) :: month_GFS_name = 'MONTH'
   character(len=25) :: day_GFS_name = 'DAY'
   character(len=25) :: hour_GFS_name = 'HOUR'
   character(len=25) :: Tincr_GFS_name = 'TIME_INCR'
   character(len=25) :: lonmin_GFS_name = 'LONMIN'
   character(len=25) :: lonmax_GFS_name = 'LONMAX'
   character(len=25) :: latmin_GFS_name = 'LATMIN'
   character(len=25) :: latmax_GFS_name = 'LATMAX'

CONTAINS
   !
   !
   !

   subroutine opengfs(nm)
      !******************************************************************************************
      !*    OPENGFS
      !*    Opens the Meteorological file to create the meteo mesh
      !*
      !*    nm       Correspond to model number readed by input file  (Input value)
      !*
      !*    This subroutine reads the netcdf file converted form grib by the Grib2nc program.
      !******************************************************************************************
      use Master
      use KindType
      use netcdf
      use Elsest_mod
      use TimeFun
      implicit none
      !
      logical               :: found
      !
      integer(ip)         :: nm
      integer(ip)         :: ielem, ipoin, ix, iy, it
      !
      integer(ip)         :: ncid                   ! Id Number of the meteo mesh
      integer(ip)         :: NXdimid                ! NX Dimension Id.
      integer(ip)         :: NYdimid                ! NY Dimension Id.
      integer(ip)         :: NPdimid                ! NP Dimension Id (pressure levels)
      integer(ip)         :: Ntdimid                ! Dt Dimension Id.
      integer(ip)         :: lonvarid               ! Dt Dimension Id.
      integer(ip)         :: latvarid               ! Dt Dimension Id.
      integer(ip)         :: hgtvarid               ! hgt Variable Id.
      integer(ip)         :: landmaskvarid          ! Landmask Variable Id.
      integer(ip)         :: tvarid                 ! Times Variable Id.
      real(rp)         :: tincr                  ! Increment Time.in seconds

      !
      real(rp)         :: tduration              ! Total Time Duration.in seconds
      real(rp)         :: lonmin_gfs, lonmax_gfs ! Longitud max and min. Attribute
      real(rp)         :: latmin_gfs, latmax_gfs ! Latitud max and min. Attribute
      real(rp)         :: dy_gfs                 ! Y resolution
      !
      !
      character(len=s_mess) :: NXname, NYname, NPname   ! Names of dimension NX, NY and PZ
      character(len=s_mess) :: Ntname                 ! Name of time step dimension.
      character(len=s_mess) :: str
      character(len=1) :: ext                    !

      !
      real(rp), allocatable :: GFS_coord(:, :)
      real(rp), allocatable :: work2d(:, :)            ! Local array used for topography
      real(rp), allocatable :: lwork2d(:, :)           ! Local array used for landmask
      real(rp), allocatable :: lon_GFS(:)
      real(rp), allocatable :: lat_GFS(:)
      !
      !*** GFS Time
      !
      integer(ip)   :: yr0_gfs, mo0_gfs, dy0_gfs, hr0_gfs
      integer(ip)   :: iyr, imo, idy, ihr, imi, ise
      !
      !
      !    Elsest variables
      !
      real(rp) :: xcoor(3), coloc(3)
      real(rp) :: shapf(4), deriv(2, 4), relse(20)
      integer(ip) :: ielse(10), ltopo, ifoun, i
      !
      !*** Elsest initializations
      !
      ielse(1) = 10                                     ! number nx boxes (search)
      ielse(2) = 10                                     ! ny
      ielse(3) = 10                                     ! nz
      ielse(4) = 1                                      ! data format (0=type,1=linked list)
      ielse(5) = 10                                     ! Maximum number of possible meshes
      ielse(6) = 1                                      ! Not used
      ielse(7) = 0                                      ! Output unit
      ielse(8) = 0                                      ! Search strategy (0=bin,1=Quad)
      !
      relse(1) = 0.01_rp                                ! Tolerance for iteration
      relse(3) = -1e8_rp                                 ! Patch bounding box
      relse(4) = 1e8_rp
      relse(5) = -1e8_rp
      relse(6) = 1e8_rp
      relse(7) = -1e8_rp
      relse(8) = 1e8_rp
      !
      ltopo = 0 ! (0 quads hexahedra , 1 for triangle tetra, 2 Pentahedra, 3 Pyramid)

      !
      !*** Open the GFS model. Gets the GFS file ID
      !
      if (NF90_OPEN(meteo(nm)%file_path, NF90_NOWRITE, ncID) /= NF90_NOERR) &
         call runend('Error in  nf90_open gfs file')
      !
      write (ext, '(i1)') nm
      str = 'METEO_DATA_'//ext
      !
      !*** Gets the Id number for the Nx Dimension. and read it
      !
      if (nf90_inq_dimid(ncID, nx_gfs_name, NXdimid) /= NF90_NOERR) &
         call runend('Error inquiring nx id dimension in gfs file')
      if (nf90_inquire_dimension(ncID, NXdimid, NXname, meteo(nm)%nx) /= NF90_NOERR) &
         call runend('Error inquiring nx dimension in gfs file')
      !
      !*** Gets the Id number for the Ny Dimension.
      !
      if (nf90_inq_dimid(ncID, ny_gfs_name, NYdimid) /= NF90_NOERR) &
         call runend('Error inquiring ny id dimension in gfs file')
      if (nf90_inquire_dimension(ncID, NYdimid, NYname, meteo(nm)%ny) /= NF90_NOERR) &
         call runend('Error inquiring NY dimension in gfs file')
      !
      !*** Gets the Id number for the number of pressure levels Dimension and read it. np
      !
      if (nf90_inq_dimid(ncID, np_gfs_name, NPdimid) /= NF90_NOERR) &
         call runend('Error inquiring np id dimension in gfs file')
      if (nf90_inquire_dimension(ncID, NPdimid, NPname, meteo(nm)%np) /= NF90_NOERR) &
         call runend('Error inquiring NP dimension in gfs file')
      !
      !*** Gets the ID number for Time dimension. nt
      !
      if (nf90_inq_dimid(ncID, nt_gfs_name, Ntdimid) /= NF90_NOERR) &
         call runend('Error inquiring nt id dimension in gfs file')
      if (nf90_inquire_dimension(ncID, Ntdimid, Ntname, meteo(nm)%nt) /= NF90_NOERR) &
         call runend('Error inquiring NT dimension in gfs file')
      !
      !*** allocate memory
      !
      allocate (meteo(nm)%timesec(meteo(nm)%nt))
      !
      !*** Allocate memory for local arrays (need to be deallocated later)
      !
      allocate (lon_gfs(meteo(nm)%nx))
      allocate (lat_gfs(meteo(nm)%ny))
      allocate (work2d(meteo(nm)%nx, meteo(nm)%ny))
      allocate (lwork2d(meteo(nm)%nx, meteo(nm)%ny))
      !
      !*** Get the day, month and year. initial and final (ATTRIBUTES)
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, year_GFS_name, yr0_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute YEAR in opengfs')
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, month_GFS_name, mo0_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute MONTH in opengfs')
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, day_GFS_name, dy0_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute DAY in opengfs')
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, hour_GFS_name, hr0_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute DAY in opengfs')
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, Tincr_GFS_name, tincr) /= NF90_NOERR) &
         call runend('Error getting attribute DAY in opengfs')
      !
      !*** Get lon lat min and max. attributes
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, lonmax_GFS_name, lonmax_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute lonmax in opengfs')
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, lonmin_GFS_name, lonmin_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute lonmin in opengfs')
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, latmax_GFS_name, latmax_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute latmax in opengfs')
      !
      if (nf90_get_att(ncID, NF90_GLOBAL, latmin_GFS_name, latmin_gfs) /= NF90_NOERR) &
         call runend('Error getting attribute latmin in opengfs')
      !
      !*** Get time variable
      !
      if (nf90_inq_varid(ncID, time_gfs_name, tvarID) /= NF90_NOERR) &
         call runend('Error getting Time variable id in opengfs')
      if (nf90_get_var(ncID, tvarID, meteo(nm)%timesec) /= NF90_NOERR) &
         call runend('Error getting variable time in opengfs')
      !
      !*** Calculate the day, month, year, hour, minute for each time step in format yyyymmddhhmm
      !
      allocate (meteo(nm)%time(meteo(nm)%nt))
      do it = 1, meteo(nm)%nt
         call addtime(yr0_gfs, mo0_gfs, dy0_gfs, hr0_gfs, iyr, imo, idy, ihr, imi, ise, meteo(nm)%timesec(it))
         meteo(nm)%time(it) = 1d8*iyr + 1d6*imo + 1d4*idy + 1d2*ihr + imi
      end do
      !
      !***  Calculates time lag by iteration (that is, the time in seconds between the GFS origin
      !***  and the DBS origin). Both origins are referred to 0000UTC, but may belong to different
      !***  days
      !
      found = .false.
      meteo(nm)%time_lag = 0.0_rp
      do while (.not. found)
         call addtime(yr0_GFS, mo0_GFS, dy0_GFS, 0, &
                      iyr, imo, idy, ihr, imi, ise, meteo(nm)%time_lag)
         !
         if ((iyr .eq. ibyr) .and. (imo .eq. ibmo) .and. (idy .eq. ibdy)) then
            found = .true.
         else
            meteo(nm)%time_lag = meteo(nm)%time_lag + 86400.0_rp
         end if
         if (meteo(nm)%time_lag .gt. meteo(nm)%timesec(meteo(nm)%nt)) &
            call runend('Time lag not found. Check interval consistency'//TRIM(str))
      end do
      !
      !*** Checks that the required time is within the bounds
      !
      if (sidt .gt. 0) then
         if ((meteo(nm)%time_lag + sist) .lt. meteo(nm)%timesec(1)) &
            call runend('Minimum datatime is after required time'//TRIM(str))
         if ((meteo(nm)%time_lag + siend) .gt. meteo(nm)%timesec(meteo(nm)%nt)) &
            call runend('Maximum datatime is before required time'//TRIM(str))
      else if (sidt .lt. 0) then
         if ((meteo(nm)%time_lag + siend) .lt. meteo(nm)%timesec(1)) &
            call runend('Minimum datatime is after  required time'//TRIM(str))
         if ((meteo(nm)%time_lag + sist) .gt. meteo(nm)%timesec(meteo(nm)%nt)) &
            call runend('Maximum datatime is before required time'//TRIM(str))
      end if
      !
      !*** Redifine the timesec adding the time_lag.
      !*** So, the meteo file and the simulation have the same origin.
      !
      meteo(nm)%timesec(:) = meteo(nm)%timesec(:) - meteo(nm)%time_lag
      !
      !*** Get GFS Lon and Lat variables
      !
      if (nf90_inq_varid(ncID, lon_gfs_name, lonvarID) /= NF90_NOERR) &
         call runend('Error getting lon variable id in opengfs')
      if (nf90_get_var(ncID, lonvarID, lon_GFS) /= NF90_NOERR) &
         call runend('Error getting variable lon in opengfs')
      !
      if (nf90_inq_varid(ncID, lat_gfs_name, latvarID) /= NF90_NOERR) &
         call runend('Error getting lat variable id in opengfs')
      if (nf90_get_var(ncID, latvarID, lat_GFS) /= NF90_NOERR) &
         call runend('Error getting variable lat in opengfs')
      !
      !*** Get gfs topography at box center 'HGT:surface' variable in GFS
      !
      if (nf90_inq_varid(ncID, hgt_gfs_name, hgtvarID) /= NF90_NOERR) &
         call runend('Error getting hgt variable id in opengfs')
      if (nf90_get_var(ncID, hgtvarID, work2d) /= NF90_NOERR) &
         call runend('Error getting variable hgt in opengfs')
      !
      !*** Get gfs landmask at box center
      !
      if (nf90_inq_varid(ncID, landmask_gfs_name, landmaskvarID) /= NF90_NOERR) &
         call runend('Error getting landmask variable id in opengfs')
      if (nf90_get_var(ncID, landmaskvarID, lwork2d) /= NF90_NOERR) &
         call runend('Error getting variable landmask in opengfs')
      work2d = work2d*lwork2d
      !
      !*** Calculate the model resolution. In degrees
      !
      meteo(nm)%res = (latmax_GFS - latmin_GFS)/(meteo(nm)%ny - 1)
      !
      !*** Determine if the meteorological model is global
      !
      if (lon_GFS(meteo(nm)%nx) + meteo(nm)%res .eq. 180 .and. lon_GFS(1) .eq. -180) meteo(nm)%global = .true.
      !
      !*** Allocate memory and save lon/lat dimension and variables according the global condition
      !
      !
      !*** No Global model case
      !
      if (meteo(nm)%global .eqv. .false.) then
         meteo(nm)%npoin2d = meteo(nm)%nx*meteo(nm)%ny
         allocate (meteo(nm)%hgt(meteo(nm)%npoin2d))
         allocate (meteo(nm)%landmask(meteo(nm)%npoin2d))
         ! Local arrays
         allocate (GFS_coord(2, meteo(nm)%npoin2d))
         !
         allocate (meteo(nm)%lon(meteo(nm)%nx, meteo(nm)%ny))
         allocate (meteo(nm)%lat(meteo(nm)%nx, meteo(nm)%ny))
         meteo(nm)%nelem2d = (meteo(nm)%nx - 1)*(meteo(nm)%ny - 1)
         allocate (meteo(nm)%lnods(meteo(nm)%nnode, meteo(nm)%nelem2d))
         !
         allocate (meteo(nm)%z(meteo(nm)%npoin2d, meteo(nm)%np))       ! stored in 2D planes to interpolate with Elsest
         allocate (meteo(nm)%u(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%v(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%w(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%T(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%p(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%qv(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%tp(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%Tv(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%ro(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%pblh(meteo(nm)%npoin2d))
         ! Create meteo(nm)%lon(nx,ny) and meteo(nm)%lat(nx,ny)
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               meteo(nm)%lon(ix, iy) = lon_GFS(ix)
               meteo(nm)%lat(ix, iy) = lat_GFS(iy)
            end do
         end do
         ! Save Time independent Variables. Topography, Landmask
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%hgt(ipoin) = work2d(ix, iy)
               meteo(nm)%landmask(ipoin) = lwork2d(ix, iy)
            end do
         end do
         ! Builds mesh structure. 2d
         ! Associate each element with 4 nods
         ielem = 0
         do iy = 1, meteo(nm)%ny - 1
            do ix = 1, meteo(nm)%nx - 1
               ielem = ielem + 1
               meteo(nm)%lnods(1, ielem) = (iy - 1)*meteo(nm)%nx + ix
               meteo(nm)%lnods(2, ielem) = (iy - 1)*meteo(nm)%nx + ix + 1
               meteo(nm)%lnods(3, ielem) = (iy)*meteo(nm)%nx + ix + 1
               meteo(nm)%lnods(4, ielem) = (iy)*meteo(nm)%nx + ix
            end do
         end do
         !
         ! Associate each point with their coordinates. 1: lon. 2: lat.
         !
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               GFS_coord(1, ipoin) = meteo(nm)%lon(ix, iy)
               GFS_coord(2, ipoin) = meteo(nm)%lat(ix, iy)
            end do
         end do
         !
         !
      else if (meteo(nm)%global .eqv. .true.) then
         !
         !*** Global model case
         !
         meteo(nm)%npoin2d = (meteo(nm)%nx + 1)*meteo(nm)%ny
         allocate (meteo(nm)%hgt(meteo(nm)%npoin2d))
         allocate (meteo(nm)%landmask(meteo(nm)%npoin2d))
         ! Local arrays
         allocate (GFS_coord(2, meteo(nm)%npoin2d))
         !
         allocate (meteo(nm)%lon(meteo(nm)%nx + 1, meteo(nm)%ny))
         allocate (meteo(nm)%lat(meteo(nm)%nx + 1, meteo(nm)%ny))
         meteo(nm)%nelem2d = (meteo(nm)%nx)*(meteo(nm)%ny - 1)
         allocate (meteo(nm)%lnods(meteo(nm)%nnode, meteo(nm)%nelem2d))
         !
         allocate (meteo(nm)%z(meteo(nm)%npoin2d, meteo(nm)%np))       ! stored in 2D planes to interpolate with Elsest
         allocate (meteo(nm)%u(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%v(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%w(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%T(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%tp(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%p(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%qv(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%Tv(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%ro(meteo(nm)%npoin2d, meteo(nm)%np))
         allocate (meteo(nm)%pblh(meteo(nm)%npoin2d))
         ! Create meteo(nm)%lon(nx,ny) and meteo(nm)%lat(nx,ny)
         do iy = 1, meteo(nm)%ny
            meteo(nm)%lon(meteo(nm)%nx + 1, iy) = 180
            meteo(nm)%lat(meteo(nm)%nx + 1, iy) = lat_GFS(iy)
            do ix = 1, meteo(nm)%nx
               meteo(nm)%lon(ix, iy) = lon_GFS(ix)
               meteo(nm)%lat(ix, iy) = lat_GFS(iy)
            end do
         end do
         ! Save Time independent Variables: Topography, Landmask
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%hgt(ipoin) = work2d(ix, iy)
               meteo(nm)%landmask(ipoin) = lwork2d(ix, iy)
            end do
            ipoin = ipoin + 1
            meteo(nm)%hgt(ipoin) = work2d(1, iy)
            meteo(nm)%landmask(ipoin) = lwork2d(1, iy)
         end do
         ! Builds mesh structure. 2d
         ! Associate each element with 4 nods
         ielem = 0
         do iy = 1, meteo(nm)%ny - 1
            do ix = 1, meteo(nm)%nx
               ielem = ielem + 1
               meteo(nm)%lnods(1, ielem) = (iy - 1)*(meteo(nm)%nx + 1) + ix
               meteo(nm)%lnods(2, ielem) = (iy - 1)*(meteo(nm)%nx + 1) + ix + 1
               meteo(nm)%lnods(3, ielem) = (iy)*(meteo(nm)%nx + 1) + ix + 1
               meteo(nm)%lnods(4, ielem) = (iy)*(meteo(nm)%nx + 1) + ix
            end do
         end do
         !
         ! Associate each point with their coordinates. 1: lon. 2: lat.
         !
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx + 1
               ipoin = ipoin + 1
               GFS_coord(1, ipoin) = meteo(nm)%lon(ix, iy)
               GFS_coord(2, ipoin) = meteo(nm)%lat(ix, iy)
            end do
         end do
         !
      end if
      !
      !*** Computes the interpolation functions in 2D for all points of the background mesh
      !*** to interpolate GFS variables (later on)
      !
      do ipoin = 1, grid%npoin2d   ! loop over background mesh points
         !
         xcoor(1) = grid%coord(1, ipoin)
         xcoor(2) = grid%coord(2, ipoin)
         !
         ! elses 2_ search element.
         !  returns:
         !    ielem > 0  element number
         !    ielem < 0  is not in the domain
         !
         call elsest( &
            2_ip, ielse, meteo(nm)%nnode, 2_ip, meteo(nm)%npoin2d, meteo(nm)%nelem2d, &
            meteo(nm)%lnods, ltopo, GFS_coord, xcoor, relse, &
            ielem, shapf, deriv, coloc)
         !
         ! Restrict the shape factor to [0,1]
         !
         do i = 1, 4
            shapf(i) = min(shapf(i), 1.0)
            shapf(i) = max(shapf(i), 0.0)
         end do
         !
         !*** Put the flag 'nm' in the background grid elements which are in this meteo domain and have better resolution than others.
         !
         if (ielem .gt. 0) then
            if (grid%model(ipoin) .eq. -1) then
               grid%model(ipoin) = nm
               grid%element(ipoin) = ielem
               grid%shape(1:4, ipoin) = shapf(1:4)
            else if (meteo(grid%model(ipoin))%res .gt. meteo(nm)%res) then
               grid%model(ipoin) = nm
               grid%element(ipoin) = ielem
               grid%shape(1:4, ipoin) = shapf(1:4)

            end if
         end if
         !
      end do
      !
      !*** Deallocates
      !
      deallocate (GFS_coord)
      deallocate (work2d)
      deallocate (lwork2d)
      deallocate (lon_gfs)
      deallocate (lat_gfs)
      !
      !*** Release elsest memory
      !
      call elsest_deallo()
      !
      !
      !*** Print Information in log file
      !
      !***  Opens and writes the log file
      !
      if (my_id .eq. 0) write (lulog, 2) 'Opened GFS file and readed necesary Time Independent variables. Created mesh structure.'
      if (out_screen) write (*, 2) 'Opened GFS file and readed necesary Time Independent variables. Created mesh structure.'
      !
2     format(/, a/)
      !
      return
      !
   end subroutine opengfs
   !
   !
   !
   subroutine readgfs(nm, itime)
      !******************************************************************************************
      !*    READGFS
      !*    Opens the Meteorological file and read variables for two times
      !*
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
      integer(ip)            :: nm, itime
      character(len=s_mess)    :: str
      character(len=1)    :: ext
      integer(ip)            :: ipoin, iz, ix, iy, i, ielem
      integer(ip)            :: jpoin, jz, kz
      logical                  :: found
      real(rp)            :: zshape
      !
      !Variables ID
      integer(ip)        :: ncID         ! Netcdf file ID
      integer(ip)        :: fiID         ! Netcdf file ID
      integer(ip)        :: uID          ! Horizontal winds ID
      integer(ip)        :: vID          ! Vertical winds ID
      integer(ip)        :: wID          ! Vertical omega winds ID
      integer(ip)        :: TID          ! Temperature
      integer(ip)        :: qID          ! Specific humidity
      integer(ip)        :: pID          ! Pressure
      integer(ip)        :: pblhid       ! planetary boundary layer heights id
      integer(ip)        :: hfxid        ! Surface heat flux id
      integer(ip)        :: ustid        ! U* id
      !
      real(rp), allocatable :: work_GFS1(:, :, :) !work_GFS1(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np)
      real(rp), allocatable :: work_GFS2(:, :, :)   !work_GFS2(meteo(nm)%nx,meteo(nm)%ny)
      real(rp), allocatable :: work_GFS(:, :) !work_GFS (npoind2d,meteo(nm)%np)
      real(rp), allocatable :: p_GFS(:)         !p(meteo(nm)%np)    pressure
      real(rp), allocatable :: p(:, :, :)         !p pressure(nx,ny,nz) levels
      !
      !*** Initializations
      !
      write (ext, '(i1)') nm
      str = 'METEO_DATA_'//ext
      !
      !*** Open the GFS model. Gets the GFS file ID
      !
      if (NF90_OPEN(meteo(nm)%file_path, NF90_NOWRITE, ncID) /= NF90_NOERR) &
         call runend('Error in  nf90_open gfs file')
      !
      !*** Gets Variables ID
      !
      if (nf90_inq_varid(ncID, fi_GFS_name, fiID) /= NF90_NOERR) &
         call runend('Error getting fiID in '//TRIM(str))
      !
      if (nf90_inq_varid(ncID, u_GFS_name, uID) /= NF90_NOERR) &
         call runend('Error getting uID in '//TRIM(str))
      !
      if (nf90_inq_varid(ncID, v_GFS_name, vID) /= NF90_NOERR) &
         call runend('Error getting vID in '//TRIM(str))
      !
      if (nf90_inq_varid(ncID, w_GFS_name, wID) /= NF90_NOERR) &
         call runend('Error getting wID in '//TRIM(str))
      !
      allocate (work_GFS1(meteo(nm)%nx, meteo(nm)%ny, meteo(nm)%np))
      allocate (work_GFS2(meteo(nm)%nx, meteo(nm)%ny, 1))
      allocate (work_GFS(meteo(nm)%npoin2d, meteo(nm)%np))
      !
      !*** Heights. This is necessary because pressure levels change with time.
      !*** GFS already gives the geopotential in m (i.e. height). The topography is substracted in order to interpolate
      !*** in terrain following
      !
      if (nf90_get_var(ncID, fiID, work_GFS1, start=(/1, 1, 1, itime/)) /= NF90_NOERR) &
         call runend('Error reading fi in '//TRIM(str))
      !
      do iz = 1, meteo(nm)%np
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%z(ipoin, iz) = work_GFS1(ix, iy, iz) - meteo(nm)%hgt(ipoin)
            end do
         end do
      end do
      if (MAXVAL(meteo(nm)%z(:, :)) .lt. (grid%ztop)) &
         call wriwar('Input max vertical layer is higher than meteo max layer. Values extrapolated in '//TRIM(str))
      !
      !*** Read U Variable, U(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np)
      !
      if (nf90_get_var(ncID, uID, work_GFS1, start=(/1, 1, 1, itime/)) /= NF90_NOERR) &
         call runend('read_GFS: Error reading u in '//TRIM(str))
      !
      do iz = 1, meteo(nm)%np
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%u(ipoin, iz) = work_GFS1(ix, iy, iz)
            end do
         end do
      end do
      !
      !*** Read V Variable, V(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np)
      !
      if (nf90_get_var(ncID, vID, work_GFS1, start=(/1, 1, 1, itime/)) /= NF90_NOERR) &
         call runend('read_GFS: Error reading v in '//TRIM(str))
      !
      do iz = 1, meteo(nm)%np
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%v(ipoin, iz) = work_GFS1(ix, iy, iz)
            end do
         end do
      end do
      !
      !*** Reads p levels for GFS mesh
      !
      allocate (p_GFS(meteo(nm)%np))
      allocate (p(meteo(nm)%nx, meteo(nm)%ny, meteo(nm)%np))
      if (nf90_inq_varid(ncID, p_GFS_name, pID) /= 0) call runend('read_GFS : Error getting pID')
      if (nf90_get_var(ncID, pid, p_GFS) /= 0) call runend('read_GFS : Error in nf90_get_var')
      p_GFS(:) = p_GFS(:)*101.3_rp  ! (converted from mb to Pa)
      !
      !*** Read W Variable, W(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np)
      !
      if (nf90_get_var(ncID, wID, work_GFS1, start=(/1, 1, 1, itime/)) /= NF90_NOERR) &
         call runend('read_GFS: Error reading w in '//TRIM(str))
      !
      do iz = 1, meteo(nm)%np
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%w(ipoin, iz) = work_GFS1(ix, iy, iz)!/p_GFS(iz)
            end do
         end do
      end do
      !
      !*** meteo(nm)%p(meteo(nm)%npoin2d,meteo(nm)%np) pressure
      !
      do iz = 1, meteo(nm)%np
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%p(ipoin, iz) = p_GFS(iz)         ! P
               p(ix, iy, iz) = p_GFS(iz)
            end do
         end do
      end do
      !
      !*** TMP(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np) --> T(meteo(nm)%nx,meteo(nm)%ny,meteo(nm)%np)  Temperature
      !
      if (nf90_inq_varid(ncID, T_GFS_name, tID) /= 0) call runend('read_GFS : Error getting tID')
      if (nf90_get_var(ncID, tID, work_GFS1, start=(/1, 1, 1, itime/)) /= 0) call runend('read_GFS: Error reading T')
      !
      do iz = 1, meteo(nm)%np
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%T(ipoin, iz) = work_GFS1(ix, iy, iz)
               meteo(nm)%tp(ipoin, iz) = work_GFS1(ix, iy, iz)*((1.01e5_rp/p(ix, iy, iz))**2)**(1/7)!(0.285_rp)
            end do
         end do
      end do
      !
      !*** RH(nx_GFS,ny_GFS,nz_GFS) Relative humidity --> qv(nx,ny,nz)  Specific humidity is computed later
      !
      if (nf90_inq_varid(ncID, Q_GFS_name, qID) /= 0) call runend('read_GFS : Error getting qID')
      if (nf90_get_var(ncID, qID, work_GFS1, start=(/1, 1, 1, itime/)) /= 0) call runend('read_GFS: Error reading QV')
      !
      do iz = 1, meteo(nm)%np
         ipoin = 0
         do iy = 1, meteo(nm)%ny
            do ix = 1, meteo(nm)%nx
               ipoin = ipoin + 1
               meteo(nm)%qv(ipoin, iz) = work_GFS1(ix, iy, iz)
            end do
         end do
      end do
      !
      !
      work_GFS(:, :) = 6.112*exp((17.67*(meteo(nm)%T(:, :) - 273.15))/(243.5 + (meteo(nm)%T(:, :) - 273.15)))  ! Saturation Pressure in mb
      work_GFS(:, :) = work_GFS(:, :)*meteo(nm)%qv(:, :)/100.0_rp                                ! Vapor      Pressure in mb
      meteo(nm)%qv(:, :) = 0.622*work_GFS(:, :)/((meteo(nm)%p(:, :)/101_rp) - 0.378*work_GFS(:, :))          ! Specific   humidity in kg/kg
      !
      !*** Tv(nx,ny,nz) virtual temperature
      !
      work_GFS = meteo(nm)%qv(:, :)/(1.0_rp - meteo(nm)%qv(:, :))      ! Mixing ratio (humidity ratio)
      !
      meteo(nm)%Tv(:, :) = meteo(nm)%T(:, :)*(1.0_rp + 1.6077_rp*work_GFS(:, :))/(1.0_rp + work_GFS(:, :))
      !
      !*** ro(nx,ny,nz) density. Gas law using virtual temperature (account for water in air)
      !
      meteo(nm)%ro(:, :) = meteo(nm)%p(:, :)/(287.06*meteo(nm)%Tv(:, :))
      !
      !*** w(npoin2d,np) vertical velocity  w = - omega/(rho*g) = - omega*R*Tv/(P*g)
      !
      meteo(nm)%w(:, :) = -meteo(nm)%w(:, :)/(meteo(nm)%ro(:, :)*9.81_rp)
      !
      !*** pblh(nx_GFS,ny_GFS) --> pblh(nx,ny) Atmospheric Boundary Layer
      !
      if (nf90_inq_varid(ncID, pblh_GFS_name, pblhID) /= 0) call runend('read_GFS : Error getting pblhID')
      if (nf90_get_var(ncID, pblhID, work_GFS2, start=(/1, 1, itime/)) /= 0) call runend('read_GFS: Error reading hpbl')
      !
      ipoin = 0
      do iy = 1, meteo(nm)%ny
         do ix = 1, meteo(nm)%nx
            ipoin = ipoin + 1
            meteo(nm)%pblh(ipoin) = work_GFS2(ix, iy, 1)
         end do
      end do
      !
      !
      !*** close file
      !
      deallocate (work_GFS1)
      deallocate (work_GFS2)
      deallocate (work_GFS)
      deallocate (p_GFS)
      if (NF90_close(ncID) /= NF90_NOERR) call runend('Error in  nf90_close gfs file')
      !
      !*** Interpolate in the grid box position for 3D variables
      !
      do iz = 1, grid%nz
         ipoin = 0
         do iy = 1, grid%ny
            do ix = 1, grid%nx
               ipoin = ipoin + 1
               !
               if (grid%model(ipoin) .eq. nm) then
                  !
                  !*** Store previous values
                  !
                  grid%u(ix, iy, iz, 1) = grid%u(ix, iy, iz, 2)
                  grid%u(ix, iy, iz, 2) = 0
                  !
                  grid%v(ix, iy, iz, 1) = grid%v(ix, iy, iz, 2)
                  grid%v(ix, iy, iz, 2) = 0
                  !
                  grid%w(ix, iy, iz, 1) = grid%w(ix, iy, iz, 2)
                  grid%w(ix, iy, iz, 2) = 0
                  !
                  grid%T(ix, iy, iz, 1) = grid%T(ix, iy, iz, 2)
                  grid%T(ix, iy, iz, 2) = 0
                  !
                  grid%ro(ix, iy, iz, 1) = grid%ro(ix, iy, iz, 2)
                  grid%ro(ix, iy, iz, 2) = 0
                  !
                  grid%drhodz(ix, iy, iz, 1) = grid%drhodz(ix, iy, iz, 2)
                  grid%drhodz(ix, iy, iz, 2) = 0.0
                  !
                  grid%qv(ix, iy, iz, 1) = grid%qv(ix, iy, iz, 2)
                  grid%qv(ix, iy, iz, 2) = 0
                  !
                  grid%tp(ix, iy, iz, 1) = grid%tp(ix, iy, iz, 2)
                  grid%tp(ix, iy, iz, 2) = 0
                  !
                  !*** Interpolate
                  !
                  ielem = grid%element(ipoin)   ! meteo model element (2d)
                  do i = 1, 4
                     jpoin = meteo(nm)%lnods(i, ielem)
                     jz = 0
                     found = .false.
                     do while (.not. found)
                        jz = jz + 1
                        !
                        if (grid%z(iz) .lt. meteo(nm)%z(jpoin, 1)) then
                           kz = 1
                           zshape = 1.0
                           found = .true.
                        else if ((grid%z(iz) .ge. meteo(nm)%z(jpoin, jz)) .and. &
                                 (grid%z(iz) .lt. meteo(nm)%z(jpoin, jz + 1))) then
                           kz = jz
                           zshape = (meteo(nm)%z(jpoin, kz + 1) - grid%z(iz))/(meteo(nm)%z(jpoin, kz + 1) - meteo(nm)%z(jpoin, kz))
                           found = .true.
                        else if (jz .eq. (meteo(nm)%np - 1)) then
                           kz = jz
                           zshape = 0.0
                           found = .true.     ! top not found. Use the highest value
                        end if
                        !
                     end do
                     !
                     ! Interpolate U
                     !
                     grid%u(ix, iy, iz, 2) = grid%u(ix, iy, iz, 2) + &
                                             zshape*(grid%shape(i, ipoin)*meteo(nm)%u(jpoin, kz)) + &
                                             (1 - zshape)*(grid%shape(i, ipoin)*meteo(nm)%u(jpoin, kz + 1))
                     !
                     ! Interpolate V
                     !
                     grid%v(ix, iy, iz, 2) = grid%v(ix, iy, iz, 2) + &
                                             zshape*(grid%shape(i, ipoin)*meteo(nm)%v(jpoin, kz)) + &
                                             (1 - zshape)*(grid%shape(i, ipoin)*meteo(nm)%v(jpoin, kz + 1))
                     !
                     ! Interpolate W
                     !
                     grid%w(ix, iy, iz, 2) = grid%w(ix, iy, iz, 2) + &
                                             zshape*(grid%shape(i, ipoin)*meteo(nm)%w(jpoin, kz)) + &
                                             (1 - zshape)*(grid%shape(i, ipoin)*meteo(nm)%w(jpoin, kz + 1))
                     !
                     ! Interpolate T
                     !
                     grid%T(ix, iy, iz, 2) = grid%T(ix, iy, iz, 2) + &
                                             zshape*(grid%shape(i, ipoin)*meteo(nm)%T(jpoin, kz)) + &
                                             (1 - zshape)*(grid%shape(i, ipoin)*meteo(nm)%T(jpoin, kz + 1))
                     !
                     ! Interpolate ro
                     !
                     grid%ro(ix, iy, iz, 2) = grid%ro(ix, iy, iz, 2) + &
                                              zshape*(grid%shape(i, ipoin)*meteo(nm)%ro(jpoin, kz)) + &
                                              (1 - zshape)*(grid%shape(i, ipoin)*meteo(nm)%ro(jpoin, kz + 1))
                     !
                     ! Interpolate qv
                     !
                     grid%qv(ix, iy, iz, 2) = grid%qv(ix, iy, iz, 2) + &
                                              zshape*(grid%shape(i, ipoin)*meteo(nm)%qv(jpoin, kz)) + &
                                              (1 - zshape)*(grid%shape(i, ipoin)*meteo(nm)%qv(jpoin, kz + 1))
                     !
                     ! Interpolate tp
                     !
                     grid%tp(ix, iy, iz, 2) = grid%tp(ix, iy, iz, 2) + &
                                              zshape*(grid%shape(i, ipoin)*meteo(nm)%tp(jpoin, kz)) + &
                                              (1 - zshape)*(grid%shape(i, ipoin)*meteo(nm)%tp(jpoin, kz + 1))
                     !
                  end do  ! i = 1,4
                  !
                  !Change the sign in backward case
                  if (sidt .lt. 0) then
                     grid%u(ix, iy, iz, 2) = grid%u(ix, iy, iz, 2)*(-1)
                     grid%v(ix, iy, iz, 2) = grid%v(ix, iy, iz, 2)*(-1)
                  end if
               end if !grid%model=nm
            end do
         end do
      end do
      !
      !***interpolate in the particle box position for 2D variables.
      !*** and Computes other variables not given by GFS
      !
      ipoin = 0
      do iy = 1, grid%ny
         do ix = 1, grid%nx
            ipoin = ipoin + 1
            if (grid%model(ipoin) .eq. nm) then
               !
               !*** Store previous values
               !
               grid%pblh(ix, iy, 1) = grid%pblh(ix, iy, 2)
               grid%pblh(ix, iy, 2) = 0.0
               !
               grid%rmonin(ix, iy, 1) = grid%rmonin(ix, iy, 2)
               grid%rmonin(ix, iy, 2) = 0
               !
               grid%ust(ix, iy, 1) = grid%ust(ix, iy, 2)
               grid%ust(ix, iy, 2) = 0
               !
               grid%hfx(ix, iy, 1) = grid%hfx(ix, iy, 2)
               grid%hfx(ix, iy, 2) = 0
               !
               grid%wst(ix, iy, 1) = grid%wst(ix, iy, 2)
               grid%wst(ix, iy, 2) = 0
               !
               !
               !*** Interpolate respect the box(ielem) position of the particle.
               !
               ielem = grid%element(ipoin)   ! meteo model element (2d)
               do i = 1, 4
                  jpoin = meteo(nm)%lnods(i, ielem)
                  !
                  ! Interpolate PBLH
                  !
                  grid%pblh(ix, iy, 2) = grid%pblh(ix, iy, 2) + &
                                         grid%shape(i, ipoin)*meteo(nm)%pblh(jpoin)
               end do
               !
               !*** ust(nx,ny), rmonin, and HFX(nx,ny). Estimated
               !
               call get_par_ABL(1.0_rp, grid%z(2), grid%T(ix, iy, 1, 2), grid%T(ix, iy, 2, 2), p(ix, iy, 1), p(ix, iy, 2), &
                                grid%u(ix, iy, 1, 2), grid%v(ix, iy, 1, 2), grid%ust(ix, iy, 2), grid%rmonin(ix, iy, 2), &
                                grid%hfx(ix, iy, 2))
               !
               !*** Calculate the convective velocity scale wst
               !
               if (grid%hfx(ix, iy, 2) .lt. 0.) then       !if surface heat flux is less than 0.
                  grid%wst(ix, iy, 2) = -grid%pblh(ix, iy, 2)*g/grid%tp(ix, iy, 1, 2)*grid%hfx(ix, iy, 2)/cpa
                  if (grid%wst(ix, iy, 2) .ge. 0) then
                     grid%wst(ix, iy, 2) = (grid%wst(ix, iy, 2))**0.333
                  else
                     grid%wst(ix, iy, 2) = -(-grid%wst(ix, iy, 2))**0.333
                  end if
                  !
               else
                  grid%wst(ix, iy, 2) = 0.
               end if
               !
               !*** Calculates drhodz: air density vertical gradient
               !
               grid%drhodz(ix, iy, 1, 2) = (grid%ro(ix, iy, 2, 2) - grid%ro(ix, iy, 1, 2))/(grid%z(2) - grid%z(1))
               do iz = 2, grid%nz - 1
                  grid%drhodz(ix, iy, iz, 2) = (grid%ro(ix, iy, iz + 1, 2) - grid%ro(ix, iy, iz - 1, 2))/ &
                                               (grid%z(iz + 1) - grid%z(iz - 1))
               end do
               grid%drhodz(ix, iy, grid%nz, 2) = grid%drhodz(ix, iy, grid%nz - 1, 2)
            end if
            !
         end do
      end do
      !
      !
      return
   end subroutine readgfs
   !
   subroutine get_par_ABL(z0, z1, Tz0, Tz1, Pz0, Pz1, u, v, ustar, rmonin, hfx)
      !*****************************************************************************
      !*
      !*     This routine estimates some parameters of the ABL
      !*     !FALL3D subroutine, added: hfx
      !*
      !*     INPUTS:
      !*      z0       Reference level 1
      !*      z1       Reference level 2 (z1 > z0)
      !*      Tz0      Temperature (K) at z0
      !*      Tz1      Temperature (K) at z1
      !*      Pz0      Pressure    (Pa) at z0
      !*      Pz1      Pressure    (Pa) at z1
      !*      u        u-velocity  (at z1)
      !*      v        v-velocity  (at z1)
      !*     OUTPUTS:
      !*      ustar    Friction velocity
      !*      rmonin   Monin-Obukhov lenght
      !*
      !*****************************************************************************
      use KindType
      implicit none
      !
      real(rp) :: z0, z1, Tz0, Tz1, Pz0, Pz1, u, v
      real(rp) :: ustar, rmonin, hfx
      !
      real(rp), parameter :: g = 9.81_rp
      real(rp), parameter :: k = 0.4_rp
      real(rp)            :: Thz0, Thz1, umod, Ri, Gm, Gh, thstar
      !
      !***   Potential temperature Theta=T*(1bar/P))**(R/cp)
      !***   R  =  287  J/kg K   Air specific gas  constant
      !***   cp = 1006  J/kg K   Air specific heat capacity
      !
      Thz0 = Tz0*(1.01d5/Pz0)**(0.285)
      Thz1 = Tz1*(1.01d5/Pz1)**(0.285)
      !
      !***   Velocity
      !
      umod = max(sqrt(u*u + v*v), 1e-4_rp)
      !
      !***   Bulk Richardson number
      !
      Ri = g*(z1 - z0)*(Thz1 - Thz0)/(umod*umod*0.5_rp*(Tz0 + Tz1))
      !
      !***   Stable/Unstable ABL
      !
      if (Ri .ge. 0.0_rp) then
         !                             Stable
         Gm = 1.0_rp + 4.7_rp*Ri
         Gm = 1.0_rp/(Gm*Gm)
         Gh = Gm
      else
         !                             Unstable
         Gm = log(z1/z0)
         Gm = 1.0_rp + (70.0_rp*k*k*sqrt(abs(Ri)*z1/z0))/(Gm*Gm)
         Gh = 1.0_rp + (50.0_rp*k*k*sqrt(abs(Ri)*z1/z0))/(Gm*Gm)
         Gm = 1.0_rp - ((9.4_rp*Ri)/Gm)
         Gh = 1.0_rp - ((9.4_rp*Ri)/Gh)
      end if
      !
      !***   ustar
      !
      ustar = log(z1/z0)
      ustar = k*umod*sqrt(Gm)/ustar
      !
      !***   thstar (Pr=1 assumed)
      !
      thstar = log(z1/z0)
      thstar = k*k*umod*(Thz1 - Thz0)*Gh/(ustar*thstar*thstar)
      thstar = max(thstar, 1e-4_rp)
      !
      !***   Monin-Obukhov
      !
      rmonin = ustar*ustar*0.5_rp*(Thz0 + Thz1)/(k*g*thstar)
      hfx = -1006*ustar*thstar
      !
      return
   end subroutine get_par_ABL
   !
END MODULE GFS
