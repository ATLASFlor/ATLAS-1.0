!***************************************************************
!*
!*              Module for master operations
!* 
!***************************************************************
MODULE Master
  use KindType
  IMPLICIT NONE
  SAVE
  !
  !*** General
  !
  real(rp), parameter   :: version = 1.5
  !
  character(len=s_name) :: problemname    
  !
  !*** MPI variables
  !
  integer (ip) ::  mpi_global
  integer (ip)    :: ierr,num_procs, my_id
  !
  !*** Constants
  !
  real(rp) :: pi     = 3.14159265358979323846_rp   ! pi
  real(rp) :: Rearth = 6356000                     ! Earth's radius in meters.
  real(rp), parameter :: g      = 9.81_rp          ! gravity
  real(rp), parameter :: cpa    = 1004.6_rp        ! specific heat for dry air
  !
  !*** Time related variables
  !
  integer(ip) :: ibyr,ibmo,ibdy,ibhr           ! Begin year, month, day and hour
  real   (rp) :: sist, siend                   ! Times for simulation start and simulation end in seconds after 00:00 UTC 
  real   ( 8) :: time_start                    ! Simulation start time in format YYYYMMDDHHMM 
  real   ( 8) :: time_end                      ! Simulation end time in format YYYYMMDDHHMM 
  integer(ip) :: fiyr,fimo,fidy,fihr,fimi,fise ! Simulation end year, month, day, hour, minute, second
  real   (rp) :: sidt                          ! Simulation time increment. is negative in backward case
  integer(ip) :: iiter                         ! Simulation iteration Number.
  real   (rp) :: time                          ! Simulation time, begin in sist, increment in sidt
  character(len=19) :: actualtime              ! In format: yyyy/mm/dd-hh:mm:ss
  real   (rp) :: sit                           ! Simulation duration time. Seconds since the initial simulation hour
  !
  !*** random numbers related variables
  !
  integer(ip) :: maxrand= 20000                  ! random numbers generated 
  real(rp) :: rannumb(20000)                     ! random numbers array
  integer(ip) :: nran                            ! random number position. Set to 1 in initialize subroutine. 
                                                 ! and back to one only when reaches the maximium.
  !
  !*** Diffussion related variables/ctes
  !
  real(rp)    :: d_trop  = 50.
  real(rp)    :: d_strat = 0.1
  real(rp)    :: tropopause    !tropopause value
  logical     :: split=.false.  !if true, make the splitting (double of particles with half of the mass, after the splitting time)
  integer(ip) :: dts = 43200   ! delta time splitting: in seconds since eruption start (24hs). After that time make a split 
  !
  !*** Define the background mesh strucutre (this mesh is used to define the computational
  !*** domain and to interpolate meteorological data to the particle position).
  !
  type background_mesh
     !
     logical :: PBC_EW = .false.               ! Boundary Condition East-West (grid%lon: (-180,180))
     logical :: PBC_S  = .false.               ! Boundary Condition South     (grid%latmin=-90)
     logical :: PBC_N  = .false.               ! Boundary Condition North     (grid%latmax= 90)
     !
     integer(ip) :: nx                        ! Number of x-grid points
     integer(ip) :: ny                        ! Number of y-grid points
     integer(ip) :: nz                        ! Number of z-grid points
     integer(ip) :: nt                        ! Number of times.
     integer(ip) :: npoin                     ! Number of points. 3d
     integer(ip) :: npoin2d                   ! Number of points. 2d 
     integer(ip) :: nelem                     ! Number of elements. Each element is a hexaedron.
     integer(ip) :: ndime = 3                 ! spatial dimension
     integer(ip) :: nnode = 8                 ! number of points of a cell
     !
     real(rp) :: dx		              ! x-resolution. degree   (from the input file)
     real(rp) :: dy                           ! y-resolution. degree   (from the input file)
     real(rp) :: dz                           ! z-resolution. meters   (from the input file)
     real(rp) :: ztop                         ! Maximum heigth. meters (from the input file)
     real(rp) :: lonmin                   
     real(rp) :: lonmax
     real(rp) :: latmin
     real(rp) :: latmax
     !
     !  Fixed variables
     !
     real(rp), allocatable :: lon(:)          ! lon(nx)
     real(rp), allocatable :: lat(:)          ! lat(ny)
     real(rp), allocatable :: z  (:)          ! z  (nz), height of each level.
     real(rp), allocatable :: hgt     (:,:)   ! hgt     (nx,ny)     Topography
     real(rp), allocatable :: landmask(:,:)   ! landmask(nx,ny)     Landmask
     !
     !  Time dependent variables
     !
     real(rp), allocatable :: u (:,:,:,:)     ! u (nx,ny,nz,2)     Horizontal winds component. 
                                              ! 1 first time. 2 second time   
     real(rp), allocatable :: v (:,:,:,:)     ! v (nx,ny,nz,2)     Vertical   winds component.  
     real(rp), allocatable :: w (:,:,:,:)     ! w (nx,ny,nz,2)     Omega Vertical Component of wind
     real(rp), allocatable :: T (:,:,:,:)     ! T (nx,ny,nz,2)     Temperature
     real(rp), allocatable :: ro(:,:,:,:)     ! ro(nx,ny,nz,2)     Air Density
     real(rp), allocatable :: qv(:,:,:,:)     ! qv(nx,ny,nz,2)     Specific humidity
     real(rp), allocatable :: tp(:,:,:,:)     ! Tp(nx,ny,nz,2)     Potential Temperature
     real(rp), allocatable :: tropopause(:,:) ! tropopause(nx,ny)  Tropopause level
     real(rp), allocatable :: pblh(:,:,:)     ! pblh(nx,ny,2)      Planetary Boudary Layer.
     real(rp), allocatable :: hfx(:,:,:)      ! hfx(nx,ny,2)       Surface heat flux
     real(rp), allocatable :: ust(:,:,:)      ! ust(nx,ny,2)       Friction Velocity U*
     real(rp), allocatable :: rmonin(:,:,:)   ! rmonin(nx,ny,2)      Monin Obukhov Lenght 
     real(rp), allocatable :: wst(:,:,:)      ! wst(nx,ny,2)       Convective Velocity scale W*
     real(rp), allocatable :: drhodz(:,:,:,:) ! drhodz(nx,ny,nz,2)   Air density vertical gradient
     !
     !
     integer(ip), allocatable :: lnods(:,:)   ! lnods(nnode,nelem)    8 nodes per element
     real(rp)   , allocatable :: coord(:,:)   ! coord(ndime,npoin)    lon lat per point
     !
     !  For interpolation of meteo data
     !
     integer(ip), allocatable :: model(:)     ! model  (npoin2d)      meteo model flag 
     integer(ip), allocatable :: element(:)   ! element(npoin2d)      meteo element
     real   (rp), allocatable :: shape(:,:)   ! shape  (4,npoin2d)    shape factor
     !
  end type background_mesh
  type(background_mesh) :: grid 
  !
  !*** Defines the structures for the meteo meshes 
  !
  type meteo_mesh
     !
     logical               :: global = .false.
     !
     character(len=s_name) :: file_path       ! path where the meteo file is
     character(len=s_name) :: modeltype       ! WRF, etc 
     logical               :: postprocess     ! Create netcdf file for this meteo file interpolated to computational grid
     !
     integer(ip) :: nx                        ! x-Dimension (from the meteo file)
     integer(ip) :: nx_stag                   ! x-Dimension staggered (from the meteo file)
     integer(ip) :: ny                        ! y-Dimension (from the meteo file) 
     integer(ip) :: ny_stag                   ! y-Dimension staggered (from the meteo file) 
     integer(ip) :: np                        ! z-Dimension (from the meteo file) nuber of pressure levels
     integer(ip) :: np_stag                   ! z-Dimension staggered (from the meteo file) 
     integer(ip) :: nt                        ! t-Dimension (from the meteo file) 
     integer(ip) :: nnode = 4                 ! number of points of a 2D cell
     integer(ip) :: nelem2d                   ! Number of elements. Each element is a square (2).
     integer(ip) :: npoin2d                   ! number of points. grid corners.
     !
     real   (8)  :: res                       ! Resolution in degrees.
     real   (rp) :: time_lag                  ! Time in seconds, between meteo start and simulation start
     integer(ip) :: itime1, itime2            ! time steps.
     real   (rp) :: time_s                    ! Interpolation factor
     !
     real   (8 ), allocatable :: time (:)     ! time(nt) In format YYYYMMDDHHMM. Stored always in real(8)
     real   (rp), allocatable :: timesec (:)  ! time(nt) In seconds. 
     integer(ip), allocatable :: lnods(:,:)   ! lnods(nnode,nelem2d) 
     !
     real(rp), allocatable :: lon(:,:)        ! lon(nx,ny)
     real(rp), allocatable :: lat(:,:)        ! lat(nx,ny)
     real(rp), allocatable :: hgt(:)          ! hgt(npoin2d)    
     real(rp), allocatable :: landmask(:)     ! landmask(npoin2d) 
     !
     real(rp), allocatable :: z (:,:)          ! z(npoin2d,np) heights   
     real(rp), allocatable :: u (:,:)          ! u(npoin2d,np) wind velocities
     real(rp), allocatable :: v (:,:)          ! v(npoin2d,np) 
     real(rp), allocatable :: w (:,:)          ! w(npoin2d,np) 
     real(rp), allocatable :: p (:,:)          ! p(npoin2d,np) pressure
     real(rp), allocatable :: T (:,:)          ! T(npoin2d,np) temperature
     real(rp), allocatable :: Tv(:,:)          ! Tv(npoin2d,np) virtual temperature
     real(rp), allocatable :: qv(:,:)          ! qv(npoin2d,np) specific humidity
     real(rp), allocatable :: tp(:,:)          ! tp(npoin2d,np) potential temperature
     real(rp), allocatable :: ro(:,:)          ! ro(npoin2d,np)  density
     real(rp), allocatable :: pblh(:)          ! pblh(npoin2d)  planetary boundary layer
     real(rp), allocatable :: hfx(:)           ! hfx (npoin2d)  surface heat flux
     real(rp), allocatable :: ust(:)           ! ust (npoin2d)  Friction velocity
     !
  end type meteo_mesh
  integer(ip)                   :: numod       ! number of Meteo models
  type(meteo_mesh), allocatable :: meteo(:)    ! meteo(1:numod)
  !
  !*** Defines the structure for particles
  !
     integer(ip)::numpart  ! Total number of particles
  type particles
     !general properties
     integer(ip):: state   ! 1:active, 0:still not out, -1: sedimented inside the domain, -2:outside the domain
     integer(ip):: ibin    ! number of bin
     integer(ip):: age     ! particle age to split
     !phisical properties
     real(rp)   :: rho     ! particle density
     real(rp)   :: diam    ! particle diameter
     real(rp)   :: mass    ! particle mass
     real(rp)   :: sphe    ! particle sphericity
     real(rp)   :: psi     ! particle shape factor(depend on the Sedimentation model)
     !geographical properties
     real(rp)   :: lon     ! particle longitude
     real(rp)   :: lat     ! particle latitude
     real(rp)   :: z       ! height above terrain
     !
     real(rp)   :: vset    ! terminal velocity 
     real(rp)   :: uvw(3)  ! Random velocities, uvw(1): ux, uvw(2): vy, uvw(3): wp
     !
     !  Air properties at the particle point
     !
     real(rp) :: ua     ! air velocity
     real(rp) :: va
     real(rp) :: wa
     real(rp) :: Ta     ! air temperature 
     real(rp) :: rhoa   ! air density
     real(rp) :: mua    ! air viscosity
     real(rp) :: qva    ! air viscosity
     real(rp) :: pblha  ! planetary boundary layer
     real(rp) :: hfxa   ! surface heat flux
     real(rp) :: usta   ! Friction velocity
     real(rp) :: rmonin ! rmonin obukhov lenght 
     real(rp) :: wst    ! convective velocity scale
     real(rp) :: drhodz ! air density vertical gradient
     !
  end type particles
  type(particles), allocatable :: part(:)      ! part(npart)
  !
  !*** Define the structure for the source term (phases).
  !
  logical               :: resuspension = .false.
  logical               :: restart      = .false.
  integer(ip)           :: sedmodel                ! Velocity Sedimentation Model
  integer(ip), parameter:: ncmax = 50              ! Maximum number of classes
  integer(ip)           :: nphases                 ! Number of phases
  !
  type source_phases
     logical              :: activated          ! is true since the simulation time enter in the phase time
     integer  (ip)        :: npart              ! number of particles of this phase
     integer  (ip)        :: totpart            ! total number of particles including the splitting
     integer  (ip)        :: npartsplit=0       ! number of particles splitted
     integer  (ip)        :: firstpart,lastpart ! particles range corresponding to each phase.
     integer  (ip)        :: actnpart           ! number of activated particles
     integer  (ip)        :: npdt               ! number of particles per simulation time interval.
     integer  (ip)        :: partsout           ! number of particles to output per output-time interval
     integer  (ip)        :: partstime          ! number of particles over number of outputs
     integer  (ip)        :: firstoutt,endoutt  ! first and end time to output new particles(used if the number 
                                                ! of particles is grater than a fixed)
     !
     character(len=s_name):: sname              ! Phase name
     character(len=s_name):: phase_type         ! Phase type: eruption/satellite/resuspension/restart
     character(len=s_name):: source_type        ! point/linear/top-hat/suzuki(eruption). / / (resuspension)
     !
     real(rp),allocatable :: begtime(:)         ! begin phase. In seconds after 00hs (simulation start day)
     integer  (ip)        :: nsp                ! number of sub-phases (e.g. for height change)
     real     (rp)        :: endtime            ! end phase
     real     (rp)        :: duration           ! phase duration
     !
     real(rp),allocatable :: colheight(:)       ! Column height
     real(rp),allocatable :: a_suzuki(:)        ! a_suzuki
     real(rp),allocatable :: l_suzuki(:)        ! l_suzuki
     real(rp),allocatable :: d_top_hat(:)       ! d_top_hat
     !
     character(len=s_name):: volname            ! Volcano name
     real     (rp)        :: lon                ! Phase's longitude
     real     (rp)        :: lat                ! Phase's latitude
     real     (rp)        :: source_elev        ! Surce elevation
     !
     !  Granulometric properties of the phase
     !
     character(len=s_name):: grnpath            ! Granulometric path
     character(len=s_name):: grndist            ! Granulometric distribution(Gaussian/Bigaussian)
     integer  (ip)        :: bins               ! Number of bins
     real     (rp)        :: fimin,fimax        !
     real     (rp)        :: fimean(2),fidisp(2)
     real     (rp)        :: rhomin,rhomax
     real     (rp)        :: rhomean            ! used for Source term (sum of rhop(ic)*fc(ic))
     real     (rp)        :: sphemin,sphemax  
     !
     integer  (ip)                          :: np,ng,mdist           ! .type distribution. np=number of particles. ng= number of modes(1=gaussian. 2=bigaussian)
                                                                     ! ng=1 Gaussian TGSD, ng=2 Bigaussian TGSD
     real     (rp),     dimension(ncmax)    :: rhop,diam,sphe,fc,psi ! particles properties, density, diammeter, 
                                                                     ! sphericity, fraction of mass, psi
     character(len=25), dimension(ncmax)    :: classname             ! Class name  
     !
     !   Aggregation
     !
     logical              :: aggregation = .false.
     character(len=s_name):: aggmodel              ! Aggregation Model
     real     (rp)        :: aggsize               ! Aggregation size
     real     (rp)        :: aggrho                ! Aggregation density
     real     (rp)        :: aggfrac               ! aggregation percentage (Only for percentage model)
     !
     !   Source term
     !
     integer  (ip)        :: ns                    ! number of source
     integer  (ip)        :: nps                   ! number of plume source
     character(len=s_name):: mer_vs_h              ! mass flow rate
     logical              :: MER_wind_coupling     ! is true for ESTIMATE-DEGRUYTER and ESTIMATE-WOODHOUSE mer_vs_h mode
     real(rp),allocatable :: M0(:)                 ! Mass flow rate from input file M0(1:nsp)
     real     (rp)        :: erumass               ! erupted mass
     !
  end type source_phases
  type(source_phases), allocatable :: phase(:)     ! phase(1:nphases)
  !
  !*** Define structure for the output (postprocess) mesh
  !
  type output_mesh
     !
     logical     :: classes                   ! If output information per classes
     logical     :: phases                    ! If output information per phases
     logical     :: track_points              ! If output information per track points
     integer(ip) :: nx                        ! number of grid points in x axis
     integer(ip) :: ny
     integer(ip) :: nz
     integer(ip) :: nt
     integer(ip) :: nelem2d                   ! Number of elements.
     integer(ip) :: nelem3d                   ! Number of elements. Each element is a hexaedron.
     !
     real(rp) :: dx		              ! x-resolution. degree (from the input file)
     real(rp) :: dy                           ! y-resolution. degree (from the input file)
     real(rp) :: dz                           ! z-resolution. meters (from the input file)
     real(rp) :: ztop                         ! Maximum heigth. meters (from the input file)
     real(rp), allocatable::zlayer(:)         ! zlayer(output%nz)
     real(rp) :: lonmin                   
     real(rp) :: lonmax
     real(rp) :: latmin
     real(rp) :: latmax
     real(rp) :: rhomean                      ! used for Source term (sum of rhop(ic)*fc(ic))
     !
     real(rp)             :: frequency        ! time frequency to output
     real(rp),allocatable :: timesec(:)       ! timesec(output%nt) in seconds since 0 of sist day
     !
     real(rp),allocatable :: area(:)          ! area(nelem2d). To calculate concentrations
     real(rp),allocatable :: volume(:)        ! volume(nelem3d)
     !
  end type output_mesh
  type(output_mesh) :: output 
  !
CONTAINS
  !
END MODULE Master
