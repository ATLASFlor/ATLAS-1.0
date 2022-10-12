!***************************************************************
!*
!*              Module for DEBUG Meteo
!* 
!***************************************************************
MODULE Debug
  use KindType
  use Master
  IMPLICIT NONE
  SAVE

CONTAINS
  subroutine opendebug(nm)
    !******************************************************************************************
    !*    OPENDEBUG
    !*    Opens a ficticial meteo files
    !*
    !*    nm       Correspond to the meteo file readed by input file  (Input value)
    !*
    !******************************************************************************************
    use Master
    use KindType
    use TimeFun
    use Elsest_mod
    implicit none
    integer (ip)  :: ix, iy, it, nm
    integer (ip)  :: ipoin, ielem
    integer (ip)  :: iyr,imo,idy,ihr,imi,ise
    integer (ip)  :: ibse, ibmi !needed to calculate the time and timesec arrays
    integer(rp)::dimnt
    real(rp)      :: den !auxiliar variable related to time to calculate dimnt
    !
    !*** dimensions like the background grid
    !
    meteo(nm)%nx=grid%nx
    meteo(nm)%ny=grid%ny
    meteo(nm)%np=grid%nz
    meteo(nm)%nt=grid%nt+1
    !
    !*** allocates memory
    !
    den=(siend-sist)/sidt !three variabes are in seconds
    dimnt=INT(meteo(nm)%nt/den )   !54=3600*3/200 y 200=dt segundos en el prueba.inp: 5760.  Denominator is the result of 3600*(hours of simulation)/(simulation time step)
    allocate (meteo(nm)%lon    (meteo(nm)%nx,meteo(nm)%ny))
    allocate (meteo(nm)%lat    (meteo(nm)%nx,meteo(nm)%ny))
    allocate (meteo(nm)%time   (meteo(nm)%nt             ))
    allocate (meteo(nm)%timesec(  dimnt        ))
    !
    ! determines the time and timesec
    !
    if(sidt.gt.0)then
       !forward case
       ibmi = INT((sist/3600-ibhr)*60)
       ibse = INT(((sist/3600-ibhr)*60-ibmi)*60)
       meteo(nm)%time   (1) = 1d8*ibyr + 1d6*ibmo + 1d4*ibdy + 1d2*ibhr + ibmi
       meteo(nm)%timesec(1) = sist! 3600.*ibhr + 60.*ibmi + ibse
       do it = 2,dimnt
          meteo(nm)%timesec(it) = meteo(nm)%timesec(1) + (it-1)*sidt*54  ! numero dividido en asignacion de dimnt
          call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,meteo(nm)%timesec(it))
          meteo(nm)%time(it) = 1d8*iyr + 1d6*imo + 1d4*idy + 1d2*ihr + imi
       end do
    else
       !backard case
       ibmi = INT((siend/3600-ibhr)*60)
       ibse = INT(((siend/3600-ibhr)*60-ibmi)*60)
       meteo(nm)%time   (1) = siend!1d8*ibyr + 1d6*ibmo + 1d4*ibdy + 1d2*ibhr + ibmi
       meteo(nm)%timesec(1) = 3600.*ibhr + 60.*ibmi + ibse
       do it = 2,meteo(nm)%nt
          meteo(nm)%timesec(it) = meteo(nm)%timesec(1) - (it-1)*sidt
          call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,meteo(nm)%timesec(it))
          meteo(nm)%time(it) = 1d8*iyr + 1d6*imo + 1d4*idy + 1d2*ihr + imi
       end do
    end if


    ! There is not time_lag and it is not necesary control that time is between the bounds.


    !
    !*** Lat
    !
    do iy = 1, meteo(nm)%ny
       meteo(nm)%lat(:,iy) = grid%lat(iy)
    end do
    !
    !*** Lon.
    !
    do ix = 1, meteo(nm)%nx
       meteo(nm)%lon(ix,:) = grid%lon(ix)
    end do
    !
    ! allocate memory
    !
    meteo(nm)%npoin2d    =  grid%nx   * grid%ny 
    meteo(nm)%nelem2d    = (grid%nx-1)*(grid%ny-1)
    allocate(meteo(nm)%hgt     (meteo(nm)%npoin2d))
    allocate(meteo(nm)%landmask(meteo(nm)%npoin2d))
    allocate(meteo(nm)%lnods(4,meteo(nm)%nelem2d))
    !
    !*** Topography
    !
    meteo(nm)%hgt(:)=1
    !
    !*** Landmask
    !
    meteo(nm)%landmask(:)=0
    !
    !*** grid%shape
    !
    grid%shape(:,:)=0.25
    grid%model(:)=1
    !
    !*** Resolution
    !
    meteo(nm)%res = grid%dy
    !
    !*** define elements. Nodal conectivity
    !
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
!    meteo(nm)%lnods(:,:)=grid%lnods(1:4,:,1)
! 
    do ielem=1,meteo(nm)%nelem2d
       ipoin=meteo(nm)%lnods(1,ielem)
       grid%element(ipoin)=ielem
    end do
! 
    do iy=1,meteo(nm)%ny-1
       ipoin=meteo(nm)%nx*iy
       grid%element(ipoin)=grid%element(ipoin-1)
    end do
! 
    do ix=1,meteo(nm)%nx
       ipoin=ix+(meteo(nm)%nx-1)*meteo(nm)%ny
       grid%element(ipoin)=grid%element(ipoin-meteo(nm)%nx)
    end do
! 
    !
    allocate (meteo(nm)%z (meteo(nm)%npoin2d,meteo(nm)%np))   ! stored in 2D planes to interpolate with Elsest
    allocate (meteo(nm)%u (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%v (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%w (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%T (meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%ro(meteo(nm)%npoin2d,meteo(nm)%np))
    allocate (meteo(nm)%qv(meteo(nm)%npoin2d,meteo(nm)%np))
! 
    !
    return
  end subroutine opendebug
  !
  !
  !
  subroutine readdebug(nm,itime)
    !******************************************************************************************
    !*    READDEBUG
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
    integer  (ip)         :: nm, itime
    integer  (ip)         :: ix, iy, iz, ipoin

    do iz=1,meteo(nm)%np
       meteo(nm)%z(:,iz)=(iz-1)*grid%dz
    end do
    if(sidt.gt.0)then
       meteo(nm)%u (:,:)=10
       meteo(nm)%v (:,:)=5.
       meteo(nm)%w (:,:)=0
       meteo(nm)%T (:,:)=10.
       meteo(nm)%ro(:,:)=1.0
       meteo(nm)%qv(:,:)=1.0
    else
       meteo(nm)%u (:,:)=-1.2
       meteo(nm)%v (:,:)=5.
       meteo(nm)%w (:,:)=-0.00001
       meteo(nm)%T (:,:)=40.
       meteo(nm)%ro(:,:)=1.0
       meteo(nm)%qv(:,:)=1.0
    end if
    !
    !*** Variables
    !
    do iz=1,meteo(nm)%np
       ipoin=0
       do iy=1,meteo(nm)%ny
          do ix=1,meteo(nm)%nx
             ipoin=ipoin+1
             grid%u (ix,iy,iz,1) = meteo(nm)%u(ipoin,iz)
             grid%u (ix,iy,iz,2) = meteo(nm)%u(ipoin,iz)
             grid%v (ix,iy,iz,1) = meteo(nm)%v(ipoin,iz)
             grid%v (ix,iy,iz,2) = meteo(nm)%v(ipoin,iz)
             grid%w (ix,iy,iz,1) = meteo(nm)%w(ipoin,iz)
             grid%w (ix,iy,iz,2) = meteo(nm)%w(ipoin,iz)
             grid%T (ix,iy,iz,1) = meteo(nm)%T(ipoin,iz)
             grid%T (ix,iy,iz,2) = meteo(nm)%T(ipoin,iz)
             grid%ro(ix,iy,iz,1) = meteo(nm)%ro(ipoin,iz)
             grid%ro(ix,iy,iz,2) = meteo(nm)%ro(ipoin,iz)
             grid%qv(ix,iy,iz,1) = meteo(nm)%qv(ipoin,iz)
             grid%qv(ix,iy,iz,2) = meteo(nm)%qv(ipoin,iz)
          end do
       end do
    end do
  end subroutine readdebug
end MODULE Debug
