subroutine set_meteo_model
  !******************************************************************************************
  !*    SET_METEO_MODEL
  !*    Opens the Meteorological file and get variables values to create the meteo mesh later.
  !*
  !*
  !******************************************************************************************
  use InpOut
  use Master
  use KindType
  use netcdf
  use Wrf 
  use GFS
  use EraInt  
  use Debug                                       
  implicit none

  integer  (ip)         :: imodel,ipoin,ix,iy,jpoin,ielem,i             
  !
  character(len=s_mess) :: str
  character(len=1     ) :: ext
  !
  logical               :: outside=.false.
  !
  !*** Open and read meteo models and time-independent variables (e.g. topography)
  !
  !*** Loop for models
  !
  do imodel = 1,numod  
     write(ext,'(i1)') imodel
     str = 'METEO_DATA_'//ext
     !                    
     if(TRIM(meteo(imodel)%modeltype)=='WRF') then
        call openwrf(imodel)   
     else if (TRIM(meteo(imodel)%modeltype)=='GFS')then
        call opengfs(imodel)
     else if (TRIM(meteo(imodel)%modeltype)=='ERAINT')then
        call openeraint(imodel)
     else if (TRIM(meteo(imodel)%modeltype)=='DEBUG')then
        call openDEBUG(imodel)
     else
        call runend('The model Type is not recognized')
     end if
     !
  end do ! loop for models
  !
  !*** Determine if any domain element is outside of all the meteo models
  !
  !*** loop for computational-grid points(2D)
  do ipoin = 1,grid%npoin2d
     if(grid%model(ipoin).eq.-1) outside=.true.
  end do
  if(outside.eqv..true.) call runend('The input domain is not inside the Meteorological coverage')
  !
  !*** Interpolate time-independent variables from meteo model to the background mesh
  !
  !    1. Topography
  !
  ipoin = 0
  do iy = 1,grid%ny
     do ix = 1,grid%nx 
        ipoin = ipoin + 1
        !
        imodel= grid%model(ipoin)     ! meteo model
        ielem = grid%element(ipoin)   ! meteo model element (2d)
        grid%hgt(ix,iy)      = 0.0  
        do i = 1,4
           jpoin = meteo(imodel)%lnods(i,ielem)  
           grid%hgt(ix,iy) = grid%hgt(ix,iy) + grid%shape(i,ipoin)*meteo(imodel)%hgt(jpoin) 
        end do
        !
     end do
  end do
  !
  !    2. Landmask
  !
  ipoin = 0
  do iy = 1,grid%ny
     do ix = 1,grid%nx 
        ipoin = ipoin + 1
        !
        imodel= grid%model(ipoin)     ! meteo model
        ielem = grid%element(ipoin)   ! meteo model element (2d)
        grid%landmask(ix,iy) = 0.0
        do i = 1,4
           jpoin = meteo(imodel)%lnods(i,ielem)  
           grid%landmask(ix,iy) = grid%landmask(ix,iy) + grid%shape(i,ipoin)*meteo(imodel)%landmask(jpoin) 
        end do
        if( grid%landmask(ix,iy) .le. 0.5_rp) then
           grid%landmask(ix,iy)= 0.0_rp
        else if(grid%landmask(ix,iy) .gt. 0.5_rp) then
           grid%landmask(ix,iy)= 1.0_rp
        end if
     end do
  end do
  !
  !*** Print Information in log file
  !
  !***  Writes the log file
  !
  if(my_id.eq.0)write(lulog,2)            'Interpolate Time Independet variables from Meteorological files'
  if(out_screen) write(*,2) 'Interpolate Time Independet variables from Meteorological files'
  !
2 format(/, a /)
  return
end subroutine set_meteo_model
