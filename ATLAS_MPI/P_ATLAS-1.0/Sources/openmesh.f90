subroutine openmesh
  !******************************************************************************************
  !*    OPENMESH
  !*    Opens the Meteorological file and get variables values to create the meteo mesh later.
  !*
  !*
  !******************************************************************************************
  use InpOut
  use Master
  use KindType
  use netcdf
  use Wrf 
  use EraInt
  use GFS                                         
  implicit none

  integer  (ip)         :: istat,nm,ipoin,ix,iy,jpoin,ielem,i             
  integer  (ip)         :: ivoid(1)               
  !
  character(len=s_mess) :: message, cvoid(2), str
  character(len=1     ) :: ext
  !
  real     (rp)         :: rvoid(1)
  !
  logical               :: outside=.false.
  !

  allocate(meteo(1:numod))
  !
  !*** Read each meteo mesh. Loop for meteo mesh
  !
  do nm = 1,numod                         ! Read the path and the model from the input file
     !        
     write(ext,'(i1)') nm
     str = 'METEO_DATA_'//ext
     !
     call get_input_cha(finp, 'METEO_DATA',TRIM(str), cvoid,2, istat, message) 
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     meteo(nm)%file_path=cvoid(1)
     call upcase(cvoid(2))
     meteo(nm)%modeltype=cvoid(2)
     !
     !***    Open and read meteo model mesh and time-independent variables (e.g. topography)
     !
     if(TRIM(meteo(nm)%modeltype)==TRIM('WRF') then
        call openwrf(nm)                
     else if (TRIM(meteo(nm)%modeltype)==TRIM('GFS') then
        call opengfs(nm)
     else if (TRIM(meteo(nm)%modeltype)==TRIM('ERAINT') then
        call openeraint(nm)
     else
        call runend('The model Type is not recognize')
     end if
     !
     !***    Determine if the model is global
     !
     if(meteo(nm)%lon(meteo(nm)%nx,1)+meteo(nm)%res.eq.180 .and. grid%lonmin.eq.-180 .and. grid%lonmax.eq.180) grid%global=.true.
     !
  end do
  !
  !*** Determine if any domain element is outside of all the meteo models
  !
  do ipoin=1,grid%npoin
     if(grid%model(ipoin).eq.-1) outside=.true.
  end do
  if(outside.eqv..true.) call runend('There are elements outside the Meteorological coverage')
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
        nm    = grid%model(ipoin)     ! meteo model
        ielem = grid%element(ipoin)   ! meteo model element (2d)
        grid%hgt(ix,iy) = 0.0  
        do i = 1,4
           jpoin = meteo(nm)%lnods(i,ielem)  
           grid%hgt(ix,iy) = grid%hgt(ix,iy) + grid%shape(i,ipoin)*meteo(nm)%hgt(jpoin) 
        end do
        !
     end do
  end do
  !
  return
end subroutine openmesh
