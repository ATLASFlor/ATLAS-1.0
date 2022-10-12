   subroutine checkdomain(state,lon,lat,z)
  !******************************************************************************************
  !*    
  !*
  !******************************************************************************************
  use Master
  use KindType                              
  implicit none
  integer      :: state
  real   (rp)  :: lon, lat, z
  !
  !***Check if the particle is outside the domain from longitude
  !
  if((lon.ge.grid%lon(grid%nx)).or.(lon.lt.grid%lon(1))) then
     if(grid%PBC_EW)then
write(*,*)'cambio de lon'
        lon=lon-sign(360.0_rp,lon)
     else
        state= -2
     end if
  end if
  !
  !*** Check if the particle is outside the domain from top latitude
  !
  if(lat .ge. grid%lat(grid%ny))then
     if(grid%PBC_N) then
write(*,*)'cambio de lon'
        lon = lon-sign(180.0_rp,lon)
        lat = 180.0_rp-lat
     else
        state= -2
     end if
  end if
  !
  !*** Check if the particle is outside the domain from bottom latitude
  !
  if(lat .lt. grid%lat(1))then
     if(grid%PBC_S) then
write(*,*)'cambio de lon'
        lon = lon-sign(180.0_rp,lon)
        lat = 180.0_rp+lat
     else
        state= -2
     end if
  end if
  !
  !*** Check if the particle is outside the domain from top heigth
  !
  if(z .ge. grid%ztop)then
     state = -2
  end if
  !
  !*** Check if the particle is settled on land
  !
  if(z.le.grid%z(1).and.state.eq.1) then
     if(resuspension.eqv. .true.) then
        z=grid%z(1)
     else
        state = -1
        z=grid%z(1)
     end if
  end if
  !      
  return     
  end subroutine checkdomain
