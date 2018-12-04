subroutine write_nc_grid
  !**************************************************************************
  !*
  !*   Writes grid data (including coordinate variables) in netCDF format
  !*
  !**************************************************************************
  use KindType
  use InpOut
  use Master
  use netcdf
  implicit none
  !
  !*** Open netCDF file (define mode) and get ncID
  !
  !     if( nf90_create(TRIM(luncname),NF90_CLOBBER, ncID) /= 0 ) &
  !         call runend('write_nc_grid : Error in nf90_create')

  if( nf90_create(TRIM(luncname),NF90_64BIT_OFFSET, ncID) /= 0 ) &
       call runend('write_dbs_grid : Error in nf90_create')
  !
  !**  Define dimensions
  !
  if( nf90_def_dim(ncID, lon_nc_name , nx, nx_nc_ID ) /= 0 ) &
       call runend('write_nc_grd : error in nf90_def_dim for lon')
  if( nf90_def_dim(ncID, lat_nc_name , ny, ny_nc_ID ) /= 0 ) &
       call runend('write_nc_grd : error in nf90_def_dim for lat')
  if( nf90_def_dim(ncID, pre_nc_name , np, np_nc_ID ) /= 0 ) &
       call runend('write_nc_grd : error in nf90_def_dim for pre')
  if( nf90_def_dim(ncID, time_nc_name, nt+1, nt_nc_ID) /= 0 ) &
       call runend('write_nc_grd : error in nf90_def_dim for time')
  !
  !*** Define coordinate variables
  !
  if( nf90_def_var(ncID, lon_nc_name ,NF90_FLOAT, (/nx_nc_ID/), lon_nc_ID) /= 0 ) &
       call runend('write_nc_grid : error in nf90_def_variable for variable '//TRIM(lon_nc_name))
  if( nf90_def_var(ncID, lat_nc_name ,NF90_FLOAT, (/ny_nc_ID/), lat_nc_ID) /= 0 ) &
       call runend('write_nc_grid : error in nf90_def_variable for variable '//TRIM(lat_nc_name))
  if( nf90_def_var(ncID, pre_nc_name ,NF90_FLOAT, (/np_nc_ID/), pre_nc_ID) /= 0 ) &
       call runend('write_nc_grid : error in nf90_def_variable for variable '//TRIM(pre_nc_name))
  if( nf90_def_var(ncID, time_nc_name,NF90_DOUBLE, (/nt_nc_ID/), time_nc_ID) /= 0 ) &
       call runend('write_nc_grid : error in nf90_def_variable for variable '//TRIM(time_nc_name))
  !
  !*** Define variables
  !
  do i2d = 1,n2d
     if( nf90_def_var(ncID, TRIM(var2d(i2d)), NF90_FLOAT, &
          (/nx_nc_ID,ny_nc_ID,nt_nc_ID/), var2d_ID(i2d)) /= 0 ) &
          call runend('write_nc_grid : error in nf90_def_var for variable '//TRIM(var2d(i2d)))
  end do
  !
  do i3d = 1,n3d
     if( nf90_def_var(ncID, TRIM(var3d(i3d)), NF90_FLOAT, &
          (/nx_nc_ID,ny_nc_ID,np_nc_ID,nt_nc_ID/), var3d_ID(i3d)) /= 0 ) &
          call runend('write_nc_grid : error in nf90_def_var for variable '//TRIM(var3d(i3d)))
  end do
  !
  !*** Define attributes
  !
  do i2d = 1,n2d
     if( nf90_put_att(ncID, var2d_ID(i2d), 'units', att2d_u(i2d)) /= 0 ) &
          call runend('write_nc_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, var2d_ID(i2d), 'description', att2d_d(i2d)) /= 0 ) &
          call runend('write_nc_grid : error in nf90_put_att')
  end do
  !
  do i3d = 1,n3d
     if( nf90_put_att(ncID, var3d_ID(i3d), 'units', att3d_u(i3d)) /= 0 ) &
          call runend('write_nc_grid : error in nf90_put_att')
     if( nf90_put_att(ncID, var3d_ID(i3d), 'description', att3d_d(i3d)) /= 0 ) &
          call runend('write_nc_grid : error in nf90_put_att')
  end do
  !
  !*** Put global attributes
  !
  if( nf90_put_att(ncID, NF90_GLOBAL, 'Grib2nc version', '7.1') /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'LONMIN', lonmin) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'LONMAX', lonmax) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'LATMIN', latmin) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'LATMAX', latmax) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  !
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NX', nx) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NY', ny) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NP', np) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'NT', nt+1) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'YEAR', ibyr) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'MONTH', ibmo) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'DAY', ibdy) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'HOUR', ibhr) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'TIME_INCR', dt*3600) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'CEN_LON', cen_lon) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'CEN_LAT', cen_lat) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  if( nf90_put_att(ncID, NF90_GLOBAL, 'missing_value', missing_value) /= 0 ) &
       call runend('write_nc_grid : error in nf90_put_att')
  !
  !*** Leave the define mode
  !
  if( nf90_enddef(ncID) /= 0 ) call runend('write_nc_grid: error in nf90_enddef')
  !
  !*** Write coordinate variables
  !
  if( nf90_put_var(ncID, nx_nc_ID, lon(1:nx)) /= 0 )  &
       call runend('write_nc_grid: error in nf90_put_var for varialbe lon')
  if( nf90_put_var(ncID, ny_nc_ID, lat(1:ny)) /= 0 )  &
       call runend('write_nc_grid: error in nf90_put_var for varialbe lat')
  if( nf90_put_var(ncID, pre_nc_ID, pres(1:np)) /= 0 )  &
       call runend('write_nc_grid: error in nf90_put_var for varialbe pre')
  if( nf90_put_var(ncID, time_nc_ID, timesec(0:nt)) /= 0 )  &
       call runend('write_nc_grid: error in nf90_put_var for varialbe time')
  !
  !*** Close the file
  !
  if( nf90_close(ncID) /= 0) call runend('write_nc_grid : Error in nf90_close')
  !
  return
end subroutine write_nc_grid
