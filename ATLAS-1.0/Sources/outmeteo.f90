subroutine outmeteo
  !*********************************************************************
  !* This subroutine creates a netcdf file with meteo information
  !* in each grid point. This assignment is done according the meteo file 
  !* that have better resolution for each point.
  !**********************************************************************
  use KindType
  use InpOut
  use Master
  use Netcdf
  implicit none

  character(len=s_mess) :: attr_desc, attr_units
  integer  (ip        ) :: ipoin
  !
  !*** Dictionary
  !
  !DIMENSIONS
  character(len=25) :: nx_met_name    = 'lon'
  character(len=25) :: ny_met_name    = 'lat'
  character(len=25) :: nz_met_name    = 'alt'
  character(len=25) :: nt_met_name    = 'time'
  !ID
  integer  (ip)     :: ncid
  integer  (ip)     :: nx_met_id
  integer  (ip)     :: ny_met_id
  integer  (ip)     :: nz_met_id
  integer  (ip)     :: nt_met_id
  !VARIABLES
  character(len=25) :: lon_met_name     = 'lon'
  character(len=25) :: lat_met_name     = 'lat'
  character(len=25) :: alt_met_name     = 'alt'
  character(len=25) :: hgt_met_name     = 'hgt'
  character(len=25) :: land_met_name    = 'land'
  character(len=25) :: model_met_name   = 'model'
  character(len=25) :: u_met_name       = 'U'
  character(len=25) :: v_met_name       = 'V'
  character(len=25) :: w_met_name       = 'W'
  character(len=25) :: T_met_name       = 'T'
  character(len=25) :: ro_met_name      = 'ro'
  character(len=25) :: qv_met_name      = 'qv'
  !ID
  integer  (ip)     :: lon_met_id
  integer  (ip)     :: lat_met_id
  integer  (ip)     :: alt_met_id
  integer  (ip)     :: hgt_met_id
  integer  (ip)     :: land_met_id
  integer  (ip)     :: model_met_id

  !
  !*** Create a netcdf file
  !
  if( nf90_create(TRIM(fnc_met), NF90_CLOBBER, ncid)/=NF90_NOERR) &
       call runend('Error creting netcdf file in outmeteo')
  !
  !*** Define Dimensions
  !
  if( nf90_def_dim(ncid, nx_met_name, grid%nx, nx_met_id)/=NF90_NOERR)&
       call runend('Error creting netcdf nx dimension in outmeteo')
  !
  if( nf90_def_dim(ncid, ny_met_name, grid%ny, ny_met_id)/=NF90_NOERR)&
       call runend('Error creting netcdf ny dimension in outmeteo')
  !
  if( nf90_def_dim(ncid, nz_met_name, grid%nz, nz_met_id)/=NF90_NOERR)&
       call runend('Error creting netcdf nz dimension in outmeteo')
  !
  if( nf90_def_dim(ncid, nt_met_name, NF90_UNLIMITED, nt_met_id)/=NF90_NOERR)&
       call runend('Error creting netcdf nt dimension in outmeteo')
  !
  !*** Define time-independent variables
  !
  if( nf90_def_var(ncID, lon_met_name ,NF90_FLOAT, (/nx_met_id/), lon_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(lon_met_name))
  !
  if( nf90_def_var(ncID, lat_met_name ,NF90_FLOAT, (/ny_met_id/), lat_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(lat_met_name))
  !
  if( nf90_def_var(ncID, alt_met_name ,NF90_FLOAT, (/nz_met_id/), alt_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(alt_met_name))
  !
  if( nf90_def_var(ncID, hgt_met_name ,NF90_FLOAT, (/nx_met_id,ny_met_id/), hgt_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(hgt_met_name))
  !
  if( nf90_def_var(ncID, land_met_name ,NF90_INT, (/nx_met_id,ny_met_id/), land_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(land_met_name))
  !
  if( nf90_def_var(ncID, model_met_name ,NF90_INT, (/nx_met_id,ny_met_id/), model_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(land_met_name))
  !
  !*** Define 4D Variables
  !
  if( nf90_def_var(ncID, u_met_name ,NF90_FLOAT, (/nx_met_id,ny_met_id,nz_met_id,nt_met_id/), u_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(u_met_name))
  !
  if( nf90_def_var(ncID, v_met_name ,NF90_FLOAT, (/nx_met_id,ny_met_id,nz_met_id,nt_met_id/), v_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(v_met_name))
  !
  if( nf90_def_var(ncID, w_met_name ,NF90_FLOAT, (/nx_met_id,ny_met_id,nz_met_id,nt_met_id/), w_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(w_met_name))
  !
  if( nf90_def_var(ncID, T_met_name ,NF90_FLOAT, (/nx_met_id,ny_met_id,nz_met_id,nt_met_id/), T_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(T_met_name))
  !
  if( nf90_def_var(ncID, ro_met_name ,NF90_FLOAT, (/nx_met_id,ny_met_id,nz_met_id,nt_met_id/), ro_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(ro_met_name))
  !
  if( nf90_def_var(ncID, qv_met_name ,NF90_FLOAT, (/nx_met_id,ny_met_id,nz_met_id,nt_met_id/), qv_met_ID) /= NF90_NOERR ) &
       call runend('outmeteo: error in nf90_def_variable for variable '//TRIM(qv_met_name))
  !
  !*** Assign attribute values
  !
  attr_desc  = 'Longitud'
  attr_units = 'deg'
  if( nf90_put_att(ncID,lon_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for lon Variable')
  if( nf90_put_att(ncID, lon_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for lon Variable')
  !
  attr_desc  = 'Latitud'
  attr_units = 'deg'
  if( nf90_put_att(ncID,lat_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for lat Variable')
  if( nf90_put_att(ncID, lat_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for lat Variable')
  !
  attr_desc  = 'Layer height'
  attr_units = 'm'
  if( nf90_put_att(ncID,alt_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for lat Variable')
  if( nf90_put_att(ncID, alt_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for lat Variable')
  !
  attr_desc  = 'Topography'
  attr_units = 'm'
  if( nf90_put_att(ncID,hgt_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for hgt Variable')
  if( nf90_put_att(ncID, hgt_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for hgt Variable')
  !
  attr_desc  = 'Landmask. 1 for land. 0 for water'
  attr_units = ''
  if( nf90_put_att(ncID,land_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for land Variable')
  if( nf90_put_att(ncID, land_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for land Variable')
  !
  attr_desc  = 'Model type'
  attr_units = ''
  if( nf90_put_att(ncID,model_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for land Variable')
  if( nf90_put_att(ncID, model_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for land Variable')
  !
  attr_desc  = 'U. Horizontal component of wind'
  attr_units = ''
  if( nf90_put_att(ncID,u_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for U Variable')
  if( nf90_put_att(ncID, u_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for U Variable')
  !
  attr_desc  = 'V. Vertical component of wind'
  attr_units = ''
  if( nf90_put_att(ncID,v_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for V Variable')
  if( nf90_put_att(ncID, v_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for V Variable')
  !
  attr_desc  = 'W. z wind component'
  attr_units = ''
  if( nf90_put_att(ncID,w_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for W Variable')
  if( nf90_put_att(ncID, w_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for W Variable')
  !
  attr_desc  = 'T. Temperature'
  attr_units = ''
  if( nf90_put_att(ncID,T_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for T Variable')
  if( nf90_put_att(ncID, T_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for T Variable')
  !
  attr_desc  = 'ro. Air density'
  attr_units = ''
  if( nf90_put_att(ncID,ro_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for ro Variable')
  if( nf90_put_att(ncID, ro_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for ro Variable')
  !
  attr_desc  = 'qv. Air Humidity'
  attr_units = ''
  if( nf90_put_att(ncID,qv_met_ID , 'units', attr_units) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for qv Variable')
  if( nf90_put_att(ncID, qv_met_ID, 'description', attr_desc) /= NF90_NOERR ) &
       call runend('outmeteo : error in nf90_put_att for qv Variable')
  !
  !*** Leave define mode
  !
  if( nf90_enddef(ncid) /= NF90_NOERR ) call runend('outmeteo: error in nf90_enddef')
  !
  !*** fill values
  !
  if( nf90_put_var(ncID, lon_met_ID,grid%lon , start=(/1/), count=(/grid%nx/) ) /= NF90_NOERR ) &
       call wriwar('outmeteo : error in nf90_put_var for variable lon')
  !
  if( nf90_put_var(ncID, lat_met_ID,grid%lat , start=(/1/), count=(/grid%ny/) ) /= NF90_NOERR ) &
       call wriwar('outmeteo : error in nf90_put_var for variable lat')
  !
  if( nf90_put_var(ncID, alt_met_ID,grid%z , start=(/1/), count=(/grid%nz/) ) /= NF90_NOERR ) &
       call wriwar('outmeteo : error in nf90_put_var for variable alt')
  !
  if( nf90_put_var(ncID, hgt_met_ID,grid%hgt , start=(/1,1/), count=(/grid%nx,grid%ny/) ) /= NF90_NOERR ) &
       call wriwar('outmeteo : error in nf90_put_var for variable hgt')
  !
  if( nf90_put_var(ncID, land_met_ID,grid%landmask , start=(/1,1/), count=(/grid%nx,grid%ny/) ) /= NF90_NOERR ) &
       call wriwar('outmeteo : error in nf90_put_var for variable land')
  !
  if( nf90_put_var(ncID, model_met_ID,grid%model , start=(/1,1/), count=(/grid%nx,grid%ny/) ) /= NF90_NOERR ) &
       call wriwar('outmeteo : error in nf90_put_var for variable model')
  !
  !*** Close the file
  !
  if( nf90_close(ncID) /= NF90_NOERR ) &
       call runend('outmeteo : error closing netcdf file')
  !
  !*** Print Information in log file
  !
  !***  Opens and writes the log file
  !
  write(lulog,2)            'Created the Meteorological Output File'
  if(out_screen) write(*,2) 'Created the Meteorological Output File'
  !
2 format(/,  a /)
  return
end subroutine outmeteo
