program Grib2nc
  !**********************************************************************
  !*
  !*    AUTHOR : A.Folch
  !*    PURPOSE: This program converts an grib file to a netCDF file
  !*             to be used later by the SetDbs utility program
  !*
  !**********************************************************************
  use Master
  use InpOut
  implicit none
  !
  version = '7.1'
  !
  !***  Gets the call arguments
  !
  iarg = 1                       ! date YYMMDDHH
  call GETARG(iarg,date)
  iarg = 2                       ! log file
  call GETARG(iarg,lulogname)
  iarg = 3                       ! config name
  call GETARG(iarg,lucnfname)
  iarg = 4                       ! netcdf name
  call GETARG(iarg,luncname)
  iarg = 5                       ! Directory for decoded grib files
  call GETARG(iarg,lugrbbase)
  iarg = 6                       ! model name
  call GETARG(iarg,model)
  !
  !***  Opens the log file
  !
  call openinp
  !
  !***  Initializations (read the model config file)
  !
  call inival
  !
  !***  Writes the netCDF in define mode
  !
  call write_nc_grid
  !
  !***  Loop over time steps
  !
  do itt = 0,nt
     it = itt*dt
     !
     !***     Loop over 2d variables
     !
     do i2d = 1,n2d
        call write_nc_var2d
     end do
     !
     !***     Loop over 3d variables
     !
     do i3d = 1,n3d
        call write_nc_var3d
     end do
     !
  end do
  !
  !***  Ends the program
  !
  call runend('OK')
  !
end program Grib2nc
