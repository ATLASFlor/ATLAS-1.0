  subroutine init_outvar
     !******************************************************************************************
     !*
     !* Initialize output variables, allocate memory
     !*
     !******************************************************************************************
     !
     use InpOut
     use Master
     use KindType
     use Netcdf
     implicit none
     integer(ip) :: iphase
     integer(ip) :: ix, iz, iy
     integer(ip) :: ielem3d, ielem2d
     real(rp), allocatable:: lon(:), lat(:)
     !
     real(rp)             :: x_rad(2), y_rad(2)
     !
     !

     !
     !*** Previous computations
     !
     output%nx = INT((output%lonmax - output%lonmin)/output%dx) + 1 !variables defined in master module.
     output%ny = INT((output%latmax - output%latmin)/output%dy) + 1
     !
     output%nelem2d = (output%nx - 1)*(output%ny - 1)
     output%nelem3d = (output%nx - 1)*(output%ny - 1)*(output%nz - 1) !output%nz is defined in readinp.inp
     !
     !*** Allocates memory
     !
     allocate (conc3d(output%nelem3d))
     allocate (conc2d(output%nelem2d))
     allocate (sedim(output%nelem2d))
     !
     allocate (conc3d_global(output%nelem3d))
     allocate (conc2d_global(output%nelem2d))
     allocate (sedim_global(output%nelem2d))
     !
     allocate (output%area(output%nelem2d))
     allocate (output%volume(output%nelem3d))
     !
     allocate (concen(output%nx, output%ny, output%nz))
     allocate (colmass(output%nx, output%ny))
     allocate (grload(output%nx, output%ny))
     !
     if (output%classes) then
        !names bin: To generate output per classes
        do iphase = 1, nphases
           allocate (sedim_bin(nphases, phase(iphase)%bins, output%nelem2d))
           allocate (conc3d_bin(nphases, phase(iphase)%bins, output%nelem3d))
           allocate (conc2d_bin(nphases, phase(iphase)%bins, output%nelem2d))
           !
           allocate (grload_bin(nphases, phase(iphase)%bins, output%nx, output%ny))
           allocate (concen_bin(nphases, phase(iphase)%bins, output%nx, output%ny, output%nz))
           allocate (colmass_bin(nphases, phase(iphase)%bins, output%nx, output%ny))
           !
           allocate (sedim_id_bin(nphases, phase(iphase)%bins))
           allocate (conc2d_id_bin(nphases, phase(iphase)%bins))
           allocate (conc3d_id_bin(nphases, phase(iphase)%bins))
        end do
        !
     end if
     if (output%phases) then
        !names bin: To generate output per classes
        allocate (sedim_phase(nphases, output%nelem2d))
        allocate (conc3d_phase(nphases, output%nelem3d))
        allocate (conc2d_phase(nphases, output%nelem2d))
        !
        allocate (grload_phase(nphases, output%nx, output%ny))
        allocate (concen_phase(nphases, output%nx, output%ny, output%nz))
        allocate (colmass_phase(nphases, output%nx, output%ny))
        !ID
        allocate (sedim_id_phase(nphases))
        allocate (conc2d_id_phase(nphases))
        allocate (conc3d_id_phase(nphases))
     end if
     if (output%track_points) then
        allocate (load_pts(npts))
        allocate (thick_pts(npts))
        allocate (loadpts_global(npts))
        allocate (thickpts_global(npts))
        loadpts_global=0.0
     end if

     !
     !*** Allocates memory
     !
     allocate (lon(output%nx))
     allocate (lat(output%ny))
     !
     !*** Creates the mesh stored as 1D vectors, with lon/lat values of the square centre.
     !
     do ix = 1, output%nx
        lon(ix) = output%lonmin + (ix - 1./2.)*output%dx      !!!IGUAL QUE EN INIT NETCDF... RESTAR 1 PARA MODIFICAR
     end do
     !
     do iy = 1, output%ny
        lat(iy) = output%latmin + (iy - 1./2.)*output%dy
     end do
     !
     !*** Calculates the surface area of each 2D grid cell.and the volume of each element(3D)
     !*** First convert degrees in radians
     !
     do iz = 1, output%nz - 1
        do iy = 1, output%ny - 1
           do ix = 1, output%nx - 1
              ielem2d = (iy - 1)*(output%nx - 1) + ix
              x_rad(1) = lon(ix)*pi/180
              x_rad(2) = lon(ix + 1)*pi/180
              y_rad(1) = lat(iy)*pi/180
              y_rad(2) = lat(iy + 1)*pi/180
              output%area(ielem2d) = (Rearth)*(Rearth)*ABS(x_rad(2) - x_rad(1))*ABS(sin(y_rad(2)) - sin(y_rad(1))) !m2
!           output%area(ielem2d)=(Rearth/1000)*(Rearth/1000)*ABS(x_rad(2)-x_rad(1))*ABS(sin(y_rad(2))-sin(y_rad(1))) !km2
              ielem3d = (iz - 1)*(output%nx - 1)*(output%ny - 1) + (iy - 1)*(output%nx - 1) + ix
              output%volume(ielem3d) = output%area(ielem2d)*ABS(output%zlayer(iz + 1) - output%zlayer(iz)) ! in m3
           end do
        end do
!
!      output%volume(ielem3d)=output%area(ielem2d)*ABS(output%zlayer(iz+1)-output%zlayer(iz))/1000 ! in km3
     end do !iz
!
!     do ix = 1,output%nx
!        do iy = 1,output%ny
!           ielem2d = (iy-1)*output%nx+ix
!           x_rad(1)=lon(ix)*pi/180
!           x_rad(2)=lon(ix+1)*pi/180
!           y_rad(1)=lat(iy)*pi/180
!           y_rad(2)=lat(iy+1)*pi/180
!           output%area(ielem2d)=(Rearth)*(Rearth)*ABS(x_rad(2)-x_rad(1))*ABS(sin(y_rad(2))-sin(y_rad(1))) !m2
!           output%area(ielem2d)=(Rearth/1000)*(Rearth/1000)*ABS(x_rad(2)-x_rad(1))*ABS(sin(y_rad(2))-sin(y_rad(1))) !km2
!           do iz=1,output%nz
!              ielem3d= (iz-1)*(output%nx)*(output%ny)+(iy-1)*(output%nx)+ix
!              output%volume(ielem3d)=output%area(ielem2d)*ABS(output%zlayer(iz+1)-output%zlayer(iz)) ! in m3
!              output%volume(ielem3d)=output%area(ielem2d)*ABS(output%zlayer(iz+1)-output%zlayer(iz))/1000 ! in km3
!           end do !iz
!        end do !iy
!     end do !ix
     deallocate (lon)
     deallocate (lat)

     return
  end subroutine init_outvar
