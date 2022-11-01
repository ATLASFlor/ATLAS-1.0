subroutine set_back_mesh
   !**************************************************************
   !*
   !*    Create the Background Mesh
   !*
   !**************************************************************
   use KindType
   use InpOut
   use Master
   implicit none
   !
   integer(ip) :: ix, iy, iz, ipoin, ielem
   !
   !*** Previous computations
   !
   grid%nx = INT((grid%lonmax - grid%lonmin)/grid%dx) + 1
   grid%ny = INT((grid%latmax - grid%latmin)/grid%dy) + 1
   grid%nz = INT(grid%ztop/grid%dz) + 1
   grid%nt = INT((siend - sist)/sidt) + 1
   !
   !*** Creates the mesh stored as 1D vectors
   !
   allocate (grid%lon(grid%nx))
   allocate (grid%lat(grid%ny))
   allocate (grid%z(grid%nz))
   !
   do ix = 1, grid%nx
      grid%lon(ix) = grid%lonmin + (ix - 1)*grid%dx
   end do
   !
   do iy = 1, grid%ny
      grid%lat(iy) = grid%latmin + (iy - 1)*grid%dy
   end do
   !
   do iz = 1, grid%nz
      grid%z(iz) = (iz - 1)*grid%dz
   end do
   !
   !*** Creates the mesh as a 3D structure (more general, allows irregular meshes in a future)
   !
   grid%npoin = grid%nx*grid%ny*grid%nz
   grid%nelem = (grid%nx - 1)*(grid%ny - 1)*(grid%nz - 1)
   grid%npoin2d = grid%nx*grid%ny
   !
   !*** Coordinates
   !
   allocate (grid%coord(grid%ndime, grid%npoin))
   ipoin = 0
   do iz = 1, grid%nz
      do iy = 1, grid%ny
         do ix = 1, grid%nx
            ipoin = ipoin + 1
            grid%coord(1, ipoin) = grid%lon(ix)
            grid%coord(2, ipoin) = grid%lat(iy)
            grid%coord(3, ipoin) = grid%z(iz)
         end do
      end do
   end do
   !
   !*** Creates nodal connectivities lnods
   !
   allocate (grid%lnods(grid%nnode, grid%nelem))
   ielem = 0
   do iz = 1, grid%nz - 1
      do iy = 1, grid%ny - 1
         do ix = 1, grid%nx - 1
            ielem = ielem + 1
            grid%lnods(1, ielem) = (iz - 1)*grid%nx*grid%ny + (iy - 1)*grid%nx + ix
            grid%lnods(2, ielem) = (iz - 1)*grid%nx*grid%ny + (iy - 1)*grid%nx + ix + 1
            grid%lnods(3, ielem) = (iz - 1)*grid%nx*grid%ny + (iy)*grid%nx + ix + 1
            grid%lnods(4, ielem) = (iz - 1)*grid%nx*grid%ny + (iy)*grid%nx + ix
            grid%lnods(5, ielem) = (iz)*grid%nx*grid%ny + (iy - 1)*grid%nx + ix
            grid%lnods(6, ielem) = (iz)*grid%nx*grid%ny + (iy - 1)*grid%nx + ix + 1
            grid%lnods(7, ielem) = (iz)*grid%nx*grid%ny + (iy)*grid%nx + ix + 1
            grid%lnods(8, ielem) = (iz)*grid%nx*grid%ny + (iy)*grid%nx + ix
         end do
      end do
   end do
   !
   !*** Allocates memory for meteo interpolation. Note that only 2D variables
   !*** are needed for the interpolation
   !
   allocate (grid%model(grid%npoin2d))
   grid%model(:) = -1                 ! default value. -1 means no model
   allocate (grid%element(grid%npoin2d))
   grid%element(:) = -1                 ! default value. -1 means no model
   allocate (grid%shape(4, grid%npoin2d))
   !
   !*** Allocates memory for meteo variables
   !
   allocate (grid%hgt(grid%nx, grid%ny))    ! Topography
   allocate (grid%landmask(grid%nx, grid%ny))    ! Landmask
   !
   allocate (grid%u(grid%nx, grid%ny, grid%nz, 2))
   allocate (grid%v(grid%nx, grid%ny, grid%nz, 2))
   allocate (grid%w(grid%nx, grid%ny, grid%nz, 2))
   allocate (grid%T(grid%nx, grid%ny, grid%nz, 2))
   allocate (grid%ro(grid%nx, grid%ny, grid%nz, 2))
   allocate (grid%qv(grid%nx, grid%ny, grid%nz, 2))
   allocate (grid%tp(grid%nx, grid%ny, grid%nz, 2))
   allocate (grid%tropopause(grid%nx, grid%ny))
   allocate (grid%pblh(grid%nx, grid%ny, 2))
   allocate (grid%hfx(grid%nx, grid%ny, 2))
   allocate (grid%ust(grid%nx, grid%ny, 2))
   allocate (grid%rmonin(grid%nx, grid%ny, 2))
   allocate (grid%wst(grid%nx, grid%ny, 2))
   allocate (grid%drhodz(grid%nx, grid%ny, grid%nz, 2))
   !
   !*** Print Information in log file
   !*** Writes the log file
   !
   if (my_id .eq. 0) write (lulog, 2) 'Created background mesh', 'Allocated memory for meteorological interpolation.'
   if (out_screen) write (*, 2) 'Created background mesh', 'Allocated memory for meteorological interpolation.'
   !
2  format(/, &
           a, /, &
           a/)
   !
   return
end subroutine set_back_mesh
