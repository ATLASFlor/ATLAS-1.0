subroutine calcdif
   !******************************************************************************************
   !*
   !*    Calculate the tropopause.
   !*
   !*
   !*
   !*
   !*
   !******************************************************************************************
   use Master
   use InpOut
   use KindType
   implicit none
   save
   !
   integer(ip):: ipart !number of particle.
   real(rp) :: z, lon, lat ! particle parameters
   real(rp) :: altmin ! minimium height for tropopause
   !
   logical :: search, goon
   integer(ip) :: ix, iy, iz, jz, zmin

   do iy = 1, grid%ny
      do ix = 1, grid%nx

         ! Calculate height of thermal tropopause (Hoinka, 1997)
         !
         !*** 1. Set minimium height for tropopause
         !
         if ((grid%lat(iy) .ge. -20.) .and. (grid%lat(iy) .le. 20.)) then
            altmin = 5000.
         else if ((grid%lat(iy) .gt. 20.) .and. (grid%lat(iy) .lt. 40.)) then
            altmin = 2500.+(40.-lat)*125.
         else if ((grid%lat(iy) .gt. -40.) .and. (grid%lat(iy) .lt. -20.)) then
            altmin = 2500.+(40.+grid%lat(iy))*125.
         else
            altmin = 2500.
         end if
         !
         !*** 2. Define a minimum level zmin, from which upward the tropopause is
         !    searched for. This is to avoid inversions in the lower troposphere
         !    to be identified as the tropopause
         !
         search = .true.
         iz = 0
         do while (search)
            iz = iz + 1
            if (grid%z(iz) .ge. altmin) then
               zmin = iz
               search = .false.
            end if
         end do
         !
         !*** 3. Search for first stable layer above minimum height that fulfills the
         !    thermal tropopause criterion
         !
         search = .true.
         iz = zmin - 1
         do while (search)
            iz = iz + 1
            jz = iz
            goon = .true.
            do while (goon)
               jz = jz + 1
               if ((grid%z(jz) - grid%z(iz)) .gt. 2000.) then
                  if ((((grid%T(ix, iy, iz, 1) + grid%T(ix, iy, iz, 2))/2 - (grid%T(ix, iy, jz, 1) + grid%T(ix, iy, jz, 2))/2)/ &
                       (grid%z(jz) - grid%z(iz))) .lt. 0.002) then
                     grid%tropopause(ix, iy) = grid%z(iz)
                     search = .false.
                  end if
                  goon = .false.
               end if
               if (jz .gt. grid%nz + 1) goon = .false.
               !      write(lulog,1)"warning, not tropopause found. Minimium is established"
!1       !      format(/,  a /)
               !      grid%tropopause(ix,iy)=(zmin) !poner grid%z(zmin)
               !  end if
            end do
            if (iz .gt. grid%nz) then
               search = .false.
               if (my_id .eq. 0) write (lulog, 2) "warning, not tropopause found. Minimium is established"
2              format(/, a/)
               grid%tropopause(ix, iy) = grid%z(zmin)
            end if
         end do

!
      end do !ix
   end do    !iy

end subroutine calcdif
