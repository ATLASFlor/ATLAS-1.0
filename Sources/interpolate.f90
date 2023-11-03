   subroutine interpolate_meteo(ipart, lon, lat, z)
      !******************************************************************************************
      !*
      !*    Interpolates meteo variables at the particle position
      !*
      !******************************************************************************************
      use Master
      use KindType
      implicit none
      !
      integer(ip)  :: i, ipart
      integer(ip)  :: ix1, ix2, iy1, iy2, iz1, iz2
      integer(ip)  :: ipoin1, ipoin2, ipoin3, ipoin4
      integer(ip)  :: nm1, nm2, nm3, nm4
      real(rp)  :: xshape, yshape, zshape
      real(rp)  :: lon, lat, z
      !
      real(rp)  :: ua(8)
      real(rp)  :: va(8)
      real(rp)  :: wa(8)
      real(rp)  :: Ta(8)
      real(rp)  :: rhoa(8)
      real(rp)  :: drhodz(8)
      real(rp)  :: qva(8)
      real(rp)  :: shapef(8)
      real(rp)  :: tropop(4)
      real(rp)  :: pblha(4)
      real(rp)  :: hfxa(4)
      real(rp)  :: usta(4)
      real(rp)  :: rmonin(4)
      real(rp)  :: wst(4)
      real(rp)  :: shapef2(4)
      real(rp)  :: mua0, Ta0
      !
      real(rp)::tshape1, tshape2, tshape3, tshape4
      !
      !*** Calculate indexes of the 8 nodes of the background mesh
      !
      ix1 = INT((lon - grid%lonmin)/grid%dx) + 1
      ix2 = ix1 + 1
      !
      iy1 = INT((lat - grid%latmin)/grid%dy) + 1
      iy2 = iy1 + 1
      !
      iz1 = INT(z/grid%dz) + 1
      iz2 = iz1 + 1
      !
      !*** Calculate the shape factor
      !
      xshape = 2.0_rp*((lon - grid%lon(ix1))/(grid%lon(ix2) - grid%lon(ix1))) - 1.0_rp   ! in [-1,1]
      yshape = 2.0_rp*((lat - grid%lat(iy1))/(grid%lat(iy2) - grid%lat(iy1))) - 1.0_rp
      zshape = 2.0_rp*((z - grid%z(iz1))/(grid%z(iz2) - grid%z(iz1))) - 1.0_rp
      call shape3(xshape, yshape, zshape, shapef)
      !
      !*** Inquire the meteo model of each point of the cell
      !
      ipoin1 = (iy1 - 1)*grid%nx + ix1
      ipoin2 = ipoin1 + 1
      ipoin3 = iy1*grid%nx + ix2
      ipoin4 = ipoin3 - 1
      nm1 = grid%model(ipoin1)
      nm2 = grid%model(ipoin2)
      nm3 = grid%model(ipoin3)
      nm4 = grid%model(ipoin4)
      tshape1 = meteo(nm1)%time_s
      tshape2 = meteo(nm2)%time_s
      tshape3 = meteo(nm3)%time_s
      tshape4 = meteo(nm4)%time_s
      !
      !*** Interpolate the Variables IN TIME at each of the 8 nodes
      !
      !U
      ua(1) = tshape1*grid%u(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%u(ix1, iy1, iz1, 1)
      ua(2) = tshape2*grid%u(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%u(ix2, iy1, iz1, 1)
      ua(3) = tshape3*grid%u(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%u(ix2, iy2, iz1, 1)
      ua(4) = tshape4*grid%u(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%u(ix1, iy2, iz1, 1)
      ua(5) = tshape1*grid%u(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%u(ix1, iy1, iz2, 1)
      ua(6) = tshape2*grid%u(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%u(ix2, iy1, iz2, 1)
      ua(7) = tshape3*grid%u(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%u(ix2, iy2, iz2, 1)
      ua(8) = tshape4*grid%u(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%u(ix1, iy2, iz2, 1)
      !V
      va(1) = tshape1*grid%v(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%v(ix1, iy1, iz1, 1)
      va(2) = tshape2*grid%v(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%v(ix2, iy1, iz1, 1)
      va(3) = tshape3*grid%v(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%v(ix2, iy2, iz1, 1)
      va(4) = tshape4*grid%v(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%v(ix1, iy2, iz1, 1)
      va(5) = tshape1*grid%v(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%v(ix1, iy1, iz2, 1)
      va(6) = tshape2*grid%v(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%v(ix2, iy1, iz2, 1)
      va(7) = tshape3*grid%v(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%v(ix2, iy2, iz2, 1)
      va(8) = tshape4*grid%v(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%v(ix1, iy2, iz2, 1)
      !W
      wa(1) = tshape1*grid%w(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%w(ix1, iy1, iz1, 1)
      wa(2) = tshape2*grid%w(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%w(ix2, iy1, iz1, 1)
      wa(3) = tshape3*grid%w(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%w(ix2, iy2, iz1, 1)
      wa(4) = tshape4*grid%w(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%w(ix1, iy2, iz1, 1)
      wa(5) = tshape1*grid%w(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%w(ix1, iy1, iz2, 1)
      wa(6) = tshape2*grid%w(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%w(ix2, iy1, iz2, 1)
      wa(7) = tshape3*grid%w(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%w(ix2, iy2, iz2, 1)
      wa(8) = tshape4*grid%w(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%w(ix1, iy2, iz2, 1)
      !T
      Ta(1) = tshape1*grid%T(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%T(ix1, iy1, iz1, 1)
      Ta(2) = tshape2*grid%T(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%T(ix2, iy1, iz1, 1)
      Ta(3) = tshape3*grid%T(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%T(ix2, iy2, iz1, 1)
      Ta(4) = tshape4*grid%T(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%T(ix1, iy2, iz1, 1)
      Ta(5) = tshape1*grid%T(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%T(ix1, iy1, iz2, 1)
      Ta(6) = tshape2*grid%T(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%T(ix2, iy1, iz2, 1)
      Ta(7) = tshape3*grid%T(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%T(ix2, iy2, iz2, 1)
      Ta(8) = tshape4*grid%T(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%T(ix1, iy2, iz2, 1)
      !ro: density
      rhoa(1) = tshape1*grid%ro(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%ro(ix1, iy1, iz1, 1)
      rhoa(2) = tshape2*grid%ro(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%ro(ix2, iy1, iz1, 1)
      rhoa(3) = tshape3*grid%ro(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%ro(ix2, iy2, iz1, 1)
      rhoa(4) = tshape4*grid%ro(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%ro(ix1, iy2, iz1, 1)
      rhoa(5) = tshape1*grid%ro(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%ro(ix1, iy1, iz2, 1)
      rhoa(6) = tshape2*grid%ro(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%ro(ix2, iy1, iz2, 1)
      rhoa(7) = tshape3*grid%ro(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%ro(ix2, iy2, iz2, 1)
      rhoa(8) = tshape4*grid%ro(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%ro(ix1, iy2, iz2, 1)
      !drhodz: air density vertical gradient
      drhodz(1) = tshape1*grid%drhodz(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%drhodz(ix1, iy1, iz1, 1)
      drhodz(2) = tshape2*grid%drhodz(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%drhodz(ix2, iy1, iz1, 1)
      drhodz(3) = tshape3*grid%drhodz(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%drhodz(ix2, iy2, iz1, 1)
      drhodz(4) = tshape4*grid%drhodz(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%drhodz(ix1, iy2, iz1, 1)
      drhodz(5) = tshape1*grid%drhodz(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%drhodz(ix1, iy1, iz2, 1)
      drhodz(6) = tshape2*grid%drhodz(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%drhodz(ix2, iy1, iz2, 1)
      drhodz(7) = tshape3*grid%drhodz(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%drhodz(ix2, iy2, iz2, 1)
      drhodz(8) = tshape4*grid%drhodz(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%drhodz(ix1, iy2, iz2, 1)
      !qv: specific humidity
      qva(1) = tshape1*grid%qv(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%qv(ix1, iy1, iz1, 1)
      qva(2) = tshape2*grid%qv(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%qv(ix2, iy1, iz1, 1)
      qva(3) = tshape3*grid%qv(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%qv(ix2, iy2, iz1, 1)
      qva(4) = tshape4*grid%qv(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%qv(ix1, iy2, iz1, 1)
      qva(5) = tshape1*grid%qv(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%qv(ix1, iy1, iz2, 1)
      qva(6) = tshape2*grid%qv(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%qv(ix2, iy1, iz2, 1)
      qva(7) = tshape3*grid%qv(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%qv(ix2, iy2, iz2, 1)
      qva(8) = tshape4*grid%qv(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%qv(ix1, iy2, iz2, 1)
      !
      !*** Interpolate the 2D Variables IN TIME at each of the 4 nodes
      !
      !PBLH
      pblha(1) = tshape1*grid%pblh(ix1, iy1, 2) + (1.0_rp - tshape1)*grid%pblh(ix1, iy1, 1)
      pblha(2) = tshape2*grid%pblh(ix2, iy1, 2) + (1.0_rp - tshape2)*grid%pblh(ix2, iy1, 1)
      pblha(3) = tshape3*grid%pblh(ix2, iy2, 2) + (1.0_rp - tshape3)*grid%pblh(ix2, iy2, 1)
      pblha(4) = tshape4*grid%pblh(ix1, iy2, 2) + (1.0_rp - tshape4)*grid%pblh(ix1, iy2, 1)
      !HFX
      hfxa(1) = tshape1*grid%hfx(ix1, iy1, 2) + (1.0_rp - tshape1)*grid%hfx(ix1, iy1, 1)
      hfxa(2) = tshape2*grid%hfx(ix2, iy1, 2) + (1.0_rp - tshape2)*grid%hfx(ix2, iy1, 1)
      hfxa(3) = tshape3*grid%hfx(ix2, iy2, 2) + (1.0_rp - tshape3)*grid%hfx(ix2, iy2, 1)
      hfxa(4) = tshape4*grid%hfx(ix1, iy2, 2) + (1.0_rp - tshape4)*grid%hfx(ix1, iy2, 1)
      !u*
      usta(1) = tshape1*grid%ust(ix1, iy1, 2) + (1.0_rp - tshape1)*grid%ust(ix1, iy1, 1)
      usta(2) = tshape2*grid%ust(ix2, iy1, 2) + (1.0_rp - tshape2)*grid%ust(ix2, iy1, 1)
      usta(3) = tshape3*grid%ust(ix2, iy2, 2) + (1.0_rp - tshape3)*grid%ust(ix2, iy2, 1)
      usta(4) = tshape4*grid%ust(ix1, iy2, 2) + (1.0_rp - tshape4)*grid%ust(ix1, iy2, 1)
      !rmonin
      rmonin(1) = tshape1*grid%rmonin(ix1, iy1, 2) + (1.0_rp - tshape1)*grid%rmonin(ix1, iy1, 1)
      rmonin(2) = tshape2*grid%rmonin(ix2, iy1, 2) + (1.0_rp - tshape2)*grid%rmonin(ix2, iy1, 1)
      rmonin(3) = tshape3*grid%rmonin(ix2, iy2, 2) + (1.0_rp - tshape3)*grid%rmonin(ix2, iy2, 1)
      rmonin(4) = tshape4*grid%rmonin(ix1, iy2, 2) + (1.0_rp - tshape4)*grid%rmonin(ix1, iy2, 1)
      !wst
      wst(1) = tshape1*grid%wst(ix1, iy1, 2) + (1.0_rp - tshape1)*grid%wst(ix1, iy1, 1)
      wst(2) = tshape2*grid%wst(ix2, iy1, 2) + (1.0_rp - tshape2)*grid%wst(ix2, iy1, 1)
      wst(3) = tshape3*grid%wst(ix2, iy2, 2) + (1.0_rp - tshape3)*grid%wst(ix2, iy2, 1)
      wst(4) = tshape4*grid%wst(ix1, iy2, 2) + (1.0_rp - tshape4)*grid%wst(ix1, iy2, 1)
      !
      !*** Interpolate variables at (x,y,z) position
      !
      part(ipart)%ua = 0.0_rp
      part(ipart)%va = 0.0_rp
      part(ipart)%wa = 0.0_rp
      part(ipart)%Ta = 0.0_rp
      part(ipart)%rhoa = 0.0_rp
      part(ipart)%drhodz = 0.0_rp
      part(ipart)%qva = 0.0_rp
      do i = 1, 8
         part(ipart)%ua = part(ipart)%ua + shapef(i)*ua(i)
         part(ipart)%va = part(ipart)%va + shapef(i)*va(i)
         part(ipart)%wa = part(ipart)%wa + shapef(i)*wa(i)
         part(ipart)%Ta = part(ipart)%Ta + shapef(i)*Ta(i)
         part(ipart)%rhoa = part(ipart)%rhoa + shapef(i)*rhoa(i)
         part(ipart)%drhodz = part(ipart)%drhodz + shapef(i)*drhodz(i)
         part(ipart)%qva = part(ipart)%qva + shapef(i)*qva(i)
      end do
      !
      !***  Other variables
      !
      ! mua
      !
      mua0 = 1.827e-5  ! reference viscosity
      Ta0 = 291.15    ! reference temperature
      part(ipart)%mua = mua0*((Ta0 + 120.0_rp)/(part(ipart)%Ta + 120.0_rp))*((part(ipart)%Ta/Ta0)**1.5_rp)
      !
      !*** Interpolate 2D variables at particle position
      !
      call shape2(xshape, yshape, shapef2)
      part(ipart)%pblha = 0.0_rp
      part(ipart)%hfxa = 0.0_rp
      part(ipart)%usta = 0.0_rp
      part(ipart)%rmonin = 0.0_rp
      part(ipart)%wst = 0.0_rp
      tropopause = 0.0_rp
      tropop(1) = grid%tropopause(ix1, iy1)
      tropop(2) = grid%tropopause(ix2, iy1)
      tropop(3) = grid%tropopause(ix2, iy2)
      tropop(4) = grid%tropopause(ix1, iy2)
      do i = 1, 4
         tropopause = tropopause + shapef2(i)*tropop(i)
         part(ipart)%pblha = part(ipart)%pblha + shapef2(i)*pblha(i)
         part(ipart)%hfxa = part(ipart)%hfxa + shapef2(i)*hfxa(i)
         part(ipart)%usta = part(ipart)%usta + shapef2(i)*usta(i)
         part(ipart)%rmonin = part(ipart)%rmonin + shapef2(i)*rmonin(i)
         part(ipart)%wst = part(ipart)%wst + shapef2(i)*wst(i)
      end do
      !
      return
   end subroutine interpolate_meteo
   !
   !
   subroutine interpolate_wind(lon, lat, z, uab, vab, wab, rhoab, muab)
      !******************************************************************************************
      !*
      !*    Interpolates meteo variables at the particle position, to estimate the advection
      !*    INPUTS:lon, lat,z firsts position approximation
      !*    OUTPUTS: uab,vab,wab,rhoab,muab calculated to estimate new advection iteration step
      !*
      !******************************************************************************************
      use Master
      use KindType
      implicit none
      !
      integer(ip)  :: i
      integer(ip)  :: ix1, ix2, iy1, iy2, iz1, iz2
      integer(ip)  :: ipoin1, ipoin2, ipoin3, ipoin4
      integer(ip)  :: nm1, nm2, nm3, nm4
      real(rp)  :: xshape, yshape, zshape
      real(rp)  :: lon, lat, z
      real(rp)  :: uab, vab, wab !winds in second position
      real(rp)  :: rhoab, muab, Tab
      real(rp)  :: Ta0, mua0
      !
      real(rp)  :: ub(8)
      real(rp)  :: vb(8)
      real(rp)  :: wb(8)
      real(rp)  :: rhoa(8)
      real(rp)  :: Ta(8)
      real(rp)  :: shapef(8)
      real(rp)  :: shapef2(4)
      real(rp)  :: tshape1, tshape2, tshape3, tshape4
      !
      !*** Calculate indexes of the 8 nodes of the background mesh
      !
      ix1 = INT((lon - grid%lonmin)/grid%dx) + 1
      ix2 = ix1 + 1
      if (ix2 .gt. grid%nx) then
         ix2 = ix1
         ix1 = ix1 - 1
      end if
      !
      iy1 = INT((lat - grid%latmin)/grid%dy) + 1
      iy2 = iy1 + 1
      if (iy2 .gt. grid%ny) then
         iy2 = iy1
         iy1 = iy1 - 1
      end if
      !
      iz1 = INT(z/grid%dz) + 1
      iz2 = iz1 + 1
      if (iz2 .gt. grid%nz) then
         iz2 = iz1
         iz1 = iz1 - 1
      end if
      !
      !*** Calculate the shape factor
      !
      xshape = 2.0_rp*((lon - grid%lon(ix1))/(grid%lon(ix2) - grid%lon(ix1))) - 1.0_rp   ! in [-1,1]
      yshape = 2.0_rp*((lat - grid%lat(iy1))/(grid%lat(iy2) - grid%lat(iy1))) - 1.0_rp
      zshape = 2.0_rp*((z - grid%z(iz1))/(grid%z(iz2) - grid%z(iz1))) - 1.0_rp
      call shape3(xshape, yshape, zshape, shapef)
      !
      !*** Inquire the meteo model of each point of the cell
      !
      ipoin1 = (iy1 - 1)*grid%nx + ix1
      ipoin2 = ipoin1 + 1
      ipoin3 = iy1*grid%nx + ix2
      ipoin4 = ipoin3 - 1
      nm1 = grid%model(ipoin1)
      nm2 = grid%model(ipoin2)
      nm3 = grid%model(ipoin3)
      nm4 = grid%model(ipoin4)
      tshape1 = meteo(nm1)%time_s
      tshape2 = meteo(nm2)%time_s
      tshape3 = meteo(nm3)%time_s
      tshape4 = meteo(nm4)%time_s
      !
      !*** Interpolate the Variables IN TIME at each of the 8 nodes
      !
      !U
      ub(1) = tshape1*grid%u(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%u(ix1, iy1, iz1, 1)
      ub(2) = tshape2*grid%u(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%u(ix2, iy1, iz1, 1)
      ub(3) = tshape3*grid%u(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%u(ix2, iy2, iz1, 1)
      ub(4) = tshape4*grid%u(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%u(ix1, iy2, iz1, 1)
      ub(5) = tshape1*grid%u(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%u(ix1, iy1, iz2, 1)
      ub(6) = tshape2*grid%u(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%u(ix2, iy1, iz2, 1)
      ub(7) = tshape3*grid%u(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%u(ix2, iy2, iz2, 1)
      ub(8) = tshape4*grid%u(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%u(ix1, iy2, iz2, 1)
      !V
      vb(1) = tshape1*grid%v(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%v(ix1, iy1, iz1, 1)
      vb(2) = tshape2*grid%v(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%v(ix2, iy1, iz1, 1)
      vb(3) = tshape3*grid%v(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%v(ix2, iy2, iz1, 1)
      vb(4) = tshape4*grid%v(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%v(ix1, iy2, iz1, 1)
      vb(5) = tshape1*grid%v(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%v(ix1, iy1, iz2, 1)
      vb(6) = tshape2*grid%v(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%v(ix2, iy1, iz2, 1)
      vb(7) = tshape3*grid%v(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%v(ix2, iy2, iz2, 1)
      vb(8) = tshape4*grid%v(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%v(ix1, iy2, iz2, 1)
      !W
      wb(1) = tshape1*grid%w(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%w(ix1, iy1, iz1, 1)
      wb(2) = tshape2*grid%w(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%w(ix2, iy1, iz1, 1)
      wb(3) = tshape3*grid%w(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%w(ix2, iy2, iz1, 1)
      wb(4) = tshape4*grid%w(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%w(ix1, iy2, iz1, 1)
      wb(5) = tshape1*grid%w(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%w(ix1, iy1, iz2, 1)
      wb(6) = tshape2*grid%w(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%w(ix2, iy1, iz2, 1)
      wb(7) = tshape3*grid%w(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%w(ix2, iy2, iz2, 1)
      wb(8) = tshape4*grid%w(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%w(ix1, iy2, iz2, 1)
      !T
      Ta(1) = tshape1*grid%T(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%T(ix1, iy1, iz1, 1)
      Ta(2) = tshape2*grid%T(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%T(ix2, iy1, iz1, 1)
      Ta(3) = tshape3*grid%T(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%T(ix2, iy2, iz1, 1)
      Ta(4) = tshape4*grid%T(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%T(ix1, iy2, iz1, 1)
      Ta(5) = tshape1*grid%T(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%T(ix1, iy1, iz2, 1)
      Ta(6) = tshape2*grid%T(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%T(ix2, iy1, iz2, 1)
      Ta(7) = tshape3*grid%T(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%T(ix2, iy2, iz2, 1)
      Ta(8) = tshape4*grid%T(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%T(ix1, iy2, iz2, 1)
      !ro: density
      rhoa(1) = tshape1*grid%ro(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%ro(ix1, iy1, iz1, 1)
      rhoa(2) = tshape2*grid%ro(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%ro(ix2, iy1, iz1, 1)
      rhoa(3) = tshape3*grid%ro(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%ro(ix2, iy2, iz1, 1)
      rhoa(4) = tshape4*grid%ro(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%ro(ix1, iy2, iz1, 1)
      rhoa(5) = tshape1*grid%ro(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%ro(ix1, iy1, iz2, 1)
      rhoa(6) = tshape2*grid%ro(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%ro(ix2, iy1, iz2, 1)
      rhoa(7) = tshape3*grid%ro(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%ro(ix2, iy2, iz2, 1)
      rhoa(8) = tshape4*grid%ro(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%ro(ix1, iy2, iz2, 1)
      !
      !*** Interpolate variables at (x,y,z) position
      !
      uab = 0.0_rp
      vab = 0.0_rp
      wab = 0.0_rp
      Tab = 0.0_rp
      rhoab = 0.0_rp
      do i = 1, 8
         uab = uab + shapef(i)*ub(i)
         vab = vab + shapef(i)*vb(i)
         wab = wab + shapef(i)*wb(i)
         Tab = Tab + shapef(i)*Ta(i)
         rhoab = rhoab + shapef(i)*rhoa(i)
      end do
      !
      !
      !***  Other variables
      !
      ! mua
      !
      mua0 = 1.827e-5  ! reference viscosity
      Ta0 = 291.15    ! reference temperature
      muab = mua0*((Ta0 + 120.0_rp)/(Tab + 120.0_rp))*((Tab/Ta0)**1.5_rp)
      return
   end subroutine interpolate_wind
   !
   !
   subroutine interpolate_sm(lon, lat, hgt)
      !******************************************************************************************
      !*
      !*    Interpolates meteo variables at the Source position
      !*    Used in source determination.
      !*
      !******************************************************************************************
      use Master
      use KindType
      use Sourcemod
      implicit none
      !
      real(rp)  :: lon, lat, hgt
      integer(ip)  :: ix1, ix2, iy1, iy2, iz1, iz2, iz, i
      integer(ip)  :: ipoin1, ipoin2, ipoin3, ipoin4
      integer(ip)  :: nm1, nm2, nm3, nm4
      real(rp)  :: xshape, yshape, zshape
      !
      real(rp)  :: shapef(8)
      real(rp)  :: ua(8)
      real(rp)  :: va(8)
      real(rp)  :: Ta(8)
      real(rp)  :: rhoa(8)
      !real   (rp)  :: qva (8)
      real(rp)  ::tshape1, tshape2, tshape3, tshape4
      !
      real(rp)  :: ux, uy, angle, z, Cao, dTdz
      !
      !*** Calculate indexes of the 8 Vent nodes of the background mesh
      !
      ix1 = INT((lon - grid%lonmin)/grid%dx) + 1
      ix2 = ix1 + 1
      !
      iy1 = INT((lat - grid%latmin)/grid%dy) + 1
      iy2 = iy1 + 1
      !
      !*** Calculate the shape factor
      !
      xshape = 2.0_rp*((lon - grid%lon(ix1))/(grid%lon(ix2) - grid%lon(ix1))) - 1.0_rp   ! in [-1,1]
      yshape = 2.0_rp*((lat - grid%lat(iy1))/(grid%lat(iy2) - grid%lat(iy1))) - 1.0_rp
      !
      !*** Inquire the meteo model of each point of the cell
      !
      ipoin1 = (iy1 - 1)*grid%nx + ix1
      ipoin2 = ipoin1 + 1
      ipoin3 = iy1*grid%nx + ix2
      ipoin4 = ipoin3 - 1
      nm1 = grid%model(ipoin1)
      nm2 = grid%model(ipoin2)
      nm3 = grid%model(ipoin3)
      nm4 = grid%model(ipoin4)

      tshape1 = meteo(nm1)%time_s
      tshape2 = meteo(nm2)%time_s
      tshape3 = meteo(nm3)%time_s
      tshape4 = meteo(nm4)%time_s
      !
      z = 0.0
      PSIa = 0.0
      do iz = 1, grid%nz
         ! z- indexes
         !
         iz1 = INT(hgt/grid%dz) + 1
         iz2 = iz1 + 1
         ! z shape factor
         zshape = 2.0_rp*((hgt - grid%z(iz1))/(grid%z(iz2) - grid%z(iz1))) - 1.0_rp
         call shape3(xshape, yshape, zshape, shapef)
         !
         !*** Interpolate the Variables in time at each of the 8 nodes
         !
         !U
         ua(1) = tshape1*grid%u(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%u(ix1, iy1, iz1, 1)
         ua(2) = tshape2*grid%u(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%u(ix2, iy1, iz1, 1)
         ua(3) = tshape3*grid%u(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%u(ix2, iy2, iz1, 1)
         ua(4) = tshape4*grid%u(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%u(ix1, iy2, iz1, 1)
         ua(5) = tshape1*grid%u(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%u(ix1, iy1, iz2, 1)
         ua(6) = tshape2*grid%u(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%u(ix2, iy1, iz2, 1)
         ua(7) = tshape3*grid%u(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%u(ix2, iy2, iz2, 1)
         ua(8) = tshape4*grid%u(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%u(ix1, iy2, iz2, 1)
         !
         !V
         va(1) = tshape1*grid%v(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%v(ix1, iy1, iz1, 1)
         va(2) = tshape2*grid%v(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%v(ix2, iy1, iz1, 1)
         va(3) = tshape3*grid%v(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%v(ix2, iy2, iz1, 1)
         va(4) = tshape4*grid%v(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%v(ix1, iy2, iz1, 1)
         va(5) = tshape1*grid%v(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%v(ix1, iy1, iz2, 1)
         va(6) = tshape2*grid%v(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%v(ix2, iy1, iz2, 1)
         va(7) = tshape3*grid%v(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%v(ix2, iy2, iz2, 1)
         va(8) = tshape4*grid%v(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%v(ix1, iy2, iz2, 1)
         !
         !T
         Ta(1) = tshape1*grid%T(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%T(ix1, iy1, iz1, 1)
         Ta(2) = tshape2*grid%T(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%T(ix2, iy1, iz1, 1)
         Ta(3) = tshape3*grid%T(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%T(ix2, iy2, iz1, 1)
         Ta(4) = tshape4*grid%T(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%T(ix1, iy2, iz1, 1)
         Ta(5) = tshape1*grid%T(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%T(ix1, iy1, iz2, 1)
         Ta(6) = tshape2*grid%T(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%T(ix2, iy1, iz2, 1)
         Ta(7) = tshape3*grid%T(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%T(ix2, iy2, iz2, 1)
         Ta(8) = tshape4*grid%T(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%T(ix1, iy2, iz2, 1)
         !
         !ro
         rhoa(1) = tshape1*grid%ro(ix1, iy1, iz1, 2) + (1.0_rp - tshape1)*grid%ro(ix1, iy1, iz1, 1)
         rhoa(2) = tshape2*grid%ro(ix2, iy1, iz1, 2) + (1.0_rp - tshape2)*grid%ro(ix2, iy1, iz1, 1)
         rhoa(3) = tshape3*grid%ro(ix2, iy2, iz1, 2) + (1.0_rp - tshape3)*grid%ro(ix2, iy2, iz1, 1)
         rhoa(4) = tshape4*grid%ro(ix1, iy2, iz1, 2) + (1.0_rp - tshape4)*grid%ro(ix1, iy2, iz1, 1)
         rhoa(5) = tshape1*grid%ro(ix1, iy1, iz2, 2) + (1.0_rp - tshape1)*grid%ro(ix1, iy1, iz2, 1)
         rhoa(6) = tshape2*grid%ro(ix2, iy1, iz2, 2) + (1.0_rp - tshape2)*grid%ro(ix2, iy1, iz2, 1)
         rhoa(7) = tshape3*grid%ro(ix2, iy2, iz2, 2) + (1.0_rp - tshape3)*grid%ro(ix2, iy2, iz2, 1)
         rhoa(8) = tshape4*grid%ro(ix1, iy2, iz2, 2) + (1.0_rp - tshape4)*grid%ro(ix1, iy2, iz2, 1)
         !
         !
         !*** Interpolate variables at (x,y,z) position.
         !
         ux = 0.0_rp
         uy = 0.0_rp
         Tair(iz) = 0.0_rp
         rair(iz) = 0.0_rp
         do i = 1, 8
            ux = ux + shapef(i)*ua(i)
            uy = uy + shapef(i)*va(i)
            Tair(iz) = Tair(iz) + shapef(i)*Ta(i)
            rair(iz) = rair(iz) + shapef(i)*rhoa(i)
         end do
         Vair(iz) = sqrt(ux*ux + uy*uy)
         !
         !***     Gets air direction in terrain following --> Aair
         !
         if (abs(ux) .gt. 1d-8) then
            angle = atan2(uy, ux)*180.0/pi
         else
            if (uy .gt. 1d-8) then
               angle = 90.0
            else if (uy .lt. -1d-8) then
               angle = 270.0
            else
               angle = 0.0
            end if
         end if
         !
         Aair(iz) = angle
         if (Aair(iz) .lt. 0.0) Aair(iz) = 360.0 + Aair(iz)    ! Angle in deg. (0 to 360)
         !
         PSIa = PSIa + angle*(grid%z(iz) - z)
         z = grid%z(iz)

      end do ! iz
      !
      !***  Averaged wind direction
      !
      PSIa = PSIa/(grid%z(grid%nz))                ! Angle in deg. (-180 to +180)
      if (PSIa .lt. 0.0) PSIa = 360.0 + PSIa          ! Angle in deg. (0 to 360)
      PSIa = PSIa*pi/180.0                         ! Angle in rad.
      !
      !***  Computes Tair0 and Rair0
      !
      Tair0 = Tair(1)
      Rair0 = Rair(1)
      !
      !*** Computes the buoyancy frequency (squared)
      !
      Cao = 998.    ! specific heat capacity at constant pressure of dry air (J kg^-1 K^-1)
      do iz = 1, grid%nz
         !
         if (iz .eq. 1) then
            dTdz = (Tair(2) - Tair(1))/(grid%z(2) - grid%z(1))
         else if (iz .eq. grid%nz) then
            dTdz = (Tair(grid%nz) - Tair(grid%nz - 1))/(grid%z(grid%nz) - grid%z(grid%nz - 1))
         else
            dTdz = (Tair(iz + 1) - Tair(iz - 1))/(grid%z(iz + 1) - grid%z(iz - 1))
         end if
         !
         Nair(iz) = g*g*(1 + Cao*dTdz/g)/(Cao*Tair0)
         !
      end do

      return
   end subroutine interpolate_sm
!
   subroutine shape3(s, t, z, shapf)
      !******************************************************************************************
      !*
      !*    Calculates the factor shape of each node
      !*
      !******************************************************************************************
      use KindType
      implicit none
      !
      real(rp)  :: s, t, z !factors shape 1D
      real(rp)  :: sm, tm, zm, sq, tp, zp
      real(rp)  :: shapf(8)
      !
      sm = 0.5_rp*(1.0_rp - s)
      tm = 0.5_rp*(1.0_rp - t)
      zm = 0.5_rp*(1.0_rp - z)
      sq = 0.5_rp*(1.0_rp + s)
      tp = 0.5_rp*(1.0_rp + t)
      zp = 0.5_rp*(1.0_rp + z)
      shapf(1) = sm*tm*zm
      shapf(2) = sq*tm*zm
      shapf(3) = sq*tp*zm
      shapf(4) = sm*tp*zm
      shapf(5) = sm*tm*zp
      shapf(6) = sq*tm*zp
      shapf(7) = sq*tp*zp
      shapf(8) = sm*tp*zp
      !
      return
   end subroutine shape3
!
   subroutine shape2(s, t, shapf)
      !******************************************************************************************
      !*
      !*    Calculates the factor shape of each node
      !*
      !******************************************************************************************
      use KindType
      implicit none
      !
      real(rp)  :: s, t !factors shape 1D
      real(rp)  :: sm, tm, sq, tp
      real(rp)  :: shapf(4)
      !
      sm = 0.5_rp*(1.0_rp - s)
      tm = 0.5_rp*(1.0_rp - t)
      sq = 0.5_rp*(1.0_rp + s)
      tp = 0.5_rp*(1.0_rp + t)
      !
      shapf(1) = sm*tm
      shapf(2) = sq*tm
      shapf(3) = sq*tp
      shapf(4) = sm*tp
      !
      return
   end subroutine shape2

