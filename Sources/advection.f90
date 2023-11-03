subroutine advection(z, lon, lat, ua, va, wa, rho, diam, sphe, psi)
   !******************************************************************************************
   !*
   !*    Advection.
   !*    INPUTS:
   !*    z : height
   !*    ua: air horizontal wind component
   !*    va: air vertical   wind component
   !*    wa: air normal     wind component, positive down
   !*    rho,diam,sphe,psi: particle properties
   !*
   !******************************************************************************************
   use Master
   use KindType
   implicit none
   save
   !
   real(rp) :: z, z0, lon, lon0, lat, lat0       ! particle position
   real(rp) :: ua, va, wa, vset                ! Velocities in initial position
   real(rp) :: ub, vb, wb                   ! velocities without turb in first approximated position.
   real(rp) :: a, b                          ! Rearth+z, angle
   real(rp) :: rho, diam, sphe, psi            ! Particle properties
   real(rp) :: rhoa, mua                    ! calculate in interpolate_winds to use in sedimentation, outputs in interpolate_wind, input in sedimentation

   !
   a = Rearth + z
   b = lat*pi/180
   !***Save the position P(t) to use in the Pettersen iteration
   lon0 = lon
   lat0 = lat
   z0 = z
   !
   !*** Calculate the first approximation: P'(t+sidt)
   lon = lon0 + ABS(sidt)*ua/(a*cos(b))*180/pi
   b = a*cos(b)
   lat = lat0 + ABS(sidt)*va/a*180/pi
   z = z - (wa)*sidt
   if (z .lt. 0) z = 0
   !
   !*** One iteration of Pettersen scheme. It is accurate to the second order
   !
   ! *** Interpolate grid-winds at new postion P'(t+sidt): V(P',t+sidt) the pettersen scheme use only grid-winds
   ! ** ub, vb, wb, rhoa and mua are outputs.
!  call interpolate_wind(lon,lat,z,ub,vb,wb,rhoa,mua)
   !
   !*** Calculate the settling velocity in this position
   ! * vset is output
!  call sedimentation(vset,rho,rhoa,diam,mua,sphe,psi,sedmodel,sidt)
   !*** Change wb
!  wb = -wb + vset  !wb is positive upward, and vset is positive down. After this line, wb is positive down.
   ! *** Calculate the wind 0.5*[V(P,t)+V(P',t+sidt)]
!  ua = (ub + ua)/2
!  va = (vb + va)/2
!  wa = (wb + wa)/2 !wb and wa are positive down
   !*** Calculate the final position P(t+sidt)
!  lon = lon0 + ABS(sidt)*ua/(a*cos(b))*180/pi
!  lat = lat0 + ABS(sidt)*va/a*180/pi
!  z   = z0   - (wa)*sidt
!  if(z.lt.0)z=0
   !
   return
end subroutine advection
