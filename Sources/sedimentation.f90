subroutine sedimentation(vsedim, rho, rhoa, diam, mua, sphe, psi, sedmodel, sidt)
   !******************************************************************************************
   !*
   !*    Calculate the terminal velocity
   !* INPUTS
   !*  sidt     : silaion time interval
   !*  sedmodel : sedimentation model choosen by the user
   !*  psi      : particle properties
   !*  sphe
   !*  mua
   !*  diam
   !*  rhoa
   !*  rho
   !* OUTPUT
   !*  vsedim   : Sedimentation velocity
   !*
   !******************************************************************************************
   use KindType
   implicit none
   save
   !
   integer(ip)    :: it
   logical        :: conve
   !
   integer(ip)    :: maxit
   real(rp)    :: vset, v1, v2, v3
   real(rp)    :: gi, eps, rey, rk1, rk2, a, b, cd100, vold, dnd, gama   ! work variables
   real(rp)    :: cd
   real(rp)    :: vsedim, rho, rhoa, diam, mua, sphe, psi
   integer(ip)    :: sedmodel ! sedimentation Model
   real(rp)    :: sidt     ! Simulation increment step (Negative in backward case)
   !
   conve = .true.
   gi = 9.81_rp              ! Gravity acceleration
   eps = 1d-3                ! Tolerance
   maxit = 100               ! Maximum number of iterations
   vold = 1d5
   cd = 1.0_rp
   !
   !***  Begin iteration to compute vset
   !
   !
   !***  0. Stokes Law
   !
   if (sedmodel == 0) then
      vsedim = (2._rp/9._rp)*(rho - rhoa)*gi*((diam/2._rp)**2._rp)/mua
      !  if(sidt.lt.0)vsedim=-vsedim
      return
   end if
   !
   !***  1. Model of Arastoopour et al. 1982
   !
   if (sedmodel == 1) then
      do it = 1, maxit
         vset = sqrt(4.0_rp*gi*diam*rho/(3.0_rp*cd*rhoa))
         rey = rhoa*vset*diam/mua
         if (rey <= 988.947_rp) then
            cd = 24.0_rp/rey*(1.0_rp + 0.15_rp*rey**0.687_rp)
         else
            cd = 0.44_rp
         end if
         if (it > 1 .and. abs(vset - vold) <= eps) goto 10
         vold = vset
      end do
   end if
   !
   !***  2. Model of Ganser 1993
   !
   if (sedmodel == 2) then
      call get_gama(diam, sphe, gama)  ! Get gama=a/c
      dnd = 0.5_rp*(1.0_rp + gama)*gama**(-2.0/3.0)
      !
      do it = 1, maxit
         vset = sqrt(4.0_rp*gi*diam*rho/(3.0_rp*cd*rhoa))
         rey = rhoa*vset*diam/mua
         !                   rk1=1.0_rp/(1.0_rp/3.0_rp+2.0_rp/(3.0_rp*sqrt(psi)))

         rk1 = 1.0_rp/(dnd/3.0_rp + 2.0_rp/(3.0_rp*sqrt(psi)))
         rk2 = 10.0_rp**(1.8148_rp*(-log10(psi))**0.5743_rp)
         cd = 24.0_rp/(rey*rk1)*(1.0_rp + 0.1118_rp* &
                                 (rey*rk1*rk2)**0.6567_rp) + 0.4305_rp*rk2/ &
              (1.0_rp + 3305.0_rp/(rey*rk1*rk2))

         if (it > 1 .and. abs(vset - vold) <= eps) goto 10
         vold = vset
      end do
   end if
   !
   !***  3. Model of Wilson & Huang 1979
   !
   if (sedmodel == 3) then
      do it = 1, maxit
         vset = sqrt(4.0_rp*gi*diam*rho/(3.0_rp*cd*rhoa))
         rey = rhoa*vset*diam/mua

         if (rey <= 100.0_rp) then
            cd = 24.0_rp/rey*psi**(-0.828_rp) + 2.0_rp*sqrt(1.0_rp - psi)
         elseif (rey > 100.0_rp .and. rey < 1000.0_rp) then
            cd100 = 0.24_rp*psi**(-0.828_rp) + 2.0_rp*sqrt(1.0_rp - psi)
            a = (1.0_rp - cd100)/900.0_rp
            b = 1.0_rp - 1000.0_rp*a
            cd = a*rey + b
         else
            cd = 1.0_rp
         end if

         if (it > 1 .and. abs(vset - vold) <= eps) goto 10
         vold = vset
      end do
   end if
   !
   !***  4. Model of dellino et al. 2005
   !
   if (sedmodel == 4) then
      vset = ((diam*diam*diam*gi*(rho - rhoa)*rhoa*(psi**1.6_rp))/ &
              (mua*mua))**0.5206_rp
      vsedim = (1.2065_rp*mua*vset)/(diam*rhoa)
      return
   end if
   !
   !*** No convergence
   !
   conve = .false.
   return
   !
   !*** Convergence
   !
10 continue
   vsedim = sqrt(4.0_rp*gi*diam*rho/(3.0_rp*cd*rhoa))
   !
   return
end subroutine sedimentation
!

