subroutine diffusion(ugrid,vgrid,wgrid,z,pblha,rmonin,usta,wst,drhodz,rhoa,lon,lat,ux,vy,wp,age)
  !******************************************************************************************
  !* ** Calculate the turbulence velocities according the height   
  !* ** Above the PBL, constant diffusion is assumed, in the PBL the Hanna scheme is used.
  !*   
  !*    INPUTS:
  !*    ugrid, vgrid,wgrid   : grid-wind velocities, wgrid is positive down
  !*    z                    : particle height
  !*    pblha                : PBL height
  !*    rmonin               : onin obukhoc lenght
  !*    usta                 : friction velocity
  !*    wst                  : convective velocity scale
  !*    drhodz               : air density vertical gradient
  !*    rhoa                 : air density in particle position
  !*    pux, pvy,pwp         : random velocities stored for previous step
  !*
  !*    CALCULATE TURBULENCE VELOCITIES AND SAVE THEM IN UGRID, VGRID,WGRID
  !*    ux,vy,wp  : random velocities due to turbulence (along wind, cross    
  !*                wind, vertical wind 
  !*    
  !******************************************************************************************
  use Master
  use KindType                              
  implicit none
  real (rp) :: uxscale
  real (rp) :: ux, vy, wp          ! random velocities due to turbulence, along, cross and vertical wind
  real (rp) :: ru, rv, rw          ! auxiliay variables
  real (rp) :: wpscale, weight
  real (rp) :: z                   ! particle height. INPUT
  real (rp) :: pblha !I
  real (rp) :: rmonin!I
  real (rp) :: usta  !I
  real (rp) :: wst   !I
  real (rp) :: drhodz!I
  real (rp) :: rhoa  !I
  real (rp) :: ugrid, vgrid, wgrid ! grid scales velocities>: I/O                      
  real (rp) :: sigu,sigv,sigw      ! Velocities standard deviation
  real (rp) :: tlu,tlv,tlw         ! lagrgian time scales
  real (rp) :: dsigw2dz            !  vertical gradient of sigw  
  integer (ip):: ipart 
  real (rp) :: delz                ! wp*sidt
  real (rp) :: lon, lat    !
  real (rp) :: a,b                 ! Rearth+z, angle
  real (rp) :: age
!
!
  a    = Rearth + z
  b    = lat*pi/180

! Assume constant, uncorrelated, turbulent perturbations
! In the stratosphere, use a small vertical diffusivity d_strat,
! in the troposphere, use a larger horizontal diffusivity d_trop.
! Turbulent velocity scales are determined based on sqrt(d_trop/dt)
! According Legras et al., 2003
!******************************************************************
!
if(z.gt.pblha)then
  if (z.lt.tropopause) then  ! in the troposphere
    uxscale=sqrt(2.*d_trop/sidt) !d_trop is defined in master, it's equal 50. 
    if (nran+1.gt.maxrand) nran=1
    ux=rannumb(nran)*uxscale
    vy=rannumb(nran+1)*uxscale
    nran=nran+2
    wp=0.
  else if (z.lt.tropopause+1000.) then     ! just above the tropopause: make transition
    weight=(z-tropopause)/1000.
    uxscale=sqrt(2.*d_trop/sidt*(1.-weight))
    if (nran+2.gt.maxrand) nran=1
    ux=rannumb(nran)*uxscale
    vy=rannumb(nran+1)*uxscale
    wpscale=sqrt(2.*d_strat/sidt*weight)
    wp=rannumb(nran+2)*wpscale+d_strat/1000.
    nran=nran+3
  else                 ! in the stratosphere
    if (nran.gt.maxrand) nran=1
    ux=0.
    vy=0.
    wpscale=sqrt(2.*d_strat/sidt)
    wp=rannumb(nran)*wpscale
    nran=nran+1
  endif
! 
!********************************************
!**** Turbulent in the ABL
!     
!********************************************
! 
else   
  !*** Calculate the turbulence disturbance
  !*** lagrangian time scales
  !*** Sigmas veocities
  ! 
  call hanna(z,rmonin,pblha,usta,wst,tlu,tlv,tlw,sigu,sigv,sigw,dsigw2dz)
  if (age.eq.sidt)then!if(time.eq.sist)then
        if (nran+2.gt.maxrand) nran=1
        ux=rannumb(nran)*sigu
        vy=rannumb(nran+1)*sigv
        wp=rannumb(nran+2)*sigw
  end if
  ! 
  !*****************************************
  ! Determine the new diffusivity velocities
  !*****************************************
  ! 
  ! Horizontal components                                                               
  !**********************
  ! 
  if (nran+1.gt.maxrand) nran=1
  if (sidt/tlu.lt..5) then
     ux=(1.-sidt/tlu)*ux+rannumb(nran)*sigu*sqrt(2.*sidt/tlu)  !eq 18 from Stohl et al., 2005
  else
     ru=exp(-sidt/tlu)
     ux=ru*ux+rannumb(nran)*sigu*sqrt(1.-ru**2)                !eq 17 from Stohl et al., 2005
  endif
  if (sidt/tlv.lt..5) then
     vy=(1.-sidt/tlv)*vy+rannumb(nran+1)*sigv*sqrt(2.*sidt/tlv)!eq 18 from Stohl et al., 2005
  else
     rv=exp(-sidt/tlv)
     vy=rv*vy+rannumb(nran+1)*sigv*sqrt(1.-rv**2)              !eq 17 from Stohl et al., 2005
  endif
  nran=nran+2
  !
  ! Vertical Component
  !********************
  ! 
  if (nran+1.gt.maxrand) nran=1
  if(sidt/tlw.lt. 0.5)then
    wp =  (1.-sidt/tlw)*wp+sidt*sigw*(dsigw2dz+sigw*drhodz/rhoa) &
          +rannumb(nran)*sigw*sqrt(2.*sidt/tlw)                 !eq 18 from Stohl et al., 2005
  else
     rw=exp(-sidt/tlw)  !autocorrelation of vertical wind
     wp=rw*wp+tlw*(1.-rw)*(dsigw2dz+sigw*sigw*drhodz/rhoa)+ &
        rannumb(nran+1)*sigw*sqrt(1.-rw**2)                     !eq 17 from Stohl et al., 2005
  end if
  delz=wp*sidt
  if(delz.gt.pblha) wp=mod(delz,pblha)/sidt
!
!
end if ! Diffussion type
  !
  !*** Add the turbulence velocities to the grid scale velocities.  
  !
ugrid=ugrid+ux        
vgrid=vgrid+vy
wgrid=wgrid+wp   !wgrid and wp are positive down
!  
return
end subroutine diffusion

!
subroutine hanna(z,rmonin,pblh,ust,wst,tlu,tlv,tlw,sigu,sigv,sigw,dsigw2dz)
!***************************************************************************************
! 										       *
!*** Inputs: 									       *
!           pbl pblh height at particle position				       *
!           rmonin obukhov lenght at particle position				       *
!           z particle height							       *
!           ust friction velocity at particle position				       *
!**Determine 									       *
!           sigu, sigv, sigw  standard deviations of turbulent velocity fluctuations   *
!           tlu [s]           Lagrangian time scale for the along wind component.      *
!           tlv [s]           Lagrangian time scale for the cross wind component.      *
!           tlw [s]           Lagrangian time scale for the vertical wind component.   *
!           wst [m/s]         convective velocity scale     	                       *
!           dsigw2dz [1/s]     vertical gradient of sigw                               *
! 										       *
!***************************************************************************************
  use Master
  use KindType                              
  implicit none
  real(rp) :: corr,z, s1, s2
  real(rp) :: ust,wst,rmonin,pblh,zeta,sigu,sigv,tlu,tlv,tlw
  real(rp) :: sigw,dsigw2dz

  ! 
  !**********************
  ! 1. Neutral conditions
  !**********************

  if (pblh/abs(rmonin).lt.1.) then
    ust=max(1.e-4,ust)
    corr=z/ust
    sigu=2.0*ust*exp(-3.e-4*corr)
    sigu=max(sigu,1.e-5)
    sigv=1.3*ust*exp(-2.e-4*corr)
    sigv=max(sigv,1.e-5)
    sigw=sigv
    dsigw2dz=-6.76e-4*ust*exp(-4.e-4*corr)
    tlu=0.5*z/sigw/(1.+1.5e-3*corr)
    tlv=tlu
    tlw=tlu


  !***********************
  ! 2. Unstable conditions
  !***********************

  else if (rmonin.lt.0.) then


  ! Determine sigmas
  !*****************
    zeta=z/pblh
    sigu=ust*(12-0.5*pblh/rmonin)**0.33333
    sigu=max(sigu,1.e-6)
    sigv=sigu
    if (zeta.lt.0.03) then
      sigw=0.96*wst*(3*zeta-rmonin/pblh)**0.33333
      dsigw2dz=1.8432*wst*wst/pblh*(3*zeta-rmonin/pblh)**(-0.33333)
    else if (zeta.lt.0.4) then
      s1=0.96*(3*zeta-rmonin/pblh)**0.33333
      s2=0.763*zeta**0.175
      if (s1.lt.s2) then
        sigw=wst*s1
        dsigw2dz=1.8432*wst*wst/pblh*(3*zeta-rmonin/pblh)**(-0.33333)
      else
        sigw=wst*s2
        dsigw2dz=0.203759*wst*wst/pblh*zeta**(-0.65)
      endif
    else if (zeta.lt.0.96) then
      sigw=0.722*wst*(1-zeta)**0.207
      dsigw2dz=-.215812*wst*wst/pblh*(1-zeta)**(-0.586)
    else if (zeta.lt.1.00) then
      sigw=0.37*wst
      dsigw2dz=0.
    endif
    sigw=max(sigw,1.e-6)

  ! Determine average Lagrangian time scale
  !****************************************

    tlu=0.15*pblh/sigu
    tlv=tlu
    if (z.lt.abs(rmonin)) then
      tlw=0.1*z/(sigw*(0.55-0.38*abs(z/rmonin)))
    else if (zeta.lt.0.1) then
      tlw=0.59*z/sigw
    else
      tlw=0.15*pblh/sigw*(1.-exp(-5*zeta))
    endif



  !*********************
  ! 3. Stable conditions
  !*********************

  else
    zeta=z/pblh
    sigu=2.*ust*(1.-zeta)
    sigv=1.3*ust*(1.-zeta)
    sigu=max(sigu,1.e-6)
    sigv=max(sigv,1.e-6)
    sigw=sigv
    dsigw2dz=3.38*ust*ust*(zeta-1.)/pblh
    tlu=0.15*pblh/sigu*(sqrt(zeta))
    tlv=0.467*tlu
    tlw=0.1*pblh/sigw*zeta**0.8
  endif


  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(30.,tlw)


end subroutine hanna
