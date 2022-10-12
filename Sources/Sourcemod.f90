!***************************************************************************
!*
!*		Module for master operations relating to SOURCE term
!*
!***************************************************************************
MODULE Sourcemod
  use KindType
  use InpOut
  IMPLICIT NONE
  SAVE
  !
  integer (ip) :: ns                    ! ns=npart2-npart1+1 : number of part liberated in current time=number of sources
  real    (rp) :: M0                    ! MFR from height using  MFR = rho*(H/2)**(1/0.241)
  real    (rp) :: HPlume                ! Plume height for current suphase
  real    (rp) :: a_suzuki              ! for current subphase
  real    (rp) :: l_suzuki              ! for current subphase
  real    (rp) :: d_top_hat             ! for current subphase
  real(rp), allocatable :: MFRpb(:)     ! MFR per bin M0/bins
  integer(ip), allocatable :: pb(:)     ! Particles per bin
  !
  real(rp), allocatable :: Vair(:)        ! Velocity
  real(rp), allocatable :: Tair(:)        ! Temperature
  real(rp), allocatable :: Aair(:)        ! Azimuth
  real(rp), allocatable :: Rair(:)        ! Density
  real(rp), allocatable :: Nair(:)        ! N2 buoyancy frequency
  real(rp)              :: PSIa,Tair0,Rair0
  character(len=30) :: ext5     !auxiliar variable to write in log file
  !
CONTAINS
  !
  !
  subroutine solvepoint(iphase,part1,part2)
    !********************************************************
    !*
    !*    Point source solution.
    !*     x : LON
    !*     y : LAT
    !*     z : height a.s.l
    !*
    !********************************************************
    use KindType  
    use Master
    implicit none
    integer (ip) :: iphase,ipart, part1,part2
    !
    !***  Source position
    !
    do ipart=part1,part2
    part(ipart)%lon = phase(iphase)%lon
    part(ipart)%lat = phase(iphase)%lat
    part(ipart)%z = Hplume !+ phase(iphase)%source_elev !+  dZ0 
    end do
    !
    !***  MER per bin
    !
       MFRpb(1:phase(iphase)%bins)= M0*phase(iphase)%fc(1:phase(iphase)%bins) 
    !
    return
  end subroutine solvepoint
  !
  !
  subroutine solvelinear(iphase,part1,part2)
    !********************************************************
    !*
    !*    Linear source solution.
    !*     x : LON
    !*     y : LAT
    !*     z : height a.s.l
    !*
    !********************************************************
    use KindType  
    use Master
    implicit none
    !
    integer (ip) :: part1
    integer (ip) :: part2
    integer (ip) :: iphase,ipart
    real    (rp) :: deltaz
    !
    !*** Horizontal Source position
    !
    do ipart=part1,part2
    part(ipart)%lon = phase(iphase)%lon
    part(ipart)%lat = phase(iphase)%lat
    end do
    !
    !***  Vertical position and total mass according to Linear distribution
    !
    deltaz = Hplume/ns   ! from the botttom to Hplum
    do ipart = 1,part2-part1+1
       part(ipart+part1-1)%z = ipart*deltaz
    end do
    !
    !***  MFR per bin
    !
       MFRpb(1:phase(iphase)%bins) = M0*phase(iphase)%fc(1:phase(iphase)%bins)
    !
    return
  end subroutine solvelinear
  !
  !
  !
  subroutine solvetophat(iphase,part1,part2)
    !********************************************************
    !*
    !*    Linear source solution.
    !*     x : LON
    !*     y : LAT
    !*     z : height a.s.l
    !*
    !********************************************************
    use KindType  
    use Master
    implicit none
    !
    integer (ip) :: part1
    integer (ip) :: part2
    integer (ip) :: iphase,ipart
    real    (rp) :: deltaz,z
    !
    !*** Horizontal Source position
    !
    do ipart=part1,part2
    part(ipart)%lon = phase(iphase)%lon
    part(ipart)%lat = phase(iphase)%lat
    end do
    !
    !***  Vertical position and total mass according to Top-Hat distribution
    !
    deltaz = d_top_hat/ns       ! from (Hplum-d_top_hat) to Hplum
    do ipart =  1,part2-part1+1
       part(ipart+part1-1)%z = Hplume-d_top_hat+(ipart*deltaz)
    end do
    !
    !***  MFR per bin
    !
       MFRpb(1:phase(iphase)%bins) = M0*phase(iphase)%fc(1:phase(iphase)%bins)
    !
    return
  end subroutine solvetophat
  !
  !
  subroutine solvesuzuki(iphase,part1,part2)
    !********************************************************
    !*
    !*    Suzuki source solution
    !*     x : LON
    !*     y : LAT
    !*     z : height a.s.l
    !*
    !********************************************************
    use KindType
    use Master
    implicit none
    !
    integer (ip) :: part1
    integer (ip) :: part2
    integer(ip) :: iphase, ipart, ibin
    real(rp)    :: deltaz,z,mfr
    real (8)    :: sum
    real(8), allocatable :: Sou(:)
    !
    !***  Horizontal source position
    !
    do ipart=part1,part2
    part(ipart)%lon = phase(iphase)%lon
    part(ipart)%lat = phase(iphase)%lat      
    end do
    !
    !***  Vertical position and total mass(fraction of mass) according to Suzuki distribution
    !
    allocate(Sou(part1:part2))
    sum = 0.0
    deltaz = Hplume/(ns)  
    do ipart = part1,part2
       z = (ipart-part1+1)*deltaz
       if(part1.eq.part2) z = Hplume-deltaz ! If in each time step is oly one particle liberated
       Sou(ipart) = ((1.0-z/Hplume)*exp(a_suzuki*(z/Hplume-1.0)))**l_suzuki
       part(ipart)%z = z
       sum = sum + Sou(ipart)
    end do
    !
    !*** Normalization to MFR (SUMz=MFR)
    !
    do ipart = part1,part2
    Sou(ipart) = Sou(ipart)/sum
    end do
    ! divido la masa total en subgrupos de particulas segun su clase y con su respectiva fraccion
    do ibin=1,phase(iphase)%bins
       MFRpb(ibin) = M0*phase(iphase)%fc(ibin)
    end do
    !
    pb(:)=0
    do ipart=part1,part2
       ibin=part(ipart)%ibin
       pb(ibin)=pb(ibin)+1
    end do
    ! 
    MFR=0.
    do ipart = part1,part2
       ibin = part(ipart)%ibin
       part(ipart)%mass=  0.5*(MFRpb(ibin)/pb(ibin)+M0*Sou(ipart))
     !  MFR = MFR+part(ipart)%mass
    end do
    !
    !
    deallocate(Sou)
    return
  end subroutine solvesuzuki
  !
  !
  !
  subroutine detvalue(iphase,part1,part2)
    !*************************************************************************
    !*
    !*    Determines the source file during a current time step (time,time+sidt)
    !*
    !*    part1,part2: particle interval liberated in the current time
    !*    for phase iphase
    !*************************************************************************
    use KindType
    use Master
    IMPLICIT NONE
    SAVE
    !
    integer (ip) :: part1
    integer (ip) :: part2
    integer (ip) :: iphase    ! Input parameters. Particle range
    integer (ip) :: ic                  
    real    (rp) :: x,y,z,MFR
    integer (ip) :: ipart                  ! particle iteration
    integer (ip) :: ibin, rest, totpart
!    integer (ip) :: ibinstart,ibinend      ! ibin corresponding to the part1 and part2 respectively

    !
    !*** MFR per bin to determine the particle mass
    !*** First, calculate pb(:)  -particles per bin-
    !
    pb(:)=0
    do ipart=part1,part2
       ibin=part(ipart)%ibin
       pb(ibin)=pb(ibin)+1
    end do
    !
    !*** Mplume(ipart): Mass flow rate per particle
    !
    SELECT CASE(phase(iphase)%source_type)
    case('POINT','LINEAR','TOP-HAT')
       do ipart=part1,part2
          part(ipart)%mass = MFRpb(part(ipart)%ibin)/pb(part(ipart)%ibin)     ! mass per particle, according to their bin.
       end do
    case('SUZUKI')
    case default 
       call runend('Incorrect source type')
    END SELECT
    !
    !***  Total mass erupted rate
    !
    MFR=0.
    do ipart=part1,part2
       MFR = MFR+part(ipart)%mass
    end do
    write(ext5,*) MFR
!    if(my_id.eq.0)write(lulog,5) 'The mass erupted, sum of all particle in the present range, is'//ext5// 'kg'
!5   format(/,  a /)
    !
    !
    !***  Writes header for the current interval
    !
    !  write(lusrc,10) isec1,isec2
    !  write(lusrc,11) ns,nc
    !  write(lusrc,12) MFR
    !10 format(i7,1x,i7)
    !11 format(i7,1x,i7)
    !12 format(e16.9)
    !
    !***  Writes the rest of file
    !
!    do ipart=part1,part2 !particles activated in the current times
!       part(ipart)%mass= Mplum(ipart) !Determine characteristics for each sub-source(is)
!       part(ipart)%lon = Xplum(ipart)
!       part(ipart)%lat = Yplum(ipart)
       !     if(terrain_following) then
       !        part(ipart)%z = Zplum(is) - phase(iphase)%source_elev !-dZ0  ! terrain following coordinates
       !     else
!       part(ipart)%z = Zplum(ipart)            ! a.s.l.
       !     end if
       !     SELECT CASE(coord_sys)
       !     case('LON-LAT')
       !        write(lusrc,20) x,y,z,(Mplum(ic,is),ic=1,nc)
       !20      format(2(1x,f11.6),2x,f9.0,2x,100(e16.9,1x))
       !     case('UTM')
       !       write(lusrc,21) x,y,z,(Mplum(ic,is),ic=1,nc)
       !21      format(2(1x,f10.1),2x,f9.0,2x,100(e16.9,1x))
       !     END SELECT
!    end do  !part
    !


    !
    return
  end subroutine detvalue
  !

   subroutine merwind(iphase)
  !********************************************************
  !*
  !*   Computes MER depending on the wind profile
  !*
  !********************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  integer(ip) :: iphase
  integer(ip) :: iz,jz
  real   (rp) :: Cao,Co,To,alfa,beta,z1,gprime
  real   (rp) :: v_mean,N_mean
  real   (rp) :: H1,V1,Ws
  character(len=30) :: ext6  
  !
  SELECT CASE(phase(iphase)%mer_vs_h)
  case('ESTIMATE-DEGRUYTER')
    !
    !*** Estimates MER as in Degruyter and Bonadonna (2012)
    !
    Cao     =  998.0_rp    ! specific heat capacity at constant pressure of dry air (J kg^-1 K^-1)
    Co      = 1250.0_rp    ! specific heat capacity at constant pressure of solids (J kg^-1 K^-1)
    To      = 1200.0_rp    ! initial plume temperature (K)
    alfa    = 0.1_rp       ! radial entrainment coefficient
    beta    = 0.5_rp       ! wind entrainment coefficient
    z1      = 2.8_rp       ! maximum non-dimensional height (Morton et al. 1956)
    !
    gprime = g*(Co*To-Cao*Tair0)/(Cao*Tair0)
    !
    v_mean = 0.0_rp
    N_mean = 0.0_rp
    jz     = 2
    do iz = 2,grid%nz
       if(grid%z(iz-1).le.HPlume) then
           v_mean = v_mean + 0.5_rp*(Vair(iz-1)+Vair(iz))*(grid%z(iz)-grid%z(iz-1))
           N_mean = N_mean + 0.5_rp*(Nair(iz-1)+Nair(iz))*(grid%z(iz)-grid%z(iz-1))
           jz = iz
       end if
    end do
    v_mean = v_mean/(grid%z(jz)-grid%z(1))
    N_mean = N_mean/(grid%z(jz)-grid%z(1))
    N_mean = sqrt(N_mean)
    !
    M0 = (2.0_rp**(5.0_rp/2.0_rp))*alfa*alfa*N_mean*N_mean*N_mean*HPlume*HPlume*HPlume*HPlume/(z1*z1*z1*z1)
    M0 = M0 + beta*beta*N_mean*N_mean*v_mean*HPlume*HPlume*HPlume/6.0_rp
    M0 = pi*Rair0*M0/gprime   !MFR is kg/s, then it is multiplied by simulation time step.
    M0 = M0/FLOAT(num_procs)
    write(ext6,*) M0
!    if(my_id.eq.0)write(lulog,4) ext6
!4   format(/, 'The Mass flow rate estimated  : ', a , 'kg/s')
    M0 = ABS(sidt)*M0
    ! 
    ! Save the MER for this phase into a counter 
    phase(iphase)%erumass=phase(iphase)%erumass+M0
    ! 
    !
  case('ESTIMATE-WOODHOUSE')
    !
    !*** Estimates MER as in Woodhouse et al. (2012)
    !
    N_mean = 0.0_rp
    jz     = 2
    do iz = 2,grid%nz
       if(grid%z(iz-1).le.HPlume) then
           N_mean = N_mean + 0.5_rp*(Nair(iz-1)+Nair(iz))*(grid%z(iz)-grid%z(iz-1))
           jz = iz
       end if
    end do
    N_mean = N_mean/(grid%z(jz)-grid%z(1))
    N_mean = sqrt(N_mean)
    !
    V1 = Vair(jz)                 ! reference velocity
    H1 = grid%z(jz) - grid%z(1)   ! reference height
    !
    Ws = 1.44_rp*V1/(N_mean*H1)
    Ws = 0.318_rp*(1.0_rp+1.373_rp*Ws)/(1.0_rp+4.266_rp*Ws+0.3527_rp*Ws*Ws)
    M0 = ABS(sidt)*(HPlume/1d3/Ws)**(1.0_rp/0.253_rp)
    M0=M0/FLOAT(num_procs)
    ! 
    ! Save the MER for this phase into a counter 
    phase(iphase)%erumass=phase(iphase)%erumass+M0
    ! 
    !
  END SELECT
  !
  return
  end subroutine merwind

END MODULE Sourcemod
