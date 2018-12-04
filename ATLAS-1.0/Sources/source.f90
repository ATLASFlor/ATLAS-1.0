   subroutine source(iphase,part2)
  !******************************************************************************************
  !*    
  !*    Calculate the terminal velocity
  !*    iphase      : Input parameter. Phase of source term
  !*    part2       : Output parameter. particle range to calculate source term
  !*
  !******************************************************************************************
  use KindType 
  use Master                             
  use Sourcemod
  implicit none
  !
  integer (ip) :: part1
  integer (ip) :: part2
  integer (ip) :: iphase
  integer (ip) :: isp          ! i-sub phase-for loop-
  integer (ip) :: subphase     ! suphase corresponding to the actual time
  ! interpolate variables
  integer(ip)  :: ix1,ix2,iy1,iy2,iz1,iz2
  integer(ip)  :: nm1,nm2,nm3,nm4
  integer(ip)  :: ipoin1, ipoin2,ipoin3,ipoin4
  real   (rp)  :: xshape, yshape, zshape
  real   (rp)  :: shapef(8)
  integer(ip),save  :: ipass=0
  character(len=30) :: ext6     !auxiliar variable to write in log file
  !
  !*** Determine in which sub-phases is the current time. nsp: number of sub-phases
  !*** subphase: number of curent subphase
  !
  if(phase(iphase)%nsp.eq.1)then
     subphase=1
  else
     do isp=1,phase(iphase)%nsp-1
        if((time.ge.phase(iphase)%begtime(isp)).and.(time.lt.phase(iphase)%begtime(isp+1))) subphase=isp
     end do
     if(time.ge.phase(iphase)%begtime(phase(iphase)%nsp).and.time.le.phase(iphase)%endtime) subphase=phase(iphase)%nsp
  end if
  ! 
  ! 
  !*** Averiguate if there is a subphase, and if the subphase have a column height greater than zero
  ! 
  if (phase(iphase)%colheight(subphase).le.0) then
  if(subphase.eq.1)call runend('The first column height of each phase must be greater than 0')!part2=phase(iphase)%actnpart
  return
  end if
  ! 
  ! 
  !*** Update number of activated particles and the range of particles (from part1 to part2)
  ! 
  part1                  = phase(iphase)%actnpart + phase(iphase)%firstpart                                         
  phase(iphase)%actnpart = phase(iphase)%actnpart + phase(iphase)%npdt
  !
  !*** Adjust the number of activated particles if this is the last time step in this phase
  !
  if((time+sidt).gt.phase(iphase)%endtime) &
  phase(iphase)%actnpart = phase(iphase)%npart
  !
  part2 = phase(iphase)%actnpart+phase(iphase)%firstpart-1
  !
  !*** Release particles in the calculated range
  !
  part(part1:part2)%state=1     
  !
  !*** Allocate memory
  !
  ns = part2-part1+1                  ! number of particles to perform the source term
  allocate(MFRpb(phase(iphase)%bins)) ! mass flow rate per bin
  allocate(pb(phase(iphase)%bins))    ! particles per bin
  !
  if(ipass.eq.0)then
     ipass=1
     !
     !*** Allocate memory for Source Variable
     !
     allocate(Tair(grid%nz))
     allocate(Vair(grid%nz))
     allocate(Aair(grid%nz))
     allocate(Rair(grid%nz))
     allocate(Nair(grid%nz))
  end if
  !
  !*** Estimate MFR from height using  MFR = rho*(H/2)**(1/0.241)
  !
  if(phase(iphase)%mer_vs_h.eq.'ESTIMATE-MASTIN')then
  !
     M0 = phase(iphase)%rhomean*((0.5_rp*phase(iphase)%colheight(subphase)/1d3)**(1.0_rp/0.241_rp)) !MFR
     M0 = M0*ABS(sidt)                                                                              !Mass erupted in this time
  ! 
  ! Save the MER for this phase into a counter 
  phase(iphase)%erumass=phase(iphase)%erumass+M0
  ! 
  else if(phase(iphase)%mer_vs_h.eq.'NONE') then
     M0=phase(iphase)%M0(subphase)
     M0=M0*ABS(sidt)
  ! 
  ! Save the MER for this phase into a counter 
  phase(iphase)%erumass=phase(iphase)%erumass+M0
  ! 
  end if
  !
  !*** Source-dependent choice
  !
  SELECT CASE(phase(iphase)%source_type)
  case('POINT','LINEAR','TOP-HAT','SUZUKI')
     !
     !***  POINT, LINEAR, TOP-HAT and SUZUKI cases
     !
     SELECT CASE(phase(iphase)%mer_vs_h)
     case('NONE','ESTIMATE-MASTIN')
        !
        !*** No plume-wind coupling
        !
        HPlume = phase(iphase)%colheight(subphase)
        !
        SELECT CASE(phase(iphase)%source_type)
        case('POINT')
           call solvepoint(iphase,part1,part2)  !in sourcemod
        case('LINEAR')
           call solvelinear(iphase,part1,part2) !in sourcemod
        case('TOP-HAT')
           d_top_hat = phase(iphase)%d_top_hat(subphase)
           call solvetophat(iphase,part1,part2) !in sourcemod
        case('SUZUKI')
           a_suzuki  = phase(iphase)%a_suzuki(subphase)
           l_suzuki  = phase(iphase)%l_suzuki(subphase)
           call solvesuzuki(iphase,part1,part2) !in sourcemod
        case default 
           call runend('Incorrect source type')
        END SELECT ! source type
        !
        !
     case('ESTIMATE-DEGRUYTER','ESTIMATE-WOODHOUSE')
        !
        !*** Plume-wind coupling
        !
        call interpolate_sm(phase(iphase)%lon,phase(iphase)%lat,phase(iphase)%source_elev)  !in interpolate.f90 subroutine
        !
        !*** M0 depending on height and wind profile
        !
        HPlume = phase(iphase)%colheight(subphase)
        call merwind(iphase)                    !in sourcemod. output M0 value
        !
        SELECT CASE(phase(iphase)%source_type)
        case('POINT')
           call solvepoint(iphase,part1,part2) !in sourcemod
        case('LINEAR')
           call solvelinear(iphase,part1,part2) !in sourcemod
        case('TOP-HAT')
           d_top_hat = phase(iphase)%d_top_hat(subphase)
           call solvetophat(iphase,part1,part2) !in sourcemod
        case('SUZUKI')
           a_suzuki  = phase(iphase)%a_suzuki(subphase)
           l_suzuki  = phase(iphase)%l_suzuki(subphase)
           call solvesuzuki(iphase,part1,part2) !in sourcemod
        case default 
           call runend('Incorrect source type')
        END SELECT
     !          
     END SELECT   ! mer_vs_h
     ! Determine the mass per particle
     call detvalue(iphase,part1,part2) ! in Sourcemod.f90
     !
  END SELECT      ! source_type
  !
  !*** Writes in log file the MER for this Dt
  !
  write(ext6,*) phase(iphase)%erumass
  write(lulog,5) iphase, ext6
5 format(/,  'The mass erupted in the phase' //i2//' until this time is'//a// 'kg' /)
  !
  !*** Deallocate memory
  !
  deallocate(pb)
  deallocate(MFRpb)

  !
  return    
end subroutine source

