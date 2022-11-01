   subroutine readinp
      !**************************************************************
      !*
      !*  Opens and reads the input file
      !*
      !**************************************************************
      use KindType
      use InpOut
      use Master
      use TimeFun
      implicit none

      integer(ip)         :: istat
      integer(ip)         :: ivoid(1)
      !integer  (ip)         :: iyr, imo, idy, ihr, imi, ise              ! Times needed for addtime
      integer(ip)         :: hrst, mist, sest, hrend, miend, seend       ! hours, minutes and seconds times for simulation
      !
      character(len=s_mess) :: message, cvoid(2)
      !
      real(rp)              :: rvoid(1)
      real(rp)              :: work(2)
      integer(ip)              :: i, iz, nz
      real(rp), allocatable :: work2(:)
      real(rp)              :: nout         ! auxiliar variable. number of output
      !
      ! Granulometry/source variables
      !
      integer(ip)         :: iphase
      character(len=1)      :: ipha
      character(len=s_mess) :: phase_str
      !
      !*** 1. Reads the SIMULATION_TIME Block. Integer values: year, month and day
      !
      call get_input_int(finp, 'SIMULATION_TIME', 'YEAR', ivoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      ibyr = ivoid(1)
      !
      call get_input_int(finp, 'SIMULATION_TIME', 'MONTH', ivoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      ibmo = ivoid(1)
      !
      call get_input_int(finp, 'SIMULATION_TIME', 'DAY', ivoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      ibdy = ivoid(1)
      !
      call get_input_rea(finp, 'SIMULATION_TIME', 'SIMULATION_START', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      sist = rvoid(1)
      ibhr = INT(sist)                     ! hours. (hh) utilized in other subroutines (readrestart, debug)
      hrst = int(sist)                     ! To calculate times with addtime
      mist = int((sist - hrst)*60)
      sest = int(((sist - hrst)*60 - mist)*60)
      sist = sist*3600                     ! Time since 0000 UTC in seconds
      !
      call get_input_rea(finp, 'SIMULATION_TIME', 'SIMULATION_END', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      siend = rvoid(1)
      siend = siend*3600    ! Time since 0000 UTC in seconds
      !
      call get_input_int(finp, 'SIMULATION_TIME', 'TIME_STEP', ivoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      sidt = ivoid(1)
      if (sist .gt. siend) sidt = (-1)*sidt  ! backwards case: sidt is negative
      !
      call get_input_cha(finp, 'SIMULATION_TIME', 'RESTART', cvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      call upcase(cvoid(1))
      if (cvoid(1) .eq. 'YES') restart = .true.
      !
      !*** 2. Calculate time_start and  time_end in the form YYYYMMDDHHMM
      !
      if (sidt .gt. 0) then
         time_start = 1d8*ibyr + 1d6*ibmo + 1d4*ibdy + 1d2*hrst + mist
         sit = siend - sist + mist + sest         ! Total time in seconds beginning in hrst. used in addtime
         call addtime(ibyr, ibmo, ibdy, hrst, fiyr, fimo, fidy, fihr, fimi, fise, sit)
         time_end = 1d8*fiyr + 1d6*fimo + 1d4*fidy + 1d2*fihr + fimi
      else if (sidt .lt. 0) then
         time_start = 1d8*ibyr + 1d6*ibmo + 1d4*ibdy + 1d2*hrst + mist
         sit = -siend + sist - mist - sest            ! Total time in seconds beginning in hrst
         call subtime(ibyr, ibmo, ibdy, hrst, fiyr, fimo, fidy, fihr, fimi, fise, sit)
         time_end = 1d8*fiyr + 1d6*fimo + 1d4*fidy + 1d2*fihr + fimi
      end if
      !
      !*** 3. Reads the COMPUTATIONAL_DOMAIN BLOCK
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'LATMAX', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%latmax = rvoid(1)
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'LATMIN', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%latmin = rvoid(1)
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'LONMAX', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%lonmax = rvoid(1)
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'LONMIN', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%lonmin = rvoid(1)
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'ZTOP', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%ztop = rvoid(1)
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'LONGITUDE_RESOLUTION', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%dx = rvoid(1)
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'LATITUDE_RESOLUTION', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%dy = rvoid(1)
      !
      call get_input_rea(finp, 'COMPUTATIONAL_DOMAIN', 'VERTICAL_RESOLUTION', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      grid%dz = rvoid(1)
      !
      !*** Check boundary condition
      !
      if ((grid%lonmax .eq. 180) .and. (grid%lonmin .eq. -180)) then  ! east-west boundary condition
         grid%PBC_EW = .true.
         if (grid%latmax .eq. 90) grid%PBC_N = .true.  ! north boundary condition
         if (grid%latmin .eq. -90) grid%PBC_S = .true.  ! south boundary condition
      end if
      !
      !*** 4. Reads the OUTPUT_GRID block the time to output meteo.
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'OUTPUT_LATMAX', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      output%latmax = rvoid(1)
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'OUTPUT_LATMIN', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      output%latmin = rvoid(1)
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'OUTPUT_LONMAX', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      output%lonmax = rvoid(1)
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'OUTPUT_LONMIN', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      output%lonmin = rvoid(1)
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'OUTPUT_FREQUENCY', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      output%frequency = rvoid(1)*3600
      nout = ABS(siend - sist)/(rvoid(1)*3600)
      if (INT(nout) .eq. nout) then
         output%nt = INT(nout) + 1
      else
         output%nt = INT(nout) + 2
      end if
      !
      allocate (output%timesec(output%nt))
      !
      ! Forward case
      if (sidt .gt. 0) then !forward case
         do i = 1, output%nt - 1
            output%timesec(i) = sist + (i - 1)*output%frequency! Determine time in seconds since 00:00
         end do
         output%timesec(output%nt) = siend                  ! final output
         !
         ! Backward case
      else
         do i = 1, output%nt - 1
            output%timesec(i) = sist - (i - 1)*output%frequency! Determine time in seconds since 00:00
         end do
         output%timesec(output%nt) = siend                  ! final output
      end if
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'LONGITUDE_RESOLUTION', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      output%dx = rvoid(1)
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'LATITUDE_RESOLUTION', rvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      output%dy = rvoid(1)
      !
      call get_input_cha(finp, 'OUTPUT_GRID', 'OUTPUT_CLASSES', cvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      call upcase(cvoid(1))
      if (TRIM(cvoid(1)) .eq. 'YES') then
         output%classes = .true.
      else
         output%classes = .false.
      end if
      !
      call get_input_cha(finp, 'OUTPUT_GRID', 'OUTPUT_PHASES', cvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      call upcase(cvoid(1))
      if (TRIM(cvoid(1)) .eq. 'YES') then
         output%phases = .true.
      else
         output%phases = .false.
      end if
      !
      call get_input_cha(finp, 'OUTPUT_GRID', 'OUTPUT_TRACK_POINTS', cvoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      call upcase(cvoid(1))
      if (TRIM(cvoid(1)) .eq. 'YES') then
         output%track_points = .true.
      else
         output%track_points = .false.
      end if
      !
      !*** Give values to output file dimensions
      !
      output%ztop = grid%ztop
      !
      call get_input_npar(finp, 'OUTPUT_GRID', 'VERTICAL_LAYERS', nz, istat, message)
      if (istat .gt. 0) call wriwar(message)
      if (istat .lt. 0) call runend(message)
      !
      allocate (work2(nz))
      !
      call get_input_rea(finp, 'OUTPUT_GRID', 'VERTICAL_LAYERS', work2, nz, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      !
      !*** Save z-layers for the output file. given as one dz value, or enumerated each z-level
      !
      if (nz == 1) then                             ! One value are readed in inp file: dz
         output%nz = INT(output%ztop/work2(1))   !=work2(1)
         allocate (output%zlayer(output%nz))
         !
         do iz = 1, output%nz
            output%zlayer(iz) = (iz)*work2(1)  ! the last level is not ztop necesarily (saco el (-1) en el iz, no incluyo h=0)
         end do
      else                                      ! Are readed each z-level
         output%nz = nz
         allocate (output%zlayer(output%nz))
         output%zlayer = work2
      end if
      deallocate (work2)
      !
      !*** 5. Reads the PHYSICS Block
      !
      call get_input_int(finp, 'PHYSICS', 'TERMINAL_VELOCITY_MODEL', ivoid, 1, istat, message)
      if (istat > 0) call wriwar(message)
      if (istat < 0) call runend(message)
      sedmodel = ivoid(1)
      !
      !*** Calculate the number of meteo models (numod)
      !
      call readinp_meteo
      !
      !*** Reads the number and characteristics of each phase
      !
      call readinp_phases
      !
      !*** Reads the tracking points if it is neccesary
      !
      if (output%track_points) then
         call readinp_trackpts
      end if
      !
   end subroutine readinp
!
!
!
   subroutine readinp_phases
      !**************************************************************
      !*
      !*   Reads the different eruption phases
      !*
      !**************************************************************
      use KindType
      use InpOut
      use Master
      implicit none
      !
      logical               :: go_on, active(nparmax), goon
      integer(ip)         :: istat
      integer(ip)         :: ivoid(1), nword, npar, iphases, ig
      character(len=s_long) :: card
      character(len=s_long) :: words(nwormax)
      character(len=s_mess) :: message, cvoid(2)
      real(rp)              :: param(nparmax)
      integer(ip)           :: isp                ! i-sub-phase
      character(len=s_name) :: path
      real(rp)              :: zerosubphase
      real(rp)              :: iphaseduration !auxiliary variables to calculate the phase duration while particles must be liberated
      !
      !*** First, averiguates the number of phases (and how many are activated)
      !
      nphases = 0    ! total number of active phases
      iphases = 0    ! total number of phases
      !
      open (90, FILE=TRIM(finp), status='unknown')
      do while (0 .eq. 0)
         !
         read (90, '(a256)', END=100) card
         call sdecode(card, words, param, nword, npar)
         !
         if (TRIM(words(1)) .eq. 'PHASE_DEFINITION') then
            !
            iphases = iphases + 1
            go_on = .true.
            do while (go_on)                              ! checks if the phase is active
               read (90, '(a256)', END=103) card            ! if the input file is incomplete witout activate nor end_phase_definition words.
               call sdecode(card, words, param, nword, npar)
               if (TRIM(words(1)) .eq. 'ACTIVATE') then
                  call upcase(words(2))
                  if (TRIM(words(2)) .eq. 'YES') then
                     nphases = nphases + 1
                     active(iphases) = .true.
                  else
                     active(iphases) = .false.
                  end if
                  go_on = .false.
               else if (TRIM(words(1)) .eq. 'END_PHASE_DEFINITION') then
                  call runend('readinp_phases: ACTIVATE block not found')
               end if
            end do
            !
         end if
      end do
      !
100   close (90)
      !
      !*** Check if it is necesary sum one phase for restart file. To allocate
      !*** phase structure must take into account the restart phase
      !
      if (restart) nphases = nphases + 1
      if (sidt .lt. 0) nphases = 1        !backward mode: only one phase is considered.
      !one input file is readed, with all the particle characteristics.
      !
      !*** Allocates memory
      !
      allocate (phase(nphases))
      !
      !*** Reads again the input file for the active phases
      !*** Open the .Phase.inp files for activated phases
      if (sidt .gt. 0) then ! Only if the file exists: only in forward case.
         iphases = 0
         nphases = 0
         open (90, FILE=TRIM(finp), status='unknown')
         do while (0 .eq. 0)
            !
            read (90, '(a256)', END=200) card
            call sdecode(card, words, param, nword, npar)
            !
            if (TRIM(words(1)) .eq. 'PHASE_DEFINITION') then
               !
               iphases = iphases + 1
               if (active(iphases)) then
                  nphases = nphases + 1
                  go_on = .true.
                  do while (go_on)                  !Read the path for the corresponding file
                     read (90, '(a256)', END=101) card
                     call sdecode(card, words, param, nword, npar)
                     if (TRIM(words(1)) .eq. 'INCLUDE') then
                        path = words(2)
                        go_on = .false.
                     end if
                  end do
                  !
                  !*** Open and read each .Phase.inp file
                  !
                  open (91, FILE=TRIM(path), status='unknown')  ! Open the .Phase.inp file
                  do while (0 .eq. 0)
                     read (91, '(a256)', END=201) card          ! Read all data. End: close 91 file and continue with .inp file
                     call sdecode(card, words, param, nword, npar)
                     !if(TRIM(words(1)).eq.'END_PHASE_DEFINITION') then
                     !   go_on = .false.
                     !else
                     if (TRIM(words(1)) .eq. 'PHASE_NAME') then
                        phase(nphases)%sname = words(2)
                     else if (TRIM(words(1)) .eq. 'NUMBER_PARTICLES') then
                        phase(nphases)%npart = INT(param(1))
                        if (my_id .eq. 0) write (lulog, 4) nphases, phase(nphases)%npart
4                       format(/, 'The total number of particles to simulate in phase        ', i2, '        is     : ', i10)
                     else if (TRIM(words(1)) .eq. 'PHASE_TYPE') then
                        phase(nphases)%phase_type = words(2)
                     else if (TRIM(words(1)) .eq. 'INITIAL_TIME') then
                        ! save each value
                        allocate (phase(nphases)%begtime(npar))
                        phase(nphases)%nsp = npar                             ! nsp: number of sub phases
                        phase(nphases)%begtime(1:npar) = param(1:npar)*3600 ! convert to second
                     else if (TRIM(words(1)) .eq. 'END_TIME') then
                        phase(nphases)%endtime = param(1)*3600              ! convert to second
                     else if (TRIM(words(1)) .eq. 'SOURCE_TYPE') then
                        call upcase(words(2))
                        phase(nphases)%source_type = words(2)
                        if (TRIM(phase(nphases)%source_type) .eq. 'POINT') then
                           phase(nphases)%source_type = 'POINT'
                        else if (TRIM(phase(nphases)%source_type) .eq. 'SUZUKI') then
                           phase(nphases)%source_type = 'SUZUKI'
                        else if (TRIM(phase(nphases)%source_type) .eq. 'LINEAR') then
                           phase(nphases)%source_type = 'LINEAR'
                        else if (TRIM(phase(nphases)%source_type) .eq. 'TOP-HAT') then
                           phase(nphases)%source_type = 'TOP-HAT'
                        else if (TRIM(phase(nphases)%source_type) .eq. 'RESUSPENSION') then
                           resuspension = .true.
                        else
                           call runend('Incorrect source type')
                        end if
                     else if (TRIM(words(1)) .eq. 'COLUMN_HEIGHT') then
                        ! save each value
                        allocate (phase(nphases)%colheight(phase(nphases)%nsp))
                        if (phase(nphases)%nsp .le. npar) then
                           phase(nphases)%colheight(1:phase(nphases)%nsp) = param(1:phase(nphases)%nsp)
                        else                              ! If there are less parameters than sub-phases
                           phase(nphases)%colheight(1:npar) = param(1:npar)
                           phase(nphases)%colheight(npar:phase(nphases)%nsp) = param(npar)
                        end if
                        !
                     else if (TRIM(words(1)) .eq. 'MASS_FLOW_RATE') then
                        !
                        call upcase(words(2))
                        if (TRIM(words(2)) .eq. 'ESTIMATE-MASTIN') then
                           !
                           phase(nphases)%mer_vs_h = 'ESTIMATE-MASTIN'
                           phase(nphases)%MER_wind_coupling = .false.
                           !
                        else if (words(2) .eq. 'ESTIMATE-DEGRUYTER') then
                           !
                           phase(nphases)%mer_vs_h = 'ESTIMATE-DEGRUYTER'
                           phase(nphases)%MER_wind_coupling = .true.
                           !
                        else if (words(2) .eq. 'ESTIMATE-WOODHOUSE') then
                           !
                           phase(nphases)%mer_vs_h = 'ESTIMATE-WOODHOUSE'
                           phase(nphases)%MER_wind_coupling = .true.
                           !
                        else
                           !
                           phase(nphases)%mer_vs_h = 'NONE'
                           phase(nphases)%MER_wind_coupling = .false.
                           ! Allocate the nsp parameters given
                           allocate (phase(nphases)%M0(phase(nphases)%nsp))
                           if (phase(nphases)%nsp .le. npar) then
                              phase(nphases)%M0(1:phase(nphases)%nsp) = param(1:phase(nphases)%nsp)
                           else                              ! If there are less parameters than sub-phases
                              phase(nphases)%M0(1:npar) = param(1:npar)
                              phase(nphases)%M0(npar:phase(nphases)%nsp) = param(npar)
                           end if
                        end if
                     else if (TRIM(words(1)) .eq. 'A_SUZUKI') then
                        allocate (phase(nphases)%a_suzuki(phase(nphases)%nsp))
                        if (phase(nphases)%nsp .le. npar) then
                           phase(nphases)%a_suzuki(1:phase(nphases)%nsp) = param(1:phase(nphases)%nsp)
                        else                              ! If there are less parameters than sub-phases
                           phase(nphases)%a_suzuki(1:npar) = param(1:npar)
                           phase(nphases)%a_suzuki(npar:phase(nphases)%nsp) = param(npar)
                        end if
                     else if (TRIM(words(1)) .eq. 'L_SUZUKI') then
                        allocate (phase(nphases)%l_suzuki(phase(nphases)%nsp))
                        if (phase(nphases)%nsp .le. npar) then
                           phase(nphases)%l_suzuki(1:phase(nphases)%nsp) = param(1:phase(nphases)%nsp)
                        else                              ! If there are less parameters than sub-phases
                           phase(nphases)%l_suzuki(1:npar) = param(1:npar)
                           phase(nphases)%l_suzuki(npar:phase(nphases)%nsp) = param(npar)
                        end if
                     else if (TRIM(words(1)) .eq. 'D_TOP_HAT') then
                        allocate (phase(nphases)%d_top_hat(phase(nphases)%nsp))
                        if (phase(nphases)%nsp .le. npar) then
                           phase(nphases)%d_top_hat(1:phase(nphases)%nsp) = param(1:phase(nphases)%nsp)
                        else                              ! If there are less parameters than sub-phases
                           phase(nphases)%d_top_hat(1:npar) = param(1:npar)
                           phase(nphases)%d_top_hat(npar:phase(nphases)%nsp) = param(npar)
                        end if
                     else if (TRIM(words(1)) .eq. 'VOLCANO_NAME') then
                        phase(nphases)%volname = words(2)
                     else if (TRIM(words(1)) .eq. 'SOURCE_LONGITUDE') then
                        phase(nphases)%lon = param(1)
                        if ((phase(nphases)%lon .lt. grid%lonmin) .or. (phase(nphases)%lon .gt. grid%lonmax)) &
                           call runend('The Source Longitude is not inside the Computational Domain')
                     else if (TRIM(words(1)) .eq. 'SOURCE_LATITUDE') then
                        phase(nphases)%lat = param(1)
                        if ((phase(nphases)%lat .lt. grid%latmin) .or. (phase(nphases)%lat .gt. grid%latmax)) &
                           call runend('The Source Latitude is not inside the Computational Domain')
                     else if (TRIM(words(1)) .eq. 'SOURCE_ELEVATION') then
                        phase(nphases)%source_elev = param(1)
                     else if (TRIM(words(1)) .eq. 'PHASE_GRANULOMETRY') then
                        phase(nphases)%grnpath = words(2)
                     else if (TRIM(words(1)) .eq. 'DISTRIBUTION') then
                        phase(nphases)%grndist = words(2)
                        if (TRIM(phase(nphases)%grndist) .eq. 'GAUSSIAN') then
                           phase(nphases)%mdist = 0
                           phase(nphases)%ng = 1
                        else if (TRIM(phase(nphases)%grndist) .eq. 'BIGAUSSIAN') then
                           phase(nphases)%mdist = 1
                           phase(nphases)%ng = 2
                        end if
                     else if (TRIM(words(1)) .eq. 'NUMBER_OF_BINS') then
                        phase(nphases)%bins = INT(param(1))
                     else if (TRIM(words(1)) .eq. 'FI_MEAN') then
                        do ig = 1, phase(nphases)%ng
                           phase(nphases)%fimean(ig) = param(ig)
                        end do
                     else if (TRIM(words(1)) .eq. 'FI_DISP') then
                        do ig = 1, phase(nphases)%ng
                           phase(nphases)%fidisp(ig) = param(ig)
                        end do
                     else if (TRIM(words(1)) .eq. 'FI_RANGE') then
                        phase(nphases)%fimin = param(1)
                        phase(nphases)%fimax = param(2)
                     else if (TRIM(words(1)) .eq. 'DENSITY_RANGE') then
                        phase(nphases)%rhomin = param(1)
                        phase(nphases)%rhomax = param(2)
                     else if (TRIM(words(1)) .eq. 'SPHERICITY_RANGE') then
                        phase(nphases)%sphemin = param(1)
                        phase(nphases)%sphemax = param(2)
                     else if (TRIM(words(1)) .eq. 'AGGREGATION_MODEL') then
                        call upcase(words(2))
                        phase(nphases)%aggmodel = words(2)
                        if (phase(nphases)%aggmodel .eq. 'PERCENTAGE') then
                           phase(nphases)%aggregation = .true.
                        else if (phase(nphases)%aggmodel .eq. 'CORNELL') then
                           phase(nphases)%aggregation = .true.
                        else if (phase(nphases)%aggmodel .eq. 'NONE') then
                           phase(nphases)%aggregation = .false.
                        else
                           call runend('readinp: Aggregation model is not valid')
                        end if
                        !
                        if (phase(nphases)%source_type .eq. 'RESUSPENSION') phase(nphases)%aggregation = .false.  ! Switch-off aggregation in resuspension mode
                        if (phase(nphases)%aggregation) then
                           goon = .true.
                           do while (goon)
                              read (91, '(a256)', END=104) card     !If there is no END_AGGREGATION statement, input file incomplete.
                              call sdecode(card, words, param, nword, npar)
                              if (TRIM(words(1)) .eq. 'END_AGGREGATION') then
                                 goon = .false.                                       ! final aggregation block
                              else if (TRIM(words(1)) .eq. 'AGGREGATE_SIZE') then
                                 phase(nphases)%aggsize = param(1)
                                 phase(nphases)%aggsize = phase(nphases)%aggsize/1d6  ! microns -->m
                              else if (TRIM(words(1)) .eq. 'AGGREGATE_DENSITY') then
                                 phase(nphases)%aggrho = param(1)
                              else if (TRIM(words(1)) .eq. 'PERCENTAGE_(%)') then
                                 phase(nphases)%aggfrac = param(1)
                                 phase(nphases)%aggfrac = phase(nphases)%aggfrac/100.0_rp   ! [0,1]
                              end if
                           end do
                        end if
                        !
                        !*** Check for inconsistencies
                        !
                        !*** Check that the initial and final time are inside the simulation times
                        !
                        if ((phase(nphases)%begtime(1) .lt. sist) .or. (phase(nphases)%endtime .gt. siend)) &
                           call runend('readinp: The times of a Source Phase is not inside the simulation times')
                        if (phase(nphases)%begtime(phase(nphases)%nsp) .gt. phase(nphases)%endtime) &
                           call runend('readinp: The times of a Source Phase are not consistents')
                        !
                        !*** Check that the column height is inside the meteorological domain
                        !
                        zerosubphase = 0.
                        do isp = 1, phase(nphases)%nsp
                           if ((grid%ztop) .lt. (phase(nphases)%colheight(isp))) &
                              call runend('The column height is not inside the Computational Domain')
                           if (phase(nphases)%colheight(isp) .le. 0) then
                              if (isp .eq. phase(nphases)%nsp) then
                                 iphaseduration = phase(nphases)%endtime - phase(nphases)%begtime(isp)
                              else
                                 iphaseduration = phase(nphases)%begtime(isp + 1) - phase(nphases)%begtime(isp)
                              end if
                              zerosubphase = zerosubphase + iphaseduration
                           end if
                        end do
                        ! Calculate phases duration
                        phase(nphases)%duration = phase(nphases)%endtime - phase(nphases)%begtime(1) - zerosubphase
                     end if
                  end do
201               close (91)
               end if      ! if(active(iphases))
            end if         ! if(TRIM(words(1)).eq.'PHASE_DEFINITION')
         end do

200      close (90)
      end if
      return
      !
      !*** List of errors
      !
101   call runend('readinp_phases: ACTIVATE path not found')
103   call runend('readinp_phases: There is a PHASE_DEFINITION block without ACTIVATE nor END_PHASE_DEFINITION')
104   call runend('readinp_phases: Is not found the END_AGGREGATION statement')
      !
      return
   end subroutine readinp_phases
!
!
!
   subroutine readinp_meteo
      !**************************************************************
      !*
      !*   Reads the different Meteo information
      !*
      !**************************************************************
      use KindType
      use InpOut
      use Master
      implicit none
      !
      logical               :: go_on, active(nparmax)
      integer(ip)         :: istat
      integer(ip)         :: ivoid(1), nword, npar, imodels
      character(len=s_long) :: card
      character(len=s_long) :: words(nwormax)
      real(rp)              :: param(nparmax)
      !
      !*** First, averiguates the number of Meteorological data (and how many are activated)
      !
      numod = 0    ! total number of active Meteo Models
      imodels = 0    ! total number of meteo models
      !
      open (90, FILE=TRIM(finp), status='unknown')
      do while (0 .eq. 0)
         !
         read (90, '(a256)', END=100) card
         call sdecode(card, words, param, nword, npar)
         !
         if (TRIM(words(1)) .eq. 'METEO_MODEL_DEFINITION') then
            !
            imodels = imodels + 1
            go_on = .true.
            do while (go_on)                              ! checks if the meteo is active
               read (90, '(a256)', END=101) card
               call sdecode(card, words, param, nword, npar)
               if (TRIM(words(1)) .eq. 'ACTIVATE') then
                  call upcase(words(2))
                  if (TRIM(words(2)) .eq. 'YES') then
                     numod = numod + 1
                     active(imodels) = .true.
                  else
                     active(imodels) = .false.
                  end if
                  go_on = .false.
               else if (TRIM(words(1)) .eq. 'END_METEO_MODEL_DEFINITON') then
                  call runend('readinp_meteo: ACTIVATE block not found')
               end if
            end do
            !
         end if
      end do
      !
100   close (90)
      !
      !*** Allocates memory
      !
      allocate (meteo(numod))
      !
      !*** Reads again the file for the active Meteo
      !
      numod = 0
      imodels = 0
      open (90, FILE=TRIM(finp), status='unknown')
      do while (0 .eq. 0)
         !
         read (90, '(a256)', END=200) card
         call sdecode(card, words, param, nword, npar)
         !
         if (TRIM(words(1)) .eq. 'METEO_MODEL_DEFINITION') then
            !
            imodels = imodels + 1
            if (active(imodels)) then
               numod = numod + 1     ! numod: number of meteo models.
               go_on = .true.
               do while (go_on)
                  read (90, '(a256)', END=101) card
                  call sdecode(card, words, param, nword, npar)
                  if (TRIM(words(1)) .eq. 'END_METEO_MODEL_DEFINITON') then
                     go_on = .false.
                  else if (TRIM(words(1)) .eq. 'MODEL_TYPE') then
                     meteo(numod)%modeltype = words(2)
                     call upcase(meteo(numod)%modeltype)
                  else if (TRIM(words(1)) .eq. 'FILE') then
                     meteo(numod)%file_path = words(2)
                  else if (TRIM(words(1)) .eq. 'POSTPROCESS') then
                     call upcase(words(2))
                     if (words(2) .eq. 'YES') then
                        meteo(numod)%postprocess = .true.
                     else
                        meteo(numod)%postprocess = .false.
                     end if
                  end if
               end do
            end if
         end if
      end do

200   close (90)
      return
      !
      !*** List of errors
      !
101   call runend('readinp_meteo: ACTIVATE path not found')
      !
      return
   end subroutine readinp_meteo
!
!
!
   subroutine readinp_trackpts
      !**************************************************************
      !*
      !*   Reads the Trackin points information
      !*
      !**************************************************************
      use KindType
      use InpOut
      use Master
      implicit none
      integer(ip)::i, istat
      character(len=s_mess) :: message

      !
      !***    Get the number of points
      !
      call get_points_npts(fpts, npts, istat, message)
      !
      !***    Allocates memory for vectors related to npts (time dependent)
      !
      allocate (name_pts(npts))
      name_pts = ' '
      !
      allocate (xpts(npts))
      xpts = 0.0_rp
      !
      allocate (ypts(npts))
      ypts = 0.0_rp
      !
      allocate (use_pts(npts))
      use_pts = .false.
      !
      allocate (ielem2dpts(npts))
      ielem2dpts = 0.0_rp
      !
      allocate (name_file_pts(npts))
      name_file_pts = ' '
      !
      !***    Loads vectors
      !
      call get_points_coordinates(fpts, name_pts, xpts, ypts, npts, istat, message)
      !
      !***    Averiguates which points lay within the computational domain
      !
      do i = 1, npts
         if ((xpts(i) >= output%lonmin) .and. (xpts(i) <= output%lonmax) .and. &
             (ypts(i) >= output%latmin) .and. (ypts(i) <= output%latmax)) use_pts(i) = .true.
      end do
      !
      !***    Set file names
      !
      call get_pts_fname
      !
      !
      return
   end subroutine readinp_trackpts
!
!
!
   subroutine get_pts_fname
      !*******************************************************************
      !*
      !*   Gets the tracking points file names
      !*
      !*******************************************************************
      use KindType
      use InpOut
      implicit none
      !
      logical     :: go_on
      integer(ip) :: ilen, i
      !
      !*** Remove extension
      !
      i = LEN(fpts) + 1
      go_on = .true.
      do while (go_on)
         i = i - 1
         if (fpts(i:i) == '.' .or. i == 2) then
            go_on = .false.
            ilen = i - 1
         end if
      end do
      !
      do i = 1, npts
         name_file_pts(i) = fpts(1:ilen)//'.tps.'//TRIM(name_pts(i))  ! no extension
!     open(15,FILE=TRIM(name_file_pts(i)),status='REPLACE')
!     write(15,10) TRIM(name_pts(i)),xpts(i),ypts(i)
!10   format(&
!     'Tracking point file for : ',a              ,/, &
!     'Coordinates             : ',f13.4,1x,f13.4 ,/, &
!     '  Time      load     thickness    conc.       cumul',/, &
!     '           ground   (compact=1)   ground     conc. ',/, &
!     '  (min)    (kg/m2)     (cm)       (g/m3)     (g/m2)',/, &
!     '------------------------------------------------------------------------------------------------')
      end do
      !
      return
   end subroutine get_pts_fname
!
