program ATLAS
   !******************************************************************************
   !*
   !*    ATmospheric LAgrangian diSpersion model (ATLAS)
   !*
   !*    Lagrangian model to calculate particle dispersion, sedimentation
   !*    and accumulation at different levels and on the ground
   !*
   !*    Paralel Version.
   !*    Author: Florencia Reckziegel & Arnau Folch & Oscar Reula
   !*    Version 1.0 (May 2014-Oct 2016)
   !*
   !******************************************************************************
   use KindType
   use InpOut
   use Master
   use TimeFun
   use MPI
   implicit none
   !
   !
   !
   logical         :: go_on
   integer(ip)   :: i, iphase, iclass, ipart
   integer(ip)   :: iyr, imo, idy, ihr, imi, ise ! needed for addtime subroutine
   character(len=4):: ciyr                    ! Date variables
   character(len=2):: cimo, cidy, cihr, cimi, cise! month day.. to calculate the actual time inside the loop over time
   integer(ip)   :: npart2                  ! Particle range to consider in the main loop. per phase.
   logical         :: finaloutput = .true.      ! is false if the end_time concur with one output time
   real(rp)   :: uold, vold, wold          ! old velocities, auxiliary variables
   !
   ! Inicializo el mpi
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
   mpi_global = MPI_COMM_WORLD
   !
   !
   !*** Gets input problemname from the call argument
   !
   call GETARG(1, problemname)
   if (TRIM(problemname) .eq. '') then
      stop 'Specify the problemname'
   end if
   !
   !*** Creates the log file
   !
   flog = TRIM(problemname)//'.log'
   if (my_id .eq. 0) call openlog
   !
   !*** Initialization
   !
   call initialize
   !
   !*** Loop over time
   !
   go_on = .true.
   iiter = 0           ! count the number of iteration
   time = sist
   !
   do while (go_on)
      !
      !*** Determine the current time
      !
      call addtime(ibyr, ibmo, ibdy, 0, iyr, imo, idy, ihr, imi, ise, time)
      write (ciyr, '(i4.4)') iyr
      write (cimo, '(i2.2)') imo
      write (cidy, '(i2.2)') idy
      write (cihr, '(i2.2)') ihr
      write (cimi, '(i2.2)') imi
      write (cise, '(i2.2)') ise
      actualtime = ciyr//'/'//cimo//'/'//cidy//'-'//cihr//':'//cimi//':'//cise
      !
      !
      !*** Read Meteo data (only if necessary) and interpolates in time for each time step
      !
      call meteos
      !
      !*** Loop over phases
      !
      do iphase = 1, nphases
         !
         !*** If necessary, activate the phase
         !
         if (time .ge. phase(iphase)%begtime(1)) phase(iphase)%activated = .true.
         !
         !*** For active phases only
         !
         if (phase(iphase)%activated) then
            !
            !*** Check if it is necesary liberate more particles at time. Only in forward mode.
            !
            if ((sidt .gt. 0) .and. (time .le. phase(iphase)%endtime)) then
               !
               !*** SOURCE term (only for particles released in the current dt)
               !*** Particles are distributed in vertical along the vent
               !
               call source(iphase, npart2)
            end if
            !
            !*** Determine particles range in restart case
            !
            !
            if (phase(iphase)%phase_type == 'restart') then
               npart2 = phase(iphase)%lastpart
            else if (phase(iphase)%phase_type == 'backward') then
               npart2 = phase(iphase)%lastpart
            end if
            !
            !*** check if it is necessary split particles
            !
            if ((split) .and. (time .gt. dts)) then
               call splitting(iphase, npart2)
            end if
            !
            !*** Loop over particles inside the current source phase
            !
            do ipart = phase(iphase)%firstpart, npart2
               !
               !*** Integrate particle movement (only if the particle is in the air and inside the domain)
               !
               if (part(ipart)%state .eq. 1) then
                  !Set the particle age
                  part(ipart)%age = part(ipart)%age + sidt
                  !
                  !*** Interpolate meteo to the particle coordinates
                  !
                  call interpolate_meteo(ipart, part(ipart)%lon, part(ipart)%lat, part(ipart)%z)
                  !
                  !*** Computes particle setling velocity
                  !
                  call sedimentation(part(ipart)%vset, part(ipart)%rho, part(ipart)%rhoa, part(ipart)%diam, &
                                     part(ipart)%mua, part(ipart)%sphe, part(ipart)%psi, sedmodel, sidt)
                  !
                  !*** Turbulence disturbance.
                  !
                  !**First, save the grid-winds and settlings velocity to compute Pettersen scheme later
                  part(ipart)%wa = -part(ipart)%wa + part(ipart)%vset !wa is positive upward, vset is positive down. now wa is positive down
                  !
                  if (sidt .gt. 0) then !! Only in forward mode
                     call diffusion(part(ipart)%ua, part(ipart)%va, part(ipart)%wa, part(ipart)%z, &
                                    part(ipart)%pblha, part(ipart)%rmonin, part(ipart)%usta, part(ipart)%wst, &
                                    part(ipart)%drhodz, part(ipart)%rhoa, part(ipart)%lon, part(ipart)%lat, &
                                    part(ipart)%uvw(1), part(ipart)%uvw(2), part(ipart)%uvw(3), part(ipart)%age)
                  end if
                  !
                  !*** Wind advection plus particle gravity sedimentation and turbulence
                  !
                  call advection(part(ipart)%z, part(ipart)%lon, part(ipart)%lat, part(ipart)%ua, &
                                 part(ipart)%va, part(ipart)%wa, &
                                 part(ipart)%rho, part(ipart)%diam, part(ipart)%sphe, part(ipart)%psi)
                  !
                  !*** Other deposition mechanisms  (dry/wet deposition)
                  !
                  !
                  !*** Check if the particle is outside the domain, then freeze it or change the latitude and/or longitude
                  !*** only if the background grid and meteo model contemplate:
                  !
                  call checkdomain(part(ipart)%state, part(ipart)%lon, part(ipart)%lat, part(ipart)%z)
                  !
                  ! si resuspension = .true.
                  !
               end if  ! Active particles
            end do     ! Particles
            !
         end if        ! if(phase(iphase)%activated)  (active phases)
      end do           ! iphase = 1,nphases           (phase loop)
      !
      !*** Check if it is necessary to generate output.
      !
      if (sidt .gt. 0) then !Forward case
         do i = 1, output%nt
            if ((time .le. output%timesec(i)) .and. (time + sidt .gt. output%timesec(i))) then
               call make_part_nc(i)
               call reduce_mass
!                        include 'reduce_mass_inc.f90'
!             call make_part_kml
               if (my_id .eq. 0) then
                  call outpart_nc(i)   !saque npart2, porque el outpart_nc no lo usa
!              call outpart_kml(i)
               end if
               if (i .eq. output%nt) finaloutput = .false.
            end if
         end do
      else              !Backward case
         do i = 1, output%nt
            if ((time .ge. output%timesec(i)) .and. (time + sidt .lt. output%timesec(i))) then
               if (my_id .eq. 0) then
                  call outpart_nc(i)        !saque npart2, porque el outpart_nc no lo usa
!              call outpart_kml(i)
               end if
               if (i .eq. output%nt) finaloutput = .false.
            end if
         end do
      end if
      !
      !*** Update counter
      !
      iiter = iiter + 1
      time = time + sidt
      if ((sidt .gt. 0 .and. time >= siend) .or. (sidt .lt. 0 .and. time <= siend)) go_on = .false.
      !
   end do  !Time loop
   !
   !*** Write final output
   !
   if (finaloutput) then
      if (my_id .eq. 0) then
         call outpart_nc(output%nt)
!        call outpart_kml(output%nt)
      end if
   end if
   !
   !*** Writes the restart file
   !
   call outrest
   !
   if (my_id .eq. 0) then
      call writerest
   end if
   ! Termino la paralelizacion
   call MPI_FINALIZE(ierr)
   !
   !*** Ends the run
   !
   call runend('OK')
   !
end program ATLAS

