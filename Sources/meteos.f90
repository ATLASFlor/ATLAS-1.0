subroutine meteos
   !******************************************************************************************
   !*    METEOS
   !*    Determines time steps
   !*    Reads a meteo file when it is necesary
   !*    Interpolates variables in time
   !*
   !******************************************************************************************
   use InpOut
   use Master
   use KindType
   use Netcdf
   implicit none
   !
   integer(ip)       :: imod
   integer(ip), save :: ipass = 0
   integer(ip), save :: ioutput = 0
   logical            :: found
   logical            :: outputmet
   integer(ip)       :: it_nm, ipoin, ix, iy, iz
   !
   character(len=s_mess) :: str
   character(len=1) :: ext
   integer(ip)           ::ncid
   !
   !***  Initializations
   !
   outputmet = .false. ! is true only when new information is readed
   !
   if (ipass .eq. 0) then
      !
      !*** For the first time only, found meteo(imod)%itime1 and meteo(imod)%itime2, according all the meteo models
      !
      do imod = 1, numod
         write (ext, '(i1)') imod
         str = 'METEO_DATA_'//ext
         !
         if (meteo(imod)%nt .gt. 1) then
            found = .false.
            it_nm = 0
            if (sidt .gt. 0) then
               do while (.not. found)
                  it_nm = it_nm + 1
                  if ((it_nm + 1) .gt. meteo(imod)%nt) call runend('Time not found in'//TRIM(str))
                  if ((sist .ge. meteo(imod)%timesec(it_nm)) .and. &
                      (sist .lt. meteo(imod)%timesec(it_nm + 1))) then
                     found = .true.
                     meteo(imod)%itime1 = it_nm
                     meteo(imod)%itime2 = it_nm + 1
                  end if
               end do
            else!sidt.lt.0
               do while (.not. found)
                  it_nm = it_nm + 1
                  if ((it_nm + 1) .gt. meteo(imod)%nt) call runend('Time not found in'//TRIM(str))
                  if ((sist .gt. meteo(imod)%timesec(it_nm)) .and. &
                      (sist .le. meteo(imod)%timesec(it_nm + 1))) then
                     found = .true.
                     meteo(imod)%itime1 = it_nm + 1
                     meteo(imod)%itime2 = it_nm
                  end if
               end do
            end if !sidt.
         else
            meteo(imod)%itime1 = 1
            meteo(imod)%itime2 = 1
         end if
         !
         !*** Read Meteo Variables for the first time
         !
         call readmeteo(imod, meteo(imod)%itime1)
         call readmeteo(imod, meteo(imod)%itime2)
         !
         !*** Check if it is necesary write information in output file.
         if (meteo(imod)%postprocess) outputmet = .true.
         !
         !*** Calculate weight factor in time (for later interpolation)
         !
         meteo(imod)%time_s = (time - meteo(imod)%timesec(meteo(imod)%itime1))/ &    ! s in [0,1]
                              (meteo(imod)%timesec(meteo(imod)%itime2) - meteo(imod)%timesec(meteo(imod)%itime1))
         !
         !*** Calculate the tropopause and another calculations neccesary in
         !    each meteorological data actualization
         !
         call calcdif
      end do
      !
      ipass = 1
      !
   else
      !
      !*** Check if meteo needs to be updated
      !
      do imod = 1, numod
         !
         if (sidt .gt. 0) then
            !
            !  Forward in time
            !
            if (time .gt. meteo(imod)%timesec(meteo(imod)%itime2)) then
               meteo(imod)%itime1 = meteo(imod)%itime2
               meteo(imod)%itime2 = meteo(imod)%itime2 + 1
               call readmeteo(imod, meteo(imod)%itime2)
               if (meteo(imod)%postprocess) outputmet = .true.
               !
               !*** Calculate the tropopause and another calculations neccesary in
               !    each meteorological data actualization
               !
               call calcdif
            end if
            !
         else if (sidt .lt. 0) then
            !
            ! Backwards in time
            !
            if (time .lt. meteo(imod)%timesec(meteo(imod)%itime2)) then
               meteo(imod)%itime1 = meteo(imod)%itime2
               meteo(imod)%itime2 = meteo(imod)%itime2 - 1
               call readmeteo(imod, meteo(imod)%itime2)
               if (meteo(imod)%postprocess) outputmet = .true.
               !
               !*** Calculate the tropopause and another calculations neccesary in
               !    each meteorological data actualization
               !
               call calcdif
            end if
         end if
         !
         !*** Calculate weight factor in time (for later interpolation)
         !
         meteo(imod)%time_s = (time - meteo(imod)%timesec(meteo(imod)%itime1))/ &    ! s entre (0,1)
                              (meteo(imod)%timesec(meteo(imod)%itime2) - meteo(imod)%timesec(meteo(imod)%itime1))
      end do
      !
   end if
   !
   !*** Write variables(u,v,w,T,ro,qv) in outmeteo
   !
   if (outputmet) then
      !
      ioutput = ioutput + 1   ! update counter
      !
      if (nf90_open(TRIM(fnc_met), NF90_WRITE, ncid) /= NF90_NOERR) &
         call runend('Error opening netcdf file in meteos')
      !
      if (nf90_put_var(ncID, u_met_ID, grid%u(:, :, :, 1), start=(/1, 1, 1, ioutput/), &
                       count=(/grid%nx, grid%ny, grid%nz, 1/)) /= NF90_NOERR) &
         call wriwar('meteos : error in nf90_put_var for variable u')
      !
      if (nf90_put_var(ncID, v_met_ID, grid%v(:, :, :, 1), start=(/1, 1, 1, ioutput/), &
                       count=(/grid%nx, grid%ny, grid%nz, 1/)) /= NF90_NOERR) &
         call wriwar('meteos : error in nf90_put_var for variable v')
      !
      if (nf90_put_var(ncID, w_met_ID, grid%w(:, :, :, 1), start=(/1, 1, 1, ioutput/), &
                       count=(/grid%nx, grid%ny, grid%nz, 1/)) /= NF90_NOERR) &
         call wriwar('meteos : error in nf90_put_var for variable w')
      !
      if (nf90_put_var(ncID, T_met_ID, grid%T(:, :, :, 1), start=(/1, 1, 1, ioutput/), &
                       count=(/grid%nx, grid%ny, grid%nz, 1/)) /= NF90_NOERR) &
         call wriwar('meteos : error in nf90_put_var for variable T')
      !
      if (nf90_put_var(ncID, ro_met_ID, grid%ro(:, :, :, 1), start=(/1, 1, 1, ioutput/), &
                       count=(/grid%nx, grid%ny, grid%nz, 1/)) /= NF90_NOERR) &
         call wriwar('meteos : error in nf90_put_var for variable ro')
      !
      if (nf90_put_var(ncID, qv_met_ID, grid%qv(:, :, :, 1), start=(/1, 1, 1, ioutput/), &
                       count=(/grid%nx, grid%ny, grid%nz, 1/)) /= NF90_NOERR) &
         call wriwar('meteos : error in nf90_put_var for variable qv')
      !
      if (nf90_close(ncid) /= NF90_NOERR) call runend('Error closing netcdf file in meteos')
      !
      !*** Print Information in log file
      !
      !***  Opens and writes the log file
      !
      if (my_id .eq. 0) write (lulog, 2) 'Time dependent meteo variables have been written in output file'
      if (out_screen) write (*, 2) 'Time dependent meteo variables have been written in output file'
      !
2     format(/, a/)
      !
   end if
   !
   return
end subroutine meteos

