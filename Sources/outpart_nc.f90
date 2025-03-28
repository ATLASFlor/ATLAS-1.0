   subroutine outpart_nc(iout)
      !******************************************************************************************
      !*
      !* write in a netcdf output file the total mass concentration and particle trajectories
      !*
      !******************************************************************************************
      !
      use InpOut
      use Master
      use KindType
      use Netcdf
      implicit none
      !
      logical     :: found = .false.
      integer(ip) :: iout
      integer(ip) :: ielem3d, ielem2d, iphase, ipart, iclass
      integer(ip) :: ix, iz, iy, izl, i
      integer(ip) :: endpart ! last activated particle
      character(len=s_mess) :: attr_desc, attr_units
      !
      integer(ip), save :: ipass = 0
      character(len=1)  :: ext
      character(len=10) :: ext2
      character(len=30) :: ext5
      !
      real(rp), allocatable:: lon(:), lat(:), z(:)
      real(rp)             :: x_rad(2), y_rad(2)
      integer(ip)         :: istat
      integer(ip)           :: timemin     ! actual time in minute
      !
      integer(ip)         :: ncid
      character(len=s_mess) :: str
      character(len=3) :: ext3
      character(len=2) :: ext4
      !
      integer(ip)  :: info
      !real(rp), save  :: rhomean = 0 !to calculate thickness on deposits
      real(rp)     :: loadpts,thickpts
      ! To calculate total number of bins
      integer(ip)  :: nbins
      integer(ip)  :: elem1, elem2, elem3, elem4, elem5, elem6, elem7, elem8
      real(rp)     :: area1, area2, area3, area4, areaT
      real(rp)    :: a, b, c, d, e, f
      !
      !
      !
      !*** OPEN netcdf output file
      !
      if (nf90_open(TRIM(fnc_part), nf90_write, ncid) /= NF90_NOERR) &
         call runend('Error opening netcdf file in outpart')
      !

      !********************************************************************************************************
      ! *************    WRITE IN OUTPUT FILE  ****************************************************************
      !********************************************************************************************************
      !
      !*** WRITE MASS FOR THE NEW OUTTIME
      !
      write (ext, '(i1)') iout
      if (nf90_put_var(ncID, conc3d_ID, conc3d_global, start=(/1, 1, 1, iout/), &
                       count=(/(output%nx - 1), (output%ny - 1), (output%nz-1), 1/)) /= NF90_NOERR) &
         call wriwar('outpart : error in nf90_put_var for variable conc3d for output time'//ext)
      !
      if (nf90_put_var(ncID, conc2d_ID, conc2d_global, start=(/1, 1, iout/), &
                       count=(/(output%nx - 1), (output%ny - 1), 1/)) /= NF90_NOERR) &
         call wriwar('outpart : error in nf90_put_var for variable conc2d for output time'//ext)
      !
      if (nf90_put_var(ncID, sedim_ID, sedim_global, start=(/1, 1, iout/), count=(/(output%nx - 1), &
         (output%ny - 1), 1/)) /= NF90_NOERR) call wriwar('outpart : error in nf90_put_var for'&
        &' variable sedim for output time'//ext)
      !
      !*** Write the mass per class if the output per class is activated
      !
      if (output%classes) then
         do iphase = 1, nphases
            do iclass = 1, phase(iphase)%bins
               if (nf90_put_var(ncID, conc3d_ID_bin(iphase, iclass), concen_bin(iphase, iclass, :, :, :), start=(/1, 1, 1, iout/), &
                                count=(/(output%nx), (output%ny), (output%nz), 1/)) /= NF90_NOERR) &
                  call wriwar('outpart : error in nf90_put_var for variable conc3d_bin for output time'//ext)
               !
               if (nf90_put_var(ncID, conc2d_ID_bin(iphase, iclass), colmass_bin(iphase, iclass, :, :), start=(/1, 1, iout/), &
                                count=(/(output%nx), (output%ny), 1/)) /= NF90_NOERR) &
                  call wriwar('outpart : error in nf90_put_var for variable conc2d_bin for output time'//ext)
               !
               if (nf90_put_var(ncID, sedim_ID_bin(iphase, iclass), grload_bin(iphase, iclass, :, :), start=(/1, 1, iout/), &
                                count=(/(output%nx), (output%ny), 1/)) /= NF90_NOERR) &
                  call wriwar('outpart : error in nf90_put_var for variable sedim_bin for output time'//ext)
            end do
         end do
      end if
      !
      !*** Write the mass per phase if the output per phase is activated
      !
      if (output%phases) then
         do iphase = 1, nphases
            if (nf90_put_var(ncID, conc3d_ID_phase(iphase), concen_phase(iphase, :, :, :), start=(/1, 1, 1, iout/), &
                             count=(/(output%nx), (output%ny), (output%nz), 1/)) /= NF90_NOERR) &
               call wriwar('outpart : error in nf90_put_var for variable conc3d_phase for output time'//ext)
            !
            if (nf90_put_var(ncID, conc2d_ID_phase(iphase), colmass_phase(iphase, :, :), start=(/1, 1, iout/), &
                             count=(/(output%nx), (output%ny), 1/)) /= NF90_NOERR) &
               call wriwar('outpart : error in nf90_put_var for variable conc2d_phase for output time'//ext)
            !
            if (nf90_put_var(ncID, sedim_ID_phase(iphase), grload_phase(iphase, :, :), start=(/1, 1, iout/), &
                             count=(/(output%nx), (output%ny), 1/)) /= NF90_NOERR) &
               call wriwar('outpart : error in nf90_put_var for variable sedim_phase for output time'//ext)
         end do
      end if
      !
      !*** CLOSE
      !
      if (nf90_close(ncid) /= NF90_NOERR) call runend('Error closing netcdf file in outpart')
      !
      !*** Write the mass per track point if there are considered. In another files.
      !
      timemin = INT(output%timesec(iout)/60)
      if (output%track_points) then
         do i = 1, npts
            if (use_pts(i)) then
               open (79, file=TRIM(name_file_pts(i))//'.res', iostat=info, status='old', position='append')
               if (info /= 0) goto 100
               !calc thickness in cm (Divides value of sedimentation load by average particle density, and multiply by 100)
               loadpts = loadpts_global(i)!sedim(ielem2dpts(i))  !kg/m2
               !thick_pts = load_pts/rhomean     !
               thickpts = thickpts_global(i)!thick_pts*1d2 !in cm
               write (79, 20)timemin,loadpts,thickpts
20             format(i7, 1x, f13.8, 1x, f13.8)
               close (79)
            end if !use_pts
         end do !npts
      end if  !track_points
      !
      !*** Print Information in log file
      !
      !***  Writes the log file
      !
      write (ext2, '(i7)') timemin!output%timesec(iout)
      !
      !
!2     format(/, a/)


   !*** Write in log file    para modificar!!! lo de arriba no estaria funcionando.. no escribe en el log|
   !
!   write (string_time, '(F13.0)') meteo(nm)%time(itime)
!   string_time = string_time(1:4)//'/'//string_time(5:6)//'/'//string_time(7:8)//'-'//string_time(9:10) &
!                 //':'//string_time(11:12)//':00'
   write (lulog, 21) ext2, massair_global,massground_global,massinside_global,massoutside_global,masstotal_global
21  format(/, &
           'Output information for ', a, ' minutes was written. ', /,&
	   'The mass in the air amounts', f16.0, 'kg', /,&
	   'The mass in the ground amounts', f16.0, 'kg', /,&
	   'The mass inside the domain amounts', f16.0, 'kg', /,&
	   'The mass outside the domain amounts', f16.0, 'kg', /,&
	   'The total mass amounts', f18.0, 'kg'/)
!   if (out_screen) write (*, 1) TRIM(meteo(nm)%modeltype), &
!      TRIM(str), &
!      string_time!meteo(nm)%time(itime) f16.0
!1  format(/, &
!           'Read ', a, ' data file : ', a, ' at time : ', a,/)
      !
      return

      !List of errors
100   call runend('Error opening file for tracked points')
   end subroutine outpart_nc
