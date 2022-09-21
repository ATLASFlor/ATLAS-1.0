  subroutine init_netcdf
  !******************************************************************************************
  !*
  !* Initialize output variables, allocate memory
  !*
  !******************************************************************************************
  !
  use InpOut
  use Master
  use KindType                                     
  use Netcdf
  implicit none
  !
  integer(ip) :: ix,iz,iy, i
  integer(ip) :: iclass,ielem3d,ielem2d,iphase
  real (rp), allocatable:: lon(:),lat(:),z(:)
  !
  real (rp)             :: x_rad(2),y_rad(2)
  ! 
  character(len=s_mess) :: attr_desc, attr_units
  ! 
  integer  (ip)         :: ncid
  character(len=s_mess) :: str
  character(len=3     ) :: ext3 
  character(len=2     ) :: ext4
  ! 
  integer(ip)  :: info
  !
     !
     !*** Allocates memory
     !
     allocate(lon   (output%nx))
     allocate(lat   (output%ny))
     !
     !*** Creates the mesh stored as 1D vectors, with lon/lat values of the square centre.
     !
     do ix = 1,output%nx
        lon(ix) = output%lonmin + (ix-1./2.)*output%dx    !RESTAR 1 en vez de un medio!!!
     end do
     !
     do iy = 1,output%ny
        lat(iy) = output%latmin + (iy-1./2.)*output%dy
     end do
     !
     !
     !
     !*** Create a netcdf file 
     !
     if( nf90_create(TRIM(fnc_part), NF90_CLOBBER, ncid)/=NF90_NOERR) &
          call runend('Error creting netcdf file in outpart_nc')
     !
     !*** Define Dimensions
     !
     if( nf90_def_dim(ncid, nx_part_name, output%nx, nx_part_id)/=NF90_NOERR)&
          call runend('Error creting netcdf nx dimension in outpart_nc')
     !   
     if( nf90_def_dim(ncid, ny_part_name, output%ny, ny_part_id)/=NF90_NOERR)&
          call runend('Error creting netcdf ny dimension in outpart_nc')
     !
     if( nf90_def_dim(ncid, nz_part_name, output%nz, nz_part_id)/=NF90_NOERR)&
          call runend('Error creting netcdf nz dimension in outpart_nc')
     !
     if( nf90_def_dim(ncid, nt_part_name, output%nt, nt_part_id)/=NF90_NOERR)&
          call runend('Error creting netcdf nt dimension in outpart_nc')
     ! 
     !*** Define time-independent variables
     !
     if( nf90_def_var(ncID, lon_part_name ,NF90_FLOAT, (/nx_part_id/), lon_part_ID) /= NF90_NOERR ) &
          call runend('outpart_nc: error in nf90_def_variable for variable '//TRIM(lon_part_name))
     !
     if( nf90_def_var(ncID, lat_part_name ,NF90_FLOAT, (/ny_part_id/), lat_part_ID) /= NF90_NOERR ) &
          call runend('outpart_nc: error in nf90_def_variable for variable '//TRIM(lat_part_name))
     !
     if( nf90_def_var(ncID, alt_part_name ,NF90_FLOAT, (/nz_part_id/), alt_part_ID) /= NF90_NOERR ) &
          call runend('outpart_nc: error in nf90_def_variable for variable '//TRIM(alt_part_name))
     !
     !*** Define Time Dependent Variables
     !
     if( nf90_def_var(ncID, conc2d_name ,NF90_FLOAT, (/nx_part_id,ny_part_id,nt_part_id/), conc2d_ID) /= NF90_NOERR ) &
          call runend('outpart: error in nf90_def_variable for variable '//TRIM(conc2d_name))
     !
     if( nf90_def_var(ncID, conc3d_name ,NF90_FLOAT, (/nx_part_id,ny_part_id,nz_part_id,nt_part_id/), &
          conc3d_ID) /= NF90_NOERR ) &
          call runend('outpart: error in nf90_def_variable for variable '//TRIM(conc3d_name))
     !
     if( nf90_def_var(ncID, sedim_name ,NF90_FLOAT, (/nx_part_id,ny_part_id,nt_part_id/), &
          sedim_ID) /= NF90_NOERR ) &
          call runend('outpart: error in nf90_def_variable for variable '//TRIM(sedim_name))
     !*** Variables per class. Only if this option is activated.
     if(output%classes)then
        do iphase=1,nphases
           do iclass=1, phase(iphase)%bins
              write(ext3,'(i3)')iclass
              write(ext4,'(i2)')iphase
              str=TRIM(ADJUSTL(ext4))//'_'//TRIM(ADJUSTL(ext3))
              sedim_name_bin ='sedimentation_'//TRIM(str) !Determine names per class
              conc2d_name_bin='concentration_2D_'//TRIM(str)
              conc3d_name_bin='concentration_3D_'//TRIM(str)
              if( nf90_def_var(ncID, conc2d_name_bin ,NF90_FLOAT, &
              (/nx_part_id,ny_part_id,nt_part_id/), conc2d_ID_bin(iphase,iclass)) /= NF90_NOERR ) &
              call runend('outpart: error in nf90_def_variable for variable conc2d_bin_'//TRIM(str))
              !
              if( nf90_def_var(ncID, conc3d_name_bin ,NF90_FLOAT, &
              (/nx_part_id,ny_part_id,nz_part_id,nt_part_id/),conc3d_ID_bin(iphase,iclass)) /= NF90_NOERR ) &
              call runend('outpart: error in nf90_def_variable for variable conc3d_bin_'//TRIM(str))
              !
              if( nf90_def_var(ncID, sedim_name_bin ,NF90_FLOAT, &
              (/nx_part_id,ny_part_id,nt_part_id/),sedim_ID_bin(iphase,iclass)) /= NF90_NOERR ) &
              call runend('outpart: error in nf90_def_variable for variable Sedim_bin_'//TRIM(str))
           end do
        end do
     end if
     !*** Variables per phases. Only if this option is activated.
     if(output%phases)then
       do iphase=1,nphases
          write(ext4,'(i2)')iphase
          str=TRIM(ADJUSTL(ext4))
          sedim_name_phase ='sedimentation_p'//TRIM(str) !Determine names per phase
          conc2d_name_phase='concentration_2D_p'//TRIM(str)
          conc3d_name_phase='concentration_3D_p'//TRIM(str)
        if( nf90_def_var(ncID, conc2d_name_phase ,NF90_FLOAT, &
            (/nx_part_id,ny_part_id,nt_part_id/), conc2d_ID_phase(iphase)) /= NF90_NOERR ) &
            call runend('outpart: error in nf90_def_variable for variable '//TRIM(conc2d_name_phase))
        !
        if( nf90_def_var(ncID, conc3d_name_phase ,NF90_FLOAT, &
          (/nx_part_id,ny_part_id,nz_part_id,nt_part_id/), &
          conc3d_ID_phase(iphase)) /= NF90_NOERR ) &
          call runend('outpart: error in nf90_def_variable for variable '//TRIM(conc3d_name_phase))
        !
        if( nf90_def_var(ncID, sedim_name_phase ,NF90_FLOAT, &
          (/nx_part_id,ny_part_id,nt_part_id/), &
          sedim_ID_phase(iphase)) /= NF90_NOERR ) &
          call runend('outpart: error in nf90_def_variable for variable '//TRIM(sedim_name_phase))
       end do
     end if
     ! 
     !*** Assign attribute values
     ! 
     attr_desc  = 'Longitude'
     attr_units = 'deg'
     if( nf90_put_att(ncID,lon_part_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for lon Variable')
     if( nf90_put_att(ncID, lon_part_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for lon Variable')
     ! 
     attr_desc  = 'Latitude'
     attr_units = 'deg'
     if( nf90_put_att(ncID,lat_part_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for lat Variable')
     if( nf90_put_att(ncID, lat_part_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for lat Variable')
     ! 
     attr_desc  = 'Layer height'
     attr_units = 'm'
     if( nf90_put_att(ncID,alt_part_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for lat Variable')
     if( nf90_put_att(ncID, alt_part_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for lat Variable')
     ! 
     !
     attr_desc  = '2D Concentration'
     attr_units = 'kg/m2'
     if( nf90_put_att(ncID,conc2d_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for concentration2d Variable')
     if( nf90_put_att(ncID, conc2d_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for concentration2d Variable')
     !
     attr_desc  = '3D Concentration'
     attr_units = 'kg/m3'
     if( nf90_put_att(ncID,conc3d_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for concentration3d Variable')
     if( nf90_put_att(ncID, conc3d_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for concentration3d Variable')
     !
     attr_desc  = 'Ground load'
     attr_units = 'kg/m2'
     if( nf90_put_att(ncID,sedim_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for Sedimentation Variable')
     if( nf90_put_att(ncID, sedim_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for Sedimentation Variable')
     !
     !*** Variables per class. Only if this option is activated.
     !
     if(output%classes)then
        do iphase=1,nphases
           do iclass=1, phase(iphase)%bins
              write(ext3,'(i3)') iclass
              write(ext4,'(i2)')iphase
              str = ext3//'_'//ext4
              attr_desc  = '2D Concentration_'//TRIM(str)
              attr_units = 'kg/m2'
              if( nf90_put_att(ncID,conc2d_ID_bin(iphase,iclass),'units', attr_units) /= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for concentration2d_bin Variable')
              if( nf90_put_att(ncID, conc2d_ID_bin(iphase,iclass), 'description', attr_desc) /= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for concentration2d_bin Variable')
              !
              attr_desc  = '3D Concentration_'//TRIM(str)
              attr_units = 'kg/m3'
              if( nf90_put_att(ncID,conc3d_ID_bin(iphase,iclass),'units',attr_units) /= NF90_NOERR) &
              call runend('outpart_nc : error in nf90_put_att for concentration3d_bin Variable')
              if( nf90_put_att(ncID, conc3d_ID_bin(iphase,iclass), 'description', attr_desc) /= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for concentration3d_bin Variable')
              !
              attr_desc  = 'Ground load_'//TRIM(str)
              attr_units = 'kg/m2'
              if( nf90_put_att(ncID,sedim_ID_bin(iphase,iclass),'units',attr_units) /= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for Sedimentation Variable')
              if( nf90_put_att(ncID,sedim_ID_bin(iphase,iclass),'description',attr_desc)/= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for Sedimentation Variable')
           end do
        end do
     end if
     !
     !*** Variables per phases. Only if this option is activated.
     !
     if(output%phases)then
        do iphase=1,nphases
           write(ext4,'(i2)')iphase
           attr_desc  = '2D Concentration_'//TRIM(ext4)
           attr_units = 'kg/m2'
           if( nf90_put_att(ncID,conc2d_ID_phase(iphase),'units', attr_units) /= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for concentration2d_phase Variable')
           if( nf90_put_att(ncID, conc2d_ID_phase(iphase), 'description', attr_desc) /= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for concentration2d_phase Variable')
           !
           attr_desc  = '3D Concentration_'//TRIM(ext4)
           attr_units = 'kg/m3'
           if( nf90_put_att(ncID,conc3d_ID_phase(iphase),'units',attr_units) /= NF90_NOERR) &
           call runend('outpart_nc : error in nf90_put_att for concentration3d_phase Variable')
           if( nf90_put_att(ncID, conc3d_ID_phase(iphase), 'description', attr_desc) /= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for concentration3d_bin Variable')
           !
           attr_desc  = 'Ground load_'//TRIM(ext4)
           attr_units = 'kg/m2'
           if( nf90_put_att(ncID,sedim_ID_phase(iphase),'units',attr_units) /= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for Sedimentation_phase Variable')
           if( nf90_put_att(ncID,sedim_ID_phase(iphase),'description',attr_desc)/= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for Sedimentation_phase Variable')
        end do
     end if
     !
     !*** Leave define mode 
     !
     if( nf90_enddef(ncID) /= NF90_NOERR ) call runend('outpart_nc: error in nf90_enddef')
     !
     !*** Fill values for time-independent variables
     !
     if( nf90_put_var(ncID, lon_part_ID,lon , start=(/1/), count=(/output%nx/) ) /= NF90_NOERR ) &
          call wriwar('outmeteo : error in nf90_put_var for variable lon')
     !
     if( nf90_put_var(ncID, lat_part_ID,lat , start=(/1/), count=(/output%ny/) ) /= NF90_NOERR ) &
          call wriwar('outmeteo : error in nf90_put_var for variable lat')
     !
     if( nf90_put_var(ncID, alt_part_ID,output%zlayer , start=(/1/), count=(/output%nz/) ) /= NF90_NOERR ) &
          call wriwar('outmeteo : error in nf90_put_var for variable alt')
     !
     !*** Close the file
     !
     if( nf90_close(ncID) /= NF90_NOERR ) &
          call runend('outmeteo : error closing netcdf file')
     !
     !
     !*** Tracking Points. Only if this option is activated.
     !
     if(output%track_points)then
        !
        !***  Loop over tracked points to determine the position inside the grid and their 2d-grid element
        !
        output%rhomean=0
        do iphase=1,nphases
           output%rhomean=output%rhomean+phase(iphase)%rhomean/nphases
        end do
        do i = 1,npts
           if(use_pts(i)) then
              !*** Get the cell where the particle is
              !
              ix = INT((xpts(i)-output%lonmin)/output%dx) + 1 
              if(ix.eq.(output%nx+1))ix=output%nx
              iy = INT((ypts(i)-output%latmin)/output%dy) + 1  
              if(iy.eq.(output%ny+1))iy=output%ny
              ielem2dpts(i)= (iy-1)*(output%nx)+ix
              ! 
              ! Creates the files
              if(.not.restart) then
                 open(79,file=TRIM(name_file_pts(i))//'.res',iostat=info,status='unknown')
                 if( info /= 0 ) goto 100
                 write(79,17) TRIM(name_pts(i)),xpts(i),ypts(i)
17               format(&
                   'Tracking point file for : ',a              ,/, &
                   'Coordinates             : ',f13.4,1x,f13.4 ,/, &
                   '  Time    DDMMMM-HH:MM    load     thickness    conc.   conc.PM5   conc.PM10  conc.PM20   cumul',/, &
                   '                         ground   (compact=1)   ground   ground     ground    ground      conc. ',/, &
                   '  (min)      (--)        (kg/m2)     (cm)       (g/m3)    (g/m3)     (g/m3)    (g/m3)     (g/m2)',/, &
                   '------------------------------------------------------------------------------------------------')
              else ! restart
                 open(79,file=TRIM(name_file_pts(i))//'.res',iostat=info,status='old',position='append')
                 if( info /= 0 ) goto 100
              end if ! not restart
           end if !use pts
        end do !pts
     end if !track_points
     !
     !*** Deallocate memory
     !
     deallocate(lon)
     deallocate(lat)
     !
   return 
  !List of errors
100 call runend('Error opening file for tracked points')
end subroutine init_netcdf
