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
  logical     :: found=.false.
  integer(ip) :: iout,ielem3d,ielem2d,iphase,ipart,iclass
  integer(ip) :: ix,iz,iy, izl, i
  integer(ip) :: endpart ! last activated particle
  character(len=s_mess) :: attr_desc, attr_units
  !
  integer (ip),save :: ipass=0
  character(len=1)  :: ext 
  character(len=10) :: ext2
  character(len=30) :: ext5
  !
  real (rp), allocatable:: lon(:),lat(:),z(:)
  real (rp)             :: x_rad(2),y_rad(2)
  integer  (ip)         :: istat
  real(rp)              :: massair     ! sum of the total mass in air
  real(rp)              :: massground  ! sum of the total mass on ground
  real(rp)              :: massinside  ! sum of the total mass inside domain
  real(rp)              :: massoutside ! sum of the total mass outside domain
  real(rp)              :: masstotal   ! sum of the total mass rate
  integer(ip)           :: timemin     ! actual time in minute 
  !
  integer  (ip)         :: ncid
  character(len=s_mess) :: str
  character(len=3     ) :: ext3 
  character(len=2     ) :: ext4
  !
  integer(ip)  :: info
  real   (rp),save  :: rhomean=0 !to calculate thickness on deposits
  real   (rp)  :: load_pts, thick_pts
  ! To calculate total number of bins
  integer(ip)  :: nbins
  !
  if(ipass==0)then
     ipass=1
     !
     !*** Previous computations
     !
     output%nx = INT((output%lonmax-output%lonmin)/output%dx) !variables defined in master module.
     output%ny = INT((output%latmax-output%latmin)/output%dy) 
     !
     output%nelem2d=output%nx*output%ny
     output%nelem3d=output%nx*output%ny*(output%nz) !output%nz is defined in readinp.inp
     !
     !*** Allocates memory
     !
     allocate(lon   (output%nx))
     allocate(lat   (output%ny))
     allocate(conc3d(output%nelem3d))
     allocate(conc2d(output%nelem2d))
     allocate(sedim (output%nelem2d))
     allocate(output%area(output%nelem2d))
     allocate(output%volume(output%nelem3d))  
     ! 
     if(output%classes)then
        !names bin: To generate output per classes
        nbins=0
        do iphase=1,nphases
	   nbins=max(nbins,phase(iphase)%bins)
        end do
           allocate(sedim_bin (nphases,nbins,output%nelem2d))
           allocate(conc3d_bin(nphases,nbins,output%nelem3d))
           allocate(conc2d_bin(nphases,nbins,output%nelem2d))
        allocate(sedim_id_bin(nphases,nbins))
        allocate(conc2d_id_bin(nphases,nbins))
        allocate(conc3d_id_bin(nphases,nbins))
     end if
     if(output%phases)then
        !names bin: To generate output per classes
           allocate(sedim_phase (nphases,output%nelem2d))
           allocate(conc3d_phase(nphases,output%nelem3d))
           allocate(conc2d_phase(nphases,output%nelem2d))
        !ID
        allocate(sedim_id_phase(nphases))
        allocate(conc2d_id_phase(nphases))
        allocate(conc3d_id_phase(nphases))
     end if
     !
     !*** Creates the mesh stored as 1D vectors, with lon/lat values of the square centre.
     !
     do ix = 1,output%nx
        lon(ix) = output%lonmin + (ix-1./2.)*output%dx
     end do
     !
     do iy = 1,output%ny
        lat(iy) = output%latmin + (iy-1./2.)*output%dy
     end do
     !
     !*** Calculates the surface area of each 2D grid cell.and the volume of each element(3D)
     !*** First convert degrees in radians
     !
     do ix = 1,output%nx
        do iy = 1,output%ny
           ielem2d = (iy-1)*output%nx+ix
           x_rad(1)=lon(ix)*pi/180
           x_rad(2)=lon(ix+1)*pi/180
           y_rad(1)=lat(iy)*pi/180
           y_rad(2)=lat(iy+1)*pi/180
           output%area(ielem2d)=(Rearth)*(Rearth)*ABS(x_rad(2)-x_rad(1))*ABS(sin(y_rad(2))-sin(y_rad(1))) !m2
!           output%area(ielem2d)=(Rearth/1000)*(Rearth/1000)*ABS(x_rad(2)-x_rad(1))*ABS(sin(y_rad(2))-sin(y_rad(1))) !km2
           do iz=1,output%nz
              ielem3d= (iz-1)*(output%nx)*(output%ny)+(iy-1)*(output%nx)+ix
              output%volume(ielem3d)=output%area(ielem2d)*ABS(output%zlayer(iz+1)-output%zlayer(iz)) ! in m3
!              output%volume(ielem3d)=output%area(ielem2d)*ABS(output%zlayer(iz+1)-output%zlayer(iz))/1000 ! in km3
           end do !iz
        end do !iy
     end do !ix
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
              sedim_name_bin ='Ground_load_'//TRIM(str) !Determine names per class
              conc2d_name_bin='Column_mass_'//TRIM(str)
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
          sedim_name_phase ='Ground_load_p'//TRIM(str) !Determine names per phase
          conc2d_name_phase='Column_mass_p'//TRIM(str)
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
     attr_desc  = 'Column Mass'
     attr_units = 'kg/m2'
     if( nf90_put_att(ncID,conc2d_ID , 'units', attr_units) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for column mass Variable')
     if( nf90_put_att(ncID, conc2d_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for column mass Variable')
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
          call runend('outpart_nc : error in nf90_put_att for ground load Variable')
     if( nf90_put_att(ncID, sedim_ID, 'description', attr_desc) /= NF90_NOERR ) &
          call runend('outpart_nc : error in nf90_put_att for ground load Variable')
     !
     !*** Variables per class. Only if this option is activated.
     !
     if(output%classes)then
        do iphase=1,nphases
           do iclass=1, phase(iphase)%bins
              write(ext3,'(i3)') iclass
              write(ext4,'(i2)')iphase
              str = ext3//'_'//ext4
              attr_desc  = 'Column Mass'//TRIM(str)
              attr_units = 'kg/m2'
              if( nf90_put_att(ncID,conc2d_ID_bin(iphase,iclass),'units', attr_units) /= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for column_masss_bin Variable')
              if( nf90_put_att(ncID, conc2d_ID_bin(iphase,iclass), 'description', attr_desc) /= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for column_mass_bin Variable')
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
              call runend('outpart_nc : error in nf90_put_att for ground_load_bin Variable')
              if( nf90_put_att(ncID,sedim_ID_bin(iphase,iclass),'description',attr_desc)/= NF90_NOERR ) &
              call runend('outpart_nc : error in nf90_put_att for ground_load_bin Variable')
           end do
        end do
     end if
     !
     !*** Variables per phases. Only if this option is activated.
     !
     if(output%phases)then
        do iphase=1,nphases
           write(ext4,'(i2)')iphase
           attr_desc  = 'Column Mass'//TRIM(ext4)
           attr_units = 'kg/m2'
           if( nf90_put_att(ncID,conc2d_ID_phase(iphase),'units', attr_units) /= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for column_masss_phase Variable')
           if( nf90_put_att(ncID, conc2d_ID_phase(iphase), 'description', attr_desc) /= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for column_mass_phase Variable')
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
           call runend('outpart_nc : error in nf90_put_att for Ground_load_phase Variable')
           if( nf90_put_att(ncID,sedim_ID_phase(iphase),'description',attr_desc)/= NF90_NOERR ) &
           call runend('outpart_nc : error in nf90_put_att for Ground_load_phase Variable')
        end do
     end if
     !
     !*** Leave define mode 
     !
     if( nf90_enddef(ncid) /= NF90_NOERR ) call runend('outpart_nc: error in nf90_enddef')
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
     !
     !*** Tracking Points. Only if this option is activated.
     !
     if(output%track_points)then
        !
        !***  Loop over tracked points to determine the position inside the grid and their 2d-grid element
        !
        do iphase=1,nphases
           rhomean=rhomean+phase(iphase)%rhomean/nphases
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
  end if
  !
  !************************************************************************************************************
  !*** Initialize
  !
  massair     = 0.0_rp
  massground  = 0.0_rp
  massinside  = 0.0_rp
  massoutside = 0.0_rp
  masstotal   = 0.0_rp
  sedim (:)=0.0_rp
  conc2d(:)=0.0_rp
  conc3d(:)=0.0_rp
  if(output%classes)then
     sedim_bin(:,:,:) =0.0_rp
     conc2d_bin(:,:,:)=0.0_rp
     conc3d_bin(:,:,:)=0.0_rp
  end if
  if(output%phases)then
     sedim_phase(:,:) =0.0_rp
     conc2d_phase(:,:)=0.0_rp
     conc3d_phase(:,:)=0.0_rp
  end if
  !
  !*** OPEN netcdf output file
  ! 
  if(nf90_open(TRIM(fnc_part), nf90_write, ncid)/=NF90_NOERR) &
       call runend('Error opening netcdf file in outpart')
  !
  do iphase=1,nphases
     !
     !*** Determine the last activated part for each phase
     !
     endpart=phase(iphase)%actnpart+phase(iphase)%firstpart-1 
     !
     !*** Loop over particles
     !
     do ipart=phase(iphase)%firstpart,endpart   
        !
        !*** Compute concentration
        !
        if(part(ipart)%state.eq.1)then 
        if(part(ipart)%lon.ge.output%lonmin .and. part(ipart)%lon.le.output%lonmax .and. &
           part(ipart)%lat.ge.output%latmin .and. part(ipart)%lat .le.output%latmax) then !agregado reciente
           !                                  
           !*** Get the cell where the particle is
           !
           ix = INT((part(ipart)%lon-output%lonmin)/output%dx) + 1 
           if(ix.eq.(output%nx+1))ix=output%nx
           iy = INT((part(ipart)%lat-output%latmin)/output%dy) + 1  
           if(iy.eq.(output%ny+1))iy=output%ny
           !
           !
           if((part(ipart)%z.gt. 0) .and. &
              (part(ipart)%z.le.output%zlayer(1)))then
              iz=1
           else if((part(ipart)%z.ge. output%zlayer(output%nz)).and. &
              (part(ipart)%z.le. output%ztop))then
              iz=output%nz
           else 
              found=.false.
              izl=2
              do while(found.eqv..false.)
                 if ((part(ipart)%z.le.output%zlayer(izl)).and. &                
                    (part(ipart)%z.gt.output%zlayer(izl-1))) then
                    iz = izl
                    found=.true.
                 else
                    izl=izl+1
                 end if
              end do
           end if
           !
           !*** Define the grid element where the particle is
           !
           ielem3d= (iz-1)*(output%nx)*(output%ny)+(iy-1)*(output%nx)+ix
           ielem2d= (iy-1)*(output%nx)+ix
           !
           !*** Sum the mass in the corresponding element (by now stored in conc)
           !
           conc2d(ielem2d) = conc2d(ielem2d) + (part(ipart)%mass/output%area(ielem2d))
           conc3d(ielem3d) = conc3d(ielem3d) + (part(ipart)%mass/output%volume(ielem3d))
           !*** Sum the mass in massair variable to know the total
           massair=massair+part(ipart)%mass
           !*** Sum the mass if the output per class is activated
           if(output%classes)then
              conc2d_bin(iphase,part(ipart)%ibin,ielem2d)= &
              conc2d_bin(iphase,part(ipart)%ibin,ielem2d)+ part(ipart)%mass/output%area(ielem2d)
              conc3d_bin(iphase,part(ipart)%ibin,ielem3d)= &
              conc3d_bin(iphase,part(ipart)%ibin,ielem3d)+ part(ipart)%mass/output%volume(ielem3d)
           end if
           if(output%phases)then
              conc2d_phase(iphase,ielem2d)= &
              conc2d_phase(iphase,ielem2d)+ part(ipart)%mass/output%area(ielem2d)
              conc3d_phase(iphase,ielem3d)= &
              conc3d_phase(iphase,ielem3d)+ part(ipart)%mass/output%volume(ielem3d)
           end if
           !
           end if! particles inside output domain
        else if((part(ipart)%state.eq.(-1)))then
        if(part(ipart)%lon.ge.output%lonmin .and. part(ipart)%lon.le.output%lonmax .and. &
           part(ipart)%lat.ge.output%latmin .and. part(ipart)%lat .le.output%latmax) then  !agregado reciente
           !
           !*** Particle ground load
           !
           !*** Define the 2d box where the particle is
           !
           ix = INT((part(ipart)%lon-output%lonmin)/output%dx) + 1 
           if(ix.eq.(output%nx+1))ix=output%nx
           iy = INT((part(ipart)%lat-output%latmin)/output%dy) + 1  
           if(iy.eq.(output%ny+1))iy=output%ny
           !
           !*** Define the grid element where the particle is
           !
           ielem2d= (iy-1)*(output%nx)+ix
           !
           !*** Sum the mass in the corresponding element
           !
           sedim(ielem2d)= sedim(ielem2d) + (part(ipart)%mass/output%area(ielem2d))
           !*** Sum the mass in massground variable to know the total
           massground=massground+part(ipart)%mass
           !*** Sum the mass if the output per class is activated
           if(output%classes)then
              sedim_bin(iphase,part(ipart)%ibin,ielem2d)= &
              sedim_bin(iphase,part(ipart)%ibin,ielem2d)+ (part(ipart)%mass/output%area(ielem2d))
           end if
           if(output%phases)then
              sedim_phase(iphase,ielem2d)= &
              sedim_phase(iphase,ielem2d)+ (part(ipart)%mass/output%area(ielem2d))
           end if
           ! 
        end if !particle inside output domain
        else if ((part(ipart)%state.eq.(-2)))then ! Must be added particles inside simulation domain and outside output domain
           massoutside=massoutside+part(ipart)%mass
        end if ! ipart state
     end do ! particle loop
  end do    ! phase loop
  !*** Calculate the total mass inside domain and total mass
  massinside=massair+massground
  masstotal=massinside+massoutside
  !
  !*** WRITE MASS FOR THE NEW OUTTIME
  !
  write(ext,'(i1)') iout
  if(nf90_put_var(ncID, conc3d_ID,conc3d,start=(/1,1,1,iout/),count=(/(output%nx),(output%ny),(output%nz),1/))/= NF90_NOERR) &
       call wriwar('outpart : error in nf90_put_var for variable conc3d for output time'//ext)
  !
  if(nf90_put_var(ncID, conc2d_ID,conc2d,start=(/1,1,iout/),count=(/(output%nx),(output%ny),1/))/= NF90_NOERR) &
       call wriwar('outpart : error in nf90_put_var for variable conc2d for output time'//ext)
  !
  if(nf90_put_var(ncID, sedim_ID,sedim,start=(/1,1,iout/),count=(/(output%nx),(output%ny),1/))/= NF90_NOERR) &
       call wriwar('outpart : error in nf90_put_var for variable sedim for output time'//ext)
  !
  !*** Write the mass in air and the mass on ground in the .log file
  !
  write(ext5,*) massground
  write(lulog,5)            'The total mass deposited on ground is'//ext5//'kg'
  write(ext5,*) massair
  write(lulog,5)            'The total mass in air inside domain is'//ext5//'kg'
  write(ext5,*) massinside
  write(lulog,5)            'The total mass inside domain is'//ext5//'kg'
  write(ext5,*) massoutside
  write(lulog,5)            'The mass outside domain is'//ext5//'kg'
  write(ext5,*) masstotal
  write(lulog,5)            'The total mass liberated is'//ext5//'kg'
5 format(/,  a /)
  !
  !*** Write the mass per class if the output per class is activated
  !
  if(output%classes)then
     do iphase=1,nphases
        do iclass=1,phase(iphase)%bins
           if(nf90_put_var(ncID, conc3d_ID_bin(iphase,iclass),conc3d_bin(iphase,iclass,:),start=(/1,1,1,iout/),&
           count=(/(output%nx),(output%ny),(output%nz),1/))/= NF90_NOERR) &
           call wriwar('outpart : error in nf90_put_var for variable conc3d_bin for output time'//ext)
           !
           if(nf90_put_var(ncID, conc2d_ID_bin(iphase,iclass),conc2d_bin(iphase,iclass,:),start=(/1,1,iout/),&
           count=(/(output%nx),(output%ny),1/))/= NF90_NOERR) &
           call wriwar('outpart : error in nf90_put_var for variable conc2d_bin for output time'//ext)
           !
           if(nf90_put_var(ncID, sedim_ID_bin(iphase,iclass),sedim_bin(iphase,iclass,:),start=(/1,1,iout/),&
           count=(/(output%nx),(output%ny),1/))/= NF90_NOERR) &
           call wriwar('outpart : error in nf90_put_var for variable sedim_bin for output time'//ext)
        end do
     end do
  end if
  !
  !*** Write the mass per phase if the output per phase is activated
  !
  if(output%phases)then
     do iphase=1,nphases
           if(nf90_put_var(ncID, conc3d_ID_phase(iphase),conc3d_phase(iphase,:),start=(/1,1,1,iout/),&
           count=(/(output%nx),(output%ny),(output%nz),1/))/= NF90_NOERR) &
           call wriwar('outpart : error in nf90_put_var for variable conc3d_phase for output time'//ext)
           !
           if(nf90_put_var(ncID, conc2d_ID_phase(iphase),conc2d_phase(iphase,:),start=(/1,1,iout/),&
           count=(/(output%nx),(output%ny),1/))/= NF90_NOERR) &
           call wriwar('outpart : error in nf90_put_var for variable conc2d_phase for output time'//ext)
           !
           if(nf90_put_var(ncID, sedim_ID_phase(iphase),sedim_phase(iphase,:),start=(/1,1,iout/),&
           count=(/(output%nx),(output%ny),1/))/= NF90_NOERR) &
           call wriwar('outpart : error in nf90_put_var for variable sedim_phase for output time'//ext)
     end do
  end if
  !
  !*** CLOSE
  !
  if( nf90_close(ncid)/=NF90_NOERR) call runend('Error closing netcdf file in outpart')
  !
  !*** Write the mass per track point if there are considered. In another files.
  !
   timemin=int(output%timesec(iout)/60)
  if(output%track_points)then
    do i=1,npts
       if(use_pts(i)) then
          open(79,file=TRIM(name_file_pts(i))//'.res',iostat=info,status='old',position='append')
          if( info /= 0 ) goto 100
          !calc thickness in cm (Divides value of sedimentation load by average particle density, and multiply by 100)
          load_pts=sedim(ielem2dpts(i))  !kg/m2
          thick_pts=load_pts/rhomean     ! 
          thick_pts=thick_pts*1d2 !in cm
          write(79,20),  timemin, load_pts,thick_pts
20        format(i7,1x,f13.6, 1x, f13.6)
          close(79)
       end if !use_pts
    end do !npts
  end if  !track_points
  !
  !*** Print Information in log file
  !
  !***  Writes the log file
  !
  write(ext2,'(f10.0)') output%timesec(iout)
  ! 
  write(lulog,2)            'The output file have been written for time'//ext2//'seconds'
  if(out_screen) write(*,2) 'The output file have been written for time'//ext2//'seconds'
  !
2 format(/,  a /)
  !
  return

  !List of errors
100 call runend('Error opening file for tracked points')
  end subroutine outpart_nc
