subroutine write_nc_var3d
  !*******************************************************************
  !*
  !*     Writes a 3d variable
  !*
  !*******************************************************************
  use KindType
  use Master
  use InpOut
  use netcdf
  implicit none
  !
  integer(ip) :: ilen,ilev,ix,iy,info,i
  !
  !***   Loop over pressure levels
  !
  do ilev = 1,np
     !
     !***      File name
     !
     lugrbname = TRIM(lugrbbase)//'/'//TRIM(var3d(i3d))//'.'
     ilen = LEN_TRIM(lugrbname)
     if(it.lt.10) then
        write(lugrbname(ilen+1:ilen+1),'(i1.1)') it
     else if(it.lt.100) then
        write(lugrbname(ilen+1:ilen+2),'(i2.2)') it
     else if(it.lt.1000) then
        write(lugrbname(ilen+1:ilen+3),'(i3.3)') it
     else if(it.lt.10000) then
        write(lugrbname(ilen+1:ilen+3),'(i4.4)') it
     end if
     lugrbname = TRIM(lugrbname)//'.'
     ilen = LEN_TRIM(lugrbname)
     if(pres(ilev).lt.10) then
        write(lugrbname(ilen+1:ilen+1),'(i1.1)') INT(pres(ilev))
     else if(pres(ilev).lt.100) then
        write(lugrbname(ilen+1:ilen+2),'(i2.2)') INT(pres(ilev))
     else if(pres(ilev).lt.1000) then
        write(lugrbname(ilen+1:ilen+3),'(i3.3)') INT(pres(ilev))
     else if(pres(ilev).lt.10000) then
        write(lugrbname(ilen+1:ilen+4),'(i4.4)') INT(pres(ilev))
     end if
     !
     !***   Reads the file (ASCI or binary)
     !
     if(.not.wgrib_bin) then
        open(99,FILE=TRIM(lugrbname),status='unknown')
        read(99,*,iostat=info) ix,iy
        if( info /= 0 ) then
           close(99)
           work(:,:,ilev) = 0_rp
        else
           do iy = 1,ny
              do ix = 1,nx
                 read(99,*) work(ix,iy,ilev)
              end do
           end do
           close(99)
        end if
        !
     else        ! binary
        !
        open(99,FILE=TRIM(lugrbname),form='unformatted',status='unknown')
        read(99,iostat=info) work4
        if( info /= 0 ) then
           close(99)
           work(:,:,ilev) = 0_rp
        else
           close(99)
           work(1:nx,1:ny,ilev)=DBLE(work4)
        end if
     end if
     !
     write(lulog,10) TRIM(lugrbname)
10   format(1x,'Converting file  : ',a)
     !
  end do
  !
  !***  Translations.
  !
  if(invert_y) then
     do ilev = 1,np
        work2d(1:nx,1:ny) = work(1:nx,1:ny,ilev)
        do iy = 1,ny
           work(1:nx,iy,ilev) = work2d(1:nx,ny-iy+1)
        end do
     end do
  end if
  !
  if(invert_x) then
     do ilev = 1,np
        work2d(1:nx,1:ny) = work(1:nx,1:ny,ilev)      
        !
        i = 0
        do ix = ipos,nx
           i = i + 1
           work(i,1:ny,ilev) = work2d(ix,1:ny)
        end do
        do ix = 1,ipos-1
           i = i + 1
           work(i,1:ny,ilev) = work2d(ix,1:ny)
        end do
     end do 
  end if  
  !
  !***  Writes the 3D variable
  !
  if( nf90_open(TRIM(luncname),NF90_WRITE, ncID) /= 0 ) &
       call runend('write_nc_var3d : Error in nf90_open')
  if( nf90_inq_varid(ncID,TRIM(var3d(i3d)),var3d_ID(i3d)) /= 0) &
       call runend('write_nc_var3d : Error getting var_ID')
  if( nf90_put_var(ncID, var3d_ID(i3d), work, start=(/1,1,1,itt+1/), count=(/nx,ny,np,1/) ) /= 0 ) &
       call runend('write_nc_var3d: error in nf90_put_var')
  if( nf90_close(ncID) /= 0) &
       call runend('write_nc_var3d : Error in nf90_close')
  !
  return
end subroutine write_nc_var3d
