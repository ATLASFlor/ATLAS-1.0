subroutine write_nc_var2d
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
  integer(ip) :: ilen,ix,iy,info,i
  !
  !***   File name
  !
  lugrbname = TRIM(lugrbbase)//'/'//TRIM(var2d(i2d))//'.'
  ilen = LEN_TRIM(lugrbname)
  if(it.lt.10) then
     write(lugrbname(ilen+1:ilen+1),'(i1.1)') it
  else if(it.lt.100) then
     write(lugrbname(ilen+1:ilen+2),'(i2.2)') it
  else if(it.lt.1000) then
     write(lugrbname(ilen+1:ilen+3),'(i3.3)') it
  end if
  !
  !***   Reads the file (ASCII or binary)
  !
  if(.not.wgrib_bin) then
     open(99,FILE=TRIM(lugrbname),status='unknown')
     read(99,*,iostat=info) ix,iy
     if( info /= 0 ) then
        close(99)
        work(:,:,1) = 0_rp
     else
        do iy = 1,ny
           do ix = 1,nx
              read(99,*) work(ix,iy,1)
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
        work(:,:,1) = 0_rp
     else
        close(99)
        work(1:nx,1:ny,1)=DBLE(work4)
     end if
  end if
  !
  write(lulog,10) TRIM(lugrbname)
10 format(1x,'Converting file  : ',a)
  !
  !***  Translations.
  !
  if(invert_y) then
     work2d(1:nx,1:ny) = work(1:nx,1:ny,1)
     do iy = 1,ny
        work(1:nx,iy,1) = work2d(1:nx,ny-iy+1)
     end do
  end if
  !
  if(invert_x) then
     work2d(1:nx,1:ny) = work(1:nx,1:ny,1)      
     !
     i = 0
     do ix = ipos,nx
        i = i + 1
        work(i,1:ny,1) = work2d(ix,1:ny)
     end do
     do ix = 1,ipos-1
        i = i + 1
        work(i,1:ny,1) = work2d(ix,1:ny)
     end do 
  end if
  !
  !***  Writes the 2D variable
  !
  if( nf90_open(TRIM(luncname),NF90_WRITE, ncID) /= 0 ) &
       call runend('write_nc_var2d : Error in nf90_open')
  if( nf90_inq_varid(ncID,TRIM(var2d(i2d)),var2d_ID(i2d)) /= 0) &
       call runend('write_nc_var2d : Error getting var_ID')
  if( nf90_put_var(ncID, var2d_ID(i2d), work, start=(/1,1,itt+1/), &
       count=(/nx,ny,1/) ) /= 0 ) &
       call runend('write_nc_var2d: error in nf90_put_var')
  if( nf90_close(ncID) /= 0) &
       call runend('write_nc_var2d : Error in nf90_close')
  !
  return
end subroutine write_nc_var2d
