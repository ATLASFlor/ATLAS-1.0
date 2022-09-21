 subroutine make_part_nc(iout)
  !******************************************************************************************
  !*
  !* makes the total mass concentration
  !*
  !******************************************************************************************
  !
  use Master
  use KindType 
  use InpOut
  implicit none
  !
  logical     :: found=.false.
  integer(ip) :: ielem3d,ielem2d,iphase,ipart,iclass
  integer(ip) :: ix,iz,iy, izl, i,iout
  integer(ip) :: endpart ! last activated particle
  ! 
  !
  !************************************************************************************************************
  !*** Initialize variables
  !
  massair     = 0.0_rp
  massground  = 0.0_rp
  massinside  = 0.0_rp
  massoutside = 0.0_rp
  masstotal   = 0.0_rp
  sedim (:)=0.0_rp             !sedim(nelem2d)
  conc2d(:)=0.0_rp
  conc3d(:)=0.0_rp             !conc(nelem3d)
  if(output%classes)then
     sedim_bin(:,:,:) =0.0_rp  !sedim_bin(phase,bin,nelem2d)
     conc2d_bin(:,:,:)=0.0_rp
     conc3d_bin(:,:,:)=0.0_rp
  end if
  if(output%phases)then
     sedim_phase(:,:) =0.0_rp  !sedim_phase(phase,nelem2d)
     conc2d_phase(:,:)=0.0_rp
     conc3d_phase(:,:)=0.0_rp
  end if
	!
  do iphase=1,nphases
     ! 
     !*** Determine the last activated part for each phase
     !
     if(split)then
        endpart=phase(iphase)%totpart
     else
        endpart=phase(iphase)%actnpart+phase(iphase)%firstpart-1 !(modificacion 05/2020) para incluir split
     end if
     !
     !*** Loop over particles
     !
     do ipart=phase(iphase)%firstpart,endpart 
     !*********************************************** PARTICLE STATE 1  *******************************************  
        !
        !*** Compute concentration
        !
        if(part(ipart)%state.eq. 1)then 
        !if(part(ipart)%lon.ge.output%lonmin .and. part(ipart)%lon.le.output%lonmax .and. &
        !   part(ipart)%lat.ge.output%latmin .and. part(ipart)%lat .le.output%latmax) then !agregado reciente
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
              do while(found.eqv. .false.)
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
           ielem3d= (iz-1)*(output%nx-1)*(output%ny-1)+(iy-1)*(output%nx-1)+ix
           ielem2d= (iy-1)*(output%nx-1)+ix
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
           !****************************************** PARTICE STATE -1 *****************************************************
           !end if! particles inside output domain in air
        else if((part(ipart)%state.eq.(-1)))then
        !if(part(ipart)%lon.ge.output%lonmin .and. part(ipart)%lon.le.output%lonmax .and. &
        !   part(ipart)%lat.ge.output%latmin .and. part(ipart)%lat .le.output%latmax) then  !agregado reciente
           !
           !*** Particle ground load
           ! 
           !*** Define the 2d box where the particle is
           ! 
           ix = INT((part(ipart)%lon-output%lonmin)/output%dx) + 1 
           if(ix.eq.(output%nx+1)) ix=output%nx
           iy = INT((part(ipart)%lat-output%latmin)/output%dy) + 1  
           if(iy.eq.(output%ny+1)) iy=output%ny
           ! 
           !*** Define the grid element where the particle is
           !
           ielem2d= (iy-1)*(output%nx-1)+ix
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
        !end if !particle inside output domain, on ground
           !****************************************** PARTICE STATE -2 *****************************************************
        else if (part(ipart)%state.eq.(-2))then ! Must be added particles inside simulation domain and outside output domain
           massoutside=massoutside+part(ipart)%mass
        end if ! ipart state
     end do ! particle loop
  end do    ! phase loop
  !*** Calculate the total mass inside domain and total mass
  massinside=massair+massground
  masstotal=massinside+massoutside  !
  !*** Write the mass per track point if there are considered. In another files.
  !
  if(output%track_points)then
    do i=1,npts
       if(use_pts(i)) then
         ! open(79,file=TRIM(name_file_pts(i))//'.res',iostat=info,status='old',position='append')
         ! if( info /= 0 ) goto 100
          !calc thickness in cm (Divides value of sedimentation load by average particle density, and multiply by 100)
          load_pts(i)=sedim(ielem2dpts(i))  !kg/m2
          thick_pts(i)=load_pts(i)/output%rhomean     ! 
          thick_pts(i)=thick_pts(i)*1d2 !in cm
         ! write(79,20),   load_pts,thick_pts
!20        format(i7,1x,f13.6, 1x, f13.6)
         ! close(79)
       end if !use_pts
    end do !npts
  end if  !track_points
  ! 
  ! 
  return
  end subroutine make_part_nc
