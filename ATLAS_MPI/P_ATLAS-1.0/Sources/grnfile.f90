subroutine grnfile(iphase)
  !************************************************************
  !*
  !*    Gets the granulometry
  !*
  !*    output: file .tgsd containing 
  !*            fc(bins),rhop(bins),diam(bins),sphe(bins)
  !*
  !************************************************************
  use KindType
  use Master
  use InpOut
  implicit none


  integer(ip) :: ic,ig
  real   (rp) :: deltafi,fi,work
  integer(ip) :: iphase
  character(len=s_name) :: phase_str
  character(len=2)      :: ipha,ext
  !
  !***  Arrays
  !
  real(rp), allocatable :: fc  (:)      ! fc(bins)
  real(rp), allocatable :: rhop(:)      ! rhop(bins)
  real(rp), allocatable :: diam(:)      ! diam(bins)
  real(rp), allocatable :: sphe(:)      ! sphe(bins)
!  integer(ip) :: memo = 0
  !
  !allocate memory
  !
  allocate(fc(phase(iphase)%bins))
  fc = 0.
!  memo = memo + rp*phase(iphase)%bins      
  !
  allocate(rhop(phase(iphase)%bins))
  rhop = 0.
!  memo = memo + rp*phase(iphase)%bins
  !
  allocate(diam(phase(iphase)%bins))
  diam = 0.
!  memo = memo + rp*phase(iphase)%bins
  !
  allocate(sphe(phase(iphase)%bins))
  sphe = 0.
!  memo = memo + rp*phase(iphase)%bins
  ! Determine the file name
  write(ipha,'(i1)')iphase
  fgrn = TRIM(problemname)//'_'//TRIM(phase(iphase)%sname)//'.tgsd'
  phase(iphase)%grnpath=TRIM(fgrn)
  !
  !***  Compute particle diameter in mm
  !
  deltafi = (phase(iphase)%fimax-phase(iphase)%fimin)/(phase(iphase)%bins-1)
  do ic = 1,phase(iphase)%bins
     fi = phase(iphase)%fimin + (ic-1)*deltafi
     diam(ic) = 2**(-fi)
  end do
  !
  !***  Density rhop(ic) and sphericity sphe(bins)
  !
  do ic = 1,phase(iphase)%bins
     fi = phase(iphase)%fimin + (ic-1)*deltafi
     rhop(ic) = (phase(iphase)%rhomax -phase(iphase)%rhomin )*(fi-phase(iphase)%fimin)/(phase(iphase)%fimax &
          -phase(iphase)%fimin) + phase(iphase)%rhomin
     sphe(ic) = (phase(iphase)%sphemax-phase(iphase)%sphemin)*(fi-phase(iphase)%fimin)/(phase(iphase)%fimax &
          -phase(iphase)%fimin) + phase(iphase)%sphemin
  end do
  !
  !***  Fraction fc(ic) (gaussian or bigaussian)
  !
  do ig = 1,phase(iphase)%ng
     do ic = 1,phase(iphase)%bins
        fi = phase(iphase)%fimin + (ic-1)*deltafi
        call M0_fi_gaussian(work,fi,phase(iphase)%fimin,phase(iphase)%fimax,deltafi,phase(iphase)%fimean(ig) &
             ,phase(iphase)%fidisp(ig))
        fc(ic) = fc(ic) + work
     end do
  end do
  fc = fc/phase(iphase)%ng
  !
  !***  Opens the file and write the granulometry inf.
  !
  open(iphase,file=TRIM(fgrn) ,status='unknown',err=100)
  !
  write(iphase,'(i5)') phase(iphase)%bins
  do ic = 1,phase(iphase)%bins
     ext = '00'
     if(ic.lt.10) then
        write(ext(2:2),'(i1)') ic
     else
        write(ext(1:2),'(i2)') ic
     end if
     phase(iphase)%classname(ic)= 'class-'//ext
     write(iphase,10) diam(ic),rhop(ic),sphe(ic),fc(ic),TRIM(phase(iphase)%classname(ic))
     !
10   format(f10.6,1x,f8.1,1x,f7.3,1x,e16.9,2x,a)
  end do
  !
  close(iphase)
  !
  !***Deallocates memory
  !
  deallocate(fc)
  deallocate(rhop)
  deallocate(diam)
  deallocate(sphe)
  !
  return
  !
100 call runend('Error opening Granulometry file '//TRIM(fgrn))

  return
end subroutine grnfile
!
!
!
subroutine M0_fi_gaussian(fc,fi,fimin,fimax,deltafi,fimean,fidisp)
  !************************************************************************
  !*
  !*   Computes the mass for a given value of fi between fimin and fimax
  !*   (mass between fi-deltafi/2 and fi+deltafi/2) for a Gaussian distribution.
  !*   The ultimate values fimin and fimax are extended to +- inifinity
  !*   to ensure mass conservation.
  !*
  !*************************************************************************
  use KindType
  implicit none
  real(rp) :: fc,fimin,fimax,deltafi,fimean,fidisp
  !
  real(rp) :: fi,fi2
  real(8) :: fer
  !
  !***  lower limit
  !
  if(fi.eq.fimin) then
     fi2 = fimin + 0.5_rp*deltafi
     fi2 = (fi2-fimean)/fidisp
     fc  = fer(fi2)
     return
  end if
  !
  !***  upper limit
  !
  if(fi.eq.fimax) then
     fi2 = fimax - 0.5_rp*deltafi
     fi2 = (fi2-fimean)/fidisp
     fc  = 1.0_rp - fer(fi2)
     return
  end if
  !
  !***  other values
  !
  fi2 = fi + 0.5_rp*deltafi
  fi2 = (fi2-fimean)/fidisp
  fc = fer(fi2)
  !
  fi2 = fi - 0.5_rp*deltafi
  fi2 = (fi2-fimean)/fidisp
  fc = abs(fc - fer(fi2))
  !
  return
end subroutine M0_fi_gaussian
!



real(kind=8) function fer(x)
  !************************************************************************
  !*
  !*    Computes the area below the typifiyed normal distribution between
  !*    -infinity and x
  !*
  !************************************************************************
  use KindType
  implicit none
  logical  :: go_on
  real(rp) ::  x
  real(rp) ::  t1,t2,dt
  !
  !***  Computes integral between t=0 and t=abs(x)
  !
  fer = 0.0_rp
  dt  = abs(x)/ 10000.0_rp
  !
  t1 = 0.0_rp
  t2 = dt
  go_on = .true.
  do while(go_on)
     fer = fer + 0.5_rp*(exp(-t1*t1/2.0_rp)+exp(-t2*t2/2.0_rp))*dt
     t1 = t2
     t2 = t2 + dt
     if(t2.ge.abs(x)) go_on = .false.
  end do
  fer = fer/sqrt(8.0_rp*atan(1.0))
  !
  if(x.ge.0.0_rp) then
     fer = 0.5_rp + fer
  else
     fer = 0.5_rp - fer
  end if
  !
  return
end function fer
