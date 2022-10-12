subroutine readbpts
  !**************************************************************
  !*
  !*    Opens and read the background (problemname.bkw) file.
  !*
  !**************************************************************
  use KindType
  use InpOut
  use Master
  implicit none  
  logical               :: go_on=.true.
  integer  (ip)         :: ipart
  integer  (ip)         :: nword,npar
  character(len=s_long) :: card
  character(len=s_long) :: words(nwormax)
  real(rp)              :: param(nparmax)
  real(rp)              :: rho,diam,mass,sphe,lon,lat,z
  !
  !*** Open the file
  !
  open(15,file=TRIM(fbkw),status='old',err=150)
  !
  !*** Read the file, allocate particle structure.
  !
  go_on=.true.
  do while(go_on)
     ! Read the total number of particles for the bkw file
     read(15,'(a256)',END=103) card
     call sdecode(card,words,param,nword,npar)
     ! 
     if(TRIM(words(1)).eq.'TOTAL_PARTICLES') then
        go_on = .false.
        phase(1)%npart=INT(param(1))
     end if
  end do
  !
  !Allocate Particle structure
  !
  allocate(part(phase(1)%npart))
  !
  !*** Read the rest of the file
  !
  ! Read the particle information for the bkw file
  do while(1.eq.1)
     read(15,*,END=104)ipart,rho,diam,mass,sphe,lon,lat,z
     part(ipart)%rho =rho
     part(ipart)%diam=diam
     part(ipart)%mass=mass
     part(ipart)%sphe=sphe
     part(ipart)%lon=lon
     part(ipart)%lat=lat
     part(ipart)%z  =z
     part(ipart)%state=1
     !
     !*** Check that the lon/lat is inside the domain
     !
     if((part(ipart)%lon.gt.grid%lonmax).or.(part(ipart)%lon.lt.grid%lonmin))&
         call runend('Particle Longitude is not inside the domain')
     if((part(ipart)%lat.gt.grid%latmax).or.(part(ipart)%lat.lt.grid%latmin))&
         call runend('Particle Latitude is not inside the domain')
     if((part(ipart)%z.gt.grid%ztop).or.(part(ipart)%z.lt.0))&
         call runend('Particle Height is not inside the domain')
     if((part(ipart)%sphe.gt.1).or.(part(ipart)%sphe.lt.0))&
         call runend('Particle Sphericity is not in [0,1] range')
     !
     !*** Calculates particle shape factor psi
     !
     call setpsi(part(ipart)%psi,part(ipart)%sphe,part(ipart)%diam,sedmodel,1)
     numpart=ipart
  end do
  ! 
  !*** Close the file
  ! 
  104 close(15)
  !
  !*** Other initilizations
  !
  phase(1)%activated = .true.
  phase(1)%firstpart = 1
  phase(1)%lastpart = phase(1)%npart
  phase(1)%phase_type = 'backward'
  allocate(phase(1)%begtime(1))
  phase(1)%begtime(1) = sist
  phase(1)%endtime    = siend
  output%classes = .false.
  output%phases = .false.
  phase(1)%actnpart =phase(1)%npart
  return
  !
  !*** List of errors
  !
  103 call runend('Total particles is not found in backward file')
  150 call runend('No backward file is available. Check that te name is problemname.bkw')


end subroutine readbpts
