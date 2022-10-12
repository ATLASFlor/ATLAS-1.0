subroutine initialize
  !*********************************************************************
  !* 
  !*   This subroutine reads input data and creates
  !*   meteorological and postprocess grids.
  !*
  !**********************************************************************
  use InpOut
  use Master
  implicit none
  save
  logical           :: out_meteo = .false.
  integer(ip)       :: imodel,ipart,iphase, ibin
  character (len=2) :: ext
  character (len=8) :: str
  ! 
  integer(ip)       :: countpart1, countpart2
  logical           :: goon=.TRUE.
  logical           :: go_on=.TRUE.
  integer(ip)       :: a
  real   (rp)       :: r=0.7
  integer(ip),allocatable :: npbin  (:)
  integer(ip),allocatable :: npbindt(:)
  integer(ip)       :: countpart,countpartbindt
  integer(ip)       :: ibin2
  !
  !*** Define the file names
  !
  finp      =         TRIM(problemname)//'.inp'
  fbkw      =         TRIM(problemname)//'.bkw'
  fpts      =         TRIM(problemname)//'.pts'
  !
  fnc_met   = 'out_'//TRIM(problemname)//'.meteo.nc'
  fnc_part  = 'out_'//TRIM(problemname)//'.part.nc'
  fkml_part = 'out_'//TRIM(problemname)//'.part.kml'
  frest     = 'out_'//TRIM(problemname)//'_'//char(my_id)//'.rest'
  if(my_id.eq.0)  fgrest = 'out_'//TRIM(problemname)//'.rest'
  fkml_traj = 'out_'//TRIM(problemname)//'.traj.kml'
  !
  !*** Generate a random number array, neccesary to compute the diffussion.
  !
  call init_random_seed()
  call RANDOM_NUMBER(rannumb)
  nran=1 !Set the first element of random array.
  !
  !*** Reads the input file
  ! 
  call readinp
  !
  !*** Create Background mesh
  !
  call set_back_mesh
  !
  !*** Reads the meteo files 
  !
  call set_meteo_model
  !
  !*** Fill meteorological data in a netcdf file
  !
  do imodel=1,numod
     if(meteo(imodel)%postprocess) out_meteo=.true. ! If the user want output meteorological inf.
  end do
  if(out_meteo) call outmeteo
  !
  !*** What follow is only for forward mode. Initialize particle structure. "Else if" is for backward
  !
  if(sidt.gt.0)then
     !   
     !*** Reads (or creates) the granulometry for each phase and, if necessary, accounts
     !*** for aggregation
     !

     call readgrn 
     !
     !*** Redefine the number of particles per phase and the number of particles released
     !*** per time interval during which the phase is active
     !
     do iphase = 1,nphases
        !
        phase(iphase)%activated = .false.  ! Initialize phases
        phase(iphase)%actnpart  = 0        ! Initialize activated particles
        phase(iphase)%erumass   = 0        ! Initialie the MER rate counter
        !
        !*** npart must be greater than the number of duration/sidt 
        !*** (in backward case this quotient is negative, then it does not change the number of particles)
        !*** Duration is considered only for subphase with column height greater than 0
        !
        if(phase(iphase)%npart.lt.(phase(iphase)%duration/sidt))then  
           phase(iphase)%npart=INT(phase(iphase)%duration/sidt)      
           write(ext,'(i2)')iphase    ! Write in log file that the number of particles are increased
           str='PHASE_'//TRIM(ext)
           call wriwar(TRIM(str)//'The number of particles was insufficient. It was automatically increased')
        end if
        ! Determine the number of particles as a multiple of the number of bins
        phase(iphase)%npart = phase(iphase)%bins*(INT(phase(iphase)%npart/phase(iphase)%bins)+1)
        ! Determine the number of particles to be released per time interval, 
        ! After, the increment of npdt particle is done only when the simulation 
        !is in a phase/subpahse with column height greater than 0.
        phase(iphase)%npdt  = INT(phase(iphase)%npart/(INT(phase(iphase)%duration/sidt)+1))
        ! If the phase%duration is shorter than simulation dt, all particles are liberated in the same time
        if(sidt.gt.phase(iphase)%duration)phase(iphase)%npdt = phase(iphase)%npart
     end do !phases
     !
     !*** Allocate Particles Structure
     !
     numpart = 0    ! total number of particles. 
     do iphase=1,nphases
        ! Determine which particle are the first for each phase
        phase(iphase)%firstpart = numpart + 1
        ! 
        if(split)then
           phase(iphase)%totpart=phase(iphase)%npart*2**INT((siend-phase(iphase)%begtime(1))/dts)
           phase(iphase)%lastpart  = phase(iphase)%firstpart + phase(iphase)%totpart - 1
           numpart = numpart+phase(iphase)%totpart
        else
           phase(iphase)%lastpart  = phase(iphase)%firstpart + phase(iphase)%npart - 1
           ! Sum the total particles  until now
           numpart = numpart+phase(iphase)%npart
        end if
     end do
     !
     !*** Read another phases like restart file and initialize the phase, 
     !    check times, sum number of particles at total
     !
     call readnumpart        ! sum the total number of particles and allocate this structure
     allocate(part(numpart)) ! allocate particle structure 
     !
     !*** Initializes particles (except source term)
     !
     if(my_id .eq. 0) then
     if(restart) then
        call readrestart
        nphases=nphases-1 ! restart phase has its particles initialized
     end if
     end if
     ! phases loop
     do iphase=1,nphases 
        !
        !Allocate memory
        allocate(npbindt(phase(iphase)%bins))
        allocate(npbin (phase(iphase)%bins))
        !* Define the particles per time interval and per bin before initialize them
        !* Distribute the particles according their size following a geometric serie
        !* Calculate a, according the number of bins.
        countpart      = 0
        countpartbindt = 0
        a = int(phase(iphase)%npart*((1-r)/(1-r**phase(iphase)%bins))) 
        !
        !*** If there is an aggregation class
        ! 
        if(phase(iphase)%aggregation)then
           ! If there is aggregation, then the last bin correspond to aggregation class. 
           ! To avoid giving more weight to aaggregation class, the order inside the series is changed.
           ibin=0
           do while(go_on)
              ibin=ibin+1
              if (phase(iphase)%diam(phase(iphase)%bins).ge.phase(iphase)%diam(ibin))then
                 ibin2=ibin
                 go_on=.false.
              end if
           end do
           do ibin=1,ibin2
              !Correspond to a geometric serie. Then, there is a rest not taked into account, added to the last bin
              !npbin(ibin)=r**(phase(iphase)%bins-ibin)
             npbin(ibin)=int(a*(r**(phase(iphase)%bins-ibin))) 
              countpart=countpart+npbin(ibin)
              npbindt(ibin)=int(npbin(ibin)/(phase(iphase)%duration/sidt))+1
              countpartbindt=countpartbindt+npbindt(ibin)
           end do
           !Aggregation class
           npbin(phase(iphase)%bins)=int(a*(r**(phase(iphase)%bins-ibin2-1)))
           countpart=countpart+npbin(phase(iphase)%bins)
           npbindt(phase(iphase)%bins)=int(npbin(phase(iphase)%bins)/(phase(iphase)%duration/sidt))+1
           countpartbindt=countpartbindt+npbindt(phase(iphase)%bins)
           !Class with diam grather than aggregation class
           do ibin=ibin2+1,phase(iphase)%bins-2 
              !Correspond to a geometric serie. Then, there is a rest not taked into account, added to the last bin
              !npbin(ibin)=r**(phase(iphase)%bins-ibin)
              npbin(ibin)=int(a*(r**(phase(iphase)%bins-ibin-1)))
              countpart=countpart+npbin(ibin)
              npbindt(ibin)=int(npbin(ibin)/(phase(iphase)%duration/sidt))+1
              countpartbindt=countpartbindt+npbindt(ibin)
           end do
!           !The rest is for the finest class
           npbin(phase(iphase)%bins-1)   = phase(iphase)%npart-countpart 
           npbindt(phase(iphase)%bins-1) = phase(iphase)%npdt-countpartbindt
        ! 
        ! *** Without aggregation
        ! 
        else
           !If there is not aggregation, then the last bin correspond to the finest particles
           do ibin=1,phase(iphase)%bins-1
              !Correspond to a geometric serie. Then, there is a rest not taked into account, added to the last bin
              !npbin(ibin)=r**(phase(iphase)%bins-ibin)
              npbin(ibin)=int(a*(r**(phase(iphase)%bins-ibin)))
              countpart=countpart+npbin(ibin)
              npbindt(ibin)=int(npbin(ibin)/(phase(iphase)%duration/sidt))
              countpartbindt=countpart+npbin(ibin)
           end do
           npbin(phase(iphase)%bins)   = phase(iphase)%npart-countpart 
           npbindt(phase(iphase)%bins) = phase(iphase)%npdt-countpartbindt
        end if
        ! 
        ! *** INITIALIZE
        ! 
        part(:)%state=0
        part(:)%age=0
        ibin=1
        countpart1=0
        do ipart = phase(iphase)%firstpart,(phase(iphase)%firstpart+phase(iphase)%npart-1)
           countpart1=countpart1+1
           if(countpart1.gt.npbindt(ibin))then
              ibin=ibin+1
              if(ibin.gt.phase(iphase)%bins)ibin=1
              countpart1=0
           end if
           !part(ipart)%state = 0
           !part(ipart)%age=0
           part(ipart)%ibin  = ibin
           ! Determine properties according the bin
           part(ipart)%rho  = phase(iphase)%rhop(ibin)
           part(ipart)%diam = phase(iphase)%diam(ibin)
           part(ipart)%sphe = phase(iphase)%sphe(ibin)
           part(ipart)%psi  = phase(iphase)%psi (ibin)
           ! 
        end do
! 
        ! Deallocate memory
        deallocate(npbindt)
        deallocate(npbin)
     end do !phases
    ! 
    if(restart)nphases=nphases+1 !Assign the real value of nphases again
     !
  else !if sidt
     call readbpts !read backward points
  end if
  ! 
  !*** Initialize Output variables
  ! 
  call init_outvar
  ! 
  !*** Create Output-Netcdf file
  ! 
  if(my_id.eq.0)call init_netcdf
  if(my_id.eq.0)call init_kml
  return 
end subroutine initialize
