subroutine readnumpart
  !**************************************************************
  !*
  !*    Opens and read the (restart or another?) file.
  !*    Check times.
  !*    numpart: input/output parameter.
  !*
  !**************************************************************
  use KindType
  use InpOut
  use Master
  implicit none    
  !
  logical               :: go_on=.true.
  real     (rp)         :: year,month,day,hour,minute
  integer  (ip)         :: nword,npar
  character(len=s_long) :: card
  character(len=s_long) :: words(nwormax)
  real(rp)              :: param(nparmax)
  !
  !Condition to file class. (To add another option)
  if(restart)then 
    nphases=nphases+1 ! restart is the last phase
  !
  !*** Open the file
  !
  open(15,file=TRIM(frest),status='old',err=150)
  !
  !*** Read time data and check for unconsistency
  !
  do while(go_on)
     read(15,'(a256)', end=102)card
     call sdecode(card,words,param,nword,npar)
     if(TRIM(words(1)).eq.'particle') then
        go_on = .false.
     else if(TRIM(words(1)).eq.'END_YEAR'  ) then
        year = param(1)
     else if(TRIM(words(1)).eq.'END_MONTH' ) then
        month = param(1)
     else if(TRIM(words(1)).eq.'END_DAY'   ) then
        day = param(1)
     else if(TRIM(words(1)).eq.'END_HOUR'  ) then
        hour= param(1)
     else if(TRIM(words(1)).eq.'END_MINUTE') then
        minute= param(1)
     end if
  end do
  !*** Check the time
  if((year.ne.ibyr).or.(month.ne.ibmo).or.(day.ne.ibdy).or.(hour.ne.ibhr))&
       call runend('Restart time is not consistent with run time')
  !
  !*** Initialize data
  !
  ! Determine the firstpart for the restart phase
  phase(nphases)%firstpart=phase(nphases-1)%lastpart+1
  phase(nphases)%activated=.true.     ! Initialize phase(start in active mode!)
  phase(nphases)%phase_type='restart' ! Determine the phase_type
  phase(nphases)%endtime=sist-1       ! endtime is before the simulation start to avoid a condition in main loop in atlas.f90
  allocate(phase(nphases)%begtime(1)) ! allocate begtime array 
  phase(nphases)%begtime(1)=sist-1    ! begtime is before simulation start 
  !
  !*** Read data in the restart file
  !
  go_on=.true.
  do while(go_on)
     ! Read the total number of particles for the restart file
     read(15,'(a256)',END=103)card
     call sdecode(card,words,param,nword,npar)
     if(TRIM(words(1)).eq.'TOTAL_PARTICLES') then
        go_on = .false.
        phase(nphases)%npart=param(1)
     end if
  end do
  close(15)
  ! Determine the lastpart for restart phase and the number of active particles
  phase(nphases)%lastpart=phase(nphases)%firstpart+ phase(nphases)%npart
  phase(nphases)%actnpart=phase(nphases)%npart
  !
  !** Determine the total number of particles to allocate the particle structure
  !
  numpart=numpart+phase(nphases)%npart
  end if
  !
  return
  !
  !*** List of errors
  !
150 call runend('No restart file is available. Check that te name is out_problemname.rest')
102 call runend('No time registered for restart file')
103 call runend('Total particles is not found')
  end subroutine readnumpart
