   subroutine readgrn
  !******************************************************************
  !*
  !*    Reads granulometry from the granulometry files.  
  !*    If it is not granulometry file for any phase, then create it
  !*
  !******************************************************************
  use KindType
  use InpOut
  use Master
  implicit none
  !
  integer(ip)           :: istat,ic,iphase,nc
  character(len=s_mess) :: message
  character(len=2)      :: ipha, ext
  real     (rp)         :: fca ! fraction of mass in the aggregate class
  !
  !***  Reads data
  !
  do iphase=1,nphases
     !
     if(phase(iphase)%grnpath.eq.'NONE') then  ! If it is not a granulometry file for this phase then create it
        !
        !*** Creates the granulometry file (without aggregates)
        ! 
     call grnfile(iphase)
     end if
     !
     !*** Reads granulometry
     !
     fgrn  = phase(iphase)%grnpath
     fgrn2 = TRIM(problemname)//'_'//TRIM(phase(iphase)%sname)//'.grn' 
     !
     !***  First number of particles and gas species
     !
     call get_granulometry_nclass(fgrn,phase(iphase)%bins,istat,message)
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     if(phase(iphase)%bins < 1    ) call runend('At least 1 particle class must exist   ')
     if(phase(iphase)%bins > ncmax) call runend('Too many particle classes, change ncmax')
     !
     !*** Reads particles
     !
     call get_granulometry_value(fgrn,'DIAMETER',phase(iphase)%diam,istat,message)
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     phase(iphase)%diam(1:phase(iphase)%bins) = phase(iphase)%diam(1:phase(iphase)%bins)/1d3  ! Convert diameter to m
     !
     call get_granulometry_value(fgrn,'DENSITY',phase(iphase)%rhop,istat,message)
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)  
     !
     call get_granulometry_value(fgrn,'SPHERICITY',phase(iphase)%sphe,istat,message)
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     phase(iphase)%sphe(1:phase(iphase)%bins) = min(phase(iphase)%sphe(1:phase(iphase)%bins),1.0_rp)
     !
     call get_granulometry_value(fgrn,'FRACTION',phase(iphase)%fc,istat,message)
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     call get_granulometry_name(fgrn,'NAME',phase(iphase)%classname,phase(iphase)%bins,istat,message)
     if(istat > 0) call wriwar(message)
     if(istat < 0) call runend(message)
     !
     !*** Computes the averaged particle density (without taking into account aggregates)
     !
     phase(iphase)%rhomean = 0.0_rp
     do ic = 1,phase(iphase)%bins
        phase(iphase)%rhomean = phase(iphase)%rhomean + phase(iphase)%rhop(ic)*phase(iphase)%fc(ic)
     end do
     !
     !*** Aggregates   
     !
     if(phase(iphase)%aggregation)then
        !
        phase(iphase)%bins=phase(iphase)%bins+1   ! add one class
        !
        phase(iphase)%classname(phase(iphase)%bins) = 'aggregate'
        phase(iphase)%sphe     (phase(iphase)%bins) = 1.0_rp
        phase(iphase)%diam     (phase(iphase)%bins) = phase(iphase)%aggsize
        phase(iphase)%rhop     (phase(iphase)%bins) = phase(iphase)%aggrho
        phase(iphase)%fc       (phase(iphase)%bins) = 0.0_rp
        !
     end if
     !
     !*** Calculates phase(iphase)%psi(bins) depending on the model
     !
     call setpsi(phase(iphase)%psi,phase(iphase)%sphe,phase(iphase)%diam,sedmodel,phase(iphase)%bins)
     !
     !*** Aggregation
     !
     SELECT CASE(phase(iphase)%aggmodel)
     case('NONE')
        !
        continue
        !
     case('PERCENTAGE')
        !
        !   Computes aggregation according to a percentage
        !   All classes below aggsize are reduced with
        !   a fixed user-defined percentage
        !
        fca = 0.0_rp            ! fraction of mass in the aggregate class
        do ic = 1,phase(iphase)%bins-1
           if(phase(iphase)%aggsize.gt.phase(iphase)%diam(ic)) then
              fca    = fca + phase(iphase)%aggfrac*phase(iphase)%fc(ic)
              phase(iphase)%fc(ic) = (1.0_rp-phase(iphase)%aggfrac)*phase(iphase)%fc(ic)
           end if
        end do
        phase(iphase)%fc(phase(iphase)%bins) = fca
        !
     case('CORNELL')
        !
        !   Computes aggregation according to the Cornell model
        !   The aggregate class is made of
        !       50% of particles with 4<phi<3
        !       75% of particles with 4<phi<5
        !       90% of particles with phi>5
        !
        do ic = 1,phase(iphase)%bins-1
           if(phase(iphase)%diam(ic).le.0.000125_rp.and.phase(iphase)%diam(ic).gt.0.0000625_rp) then            ! 4<phi<3
              phase(iphase)%fc(phase(iphase)%bins) = phase(iphase)%fc(phase(iphase)%bins) + 0.5_rp*phase(iphase)%fc(ic)
              phase(iphase)%fc(ic   )              = 0.5_rp*phase(iphase)%fc(ic)
           else if(phase(iphase)%diam(ic).le.0.0000625_rp.and.phase(iphase)%diam(ic).ge.0.00003125_rp) then    ! 4<phi<5
              phase(iphase)%fc(phase(iphase)%bins) = phase(iphase)%fc(phase(iphase)%bins) + 0.75_rp*phase(iphase)%fc(ic)
              phase(iphase)%fc(ic   )              = 0.25_rp*phase(iphase)%fc(ic)
           else if(phase(iphase)%diam(ic).lt.0.00003125_rp) then                                 ! phi>5
              phase(iphase)%fc(phase(iphase)%bins) = phase(iphase)%fc(phase(iphase)%bins) + 0.9_rp*phase(iphase)%fc(ic)
              phase(iphase)%fc(ic   )              = 0.1_rp*phase(iphase)%fc(ic)
           end if
        end do
     END SELECT
     !
     !*** Writes the file with the modified granulometry (aggregation + volatiles)
     !
     if(my_id.eq.0)then
     	open(98,file=TRIM(fgrn2) ,status='unknown',err=100)
     	write(98,'(i5)') phase(iphase)%bins
     	do ic = 1,phase(iphase)%bins
     	   write(98,10) 1d3*phase(iphase)%diam(ic),phase(iphase)%rhop(ic), &
     	        phase(iphase)%sphe(ic),phase(iphase)%fc(ic),TRIM(phase(iphase)%classname(ic))
10   	   format(f10.6,1x,f8.1,1x,f7.3,1x,e16.9,2x,a)
     	end do
     	close(98)
     end if

     if(my_id.ne.0)then
	open(90,file=TRIM(fgrn), status='old',err=110)
	close(90, status='delete')
     end if
     !
  end do
  !
  return
  !
  !*** Errors
  !
100 call runend('Error opening Granulometry file '//TRIM(fgrn2))
110 call runend('Error opening TGSD file for removing'//TRIM(fgrn))
end subroutine readgrn
!
!
!
subroutine get_granulometry_nclass(fname,nc,istat,message)
  !**************************************************************************
  !*
  !*    Gets the number of granulometric classes
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*
  !*    OUTPUT:
  !*    integer         nc       Number of granulometric classes
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !**************************************************************************
  use kindtype
  implicit none
  !
  character(len=*)  :: fname,message
  integer(ip)       :: nc,istat
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Opens the file
  !
  open(90,FILE=fname(1:LEN_TRIM(fname)),STATUS='old',ERR=100)
  !
  !***  Reads
  !
  read(90,*,ERR=101) nc
  close(90)
  !
  !***  Successful end
  !
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_granulometry_nclass: error opening the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_granulometry_nclass: error reading the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_granulometry_nclass
!
!
!
subroutine get_granulometry_value(fname,word,val,istat,message)
  !**************************************************************************
  !*
  !*    Gets the a granulometric property
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*    character*(*)   word     Property to extract. Possible values are
  !*                             'DIAMETER' (in mm)
  !*                             'DENSITY'
  !*                             'FRACTION'
  !*                             'SPHERICITY'
  !*
  !*    OUTPUT:
  !*    real            val(nc)  Values of the property
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !*************************************************************************
  use kindtype
  implicit none
  !
  character(len=*)  :: fname,message,word
  integer(ip)       :: nc,istat
  real   (rp)       :: val(*)
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen,ipos,ic
  real   (rp)           :: rvoid(4)
  


integer(ip)::iostat
  !
  !***  Initializations
  !
  ilen = LEN(message)
  message(:)  = ' '
  istat = 0
  !
  !***  Checks that word exists in the dictionary
  !
  if(word(1:LEN_TRIM(word)).eq.'DIAMETER') then
     ipos = 1
  else if(word(1:LEN_TRIM(word)).eq.'DENSITY') then
     ipos = 2
  else if(word(1:LEN_TRIM(word)).eq.'SPHERICITY') then
     ipos = 3
  else if(word(1:LEN_TRIM(word)).eq.'FRACTION') then
     ipos = 4
  else
     istat = -1
     mymessage = 'get_granulometry_value: word not found in the dictionary'
     message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
     return
  end if
  !
  !***  Opens the file
  !
  open(90,FILE=TRIM(fname),STATUS='old',ERR=100)
  !
  !***  Reads
  !
  read(90,*,ERR=101) nc
  do ic = 1,nc
     read(90,*,IOSTAT=iostat)rvoid(1),rvoid(2),rvoid(3),rvoid(4)!desde iteger>  ERR=101) rvoid(1),rvoid(2),rvoid(3),rvoid(4)
if( iostat < 0 )then
   write(6,'(A)') 'Warning: File containts less than 10 entries', ic
   exit
  else if( iostat > 0 )then
   write(6,'(A)') 'Error: error reading file'
   stop
  end if
     val(ic) = rvoid(ipos)
  end do
  close(90)
  !
  !***  Successful end
  !
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_granulometry_value: error opening the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_granulometry_value: error reading the '// &
       'granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_granulometry_value
!
!
!
subroutine get_granulometry_name(fname,word,val,nc,istat,message)
  !**************************************************************************
  !*
  !*    Gets a string associated with a class
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*    character*(*)   word     Property to extract. Possible values are
  !*                             'NAME'
  !*
  !*    OUTPUT:
  !*    character*(*)   val(nc)  Values of the property
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message 
  !*
  !*****************************************************************************
  use kindtype
  implicit none
  !
  integer(ip)       :: nc,istat 
  character(len=*)  :: fname,message,word,val(nc)
  !
  character(len=s_mess) :: mymessage
  integer(ip)           :: ilen,ipos,ic,ncc
  real   (rp)           :: rvoid(4)
  !
  !***  Initializations
  !
  ilen = LEN(message)  
  message(:)  = ' '
  istat = 0
  !
  !***  Checks that word exists in the dictionary
  !
  if(word(1:LEN_TRIM(word)).eq.'NAME') then
     ipos = 1
  else
     istat = -1
     mymessage = 'get_granulometry_name: word not found in the dictionary'
     message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
     return
  end if
  !
  !***  Opens the file
  !
  open(90,FILE=TRIM(fname),STATUS='old',ERR=100)
  !
  !***  Reads
  !
  read(90,*,ERR=101) nc
  do ic = 1,nc
     read(90,*,ERR=101) rvoid(1),rvoid(2),rvoid(3),rvoid(4),val(ic)
  end do
  close(90)
  !
  !***  Successful end
  ! 
  return
  !
  !***  List of errors
  !
100 istat = -1
  mymessage = 'get_granulometry_name: error opening the granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
101 istat = -1
  mymessage = 'get_granulometry_name: error reading the granulometry file '//TRIM(fname)
  message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
  return
  !
end subroutine get_granulometry_name
