 subroutine writerest
  !******************************************************************************************
  !*
  !* makes the total mass concentration for restart
  !*
  !******************************************************************************************
  !
  use Master
  use KindType 
  use InpOut
  implicit none
  !
  logical     :: goon,go_on=.false.
  integer(ip) :: ipart, iproc
  character(len=s_name) :: filei
  integer (ip)           :: parti
  integer(ip)		:: countpart !counter, to write in restart file
!  integer (ip)		:: restpart !counter to total particles
  character(len=s_long) :: card
  integer  (ip)         :: nword,npar

  character(len=s_long) :: words(nwormax)
  real(rp)              :: param(nparmax)
  integer (ip)          :: year,month,day,hour,minute
  integer (ip)          :: state,ibin, age 
  real    (rp)          :: rho,diam,mass,sphe,psi,lon,lat,z
!**********************************
	open(88,FILE=TRIM(fgrest),status='REPLACE')
	!*data time
  	filei     = 'out_'//TRIM(problemname)//'_'//char(my_id)//'.rest'
	open(80,FILE=TRIM(filei),status='OLD')
	do while(go_on)
     		read(80,'(a256)', end=102)card
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
	close(80)
	! Write the information in the restart file
	write(88,564)'END_YEAR', year, 'END_MONTH',month, 'END_DAY',day,'END_HOUR',hour,'END_MINUTE',minute
  	write(88,565)'particle','state','bin','age','rho','diam','mass','sphe','psi','lon','lat','z'
	564 FORMAT(a,1X,i5,/,a,1X,i4,/,a,1X,i4,/,a,1X,i4,/,a,1X,i4)
	565 FORMAT(a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a,1X,a)
	! 
	!* Initialice restart particle counter
	countpart=0
        !*** Loop over procs to open and read each proc. restart file and write the information over one.
		iproc=0
	do iproc=0,num_procs-1
		filei     = 'out_'//TRIM(problemname)//'_'//char(iproc)//'.rest'
		open(iproc,FILE=TRIM(filei),status='OLD')
  		!*** read the total number of particles in this file.
  		go_on=.true.
  		do while(go_on)
     			! Read the total number of particles for the restart file iproc
     			read(iproc,'(a256)',END=103)card
     			call sdecode(card,words,param,nword,npar)
     			if(TRIM(words(1)).eq.'TOTAL_PARTICLES') then
        			go_on = .false.
        			ipart=param(1)
				close(iproc)
        			!restpart=restpart+param(1)
     			end if
  		end do
  		!*** read the particle information.
		!filei     = 'out_'//TRIM(problemname)//'_'//char(iproc)//'.rest'
		open(iproc,FILE=TRIM(filei),status='OLD')
		go_on=.true.
  		do while(go_on)
    			read(iproc,'(a256)', end=102)card
     			call sdecode(card,words,param,nword,npar)
     			if(TRIM(words(1)).eq.'particle') then
       				go_on = .false.
			        goon=.true.
      				do while(goon)
        				read(iproc,*)parti,state,ibin,age,rho,diam,mass,sphe,psi,lon,lat,z 
					countpart=countpart+1
            				write(88,570)countpart,state,ibin,age,rho,&
                			diam,mass,sphe,psi,&
                			lon,lat,z
570        				FORMAT(i9,1X,i2,1X,i4,1X,i9,1X,f11.6,1X,f11.6,1X,f24.6,1X,f11.6,1X,f11.6,1X,f11.6,1X,f11.6,1X,f15.6)
					if(parti.eq. ipart)then  
              					goon = .false.
 						close(iproc,status='delete') 
           				end if
        			end do
     			end if
  		end do
	end do !procs
	!write total number of particles
  	write(88,571)'TOTAL_PARTICLES',countpart
571 	FORMAT(a,1X,i9)
  	!
  	!*** Close the file
  	!
  	close(88)
  	return
102 call runend('No time registered for restart file')
103 call runend('Total particles is not found')
end subroutine writerest









