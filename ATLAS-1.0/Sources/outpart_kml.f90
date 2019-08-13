subroutine outpart_kml(iout)
  !******************************************************************************************
  !* write in a kmz output file 
  !******************************************************************************************
  use InpOut
  use Master
  use KindType
  use TimeFun                                     
  implicit none
  !
  integer (ip)     :: iout
  integer (ip)     :: iphase,ipart, i
  integer (ip)     :: iyr,imo,idy,ihr,imi,ise
  integer (ip)     :: firstpart, lastpart ! first and last particle to output in each iteration
  real    (rp)     :: partfraction,a
  integer (ip)     :: outnpart
  integer (ip)     :: ntimes
  character(len=4) :: ciyr
  character(len=2) :: cimo,cidy,cihr,cimi,cise
  character(len=25):: outtime1,outtime2
  character(len=2) :: timei
  character(len=6) :: strtime
  character(len=8) :: height
  !
  !*** Write the actual time in format yyyy/mm/dd_hh:mm:ss
  !
  write(timei,'(i2)')iout
  strtime='Time'//timei
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,output%timesec(iout))
  write(ciyr,'(i4.4)') iyr
  write(cimo,'(i2.2)') imo
  write(cidy,'(i2.2)') idy
  write(cihr,'(i2.2)') ihr
  write(cimi,'(i2.2)') imi
  write(cise,'(i2.2)') ise
  outtime1 = ciyr//'-'//cimo//'-'//cidy//'T'//cihr//':'//cimi//':'//cise//'Z'
  !
  call addtime(ibyr,ibmo,ibdy,0,iyr,imo,idy,ihr,imi,ise,output%timesec(iout)+SIGN(output%frequency,sidt))
  write(ciyr,'(i4.4)') iyr
  write(cimo,'(i2.2)') imo
  write(cidy,'(i2.2)') idy
  write(cihr,'(i2.2)') ihr
  write(cimi,'(i2.2)') imi
  write(cise,'(i2.2)') ise
  outtime2 = ciyr//'-'//cimo//'-'//cidy//'T'//cihr//':'//cimi//':'//cise//'Z'
  !
  !*** First time opens the file and writes the header
  !
  if(iout.eq.1) then
     open(90,FILE=TRIM(fkml_part),status='unknown') 
     !
     write(90,10) TRIM(problemname)
10   format(&     
          '<?xml version="1.0" encoding="UTF-8"?>'            ,/,&
          '<kml xmlns="http://www.opengis.net/kml/2.2">'      ,/,&  
          '<Document>'                                        ,/,&
          '<name>ATLAS results for ',a,'</name>'          ,/,&
          '<description>Select the time instant</description>', /,&
          ' ',/,&
          '<Style id="altura_1">',/,&
          '<IconStyle> ',/,&
          '<Icon> <href>../../ATLAS-1.0/Resources/b-red.png</href> </Icon>',/,&
          '<scale>0.2</scale>',/,&
          '</IconStyle> ',/,&
          '</Style>',/,& 
          ' ',/,&
          '<Style id="altura_2">',/,&
          '<IconStyle> ',/,&
          '<Icon> <href>../../ATLAS-1.0/Resources/b-orange.png</href> </Icon>',/,&
          '<scale>0.2</scale>',/,&
          '</IconStyle> ',/,&
          '</Style>',/,&
          ' ',/,&
          '<Style id="altura_3">',/,&
          '<IconStyle> ',/,&
          '<Icon> <href>../../ATLAS-1.0/Resources/b-yellow.png</href> </Icon>',/,&
          '<scale>0.2</scale>',/,&
          '</IconStyle> ',/,&
          '</Style>',/,&
          ' ',/,&
          '<Style id="altura_0">',/,&
          '<IconStyle> ',/,&
          '<Icon> <href>../../ATLAS-1.0/Resources/b-green.png</href> </Icon>',/,&
          '<scale>0.2</scale>',/,&
          '</IconStyle> ',/,&
          '</Style>',/,&
          ' ') 
     !
     close(90)
  !
  !*** Loop over particles
  !
    if(numpart.gt.2500)then   ! Determines maximum number of particles to output.
    do iphase=1,nphases
      partfraction=phase(iphase)%npart
      partfraction=partfraction/numpart  ! determines fraction of part per phase to output
      outnpart=INT(2500*partfraction) !number of particles to output in this phase
     !if(sidt.gt.0)then
!      do i=1,output%nt
!         if((sidt*phase(iphase)%begtime(1).le.sidt*output%timesec(i))&       !sidt change the sign and the inequality
!            .and.(sidt*(phase(iphase)%begtime(1)+SIGN(output%frequency,sidt))&
!            .gt.sidt*output%timesec(i)))then
!            phase(iphase)%firstoutt=i
!         else if((sidt*phase(iphase)%endtime.ge.sidt*output%timesec(i)) &
!            .and.(sidt*(phase(iphase)%endtime-SIGN(output%frequency,sidt))&
!            .lt.sidt*output%timesec(i)))then
!            phase(iphase)%endoutt=i
!         end if
!      end do !i
a=(phase(iphase)%begtime(1)-sist)/output%frequency
      if(a.eq.INT(a))then
         phase(iphase)%firstoutt=(phase(iphase)%begtime(1)-sist)/output%frequency
      else
         phase(iphase)%firstoutt=abs((phase(iphase)%begtime(1)-sist)/output%frequency)
      end if
      ! the same with endoutt
a=(phase(iphase)%endtime-sist)/output%frequency
      if(a.eq.INT(a))then
         phase(iphase)%endoutt=(phase(iphase)%endtime-sist)/output%frequency
      else
         phase(iphase)%endoutt=abs((phase(iphase)%endtime-sist)/output%frequency)+1
      end if
      ntimes=phase(iphase)%endoutt-phase(iphase)%firstoutt+1!output%nt-phase(iphase)%firstoutpart+1!
      phase(iphase)%partsout  = INT(outnpart/ntimes)
      phase(iphase)%partstime = INT(phase(iphase)%npart/ntimes)
    end do ! iphases
    else if(numpart.le.2500)then
    do iphase=1,nphases
      phase(iphase)%firstoutt=0
      phase(iphase)%endoutt=output%nt
      if(sidt.gt.0)then
         phase(iphase)%partsout  = phase(iphase)%npart
         phase(iphase)%partstime = 0!INT(phase(iphase)%npart/output%nt)
      else
         phase(iphase)%partsout  = phase(iphase)%npart
         phase(iphase)%partstime = 0
      end if
    end do
    !
    end if !numpart
  end if !iout=1
  !
  open(90,FILE=TRIM(fkml_part),access = 'append',status='old') 
  !
  write(90,100) TRIM(outtime1),TRIM(outtime1),TRIM(outtime2)
100 format(&
       '<Folder>'             , /,&
       '<name>Time:',a,'</name>'     , /,&
       '<visibility>0</visibility>',/,&
       '  <TimeSpan>'          , /,&
       '    <begin>',a,'</begin>'  , /,&
       '    <end>',a,'</end>'       , /,&
       '  </TimeSpan>' )
  !
  !Loop over phases
  do iphase=1,nphases
     !Determine the last activated part
     ! firstpart=phase(iphase)%firstpart+(i-1)*phase(iphase)%partstime
     ! endpart=phase(iphase)%actnpart+phase(iphase)%firstpart-1 
     !Loop over particles
!    do i=1,iout  !loop over the number of output to the current
    if(iout.ge.phase(iphase)%firstoutt.and.iout.le.phase(iphase)%endoutt)then   
       firstpart = phase(iphase)%firstpart!+(i-1)*phase(iphase)%partstime !determines the first particle of each output
       lastpart   = firstpart+phase(iphase)%partsout-1
       !
     do ipart=firstpart,lastpart
        !Output particles inside the domain
        if(part(ipart)%state.eq.1 .or. part(ipart)%state.eq.(-1) )then
           !Determine the style according the particle height
           if(part(ipart)%z.eq.grid%z(1))then
              height='altura_0'
           else if(part(ipart)%z.lt.INT(grid%ztop/3))then
              height='altura_1'
           else if(part(ipart)%z.lt.INT(2*grid%ztop/3))then
              height='altura_2'
           else if(part(ipart)%z.le.grid%ztop)then
              height='altura_3'
           end if
           !
           write(90,500) part(ipart)%lon,part(ipart)%lat,part(ipart)%z,TRIM(height) !,desc
500        format(&
                ' <Placemark>'          , /,&
                '   <Point> ',/,&
                '     <coordinates>', f11.6,',',f11.6,',',f15.6',','</coordinates>' , /,&
                '     <altitudeMode>relativeToGround</altitudeMode>',/, &
                '   </Point>',/,&
                '   <styleUrl>#',a,'</styleUrl>' , /,&
                ' </Placemark> ' )
        end if !state
     end do  ! ipart
   end if    ! iout between firstoutt and endoutt
  end do     ! iphase
  !
  write(90,200) 
200 format(&
       '</Folder>',/,&
       ' ' )
  close(90)
  !
  !*** Last time only
  !
  if(iout.eq.output%nt) then
     !
     open(90,FILE=TRIM(fkml_part),access = 'append',status='old') 
     write(90,1000)
1000 format(&
          '</Document>',/,&
          '</kml>')
     close(90)
     !
  end if

  return
end subroutine outpart_kml
