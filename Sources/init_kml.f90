subroutine init_kml
   !******************************************************************************************
   !* write in a kmz output file
   !******************************************************************************************
   use InpOut
   use Master
   use KindType
   implicit none
   integer(ip)     :: iphase, i
   real(rp)     :: partfraction
   integer(ip)     :: outnpart
   integer(ip)     :: ntimes
   !
   !*** First time opens the file and writes the header
   !
   open (90, FILE=TRIM(fkml_part), status='unknown')
   !
   write (90, 10) TRIM(problemname)
10 format( &
      '<?xml version="1.0" encoding="UTF-8"?>', /, &
      '<kml xmlns="http://www.opengis.net/kml/2.2">', /, &
      '<Document>', /, &
      '<name>ATLAS results for ', a, '</name>', /, &
      '<description>Select the time instant</description>', /, &
      ' ', /, &
      '<Style id="altura_1">', /, &
      '<IconStyle> ', /, &
      '<Icon> <href>../../ATLAS-1.0/Resources/b-red.png</href> </Icon>', /, &
      '<scale>0.2</scale>', /, &
      '</IconStyle> ', /, &
      '</Style>', /, &
      ' ', /, &
      '<Style id="altura_2">', /, &
      '<IconStyle> ', /, &
      '<Icon> <href>../../ATLAS-1.0/Resources/b-orange.png</href> </Icon>', /, &
      '<scale>0.2</scale>', /, &
      '</IconStyle> ', /, &
      '</Style>', /, &
      ' ', /, &
      '<Style id="altura_3">', /, &
      '<IconStyle> ', /, &
      '<Icon> <href>../../ATLAS-1.0/Resources/b-yellow.png</href> </Icon>', /, &
      '<scale>0.2</scale>', /, &
      '</IconStyle> ', /, &
      '</Style>', /, &
      ' ', /, &
      '<Style id="altura_0">', /, &
      '<IconStyle> ', /, &
      '<Icon> <href>../../ATLAS-1.0/Resources/b-green.png</href> </Icon>', /, &
      '<scale>0.2</scale>', /, &
      '</IconStyle> ', /, &
      '</Style>', /, &
      ' ')
   !
   close (90)
   !
   !*** Loop over particles
   !
   if (numpart .gt. 2500) then   ! Determines maximum number of particles to output.
   do iphase = 1, nphases
      partfraction = phase(iphase)%npart
      partfraction = partfraction/numpart  ! determines fraction of part per phase to output
      outnpart = INT(2500*partfraction) !number of particles to output in this phase
      !ARNAU: cómo elijo entre cuales.. con un incremento no conviene porque podria estar
      !sacando solo algunas clases y no todas. las
      !primeras de cada intervalo de tiempo: entonces hay que contar cuántos intervalos de tiempo hay dentro de
      !la fase para hacer esa division. qué opinas de esta forma de sacar algunas?
      ! A MI NO ME GUSTA ESTA FORMA PORQUE ES UN PROBLEMA PARA LA SIMULACION HACIA ATRAS EN CASO DE TENER EN UN ARCHIVO .BKW LA INFO DE VARIAS FUENTES
      !if(sidt.gt.0)then
      do i = 1, output%nt
         if ((sidt*phase(iphase)%begtime(1) .le. sidt*output%timesec(i)) &       !sidt change the sign and the inequality
             .and. (sidt*(phase(iphase)%begtime(1) + SIGN(output%frequency, sidt)) &
                    .gt. sidt*output%timesec(i))) then
            phase(iphase)%firstoutt = i
         else if ((sidt*phase(iphase)%endtime .ge. sidt*output%timesec(i)) &
                  .and. (sidt*(phase(iphase)%endtime - SIGN(output%frequency, sidt)) &
                         .lt. sidt*output%timesec(i))) then
            phase(iphase)%endoutt = i
         end if
      end do !i
      ntimes = phase(iphase)%endoutt - phase(iphase)%firstoutt + 1!output%nt-phase(iphase)%firstoutpart+1!
      phase(iphase)%partsout = INT(outnpart/ntimes)
      phase(iphase)%partstime = INT(phase(iphase)%npart/ntimes)
   end do ! iphases
   else if (numpart .le. 2500) then
   do iphase = 1, nphases
      phase(iphase)%firstoutt = 0
      phase(iphase)%endoutt = output%nt
      if (sidt .gt. 0) then
         phase(iphase)%partsout = phase(iphase)%npart
         phase(iphase)%partstime = 0!INT(phase(iphase)%npart/output%nt)
      else
         phase(iphase)%partsout = phase(iphase)%npart
         phase(iphase)%partstime = 0
      end if
   end do
   !
   end if !numpart

   return
end subroutine init_kml
