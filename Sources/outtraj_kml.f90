subroutine outtraj_kml(iout)
!******************************************************************************************
!* write in a kml output file
!******************************************************************************************
   use InpOut
   use Master
   use KindType
   use TimeFun
   implicit none

   integer(ip)   :: iout
   integer(ip)   :: ipart, itime
   real(rp)::trajlon(output%nt, npart)
   real(rp)::trajlat(output%nt, npart)
   real(rp)::trajz(output%nt, npart)
   character(len=1)::ext
   character(len=11)::traj_number

!
!*** First time opens the file and writes the header
!
   if (iout .eq. 1) then
      open (92, FILE=TRIM(fkml_traj), status='unknown')
!
      write (92, 10) TRIM(problemname)
10    format( &
         '<?xml version="1.0" encoding="UTF-8"?>', /, &
         '<kml xmlns="http://www.opengis.net/kml/2.2">', /, &
         '<Document>', /, &
         '', /, &
         '<name>Trajectories</name>', /, &
         '', /, &
         '<Style id="orange">', /, &
         '<LineStyle>', /, &
         '<color>7f00ffff</color>', /, &
         '<width>4</width>', /, &
         '</LineStyle>', /, &
         '<PolyStyle>', /, &
         '<color>7f00ff00</color>', /, &
         '</PolyStyle>', /, &
         '</Style>')
      close (90)
!
   end if
!
!*** Loop over particles
!
   do ipart = 1, npart   !ver criterio para elegir algunas particulas
      trajlon(iout, ipart) = part(ipart)%lon
      trajlat(iout, ipart) = part(ipart)%lat
      trajz(iout, ipart) = part(ipart)%z
   end do
!
!*** Last time only
!
   if (iout .eq. output%nt) then
!
      do ipart = 1, npart   !ver criterio para elegir sólo algunas particulas
         ! Define trajectory number:
         write (ext, '(i1)') ipart
         traj_number = 'Trajectory'//ext
         !
         open (92, FILE=TRIM(fkml_traj), access='append', status='old')
         write (92, 11) traj_number
11       format( &
            '<Placemark>', /, &
            '<name>', a, '</name>', /, &
            '<visibility>1</visibility>', /, &
            '<description>At a given h</description>', /, &
            '<styleUrl>#orange</styleUrl>', /, &
            '<LineString>', /, &
            '<extrude>0</extrude>', /, &
            '<tessellate>1</tessellate>', /, &
            '<altitudeMode>relativeToGround</altitudeMode>', /, &
            '<coordinates>')
         do itime = 1, output%nt
            write (92, 12) trajlon(itime, ipart), trajlat(itime, ipart), trajz(itime, ipart)
12          format( &
               f11.6, f11.6, f10.4)
         end do
         write (92, 13)
13       format( &
            '</coordinates>', /, &
            '</LineString>', /, &
            '</Placemark>')
      end do
      write (92, 14)
14    format( &
         '</Document>', /, &
         '</kml>')
      close (92)
!
   end if

   return
end subroutine outtraj_kml
