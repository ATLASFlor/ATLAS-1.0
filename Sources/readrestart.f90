subroutine readrestart
   !**************************************************************
   !*
   !*    Opens and read the restart file.
   !*
   !**************************************************************
   use KindType
   use InpOut
   use Master
   implicit none
   logical               :: goon, go_on = .true.
   integer(ip)         :: nword, npar
   character(len=s_long) :: card
   character(len=s_long) :: words(nwormax)
   real(rp)              :: param(nparmax)
   integer(ip)          :: state, ibin, age           !to read in restart file
   integer(ip)           :: parti
   real(rp)          :: rho, diam, mass, sphe, psi, lon, lat, z
   !
   !*** Open the file
   !
   open (15, file=TRIM(frest), status='old', err=150)
   !
   !*** read and save values for each particle
   !
   go_on = .true.
   do while (go_on)
      read (15, '(a256)', end=102) card
      call sdecode(card, words, param, nword, npar)
      if (TRIM(words(1)) .eq. 'particle') then
         go_on = .false.
         goon = .true.
         do while (goon)
            read (15, *) parti, state, ibin, age, rho, diam, mass, sphe, psi, lon, lat, z
            if (parti .eq. phase(nphases)%npart) then
               goon = .false.
            end if
            parti = phase(nphases)%firstpart - 1 + parti
            part(parti)%state = state
            part(parti)%ibin = ibin
            part(parti)%age = age
            part(parti)%rho = rho
            part(parti)%diam = diam
            part(parti)%mass = mass
            part(parti)%sphe = sphe
            part(parti)%psi = psi
            part(parti)%lon = lon
            part(parti)%lat = lat
            part(parti)%z = z
         end do
      end if
   end do
   !
   !*** Close the file
   !
   close (15)
   return
   !
   !*** List of errors
   !

150 call runend('No restart file is available. Check that te name is out_problemname.rest')
102 call runend('No time registered for restart file')
103 call runend('Total particles is not found')
end subroutine readrestart
