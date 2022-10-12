     subroutine splitting(iphase,npart2)
  !******************************************************************************************
  !*    
  !*    Generate the splitting only for in conditions particles
  !*    iphase      : Input parameter. Phase of source term
  !*    npart2      : Input parameter. particle range to calculate ADS equation
  !*    
  !*
  !******************************************************************************************
  use KindType 
  use Master                             
  implicit none
  !
  integer (ip) :: ipart
  integer (ip) :: iphase,npart2
  !integer(ip),save::npartsplit=0
! 

  do ipart=phase(iphase)%firstpart,npart2
     if ((part(ipart)%age.gt.dts).and.(part(ipart)%state.eq.1))then
         part(ipart)%age=0
         phase(iphase)%npartsplit=phase(iphase)%npartsplit+1
         npart2=phase(iphase)%npart+phase(iphase)%npartsplit
         part(ipart)%mass=part(ipart)%mass/2
         part(npart2)=part(ipart)
     end if
  end do
!       If there is not particles in that conditions just in this time, then change npart2 according the actual value of npartsplit
         npart2=phase(iphase)%npart+phase(iphase)%npartsplit
! 
  return    
end subroutine splitting
