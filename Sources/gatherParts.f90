subroutine gather_parts
  !******************************************************************************************
  !*
  !* makes the total particle information to write a complete restart file
  !*
  !******************************************************************************************
  !
  use Master
  use KindType 
  use InpOut
  use MPI
  implicit none
  ! 
  integer(ip):: vector_size
  integer(ip):: totpart_global(num_procs),displs(num_procs)
  integer(ip):: totalpcount
  !
displs(:)=0
	!habiendo hecho el contaador de particulas en cada procesador. el total de cada uno es totpart
	call mpi_gather(totpart, 1, mpi_int, totpart_global,1,mpi_int,  0, mpi_global, ierr)
	!Sumo el total de particulas que quedaron activas en aire o depositadas (dentro del dominio)
	call MPI_REDUCE(totpart,totalpcount ,1 , MPI_INT, MPI_SUM, 0, mpi_global, ierr)
  !
  !*** Defines the structure for particles in restart
  !
  type particlesrest
     !general properties
     integer(ip):: state   ! 1:active, 0:still not out, -1: sedimented inside the domain, -2:outside the domain
     integer(ip):: ibin    ! number of bin
     integer(ip):: age     ! particle age to split
     !phisical properties
     real(rp)   :: rho     ! particle density
     real(rp)   :: diam    ! particle diameter
     real(rp)   :: mass    ! particle mass
     real(rp)   :: sphe    ! particle sphericity
     real(rp)   :: psi     ! particle shape factor(depend on the Sedimentation model)
     !geographical properties
     real(rp)   :: lon     ! particle longitude
     real(rp)   :: lat     ! particle latitude
     real(rp)   :: z       ! height above terrain
     !
     real(rp)   :: vset    ! terminal velocity 
     real(rp)   :: uvw(3)  ! Random velocities, uvw(1): ux, uvw(2): vy, uvw(3): wp
     !
     !  Air properties at the particle point
     !
     real(rp) :: ua     ! air velocity
     real(rp) :: va
     real(rp) :: wa
     real(rp) :: Ta     ! air temperature 
     real(rp) :: rhoa   ! air density
     real(rp) :: mua    ! air viscosity
     real(rp) :: qva    ! air viscosity
     real(rp) :: pblha  ! planetary boundary layer
     real(rp) :: hfxa   ! surface heat flux
     real(rp) :: usta   ! Friction velocity
     real(rp) :: rmonin ! rmonin obukhov lenght 
     real(rp) :: wst    ! convective velocity scale
     real(rp) :: drhodz ! air density vertical gradient
     !
  end type particlesres
  type(particlesres), allocatable :: partrest(:)      ! part(totalpcount)
  allocate(partrest(totalpcount))
	! 
	! Reuno la info
	call MPI_GATHERV(part,totpart , MPI_float, partrest, totpart_global, displs,MPI_float, 0, mpi_global, ierr)

!pruebo
!write(*,*) partrest(5)%mass
stop
!  partcount=0
!  do iphase=1,nphases
!     do ipart=phase(iphase)%firstpart,phase(iphase)%lastpart
!        if ((part(ipart)%state==1).or.(part(ipart)%state==-1)) then
!           partcount=partcount+1


!        restart part(ipart)%state,part(ipart)%ibin,part(ipart)%rho,&
!                part(ipart)%diam,part(ipart)%mass,part(ipart)%sphe,part(ipart)%psi,&
!                part(ipart)%lon,part(ipart)%lat,part(ipart)%z
!	call MPI_GATHER(, restart(,), 1, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)

!    end if
!     end do
!  end do
!	call MPI_REDUCE(partcount, partcount_global, 1, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)

	!
	return
	end subroutine
