subroutine reduce_mass
  !******************************************************************************************
  !*
  !* makes the total mass concentration
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
  !
	call MPI_REDUCE(massair, massair_global, 1, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
	call MPI_REDUCE(massground, massground_global, 1, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
	call MPI_REDUCE(massinside, massinside_global, 1, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
	call MPI_REDUCE(massoutside, massoutside_global, 1, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
	call MPI_REDUCE(masstotal, masstotal_global, 1, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
        ! 
	call MPI_REDUCE(conc2d, conc2d_global, output%nelem2d, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
	call MPI_REDUCE(conc3d, conc3d_global, output%nelem3d, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
	call MPI_REDUCE(sedim, sedim_global, output%nelem2d, MPI_REAL, MPI_SUM, 0, mpi_global, ierr)
	!
  return
end subroutine
	
