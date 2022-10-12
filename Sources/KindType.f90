!***************************************************************
!*
!*    Module for variable kind definition
!*
!***************************************************************
MODULE KindType
  implicit none
  save
  !
  integer    , parameter :: ip = 4               ! integer lenght
  integer(ip), parameter :: rp = 4               ! real    lenght 
  integer(ip), parameter :: s_name = 256         ! string lenght
  integer(ip), parameter :: s_long = 256         ! string lenght
  integer(ip), parameter :: s_mess = 256         ! string lenght

  !
  !*** Number of words and parameters
  ! 
  integer(ip), parameter :: nwormax = 128 
  integer(ip), parameter :: nparmax = 128

END MODULE KindType
