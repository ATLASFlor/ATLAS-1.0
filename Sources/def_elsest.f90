module def_elsest

   use kindtype
   integer, parameter  :: lg = kind(.true.)
   type i1p
      integer(ip), pointer :: l(:)
   end type i1p
   type poiarr
      type(octbox), pointer :: o
   end type poiarr

   type octbox
      integer(ip)               :: id          ! My global ID
      integer(ip)               :: level       ! Generation
      integer(ip)               :: npoinbox    ! Number of nodes
      integer(ip)               :: nelembox    ! Number of elements
      integer(ip)               :: childid     ! Child ID (1->4 or 1->8)
      integer(ip)               :: whoiam      ! Father or have nodes
      integer(ip), allocatable :: nodes(:)
      integer(ip), allocatable :: elems(:)    ! List of elements
      real(rp)                  :: minc(3)     ! Min coordinates
      real(rp)                  :: maxc(3)     ! Max coordinates
      type(octbox), pointer     :: parent      ! Pointer to parent
      type(octbox), pointer     :: children(:) ! Pointer to children
   end type octbox
   !
   ! General data
   !
   integer(ip)                  :: iunit(3)
   integer(ip)                  :: nthre, ndime_els, nelem_els, npoin_els
   integer(ip)                  :: kfl_memor(3)
   integer(ip), pointer     :: pelpo_els(:)
   integer(ip), pointer     :: lelpo_els(:)
   real(rp), pointer     :: coord_els(:, :)
   !
   ! Bin structures ---------------------------------------------------------
   !
   type bintype
      integer(ip)           :: nboxe
      integer(ip)           :: mboel
      integer(ip)           :: nboxx(3)
      integer(ip)           :: dataf
      integer(ip)           :: iallo
      integer(ip)           :: iwhat
      integer(ip), pointer  :: lnobo(:)
      integer(ip), pointer  :: lboel(:)
      integer(ip), pointer  :: pboel(:)
      type(i1p), pointer  :: tboel(:)
      integer(ip), pointer  :: memor(:, :)
      integer(ip)           :: memax
      integer(ip), pointer  :: kstat(:, :)
      integer(ip), pointer  :: ksear(:)
      integer(ip), pointer  :: kfirs(:)
      integer(ip), pointer  :: kseco(:)
      real(rp)              :: delta(3)
      real(rp)              :: comin(3)
      real(rp)              :: comax(3)
      real(rp)              :: lmini
      real(rp)              :: lmaxi
      real(rp), pointer  :: elcod(:, :, :)
      real(rp), pointer  :: cputi(:, :)
   end type bintype
   type octtype
      integer(ip)           :: mboel
      integer(ip)           :: iallo
      integer(ip)           :: iwhat
      integer(ip), pointer  :: lnobo(:)
      integer(ip), pointer  :: nbono(:)
      integer(ip), pointer  :: pbono(:)
      type(poiarr), pointer  :: current(:)
      type(octbox), pointer  :: tree_root
      integer(ip), pointer  :: memor(:, :)
      integer(ip)           :: memax
      integer(ip), pointer  :: kstat(:, :)
      integer(ip), pointer  :: ksear(:)
      integer(ip), pointer  :: kfirs(:)
      integer(ip), pointer  :: kseco(:)
      integer(ip)           :: limit
      integer(ip)           :: divmax
      real(rp)              :: comin(3)
      real(rp)              :: comax(3)
      real(rp)              :: lmini
      real(rp)              :: lmaxi
      real(rp), pointer  :: elcod(:, :, :)
      integer(ip), pointer  :: lboel(:)
      real(rp), pointer  :: cputi(:, :)
   end type octtype

   type(bintype), pointer  :: bin_struc(:)
   type(octtype), pointer  :: oct_struc(:)

   integer(ip), pointer  :: nboxe
   integer(ip), pointer  :: mboel
   integer(ip), pointer  :: nboxx(:)
   integer(ip), pointer  :: dataf
   integer(ip), pointer  :: iallo
   integer(ip), pointer  :: iwhat
   integer(ip), pointer  :: lboel(:)
   integer(ip), pointer  :: nbono(:)
   integer(ip), pointer  :: pboel(:)
   type(i1p), pointer  :: tboel(:)
   type(poiarr), pointer  :: current(:)
   type(octbox), pointer  :: tree_root
   integer(ip), pointer  :: chbox(:, :)
   integer(ip), pointer  :: limit
   integer(ip), pointer  :: divmax
   integer(ip), pointer  :: chdel(:, :)
   integer(ip), pointer  :: nexbo(:, :)
   integer(ip), pointer  :: memor(:, :)
   integer(ip), pointer  :: memax
   integer(ip), pointer  :: kstat(:, :)
   integer(ip), pointer  :: ksear(:)
   integer(ip), pointer  :: kfirs(:)
   integer(ip), pointer  :: kseco(:)
   real(rp), pointer  :: delta(:)
   real(rp), pointer  :: comin(:)
   real(rp), pointer  :: comax(:)
   real(rp)                 :: lmini
   real(rp)                 :: lmaxi
   real(rp), pointer  :: elcod(:, :, :)
   real(rp), pointer  :: cputi(:, :)

contains

   function elsest_intost(integ)
      !-------------------------------------
      !
      !  Convert an integer(ip) to a string
      !
      !-------------------------------------
      implicit none
      integer(ip)   :: integ
      character(20) :: elsest_intost
      character(20) :: intaux

      write (intaux, *) integ
      elsest_intost = adjustl(intaux)

   end function elsest_intost

end module def_elsest
