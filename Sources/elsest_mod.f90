
MODULE Elsest_mod
!***************************************************************
!*
!*  Module for interpolation (q1 elements only)
!*
!***************************************************************
!subroutine elsest(&
!     itask,ipara,mnode,ndime,npoin,nelem,&
!     lnods,ltopo,coord,xcoor,rpara,ifoun,&
!     shapt,derit,coloc)
! INPUT
!    ITASK (2).... =0: Allocate memory for the bin
!                      structures. Should be called only once
!                  =1: Preprocess
!                  =2: Element Search
!                  =4: Deallocate memory for structures
!                  Structures would be automatically allocated in case they are not           
!    MNODE ....... number of nodes per elements
!    NDIME ....... Dimension
!    NPOIN ....... Number of nodes of the mesh
!    NELEM ....... Number of elements of the mesh
!    LNODS ....... Conectivity array
!    LTOPO ......  = 0 quads (2d) or hexahedra (3d)
!                  = 1 for triangle (2d) or tetra (3d)
!                  = 2 Pentahedra 
!                  = 3 Pyramid
!    IPARA(*) .... Integer input parameters
!                  IPARA( 1) = # of bins in x               
!                  IPARA( 2) = # of bins in y                
!                  IPARA( 3) = # of bins in z                
!                  IPARA( 4) = Data format (0=type,1=list)   
!                  IPARA( 5) = Max # of background meshes     
!                  IPARA( 6) = Not used            
!                  IPARA( 7) = output unit (0 dfor no output) 
!                  IPARA( 8) = Search strategy (0=bin)
!    RPARA(*) .... Real input parameters
!                  RPARA( 1) = Tolerance                      (BOTH)
!                  RPARA( 3) = Patch bounding box X min
!                  RPARA( 4) = Patch bounding box X max
!                  RPARA( 5) = Patch bounding box Y min
!                  RPARA( 6) = Patch bounding box Y max
!                  RPARA( 7) = Patch bounding box Z min
!                  RPARA( 8) = Patch bounding box Z max
! OUTPUT
!    IFOUN ....... =-1: Point out of bounding box
!                  = 0: Element not found
!                  > 0: Host element
!    SHAPT ....... = shape function values inside the element (0<=shapt<=1)
!    DERIT ....... = derivative values in local mesh
!    COLOC ....... = coordinates in reference mesh (-1<= COLOC<= 1 for ltopo=0 )
!                                                  (0 <= COLOC<= 1 when ltopo =1)
!-----------------------------------------------------------------------
!   EXAMPLE OF USE:
!    Elsest variables
!
! 
!     real   (rp) :: xcoor(3), coloc(3)
!     real   (rp) :: shapf(nnode), deriv(3,nnode), relse(20)
!     integer(ip) :: ielse(10), ltopo, ifoun
     
!
!   Elsest initialization
!    

!     ielse(1)      = 80                                     ! nx
!     ielse(2)      = 80                                     ! ny
!     ielse(3)      = 10                                     ! nz
!     ielse(4)      = 1                                      ! data format (0=type,1=linked list)
!     ielse(5)      = 10                                     ! Maximum number of possible meshes
!     ielse(6)      = 1                                      ! Not used
!     ielse(7)      = 0                                      ! Output unit
!     ielse(8)      = 0                                      ! Search strategy (0=bin,1=Quad)


!     relse(1)      = 0.001_rp                                ! Tolerance for iteration
!     relse(3)      =-1e12_rp                                 ! Patch bounding boxs
!     relse(4)      = 1e12_rp                                   
!     relse(5)      =-1e12_rp
!     relse(6)      = 1e12_rp
!     relse(7)      =-1e12_rp
!     relse(8)      = 1e12_rp
!
!    .....
!
!           call elsest(&
!                2,ielse,nnode,ndime,npoin,nelem,&
!                lnods,ltopo,coord,xcoor,relse,&
!                ielem,shapf,deriv,coloc)
!           ifoun = ielem 
!           if (ifoun <1)  call runend('ERROR:ELEMENT NOT FOUND')
              ! 
              !*** Interpolates values
              !
!             do inode = 1,nnode
!                kpoin = lnods(inode, ielem)
!                u(ipoin) = u(ipoin) +  U_Backg(kpoin)*shapf(inode)
!             end do
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

   use KindType

   IMPLICIT NONE
   SAVE

  interface elsest_memchk
     module procedure memrp1,memrp2,memrp3,&   ! reals
          &           memrp4,&
          &           memip1,memip2,memip3,&   ! 4-byte integers 
          &           memi81,memi82,memi83,&   ! 8-byte integers 
          &           mei1p1,mei1p2,&          ! types
          &           memlg1                   ! logical
  end interface 

!  
!
!
CONTAINS
!
!
!
subroutine elsest(&
     itask,ipara,mnode,ndime,npoin,nelem,&
     lnods,ltopo,coord,xcoor,rpara,ifoun,&
     shapt,derit,coloc)
  !-----------------------------------------------------------------------
  !****f* Elsest
  ! NAME
  !    Elsest
  ! DESCRIPTION
  !    EL      SE     ST
  !      ement   arch   rategies   L I B R A R Y
  !-----------------------------------------------------------------------
  use def_elsest, only     :  iunit,nthre,kfl_memor,lmini,lmaxi
  implicit none
  integer(ip), intent(in)  :: itask
  integer(ip), intent(in)  :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)  :: lnods(*)
  integer(ip), intent(in)  :: ltopo
  real(rp),    intent(in)  :: coord(*),xcoor(*)
  real(rp),    intent(in)  :: rpara(*)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: shapt(*),derit(*),coloc(*)
  integer(ip)              :: ithre,recmeth,ielem,idime,itotn,ipara(*)
  integer(ip)              :: ipoin,inode,itote, nnode, imesh
  real(rp)                 :: elcod(ndime,mnode)
  integer(ip), save        :: lmesh=0,ipass=0

  nnode    = mnode
  imesh    = 1
  ipara(8) = 0
  !
  ! Errors
  !
  if( imesh > ipara(5) ) call runend('ELSEST: WRONG MESH NUMBER')
  !
  ! Initialization
  !
  if( ipass == 0 ) then
     ipass = 1
     kfl_memor(1) = 0
     kfl_memor(2) = 0
     kfl_memor(3) = 0
  end if


  ithre=1

  select case(itask)

  case(-1_ip)
     !
     ! Automatically find a method
     !
!     call elsest_recomm(&
!          mnode,ndime,npoin,nelem,nnode,lnods,ltype,&
!          coord,recmeth)

  case(0_ip)
     !
     ! Allocate memory for structures
     !
     call elsest_alloca(0_ip,ipara)

  case(1_ip)
     !
     ! Preprocess
     !
     call elsest_binpre(&
          ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&         ! Bin strategy
          lnods,coord,rpara)


  case(2_ip)
     !
     ! Element search
     !
     lmini = -rpara(1)
     lmaxi = 1.0_rp+rpara(1)


     
     call elsest_binpro(&                                           ! Bin strategy
          imesh,lmesh,ipara,ithre,mnode,ndime,npoin,nelem,nnode,&
          lnods,ltopo,coord,xcoor,rpara,&
          ifoun,shapt,derit,coloc)

  case(3_ip)
     !
     ! Deallocate memory
     !
     call elsest_statis(2_ip,imesh,ipara,1)

     call elsest_bindea(ithre,imesh,lmesh)                          ! Bin strategy


  case(4_ip)
     !
     ! Deallocate memory for structures
     !
     call elsest_deallo()

  case(5_ip)
     !
     ! Check if test point is in element IFOUN
     !
     call elsest_inelem(&
          mnode,ndime,nnode,lnods,ltopo,coord,xcoor,rpara,&
          ifoun,shapt,derit,coloc)

  case(6_ip)
     !
     ! Check if test point is in element IFOUN: stop if not
     !
     ielem = ifoun
     call elsest_inelem(&
          mnode,ndime,nnode,lnods,ltopo,coord,xcoor,rpara,&
          ifoun,shapt,derit,coloc)
     if( ifoun == 0 ) then
        write(*,*) ielem,coloc(1:ndime)
        call runend('ELSEST: TEST POINT NOT IN ELEMENT')
     end if


  end select
  !
  ! Save last mesh number
  !
  lmesh=imesh

end subroutine elsest



subroutine elsest_binpro(&
     imesh,lmesh,ipara,ithre,mnode,ndime,npoin,nelem,nnode,&
     lnods,ltopo,coord,point_x,rpara,&
     ifoun,shapt,derit,coloc)
  !
  ! Bin search: look for host element
  !
  use def_elsest
!  use mod_elsest
  implicit none
  integer(ip), intent(in)  :: imesh,lmesh,ithre
  integer(ip), intent(in)  :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)  :: nnode
  integer(ip), intent(in)  :: lnods(mnode,nelem)
  integer(ip), intent(in)  :: ltopo
  real(rp),    intent(in)  :: coord(ndime,npoin),point_x(*),rpara(*)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: shapt(*),derit(*),coloc(*)
  integer(ip)              :: curr_box_coor(3),box_nr,ipara(*)
  integer(ip)              :: array_size,ielem,ii,inode,ipoin,ilook,kk
  integer(ip)              :: pnode,pelty,idime
  real(rp)                 :: time1,time2,time3
  !
  ! If not allocated, create structure of mesh IMESH
  !
  if( kfl_memor(1) == 0 ) call elsest_alloca(1_ip,ipara)
  if( bin_struc(imesh)%iallo == 0 ) then
     call elsest_binpre(&
          ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&
          lnods,coord,rpara)    
  end if

  if( imesh /= lmesh ) call elsest_binpoi(imesh)

  ksear(ithre) = ksear(ithre) + 1
  call elsest_cputim(time1)
  !
  ! Check if point is outside the bounding box
  !
  do idime = 1,ndime
     if( point_x(idime) < comin(idime) .or. point_x(idime) > comax(idime) ) then
        ifoun = -1
        return
     end if
  end do
  !
  ! Determine in which box (i,j,k) the point lies: curr_box_coor
  !
  call elsest_boxijk(ndime,nboxx,point_x,curr_box_coor,comin,comax)
  call elsest_boxnum(ndime,nboxx,curr_box_coor,box_nr)
  !
  ! Geometric CHECK
  !
  if( dataf == 0 ) then 
     array_size = size(bin_struc(imesh) % tboel(box_nr)%l,1)
  else if( dataf == 1 ) then
     array_size = bin_struc(imesh) % pboel(box_nr+1)-bin_struc(imesh)%pboel(box_nr)
     kk         = bin_struc(imesh) % pboel(box_nr) - 1
  end if

  !
  ! First try: Loop over elements in box
  !
  ii    = 0
  ifoun = 0

  do while( ifoun == 0 .and. ii < array_size )

     ii = ii + 1
     if( dataf == 0 ) then
        ielem = bin_struc(imesh) % tboel(box_nr) % l(ii)
     else if( dataf == 1 ) then
        ielem = bin_struc(imesh) % lboel(ii+kk)
     end if

     ilook = 1

     if( ilook == 1 ) then
        pnode = nnode
        do inode = 1,pnode
           ipoin= lnods(inode,ielem)
           do idime = 1,ndime
              bin_struc(imesh) % elcod(idime,inode,ithre) = coord(idime,ipoin)
           end do
        end do

        call elsest_chkelm(&
             ndime,ltopo,pnode,bin_struc(imesh)%elcod(:,:,ithre),shapt,&
             derit,point_x,coloc,ifoun,lmini,lmaxi)

        if( ifoun > 0 ) then
           ifoun = ielem
           kfirs(ithre) = kfirs(ithre)+1
        end if
     end if

  end do

  call elsest_cputim(time2)
  cputi(6,ithre)=cputi(6,ithre)+time2-time1

  call elsest_cputim(time3)
  cputi(7,ithre)=cputi(7,ithre)+time3-time2


end subroutine elsest_binpro


subroutine elsest_binpre(&
     ipara,imesh,ithre,mnode,ndime,npoin,nelem,nnode,&
     lnods,coord,rpara)
  !
  ! Construct the bin structure
  !
  use def_elsest
!  use mod_elsest
!  use def_master, only : kfl_paral
  implicit none
  integer(ip), intent(in)  :: ipara(*),imesh,ithre
  integer(ip), intent(in)  :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)  :: nnode
  integer(ip), intent(in)  :: lnods(mnode,*)
  real(rp),    intent(in)  :: coord(ndime,npoin),rpara(*)
  integer(ip)              :: idime,iboxe,istat,i,j,kk,jj,ii,ll
  integer(ip)              :: imin,imax,jmin,jmax,kmin,kmax,ni,nj,nk
  integer(ip)              :: ielem,box_nr,box_nr1,box_nr2,ninj
  real(rp),    pointer     :: xmima(:,:,:)
  real(rp)                 :: time1,time2,time3,time4
  real(rp)                 :: deltx,delty,deltz,dni,dnj,dnk

  call elsest_cputim(time1)
  bin_struc(imesh)%iallo = 1

  !----------------------------------------------------------------------
  !
  ! Initialize parameters
  !
  !----------------------------------------------------------------------

  !$OMP PARALLEL
  !$OMP SECTIONS
  !$OMP SECTION
  allocate(bin_struc(imesh)%cputi(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'CPUTI','elsest_alloc',0_ip)
  !$OMP SECTION
  allocate(bin_struc(imesh)%memor(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'MEMOR','elsest_alloc',0_ip)
  !$OMP SECTION
  allocate(bin_struc(imesh)%kstat(10,nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSTAT','elsest_alloc',0_ip)
  !$OMP SECTION
  allocate(bin_struc(imesh)%ksear(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSEAR','elsest_alloc',0_ip)
  !$OMP SECTION
  allocate(bin_struc(imesh)%kfirs(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KFIRS','elsest_alloc',0_ip)
  !$OMP SECTION
  allocate(bin_struc(imesh)%kseco(nthre),stat=istat)
  if(istat/=0) call elsest_memerr(0_ip,'KSECO','elsest_alloc',0_ip)
  !$OMP END SECTIONS
  !$OMP END PARALLEL

  !----------------------------------------------------------------------
  !
  ! Point to current mesh (IMESH) structure
  !
  !----------------------------------------------------------------------

  call elsest_binpoi(imesh)
  !$OMP PARALLEL DO PRIVATE(i,j)
  do i=1,nthre
     do j=1,10
        memor(j,i) = 0_ip
        cputi(j,i) = 0.0_rp
        kstat(j,i) = 0_ip
     end do
     kstat(1,i)    = 1e6_rp
     kstat(3,i)    = 1e6_rp
     ksear(i)      = 0_ip
     kfirs(i)      = 0_ip
     kseco(i)      = 0_ip
  end do
  !$OMP END PARALLEL DO

  memax      = 0_ip
  nboxx(1)   = ipara(1)
  nboxx(2)   = ipara(2)
  nboxx(3)   = ipara(3)  
  dataf      = ipara(4)
  comin(1)   = 0.0_rp
  comin(2)   = 0.0_rp
  comin(3)   = 0.0_rp
  comax(1)   = 0.0_rp
  comax(2)   = 0.0_rp
  comax(3)   = 0.0_rp

  !----------------------------------------------------------------------
  !
  ! Total number of boxes
  !
  !----------------------------------------------------------------------

  if( ndime == 1 ) then
     nboxe = nboxx(1)
  else if( ndime == 2 ) then
     nboxe = nboxx(1) * nboxx(2)
  else
     nboxe = nboxx(1) * nboxx(2) * nboxx(3)
  end if
  ni   = nboxx(1)
  nj   = nboxx(2)
  nk   = nboxx(3)
  bin_struc(imesh) % nboxx(1) = ni
  bin_struc(imesh) % nboxx(2) = nj
  bin_struc(imesh) % nboxx(3) = nk

  !----------------------------------------------------------------------
  !
  ! Allocate memory for bin structure
  !
  !----------------------------------------------------------------------

  allocate(bin_struc(imesh)%elcod(ndime,mnode,nthre),stat=istat)
  call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'ELCOD','elsest_binpre',bin_struc(imesh)%elcod)
  if( dataf == 0 ) then
     allocate(bin_struc(imesh)%tboel(nboxe),stat=istat)
     call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'TBOEL','elsest_binpre',bin_struc(imesh)%tboel)
  else if( dataf == 1 ) then
     allocate(bin_struc(imesh)%pboel(nboxe+1),stat=istat)
     call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'PBOEL','elsest_binpre',bin_struc(imesh)%pboel)
  end if

  elcod => bin_struc(imesh) % elcod
  tboel => bin_struc(imesh) % tboel
  pboel => bin_struc(imesh) % pboel
  lboel => bin_struc(imesh) % lboel


  !----------------------------------------------------------------------
  !
  ! Compute bounding box and bin spacing delta
  !
  !----------------------------------------------------------------------

  call elsest_boubox(ndime,npoin,coord,comin,comax)
  !
  ! Bin spacing DELTA
  !
  do idime = 1,ndime
     delta(idime) = (comax(idime)-comin(idime))/real(nboxx(idime))
  end do
  deltx = 1.0_rp / ( comax(1)-comin(1) )
  delty = 1.0_rp / ( comax(2)-comin(2) )
  dni   = real(ni,rp) * deltx
  dnj   = real(nj,rp) * delty
  if( ndime == 3 ) then
     deltz = 1.0_rp / ( comax(3)-comin(3) )
     dnk   = real(nk,rp) * deltz
  end if
  ninj  = nboxx(1) * nboxx(2)

  call elsest_cputim(time2)
  cputi(1,ithre)=time2-time1

  !----------------------------------------------------------------------
  !
  ! Compute element bounding boxes
  !
  !----------------------------------------------------------------------

  allocate( xmima(3,2,nelem), stat = istat )
  call elsest_elbbox(&
       mnode,ndime,npoin,nelem,nnode,lnods,coord,xmima)

  !----------------------------------------------------------------------
  !
  ! Number of elements per box
  !
  !----------------------------------------------------------------------

  allocate(nbono(nboxe),stat=istat)
  call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'NBONO','elsest_binpre',nbono)

  if( ndime == 1 ) then

     do ielem = 1,nelem

        imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * ni , ip ) + 1
        imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * ni , ip ) + 1  

        imin = max(imin,1_ip)
        imax = min(imax,nboxx(1))

        do ii = imin,imax
           box_nr = ii
           nbono(box_nr) = nbono(box_nr) + 1
        end do
     end do

  else if( ndime == 2 ) then

     do ielem = 1,nelem

        imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * ni , ip ) + 1
        imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * ni , ip ) + 1      

        jmin = int(( ( xmima(2,1,ielem) - comin(2)) * delty ) * nj , ip ) + 1
        jmax = int(( ( xmima(2,2,ielem) - comin(2)) * delty ) * nj , ip ) + 1      

        imin = max(imin,1_ip)
        imax = min(imax,nboxx(1))
        jmin = max(jmin,1_ip)
        jmax = min(jmax,nboxx(2))

        do ii = imin,imax
           do jj = jmin,jmax
              box_nr = (jj-1_ip) * ni + ii
              nbono(box_nr) = nbono(box_nr) + 1
           end do
        end do
     end do

  else if( ndime == 3 ) then

     do ielem = 1,nelem

        imin = int( ( xmima(1,1,ielem) - comin(1)) * dni , ip )  + 1
        imax = int( ( xmima(1,2,ielem) - comin(1)) * dni , ip )  + 1      
                                                                 
        jmin = int( ( xmima(2,1,ielem) - comin(2)) * dnj , ip )  + 1
        jmax = int( ( xmima(2,2,ielem) - comin(2)) * dnj , ip )  + 1      
                                                                 
        kmin = int( ( xmima(3,1,ielem) - comin(3)) * dnk , ip )  + 1
        kmax = int( ( xmima(3,2,ielem) - comin(3)) * dnk , ip )  + 1      

        imin = max(imin,1_ip)
        imax = min(imax,nboxx(1))
        jmin = max(jmin,1_ip)
        jmax = min(jmax,nboxx(2))
        kmin = max(kmin,1_ip)
        kmax = min(kmax,nboxx(3))

        do kk = kmin,kmax
           box_nr2 = ninj * (kk-1)
           do jj = jmin,jmax
              box_nr1 = box_nr2 + nboxx(1) * (jj-1)
              do ii = imin,imax
                 box_nr        = box_nr1 + ii
                 nbono(box_nr) = nbono(box_nr) + 1
              end do
           end do
        end do

     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Fill in box element list
  !
  !----------------------------------------------------------------------

  if( dataf == 0 ) then
     !
     ! Type
     !
     do iboxe = 1,bin_struc(imesh) % nboxe
        allocate( bin_struc(imesh) % tboel(iboxe)%l(nbono(iboxe)) , stat = istat )
        nbono(iboxe) = 0
     end do

     if( ndime == 1 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * ni , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * ni , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))

           do ii = imin,imax
              box_nr        = ii
              nbono(box_nr) = nbono(box_nr) + 1
              bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
           end do

        end do

     else if( ndime == 2 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * ni , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * ni , ip ) + 1      

           jmin = int(( ( xmima(2,1,ielem) - comin(2)) * delty ) * nj , ip ) + 1
           jmax = int(( ( xmima(2,2,ielem) - comin(2)) * delty ) * nj , ip ) + 1    

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))

           do ii = imin,imax
              do jj = jmin,jmax
                 box_nr        = (jj-1_ip) * ni + ii
                 nbono(box_nr) = nbono(box_nr) + 1
                 bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
              end do

           end do

        end do

     else

        do ielem = 1,nelem

           imin = int( ( xmima(1,1,ielem) - comin(1) ) * dni , ip ) + 1
           imax = int( ( xmima(1,2,ielem) - comin(1) ) * dni , ip ) + 1      

           jmin = int( ( xmima(2,1,ielem) - comin(2) ) * dnj , ip ) + 1
           jmax = int( ( xmima(2,2,ielem) - comin(2) ) * dnj , ip ) + 1      

           kmin = int( ( xmima(3,1,ielem) - comin(3) ) * dnk , ip ) + 1
           kmax = int( ( xmima(3,2,ielem) - comin(3) ) * dnk , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))
           kmin = max(kmin,1_ip)
           kmax = min(kmax,nboxx(3))

           do kk = kmin,kmax
              box_nr2 = ninj * (kk-1)
              do jj = jmin,jmax
                 box_nr1 = box_nr2 + nboxx(1) * (jj-1)
                 do ii = imin,imax
                    box_nr        = box_nr1 + ii
                    nbono(box_nr) = nbono(box_nr) + 1
                    bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
                 end do
              end do

           end do

        end do

     end if

  else
     !
     ! Linked list
     !
     bin_struc(imesh) % pboel(1) = 1

     do iboxe = 1,bin_struc(imesh) % nboxe
        bin_struc(imesh) % pboel(iboxe+1) = bin_struc(imesh) % pboel(iboxe) + nbono(iboxe) 
        nbono(iboxe) = 0
     end do

     allocate(bin_struc(imesh)%lboel(bin_struc(imesh) % pboel(bin_struc(imesh) % nboxe+1)),stat=istat)
     call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'LBOEL','elsest_binlis',bin_struc(imesh)%lboel)
     lboel => bin_struc(imesh) % lboel

     if( ndime == 1 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * ni , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * ni , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))

           do ii = imin,imax
              box_nr        = ii
              ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
              nbono(box_nr) = nbono(box_nr) + 1
              lboel(ll)     = ielem
           end do
        end do

     else if( ndime == 2 ) then

        do ielem = 1,nelem

           imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * ni , ip ) + 1
           imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * ni , ip ) + 1      

           jmin = int(( ( xmima(2,1,ielem) - comin(2)) * delty ) * nj , ip ) + 1
           jmax = int(( ( xmima(2,2,ielem) - comin(2)) * delty ) * nj , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))

           do ii = imin,imax
              do jj = jmin,jmax
                 box_nr        = (jj-1_ip) * ni + ii
                 ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
                 nbono(box_nr) = nbono(box_nr) + 1
                 lboel(ll)     = ielem
              end do
           end do
        end do

     else if( ndime == 3 ) then

        do ielem = 1,nelem

           imin = int( ( xmima(1,1,ielem) - comin(1)) * dni , ip ) + 1
           imax = int( ( xmima(1,2,ielem) - comin(1)) * dni , ip ) + 1      

           jmin = int( ( xmima(2,1,ielem) - comin(2)) * dnj , ip ) + 1
           jmax = int( ( xmima(2,2,ielem) - comin(2)) * dnj , ip ) + 1      

           kmin = int( ( xmima(3,1,ielem) - comin(3)) * dnk , ip ) + 1
           kmax = int( ( xmima(3,2,ielem) - comin(3)) * dnk , ip ) + 1      

           imin = max(imin,1_ip)
           imax = min(imax,nboxx(1))
           jmin = max(jmin,1_ip)
           jmax = min(jmax,nboxx(2))
           kmin = max(kmin,1_ip)
           kmax = min(kmax,nboxx(3))

           do kk = kmin,kmax
              box_nr2 = ninj * (kk-1)
              do jj = jmin,jmax
                 box_nr1 = box_nr2 + nboxx(1) * (jj-1)
                 do ii = imin,imax
                    box_nr        = box_nr1 + ii
                    ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
                    nbono(box_nr) = nbono(box_nr) + 1
                    lboel(ll)     = ielem
                 end do
              end do
           end do
        end do

     end if

  end if

  !
  ! Deallocate memory
  !
  call elsest_cputim(time3)
  deallocate( xmima, stat = istat )
  call elsest_memchk(2_ip,ithre,istat,memor(1,ithre),'NBONO','elsest_binpre',nbono)
  deallocate(nbono,stat=istat)
  if(istat/=0) call elsest_memerr(2_ip,'NBONO','elsest_binpre',0_ip)
  call elsest_cputim(time4)
  cputi(3,ithre)=time4-time3

  !
  ! Output statistics
  !
  if( ipara(7) /= 0 ) call elsest_statis(1_ip,imesh,ipara,ithre)

  call elsest_binpoi(imesh)

end subroutine elsest_binpre




subroutine elsest_memerr(itask,vanam,vacal,istat)
  !-----------------------------------------------------------------------
  !****f* elsest_memerr
  ! NAME
  !    elsest_memerr
  ! DESCRIPTION
  !    This routine ends Alya when an error has been found
  !    allocating or deallocating memory.
  ! USES
  ! USED BY
  !    mod_elsest
  !***
  !-----------------------------------------------------------------------
  use def_elsest, only       : elsest_intost
  implicit none
  integer(ip),   intent(in) :: itask,istat
  character*(*), intent(in) :: vanam,vacal

  if(itask==0) then
     call elsest_runend(&
          trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE ALLOCATED.'&
          //' RUN TIME ERROR: '//elsest_intost(istat))
  else if(itask==1) then
     call elsest_runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE REALLOCATED')
  else
     call elsest_runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE DEALLOCATED')
  end if

end subroutine elsest_memerr

subroutine elsest_binpoi(imesh)
  !
  ! Points to current mesh structure
  !
  use def_elsest
  implicit none
  integer(ip), intent(in) :: imesh

  nboxe => bin_struc(imesh)%nboxe
  nboxx => bin_struc(imesh)%nboxx
  dataf => bin_struc(imesh)%dataf
  iallo => bin_struc(imesh)%iallo
  lboel => bin_struc(imesh)%lboel
  pboel => bin_struc(imesh)%pboel
  tboel => bin_struc(imesh)%tboel
  memor => bin_struc(imesh)%memor
  memax => bin_struc(imesh)%memax
  kstat => bin_struc(imesh)%kstat
  ksear => bin_struc(imesh)%ksear
  kfirs => bin_struc(imesh)%kfirs
  kseco => bin_struc(imesh)%kseco
  comin => bin_struc(imesh)%comin
  comax => bin_struc(imesh)%comax
  delta => bin_struc(imesh)%delta
  elcod => bin_struc(imesh)%elcod
  cputi => bin_struc(imesh)%cputi

end subroutine elsest_binpoi



subroutine elsest_runend(message)
  !-----------------------------------------------------------------------
  !
  ! This routine stops the run and writes the summary of CPU time.
  !
  !-----------------------------------------------------------------------
  use def_elsest
  implicit none
  character(*) :: message
  !
  ! Write message and stop the run
  !
  if(iunit(1)/=0) then
!$OMP CRITICAL(write)
     write(iunit(1),1) 'ELSEST ERROR WAS FOUND: '//trim(message)
!$OMP END CRITICAL(write)
  else
!$OMP CRITICAL(write)
     write(6,1)     'ELSEST ERROR WAS FOUND: '//trim(message)
!$OMP END CRITICAL(write)
  end if

  stop

1 format(//,5x,'>>> ',a)

end subroutine elsest_runend

subroutine elsest_alloca(itask,ipara)
  !
  ! Allocate memory for bin and/or quad/oct structures
  !
  use def_elsest
!  use mod_elsest
  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip), intent(inout) :: ipara(*)
  integer(ip)                :: nmesh,jmesh,istat


  nthre=1

  !
  ! Allocate bin structure
  !
  if( (itask==0 .or. itask==1 ) .and. kfl_memor(1)==0 ) then
     nmesh        = ipara(5)
     kfl_memor(1) = 1
     allocate(bin_struc(nmesh),stat=istat)
     if(istat/=0) call elsest_memerr(0_ip,'BIN_STRUC','elsest_alloc',0_ip)
     do jmesh=1,nmesh
        bin_struc(jmesh)%iallo=0
     end do
  end if
  !
  ! Allocate oct structure
  !
  if( (itask==0 .or. itask==2 ) .and. kfl_memor(2)==0 ) then
     nmesh        = ipara(5)
     kfl_memor(2) = 1
     allocate(oct_struc(nmesh),stat=istat)
     if(istat/=0) call elsest_memerr(0_ip,'OCT_STRUC','elsest_alloc',0_ip)
     do jmesh=1,nmesh
        oct_struc(jmesh)%iallo=0
     end do
  end if
  !
  ! Current
  !
  if( kfl_memor(3)==0 ) then
     kfl_memor(3) = 1
     allocate(current(nthre),stat=istat)
  end if

end subroutine elsest_alloca


subroutine elsest_deallo()
  !
  ! Deallocate memory for bin and quad/oct structures
  !
  use def_elsest
!  use mod_elsest
  implicit none
  integer(ip) :: istat
  !
  ! Deallocate Bin structure
  !
  if( kfl_memor(1)==1 ) then
     kfl_memor(1)=0
     deallocate(bin_struc,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'BIN_STRUC','elsest_deallo',0_ip)
  end if
  !
  ! Deallocate Quad/Oct structure
  !
  if( kfl_memor(2)==1 ) then
     kfl_memor(2)=0
     deallocate(oct_struc,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'OCT_STRUC','elsest_deallo',0_ip)
  end if
  !
  ! Cuurent
  !
  if( kfl_memor(3)==1 ) then
     kfl_memor(3)=0
     deallocate(current,stat=istat)
  end if


end subroutine elsest_deallo

subroutine elsest_bindea(ithre,imesh,lmesh)
  !
  ! Bin deallocation
  !
  use def_elsest
!  use mod_elsest
  implicit none
  integer(ip), intent(in) :: ithre,imesh,lmesh
  integer(ip)             :: istat

  if(imesh/=lmesh) call elsest_binpoi(imesh)
  !
  ! It has not been allocated
  !
  if(  bin_struc(imesh)%iallo == 0 ) return

  bin_struc(imesh)%iallo=0
  if(bin_struc(imesh)%dataf==0) then
     call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'TBOEL','elsest_bindea',bin_struc(imesh)%tboel)
     deallocate(bin_struc(imesh)%tboel,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'TBOEL','elsest_bindea',0_ip)
  else if(bin_struc(imesh)%dataf==1) then
     call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'PBOEL','elsest_bindea',bin_struc(imesh)%pboel)
     deallocate(bin_struc(imesh)%pboel,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'PBOEL','elsest_bindea',0_ip)
  end if
  call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'ELCOD','elsest_bindea',bin_struc(imesh)%elcod)
  deallocate(bin_struc(imesh)%elcod,stat=istat)
  if(istat/=0) call elsest_memerr(2_ip,'ELCOD','elsest_bindea',0_ip)
  !
  ! Deallocate arrays for memory and cpu statistics
  !
  if(associated(bin_struc(imesh)%cputi)) then
     deallocate(bin_struc(imesh)%cputi,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'CPUTI','elsest_deallo',0_ip)
  end if

  if(associated(bin_struc(imesh)%memor)) then
     deallocate(bin_struc(imesh)%memor,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'MEMOR','elsest_deallo',0_ip)
  end if

  if(associated(bin_struc(imesh)%kstat)) then
     deallocate(bin_struc(imesh)%kstat,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSTAT','elsest_alloc',0_ip)
  end if

  if(associated(bin_struc(imesh)%ksear)) then
     deallocate(bin_struc(imesh)%ksear,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSEAR','elsest_alloc',0_ip)
  end if

  if(associated(bin_struc(imesh)%kfirs)) then
     deallocate(bin_struc(imesh)%kfirs,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KFIRS','elsest_alloc',0_ip)
  end if

  if(associated(bin_struc(imesh)%kseco)) then
     deallocate(bin_struc(imesh)%kseco,stat=istat)
     if(istat/=0) call elsest_memerr(2_ip,'KSECO','elsest_alloc',0_ip)
  end if

end subroutine elsest_bindea


subroutine elsest_boxnum(ndime,nboxx,box_coord,box_nr)
  !
  ! Returns the box number box_nr of a given box with coordinates in i j k
  !
!  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: nboxx(ndime),box_coord(ndime)
  integer(ip), intent(out) :: box_nr

  if(ndime==2) then
     box_nr=(box_coord(2)-1)*nboxx(1) + box_coord(1)
  else
     box_nr=(box_coord(3)-1)*(nboxx(1)*nboxx(2)) + (box_coord(2)-1)*nboxx(1) + box_coord(1)
  end if

end subroutine elsest_boxnum


subroutine elsest_boxijk(ndime,nboxx,point_x,curr_box_coor,comin,comax)
  !
  ! Look for the box curr_box_coor)(i,j,k) containing point point_x
  !
!  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: nboxx(ndime)
  real(rp),    intent(in)  :: point_x(ndime),comin(ndime),comax(ndime)
  integer(ip), intent(out) :: curr_box_coor(ndime)
  integer(ip)              :: idime

  do idime=1,ndime
     curr_box_coor(idime) = int((  (point_x(idime)-comin(idime))&
          &                      /(comax(idime)-comin(idime)) )*nboxx(idime))+1
  end do

end subroutine elsest_boxijk


subroutine elsest_boubox(ndime,npoin,coord,comin,comax)
  !
  ! Compute the bounding box of the domain
  !
 
  implicit none
  integer(ip), intent(in)  :: ndime,npoin
  real(rp),    intent(in)  :: coord(ndime,npoin)
  real(rp),    intent(out) :: comin(ndime),comax(ndime)
  integer(ip)              :: idime,ipoin

  do idime=1,ndime
     comax(idime)=-1e30
     comin(idime)= 1e30
  end do
  do ipoin=1,npoin
     do idime=1,ndime
        if(coord(idime,ipoin)>comax(idime)) comax(idime)=coord(idime,ipoin)
        if(coord(idime,ipoin)<comin(idime)) comin(idime)=coord(idime,ipoin)
     end do
  end do
  do idime=1,ndime
     comax(idime)= comax(idime)+0.000001_rp*(comax(idime)-comin(idime))
     comin(idime)= comin(idime)-0.000001_rp*(comax(idime)-comin(idime))
  end do

end subroutine elsest_boubox


subroutine elsest_cputim(rtime)
  !-----------------------------------------------------------------------
  !****f* elsest_cputim
  ! NAME
  !    nsi_elmope
  ! DESCRIPTION
  !    Returns the CPU time in seconds
  ! OUTPUT
  !    rtime
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------


  implicit none
  real(rp), intent(out) :: rtime
  real(4)               :: rtim4,elap4(2),etime


  rtim4 = etime(elap4)
  rtime = real(rtim4)

end subroutine elsest_cputim



subroutine elsest_elbbox(&
     mnode,ndime,npoin,nelem,nnode,lnods,coord,xmima)
  !
  ! Construct the bin structure
  !
  use def_elsest
  implicit none
  integer(ip), intent(in)  :: mnode,ndime,npoin,nelem
  integer(ip), intent(in)  :: nnode
  integer(ip), intent(in)  :: lnods(mnode,*)
  real(rp),    intent(in)  :: coord(ndime,npoin)
  real(rp),    intent(out) :: xmima(3,2,nelem)
  integer(ip)              :: ielem,inode,ipoin
  real(rp)                 :: zeror

  zeror = epsilon(1.0_rp)

  if( ndime == 1 ) then

     do ielem = 1,nelem
        xmima(1,1,ielem) = coord(1,lnods(1,ielem))
        xmima(1,2,ielem) = coord(1,lnods(1,ielem))
        do inode = 2,nnode
           ipoin            = lnods(inode,ielem)
           xmima(1,1,ielem) = min( xmima(1,1,ielem) , coord(1,ipoin) )
           xmima(1,2,ielem) = max( xmima(1,2,ielem) , coord(1,ipoin) )
        end do
        xmima(1,1,ielem) = xmima(1,1,ielem) - zeror
        xmima(1,2,ielem) = xmima(1,2,ielem) + zeror
     end do

  else if( ndime == 2 ) then

     do ielem = 1,nelem
        xmima(1,1,ielem) = coord(1,lnods(1,ielem)) 
        xmima(2,1,ielem) = coord(2,lnods(1,ielem))
        xmima(1,2,ielem) = coord(1,lnods(1,ielem))
        xmima(2,2,ielem) = coord(2,lnods(1,ielem))
        do inode = 2,nnode
           ipoin            = lnods(inode,ielem)
           xmima(1,1,ielem) = min( xmima(1,1,ielem) , coord(1,ipoin) )
           xmima(2,1,ielem) = min( xmima(2,1,ielem) , coord(2,ipoin) )
           xmima(1,2,ielem) = max( xmima(1,2,ielem) , coord(1,ipoin) )
           xmima(2,2,ielem) = max( xmima(2,2,ielem) , coord(2,ipoin) )
        end do
        xmima(1,1,ielem) = xmima(1,1,ielem) - zeror
        xmima(1,2,ielem) = xmima(1,2,ielem) + zeror
        xmima(2,1,ielem) = xmima(2,1,ielem) - zeror
        xmima(2,2,ielem) = xmima(2,2,ielem) + zeror
     end do

  else

     do ielem = 1,nelem
        xmima(1,1,ielem) = coord(1,lnods(1,ielem)) 
        xmima(2,1,ielem) = coord(2,lnods(1,ielem))  
        xmima(3,1,ielem) = coord(3,lnods(1,ielem))  
        xmima(1,2,ielem) = coord(1,lnods(1,ielem)) 
        xmima(2,2,ielem) = coord(2,lnods(1,ielem)) 
        xmima(3,2,ielem) = coord(3,lnods(1,ielem)) 
        do inode = 2,nnode
           ipoin            = lnods(inode,ielem)
           xmima(1,1,ielem) = min( xmima(1,1,ielem) , coord(1,ipoin) )
           xmima(2,1,ielem) = min( xmima(2,1,ielem) , coord(2,ipoin) )
           xmima(3,1,ielem) = min( xmima(3,1,ielem) , coord(3,ipoin) )
           xmima(1,2,ielem) = max( xmima(1,2,ielem) , coord(1,ipoin) )
           xmima(2,2,ielem) = max( xmima(2,2,ielem) , coord(2,ipoin) )
           xmima(3,2,ielem) = max( xmima(3,2,ielem) , coord(3,ipoin) )
        end do
        xmima(1,1,ielem) = xmima(1,1,ielem) - zeror
        xmima(1,2,ielem) = xmima(1,2,ielem) + zeror
        xmima(2,1,ielem) = xmima(2,1,ielem) - zeror
        xmima(2,2,ielem) = xmima(2,2,ielem) + zeror
        xmima(3,1,ielem) = xmima(3,1,ielem) - zeror
        xmima(3,2,ielem) = xmima(3,2,ielem) + zeror
     end do
    
  end if
end subroutine elsest_elbbox




subroutine elsest_statis(itask,imesh,ipara,ithre)
  !------------------------------------------------------------------------
  !****f* elsest/elsest_statis
  ! NAME 
  !    elsest_statis
  ! DESCRIPTION
  !    Output statistics
  ! USES
  ! USED BY
  !    elsest
  !***
  !------------------------------------------------------------------------
  use def_elsest
  implicit none
  integer(ip), intent(in) :: itask,imesh,ipara(*),ithre
  real(rp)                :: cputo,cputo_master
  integer(ip)             :: tomem,ii,dummi,jj
  real(rp)                :: rbyte,rbyt2,rbyt3,dummr
  character(2)            :: lbyte,lbyt2,lbyt3

  if(iunit(1)/=0) then

     select case(itask)

     case(1)
        !
        ! CPU time
        !
        if(ipara(8)==0) then
           !$OMP CRITICAL(write)
           write(iunit(1),1) imesh,'(BIN STRATEGY)'
           cputo=cputi(1,ithre)+cputi(2,ithre)+cputi(3,ithre)+cputi(4,ithre)
           !$OMP END CRITICAL(write)
           if(cputo==0.0_rp) cputo=1.0_rp
           !$OMP CRITICAL(write)
           write(iunit(1),100) cputo,&
                cputi(1,ithre),100.0_rp*cputi(1,ithre)/cputo,&
                cputi(2,ithre),100.0_rp*cputi(2,ithre)/cputo,&
                cputi(3,ithre),100.0_rp*cputi(3,ithre)/cputo,&
                cputi(4,ithre),100.0_rp*cputi(4,ithre)/cputo
           !$OMP END CRITICAL(write)
        else
           !$OMP CRITICAL(write)
           write(iunit(1),1) imesh,'(QUAD/OCT TREE STRATEGY)'
           !$OMP END CRITICAL(write)
           cputo=cputi(1,ithre)+cputi(2,ithre)+cputi(3,ithre)+cputi(4,ithre)
           if(cputo==0.0_rp) cputo=1.0_rp
           !$OMP CRITICAL(write)
           write(iunit(1),101) cputo,&
                cputi(1,ithre),100.0_rp*cputi(1,ithre)/cputo,&
                cputi(2,ithre),100.0_rp*cputi(2,ithre)/cputo,&
                cputi(3,ithre),100.0_rp*cputi(3,ithre)/cputo,&
                cputi(4,ithre),100.0_rp*cputi(4,ithre)/cputo
           !$OMP END CRITICAL(write)
        end if
        !
        ! Memory
        !
        tomem=0
        if(nthre>1) then
           do ii=2,5
              do jj=1,nthre
                 tomem=tomem+memor(ii,jj)
              end do
           end do
           do jj=2,nthre
              memor(2,1)=memor(2,1)+memor(2,jj)
              memor(3,1)=memor(3,1)+memor(3,jj)
              memor(5,1)=memor(5,1)+memor(5,jj)
           end do
        else
           do ii=2,5
              tomem=tomem+memor(ii,ithre)
           end do
        end if
        if(tomem==0) tomem=1

        call elsest_inbyte(tomem,rbyte,lbyte)
        call elsest_inbyte(memor(5,ithre),rbyt2,lbyt2)
        call elsest_inbyte(memax,rbyt3,lbyt3)

        if(ipara(8)==0) then
           !$OMP CRITICAL(write)
           write(iunit(1),200) real(tomem)/rbyte,lbyte,&
                real(memor(2,ithre))/rbyte,lbyte,100.0_rp*real(memor(2,ithre))/real(tomem),&
                real(memor(3,ithre))/rbyte,lbyte,100.0_rp*real(memor(3,ithre))/real(tomem),&
                real(memax)/rbyt3,lbyt3
           !$OMP END CRITICAL(write)
        else if(ipara(8)==1) then
           !$OMP CRITICAL(write)
           write(iunit(1),201) real(tomem)/rbyte,lbyte,&
                real(memor(2,ithre))/rbyte,lbyte,100.0_rp*real(memor(2,ithre))/real(tomem),&
                real(memor(3,ithre))/rbyte,lbyte,100.0_rp*real(memor(3,ithre))/real(tomem),&
                real(memax)/rbyt3,lbyt3
           !$OMP END CRITICAL(write)
        end if
        !
        ! Structure statistics
        !
        if(ipara(8)==0) then
           if(kstat(6,ithre)/=0) then
              dummr=real(kstat(6,ithre))
           else
              dummr=1.0_rp
           end if
           !$OMP CRITICAL(write)
           write(iunit(1),300)&
                kstat(6,ithre),kstat(5,ithre),100.0_rp*real(kstat(5,ithre))/dummr,kstat(1,ithre),kstat(2,ithre),&
                int(real(kstat(8,ithre))/dummr),kstat(3,ithre),kstat(4,ithre),int(real(kstat(7,ithre))/dummr)
           !$OMP END CRITICAL(write)
        else if(ipara(8)==1) then
           if(kstat(6,ithre)/=0) then
              dummr=real(kstat(6,ithre))
           else
              dummr=1.0_rp
           end if
           !$OMP CRITICAL(write)
           write(iunit(1),301)&
                kstat(6,ithre),kstat(1,ithre),kstat(2,ithre),int(real(kstat(8,ithre))/dummr),&
                kstat(3,ithre),kstat(4,ithre),int(real(kstat(7,ithre))/dummr)
           !$OMP END CRITICAL(write)
        end if

     case(2)
        !
        ! Process
        !
        if(ipara(8)==0) then
           !$OMP CRITICAL(write)
           write(iunit(1),2) imesh,'(BIN STRATEGY)',ithre
           !$OMP END CRITICAL(write)
        else if(ipara(8)==1) then
           !$OMP CRITICAL(write)
           write(iunit(1),2) imesh,'(QUAD/OCT TREE STRATEGY)'
           !$OMP END CRITICAL(write)
        end if

        !
        ! USING OMP
        !
        if(ksear(ithre)/=0) then
           if (nthre>1) then
              dummr=0.
              do ii=1,nthre
                 dummr=dummr+real(ksear(ii))
              end do
           else
              dummr=real(ksear(ithre))
           end if
        else
           dummr=1.0_rp
        end if
        if(nthre>1) then
           cputo=0.
           cputo_master=cputi(6,ithre)+cputi(7,ithre)
           do ii=1,nthre
              cputo=cputo+cputi(6,ii)+cputi(7,ii)
           end do
           do ii=2,nthre
              ksear(1)   = ksear(1)+ksear(ii)
              kfirs(1)   = kfirs(1)+kfirs(ii)
              kseco(1)   = kseco(1)+kseco(ii)
              cputi(1,1) = cputi(1,1)+cputi(1,ii)
              cputi(2,1) = cputi(2,1)+cputi(2,ii)
              cputi(3,1) = cputi(3,1)+cputi(3,ii)
              cputi(4,1) = cputi(4,1)+cputi(4,ii)
              cputi(5,1) = cputi(5,1)+cputi(5,ii)
              cputi(6,1) = cputi(6,1)+cputi(6,ii)
              cputi(7,1) = cputi(7,1)+cputi(7,ii)
           end do
        else
           cputo=cputi(6,ithre)+cputi(7,ithre)
           cputo_master=cputi(6,ithre)+cputi(7,ithre)
        end if
        if(cputo==0.0_rp) cputo=1.0_rp
        dummi=ksear(ithre)-kfirs(ithre)-kseco(ithre)
        !$OMP CRITICAL(write)
        write(iunit(1),402)&
             ksear(ithre),&
             nthre,&
             kfirs(ithre),&
             100.0_rp*real(kfirs(ithre))/real(ksear(ithre)),&
             kseco(ithre),100.0_rp*real(kseco(ithre))/real(ksear(ithre)),&
             dummi,100.0_rp*real(dummi)/real(ksear(ithre)),&
             cputo,&
             cputo_master,&
             cputi(6,ithre),100.0_rp*cputi(6,ithre)/cputo,&
             cputi(7,ithre),100.0_rp*cputi(7,ithre)/cputo
        !$OMP END CRITICAL(write)
     end select

  end if

1 format(/,&
       & 5x,'   ----------------','---',/&
       & 5x,'|- PREPROCESS MESH ',i3,' ',a,/,&
       & 5x,'   ----------------','---')
2 format(/,&
       & 5x,'   ----------------','---',/&
       & 5x,'|- PROCESS MESH    ',i3,' ',a,/,&
       & 5x,'   ----------------','---')
100 format(/,&
       & 5x,'   |- SUMMARY OF COMPUTING TIMES  ',//,&
       & 5x,'      TOTAL CPUT TIME:            ',f11.2,3x,/,&
       & 5x,'      INTIALIZATION:              ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      FILL BIN WITH NODES:        ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      FILL BIN WITH ELEMENTS:     ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,3x,' (',f6.2,' % )')
101 format(/,&
       & 5x,'   |- SUMMARY OF COMPUTING TIMES  ',//,&
       & 5x,'      TOTAL CPUT TIME:            ',f11.2,3x,/,&
       & 5x,'      INITIALIZATION:             ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      QUAD/OCT STRUCTURE:         ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      FILL BIN WITH ELEMENTS:     ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,3x,' (',f6.2,' % )')
200 format(/,&
       & 5x,'   |- SUMMARY OF MEMORY           ',//,&
       & 5x,'      TOTAL MEMORY:               ',f11.2,1x,a2,/,&
       & 5x,'      BIN STRUCTURE:              ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,/,&
       & 5x,'      MAXIMUM MEMORY REQUIRED:    ',f11.2,1x,a2)
201 format(/,&
       & 5x,'   |- SUMMARY OF MEMORY           ',//,&
       & 5x,'      TOTAL MEMORY:               ',f11.2,1x,a2,/,&
       & 5x,'      QUAD/OCT STRUCTURE:         ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,/,&
       & 5x,'      MAXIMUM MEMORY REQUIRED:    ',f11.2,1x,a2)
300 format(/,&
       & 5x,'   |- SUMMARY OF BIN STRUCTURE    ',//,&
       & 5x,'      # OF BINS:                  ',i11,/,&
       & 5x,'      # OF BINS WITH ELEMENTS:    ',i11,3x,' (',f6.2,' % )',/,&
       & 5x,'      MIN # NODES IN A BIN:       ',i11,/,&
       & 5x,'      MAX # NODES IN A BIN:       ',i11,/,&
       & 5x,'      AVERAGE # NODES IN BIN:     ',i11,/,&
       & 5x,'      MIN # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      MAX # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      AVERAGE # ELEMENTS IN BIN:  ',i11)
301 format(/,&
       & 5x,'   |- SUMMARY OF QUAD/OCT STRUCTURE',//,&
       & 5x,'      # OF BINS WITH ELEMENTS:    ',i11,/,&
       & 5x,'      MIN # NODES IN A BIN:       ',i11,/,&
       & 5x,'      MAX # NODES IN A BIN:       ',i11,/,&
       & 5x,'      AVERAGE # NODES IN BIN:     ',i11,/,&
       & 5x,'      MIN # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      MAX # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      AVERAGE # ELEMENTS IN BIN:  ',i11)
400 format(/,&
       & 5x,'   |- ELEMENT SEARCH:             ',//,&
       & 5x,'      # OF SEARCHES:              ',i11,/,&
       & 5x,'      # ELEMS FOUND 1ST STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS FOUND 2ND STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS NOT FOUND:          ',i11  ,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      TOTAL CPU TIME :            ',f11.2,3x,                      /,&
       & 5x,'      CPU TIME 1ST STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      CPU TIME 2ND STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )')
402 format(/,&
       & 5x,'   |- ELEMENT SEARCH:             ',//,&
       & 5x,'      # OF SEARCHES:              ',i11,/,&
       & 5x,'      # OF THREADS:               ',i11,/,&
       & 5x,'      # ELEMS FOUND 1ST STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS FOUND 2ND STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS NOT FOUND:          ',i11  ,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      TOTAL CPU TIME :            ',f11.2,3x,                      /,&
       & 5x,'      CPU TIME MASTER THREAD:     ',f11.2,3x,                      /,&
       & 5x,'      CPU TIME 1ST STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      CPU TIME 2ND STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )')

end subroutine elsest_statis

subroutine elsest_inbyte(tomem,rbyte,lbyte)
  !
  ! Give unit in bytes and conversion factor
  !
  use def_elsest
  implicit none
  integer(ip),  intent(in)  :: tomem
  real(rp),     intent(out) :: rbyte
  character(2), intent(out) :: lbyte

  if(tomem>=1024*1024*1024) then
     rbyte=1024.0_rp*1024.0_rp*1024.0_rp
     lbyte='Gb'
  else if(tomem>=1024*1024) then
     rbyte=1024.0_rp*1024.0_rp
     lbyte='Mb'
  else if(tomem>=1024) then
     rbyte=1024.0_rp
     lbyte='kb'
  else
     rbyte=1.0_rp
     lbyte=' b'
  end if

end subroutine elsest_inbyte

subroutine elsest_inelem(&
     mnode,ndime,nnode,lnods,ltopo,coord,xcoor,rpara,&
     ifoun,shapt,derit,coloc)
  !-----------------------------------------------------------------------
  !****f* Elsest
  ! NAME
  !    Inelem
  ! DESCRIPTION
  !    Identify shape function and derivatives in element
  ! INPUT
  ! OUTPUT
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_elsest, only     :  iunit,nthre,kfl_memor,lmini,lmaxi
  implicit none
  integer(ip), intent(in)  :: mnode,ndime
  integer(ip), intent(in)  :: nnode
  integer(ip), intent(in)  :: lnods(mnode,*)
  integer(ip), intent(in)  :: ltopo
  real(rp),    intent(in)  :: coord(ndime,*),xcoor(*)
  real(rp),    intent(in)  :: rpara(*)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: shapt(*),derit(*),coloc(*)
  integer(ip)              :: ptopo,pnode,inode,idime,ipoin
  real(rp)                 :: elcod(ndime,mnode)

  lmini = -rpara(1)
  lmaxi = 1.0_rp + rpara(1)
  ptopo = ltopo
  pnode = nnode
  do inode = 1,pnode
     ipoin = lnods(inode,ifoun)
     do idime = 1,ndime
        elcod(idime,inode) = coord(idime,ipoin)
     end do
  end do
  call elsest_chkelm(&
       ndime,ptopo,pnode,elcod,shapt,derit,&
       xcoor,coloc,ifoun,lmini,lmaxi)

end subroutine elsest_inelem


subroutine elsest_chkelm(&
     ndime,ptopo,pnode,elcod,shapf,deriv,&
     coglo,coloc,ifoun,lmini,lmaxi)
  !-----------------------------------------------------------------------
  !****f* Domain/elsest_chkelm
  ! NAME
  !    elsest_chkelm
  ! DESCRIPTION
  !    Check if a point belongs to an element
  ! INPUT
  !    NDIME ... Dimension
  !    PTOPO ... Element topology
  !              0 Quadrilateral  Hexahedra
  !              1 Triangle       Tetrahedra
  !              2    -           Pentahedra (wedge-shaped)
  !              3    -           Pyramid
  !    PNODE ... Number of element nodes
  !    ELCOD ... Element node coordinates
  !    COGLO ... Global coordinates of test point
  !    LMINI ... Minimum local coordinate (0.0)
  !    LMAXI ... Maximum local coordinate (1.0)
  ! OUTPUT
  !    IFOUN ... 1 if point is in element
  !              0 otherwise
  !    COLOC ... Local coordinates of test point
  !    SHAPF ... Shape function of test point in element
  !    DERIV ... Shape function derivatives of test point in element
  ! USED BY
  !***
  !-----------------------------------------------------------------------
!  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime,ptopo,pnode
  real(rp),    intent(in)  :: coglo(ndime),elcod(ndime,pnode)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: coloc(3),deriv(ndime,pnode),shapf(pnode)
  real(rp)                 :: ezzzt,lmaxi,lmini
  real(rp)                 :: xjacm(9),xjaci(9)

  if(  (ptopo==0.and.ndime==2.and.pnode==4).or.&
       (ptopo==1.and.ndime==2.and.pnode==3) ) then
     !
     ! Specific treatment: P1 and Q1 elements
     !
     call elsest_elq1p1(&
          pnode,lmini,lmaxi,elcod,coglo,coloc,&
          shapf,deriv,ifoun,xjacm,xjaci)
  else
     !
     ! Compute local from global coordinates
     !
     call elsest_newrap(&
          coglo,coloc,ndime,pnode,elcod,&
          xjacm,xjaci,shapf,deriv)
     ifoun=0

     if(ptopo==0.or.ptopo==-1) then
        !
        ! Quadrilateral and Hexahedra
        !
        if((coloc(1)>=-lmaxi).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=-lmaxi).and.(coloc(2)<=lmaxi)) then
              if((coloc(3)>=-lmaxi).and.(coloc(3)<=lmaxi)) then
                 ifoun=1
              end if
           end if
        end if

     else if(ptopo==1) then
        !
        ! Triangle and Tetrahedra
        !
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              if(ndime==2) then
                 coloc(3) = 1.0_rp-coloc(1)-coloc(2)
                 ezzzt    = 0.0_rp
              else
                 ezzzt    = 1.0_rp-coloc(1)-coloc(2)-coloc(3)
              end if
              if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
                 if((ezzzt>=lmini).and.(ezzzt<=lmaxi)) then
                    ifoun=1
                 end if
              end if
           end if
        end if

     else if(ptopo==2) then
        !
        ! Pentahedra
        !
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              if((coloc(3)>=lmini).and.(coloc(3)<=lmaxi)) then
                 ifoun=1
              end if
           end if
        end if

     end if

  end if

end subroutine elsest_chkelm




subroutine elsest_elq1p1(&
     pnode,lmini,lmaxi,elcod,coglo,coloc,&
     shapf,deriv,ifoun,xjacm,xjaci)
  !-----------------------------------------------------------------------
  !****f* domain/elmder
  ! NAME
  !    elmder
  ! DESCRIPTION
  !    Check if point with global coordinates (x,y)=COGLO is inside
  !    a triangle P1 or a quadrilateral Q1. The Q1 element is
  !    divided into two P1 elements. Returns the local coordinates
  !    (s,t)=COLOC
  !
  !    For P1 triangles we have:
  !    x = (1-s-t)*x1 + s*x2 + t*x3
  !    y = (1-s-t)*y1 + s*y2 + t*y3
  !
  !    This linear problem is solved for (s,t):
  !         (x3-x1)(y-y1) -(y3-y1)(x-x1)
  !    s =  ----------------------------
  !         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
  !
  !         (x-x1)(y2-y1) -(y-y1)(x2-x1)
  !    t =  ----------------------------
  !         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
  ! USES
  !    invmtx
  ! USED BY
  !    ***_elmope
  !    extnor
  ! SOURCE
  !-----------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: lmini,lmaxi,elcod(2,pnode),coglo(*)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: coloc(*),shapf(pnode),deriv(2,pnode)
  real(rp),    intent(out) :: xjacm(*),xjaci(*)
  real(rp)                 :: deter,colo3,x2x1,y2y1,x3x1,y3y1,xx1,yy1

  ifoun=0
  !
  ! P1 and Q1: Check if point is in first triangle: nodes 1-2-3
  !
  x2x1     = elcod(1,2)-elcod(1,1)
  y2y1     = elcod(2,2)-elcod(2,1)
  x3x1     = elcod(1,3)-elcod(1,1)
  y3y1     = elcod(2,3)-elcod(2,1)
  xx1      = coglo(1)  -elcod(1,1)
  yy1      = coglo(2)  -elcod(2,1)
  deter    = 1.0_8/(x3x1*y2y1-y3y1*x2x1)
  coloc(1) = deter*(x3x1*yy1-y3y1*xx1)
  coloc(2) = deter*(y2y1*xx1-x2x1*yy1)
  if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
     if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
        colo3 = 1.0_rp-coloc(1)-coloc(2)
        if(colo3>=lmini.and.colo3<=lmaxi) ifoun=1
     end if
  end if

  if(pnode==4) then
     !
     ! Q1: Check if point is in second triangle: nodes 1-3-4
     !
     if(ifoun==0) then
        x2x1     = elcod(1,3)-elcod(1,1)
        y2y1     = elcod(2,3)-elcod(2,1)
        x3x1     = elcod(1,4)-elcod(1,1)
        y3y1     = elcod(2,4)-elcod(2,1)
        xx1      = coglo(1)  -elcod(1,1)
        yy1      = coglo(2)  -elcod(2,1)
        deter    = 1.0_8/(x3x1*y2y1-y3y1*x2x1)
        coloc(1) = deter*(x3x1*yy1-y3y1*xx1)
        coloc(2) = deter*(y2y1*xx1-x2x1*yy1)
        if((coloc(1)>=lmini).and.(coloc(1)<=lmaxi)) then
           if((coloc(2)>=lmini).and.(coloc(2)<=lmaxi)) then
              colo3 = 1.0_rp-coloc(1)-coloc(2)
              if(colo3>=lmini.and.colo3<=lmaxi) ifoun=1
           end if
        end if
     end if
     if(ifoun==1) then
        call elsest_newrap(&
             coglo,coloc,2_ip,4_ip,elcod,xjacm,xjaci,shapf,deriv)
     end if

  else if(pnode==3.and.ifoun==1) then
     !
     ! P1: Compute shape function and derivatives
     !
     shapf(1)  = colo3
     shapf(2)  = coloc(1)
     shapf(3)  = coloc(2)
     deriv(1,1)=-1.0_rp
     deriv(1,2)= 1.0_rp
     deriv(1,3)= 0.0_rp
     deriv(2,1)=-1.0_rp
     deriv(2,2)= 0.0_rp
     deriv(2,3)= 1.0_rp

  end if

end subroutine elsest_elq1p1
!
!
!
subroutine elsest_newrap(&
     coglo,coloc,ndime,pnode,elcod,&
     xjacm,xjaci,shapf,deriv)
  !------------------------------------------------------------------------
  !
  !    Calculate the inverse transformation (x,y,z)-->(s,t,r)
  !
  !    Iterate for f(s_i)=x_i: ds = J^{-1}.xpoin - J^{-1}.f(s_i)
  !                               = J^{-1}dx
  !                            ds = s_{i+1}-s_i  (deltas)
  !                            dx = xpoin-f(s_i) (deltax)
  !    where the s_i's are the coordinates in the local basis and
  !    xpoin(idime)'s the real ones.
  !
  !------------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: coglo(3),elcod(ndime,pnode)
  real(rp),    intent(out) :: coloc(3)
  integer(ip)              :: inode,iiter,jdime,maxit,idime,ierro
  integer(ip)              :: ntens,jnode
  real(rp)                 :: shapf(pnode),deriv(ndime,pnode)
  real(rp)                 :: xjacm(ndime,ndime),xjaci(ndime,ndime)
  real(rp)                 :: deltx(3),delts(3),xnorm,detja
  real(rp)                 :: rnode,diame,coocg(3),dercg(ndime,pnode)
  real(rp)                 :: shacg(pnode)
  !
  ! Initial condition
  !
  coloc(1)=0.0_rp
  coloc(2)=0.0_rp
  coloc(3)=0.0_rp
  if(ndime==1) then
     ntens=1
  else if(ndime==2) then
     ntens=3
  else
     ntens=6
  end if
  call shafun(coloc,ndime,pnode,ntens,shapf,deriv,ierro)
  !
  ! Element diameter
  !
  do idime=1,ndime
     coocg(idime)=0.0_rp
  end do
  do inode=1,pnode
     do idime=1,ndime
        coocg(idime)=coocg(idime)+elcod(idime,inode)
     end do
  end do
  rnode=1.0_rp/real(pnode)
  do idime=1,ndime
     coocg(idime)=rnode*coocg(idime)
  end do
  call shafun(coocg,ndime,pnode,ntens,shacg,dercg,ierro)
  call mbmabt(xjacm,elcod,dercg,ndime,ndime,pnode)
  call invmtx(xjacm,xjaci,diame,ndime)
  if(diame<=0.0_rp) then
     diame=-1.0_rp
     do inode=1,pnode
        do jnode=inode+1,pnode
           do idime=1,ndime
              diame=max(diame,abs(elcod(idime,inode)-elcod(idime,jnode)))
           end do
        end do
     end do
  else
     diame=1.0_rp/(diame**(2.0_rp/real(ndime)))
  end if
  !
  ! Initialize dx=coglo-f(s_1) with s_1=(0,0,0)
  !  
  do idime=1,ndime
     deltx(idime)=0.0_rp
     do inode=1,pnode
        deltx(idime)=deltx(idime)&
             +shapf(inode)*elcod(idime,inode)
     end do
  end do
  xnorm=0.0_rp
  do idime=1,ndime
     deltx(idime) = coglo(idime)-deltx(idime)
     xnorm        = xnorm+deltx(idime)*deltx(idime)
  end do
  xnorm=xnorm*diame
  iiter=0
  maxit=10
  !
  ! Iterate for f(s_i)=x_i
  !
  do while( (xnorm>1e-8).and.(iiter<=maxit) )
     iiter=iiter+1
     !
     ! Compute J
     !
     do jdime=1,ndime
        do idime=1,ndime
           xjacm(idime,jdime)=0.0_rp
           do inode=1,pnode
              xjacm(idime,jdime)=xjacm(idime,jdime)&
                   +deriv(idime,inode)*elcod(jdime,inode)
           end do
        end do
     end do
     !
     ! Compute J^{-1}
     !
     call elsest_invmtx(xjacm,xjaci,detja,ndime)

     do idime=1,ndime
        delts(idime)=0.0_rp
        !
        ! Compute J^{-1}.dx
        !
        do jdime=1,ndime
           delts(idime)=delts(idime)+deltx(jdime)*xjaci(jdime,idime)
        end do
     end do
     do idime=1,ndime
        coloc(idime)=coloc(idime)+delts(idime)
     end do
     if ((coloc(1)>1d99).or.(coloc(2)>1d99).or.(coloc(3)>1d99)) then
        iiter=maxit+1
     else
        call shafun(coloc,ndime,pnode,ntens,shapf,deriv,ierro)
     end if
     !
     ! Compute f_i
     !
     do idime=1,ndime
        deltx(idime)=0.0_rp
        do inode=1,pnode
           deltx(idime)=deltx(idime)&
                +shapf(inode)*elcod(idime,inode)
        end do
     end do
     !
     ! Compute dx=coglo-f
     !         xnorm=sum ds^2
     !
     xnorm=0.0_rp
     do idime=1,ndime
        deltx(idime)=coglo(idime)-deltx(idime)
        xnorm=xnorm+delts(idime)*delts(idime)
     end do
     xnorm=xnorm*diame
  end do
  if(xnorm>1e-8) coloc(1)=2.0_rp

end subroutine elsest_newrap
!
!
!
subroutine mbmabt(a,b,c,n1,n2,n3)

!-----------------------------------------------------------------------
!
! This routine evaluates the matrix product A = B Ct, where
! A -> Mat(n1,n2), B -> Mat(n1,n3), C -> Mat(n2,n3)
!
!-----------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)  :: n1,n2,n3
  real(rp),    intent(in)  :: b(n1,n3), c(n2,n3)
  real(rp),    intent(out) :: a(n1,n2)
  integer(ip)              :: i,j,k
    
  do i=1,n1
     do j=1,n2
        a(i,j)=0.0_rp
        do k=1,n3
           a(i,j)=a(i,j)+b(i,k)*c(j,k)
        end do
     end do
  end do

end subroutine mbmabt
!
!
!
subroutine invmtx(a,b,deter,nsize)

!-----------------------------------------------------------------------
!
! This routine inverts a square matrix A -> Mat(nsize,nsize). The
! inverse is stored in B. Its determinant is DETER
!
!    
!-----------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)  :: nsize
  real(rp),    intent(in)  :: a(nsize,nsize)
  real(rp),    intent(out) :: b(nsize,nsize),deter
  real(rp)                 :: denom,t1,t2,t3,t4

  select case (nsize)
     
  case(1)
     deter=a(1,1)
     if(deter==0.0_rp) return
     b(1,1) = 1.0_rp/a(1,1)

  case(2)
     deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
     if(deter==0.0_rp) return
     denom=1.0_rp/deter
     b(1,1) = a(2,2)*denom
     b(2,2) = a(1,1)*denom
     b(2,1) =-a(2,1)*denom
     b(1,2) =-a(1,2)*denom  

  case(3)
     t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
     t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
     t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
     deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
     if(deter==0.0_rp) return
     denom = 1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
     b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
     b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
     b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
     b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
     b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

  case(4)
     t1= a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
          + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
          - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
     t2=-a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
          - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
          + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
     t3=+a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
          + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
          - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
     t4=-a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
          - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
          + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
     deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
     if(deter==0.0_rp) return
     denom=1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(4,1) = t4*denom
     b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
          - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
          + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
     b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
          + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
          - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
     b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
          - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
          + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
     b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
          + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
          - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
     b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
          + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
          - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
     b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
          - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
          + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
     b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
          + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
          - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
     b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
          - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
          + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
     b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
          - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
          + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
     b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
          + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
          - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
     b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
          - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
          + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
     b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
          + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
          - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom

     
  case default
     b=a
     call invert(b,nsize,nsize)

  end select

end subroutine invmtx
!
!
!
subroutine invert(a,nmax,ndm)

!-----------------------------------------------------------------------
!
! This routine performs the inversion of a ndm*ndm square matrix 
! or just part of it (nmax*nmax)
!
!-----------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)    :: ndm,nmax
  real(rp),    intent(inout) :: a(ndm,ndm)
  real(rp)                   :: d
  integer(ip)                :: n,j,i

  do n = 1,nmax
     d = a(n,n)
     do j = 1,nmax
        a(n,j) = -a(n,j)/d
     end do
     do i = 1,nmax
        if(n/=i) then
           do j = 1,nmax
              if(n/=j) a(i,j) = a(i,j) +a(i,n)*a(n,j)
           end do
        end if
        a(i,n) = a(i,n)/d
     end do
     a(n,n) = 1.0_rp/d
  end do
  
end subroutine invert
!
!
!
subroutine shafun(&
     posgp,ndime,nnode,ntens,shapf,deriv,ierro)

  !-----------------------------------------------------------------------
  !
  !    This routine evaluates shapf functions and their derivatives
  !    for linear and quadratic isoparametric elements
  !
  !-----------------------------------------------------------------------
  implicit none

  integer(ip), intent(in)  :: ndime,nnode,ntens
  real(rp),    intent(in)  :: posgp(ndime)
  integer(ip), intent(out) :: ierro
  real(rp),    intent(out) :: shapf(nnode),deriv(max(1_ip,ndime),nnode)

  !
  ! Initializations
  !
  shapf = 0.0_rp
  deriv = 0.0_rp

  !
  ! Evaluation of the shapf functions
  !
  if(ndime==0) then
!     call shape0(nnode,shapf,ierro)
  else if(ndime==1) then
!     call shape1(posgp(1),nnode,shapf,deriv,ierro)
  else if(ndime==2) then
     call shape2(posgp(1),posgp(2),nnode,shapf,deriv,ierro)
  else if(ndime==3) then
     call shape3(posgp(1),posgp(2),posgp(3),nnode,shapf,deriv,ierro)
  end if

end subroutine shafun
!
!
!
subroutine shape2(s,t,nnode,shapf,deriv,ierro)

!-----------------------------------------------------------------------
!
!    This routine evaluates shape functions and their first and
!    second derivatives for 2-d continuos standar interpolation 
!    elements.
!
!    TRIANGLES       3   6  &  10  nodes
!    QUADRILATERALS  4   9  &  16  nodes
!
!-----------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)  :: nnode
  integer(ip), intent(out) :: ierro
  real(rp),    intent(in)  :: s,t
  real(rp),    intent(out) :: deriv(2,nnode),shapf(nnode)
  real(rp)                 :: st,a1,a2,a3,ss,tt,s1,s2,s3,s4
  real(rp)                 :: t1,t2,t3,t4,s9,t9,c,a

  ierro=0
  if(nnode==3) then     
     shapf(1)=1.0_rp-s-t                                
     shapf(2)=s                                  
     shapf(3)=t                                           !  3
     deriv(1,1)=-1.0_rp                                   !   
     deriv(1,2)= 1.0_rp                                   !
     deriv(1,3)= 0.0_rp                                   !
     deriv(2,1)=-1.0_rp                                   !  1       2
     deriv(2,2)= 0.0_rp
     deriv(2,3)= 1.0_rp
  else if(nnode==4) then
     st=s*t                                           
     shapf(1)=(1.0_rp-t-s+st)*0.25_rp                     !  4         3
     shapf(2)=(1.0_rp-t+s-st)*0.25_rp                     !
     shapf(3)=(1.0_rp+t+s+st)*0.25_rp                     !      
     shapf(4)=(1.0_rp+t-s-st)*0.25_rp                     !
     deriv(1,1)=(-1.0_rp+t)*0.25_rp                       !  1         2
     deriv(1,2)=(+1.0_rp-t)*0.25_rp
     deriv(1,3)=(+1.0_rp+t)*0.25_rp
     deriv(1,4)=(-1.0_rp-t)*0.25_rp
     deriv(2,1)=(-1.0_rp+s)*0.25_rp
     deriv(2,2)=(-1.0_rp-s)*0.25_rp
     deriv(2,3)=(+1.0_rp+s)*0.25_rp
     deriv(2,4)=(+1.0_rp-s)*0.25_rp  
  else if(nnode==6) then
     a1=1.0_rp-s-t
     a2=s                                             
     a3=t
     shapf( 1)=(2.0_rp*a1-1.0_rp)*a1                      !  3
     shapf( 2)=(2.0_rp*a2-1.0_rp)*a2                      !   
     shapf( 3)=(2.0_rp*a3-1.0_rp)*a3                      !   
     shapf( 4)=4.0_rp*a1*a2                               !  6      5
     shapf( 5)=4.0_rp*a2*a3                               !     
     shapf( 6)=4.0_rp*a1*a3                               !   
     deriv(1,1)= 1.0_rp-4.0_rp*a1                         !  1     4     2
     deriv(1,2)= 4.0_rp*a2-1.0_rp    
     deriv(1,3)= 0.0_rp           
     deriv(1,4)= 4.0_rp*(a1-a2)   
     deriv(1,5)= 4.0_rp*a3        
     deriv(1,6)=-4.0_rp*a3       
     deriv(2,1)= 1.0_rp-4.0_rp*a1    
     deriv(2,2)= 0.0_rp           
     deriv(2,3)= 4.0_rp*a3-1.0_rp    
     deriv(2,4)=-4.0_rp*a2       
     deriv(2,5)= 4.0_rp*a2        
     deriv(2,6)= 4.0_rp*(a1-a3)
  else if(nnode==9) then
     ss=s*s
     st=s*t
     tt=t*t
     s1=s+1.0_rp
     t1=t+1.0_rp
     s2=s*2.0_rp
     t2=t*2.0_rp
     s9=s-1.0_rp                               
     t9=t-1.0_rp                                          !  4      7      3
     shapf( 1)=0.25_rp*s9*st*t9                           !
     shapf( 2)=0.25_rp*s1*st*t9                           !        
     shapf( 3)=0.25_rp*s1*st*t1                           !      
     shapf( 4)=0.25_rp*s9*st*t1                           !  8      9      6
     shapf( 5)=0.5_rp*(1.0_rp-ss)*t*t9                    !    
     shapf( 6)=0.5_rp*s*s1*(1.0_rp-tt)                    !   
     shapf( 7)=0.5_rp*(1.0_rp-ss)*t*t1                    !  
     shapf( 8)=0.5_rp*s*s9*(1.0_rp-tt)                    !  1      5      2
     shapf( 9)=(1.0_rp-ss)*(1.0_rp-tt)
     deriv(1,1)= 0.25_rp*t*t9*(-1.0_rp+s2)
     deriv(1,2)= 0.25_rp*(1.0_rp+s2)*t*t9
     deriv(1,3)= 0.25_rp*(1.0_rp+s2)*t*t1
     deriv(1,4)= 0.25_rp*(-1.0_rp+s2)*t*t1
     deriv(1,5)=-st*t9
     deriv(1,6)= 0.5_rp*(1.0_rp+s2)*(1.0_rp-tt)
     deriv(1,7)=-st*t1
     deriv(1,8)= 0.5_rp*(-1.0_rp+s2)*(1.0_rp-tt)
     deriv(1,9)=-s2*(1.0_rp-tt)
     deriv(2,1)= 0.25_rp*(-1.0_rp+t2)*s*s9
     deriv(2,2)= 0.25_rp*s*s1*(-1.0_rp+t2)
     deriv(2,3)= 0.25_rp*s*s1*(1.0_rp+t2)
     deriv(2,4)= 0.25_rp*s*s9*(1.0_rp+t2)
     deriv(2,5)= 0.5_rp*(1.0_rp-ss)*(-1.0_rp+t2)
     deriv(2,6)=-st*s1
     deriv(2,7)= 0.5_rp*(1.0_rp-ss)*(1.0_rp+t2)
     deriv(2,8)=-st*s9
     deriv(2,9)=-t2*(1.0_rp-ss)
  else if(nnode==10) then
     c=9.0_rp/2.0_rp
     a1=1.0_rp-s-t
     a2=2.0_rp/3.0_rp-s-t
     a3=1.0_rp/3.0_rp-s-t
     shapf( 1)=c*a1*a2*a3                                 !  3            
     shapf( 2)=c*(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s)*s    !               
     shapf( 3)=c*(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t)*t    !               
     shapf( 4)= 3.0_rp*c*a1*a2*s                          !  8    7  
     shapf( 5)=-3.0_rp*c*a1*(1.0_rp/3.0_rp-s)             !               
     shapf( 6)=-3.0_rp*c*(1.0_rp/3.0_rp-s)*s*t            !               
     shapf( 7)=-3.0_rp*c*s*(1.0_rp/3.0_rp-t)*t            !  9   10    6
     shapf( 8)=-3.0_rp*c*a1*(1.0_rp/3.0_rp-t)*t           !
     shapf( 9)= 3.0_rp*c*a1*a2*t                          !
     shapf(10)= 6.0_rp*c*a1*s*t                           !  1    4    5    2
     deriv(1, 1)=-c*(a1*a2+a1*a3+a2*a3)       
     deriv(1, 2)=-c*((2.0_rp/3.0_rp-s)*s&
          + (1.0_rp/3.0_rp-s)*s-(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s))
     deriv(1, 3)=0.0_rp
     deriv(1, 4)= 3.0_rp*c*(a1*a2-a1*s-a2*s)
     deriv(1, 5)=-3.0_rp*c*(a1*(1.0_rp/3.0_rp-s)&
          - a1*s-(1.0_rp/3.0_rp-s)*s)
     deriv(1, 6)=-3.0_rp*c*((1.0_rp/3.0_rp-s)*t-s*t)
     deriv(1, 7)=-3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
     deriv(1, 8)= 3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
     deriv(1, 9)= 3.0_rp*c*(-a1*t-a2*t)
     deriv(1,10)= 6.0_rp*c*(a1*t-s*t)
     deriv(2, 1)=-c*(a1*a2+a1*a3+a2*a3)
     deriv(2, 2)= 0.0_rp
     deriv(2, 3)=-c*((2.0_rp/3.0_rp-t)*t&
          + (1.0_rp/3.0_rp-t)*t-(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t))
     deriv(2, 4)= 3.0_rp*c*(-a1*s-a2*s)
     deriv(2, 5)=-3.0_rp*c*(-(1.0_rp/3.0_rp-s)*s)
     deriv(2, 6)=-3.0_rp*c*((1.0_rp/3.0_rp-s)*s)
     deriv(2, 7)=-3.0_rp*c*((1.0_rp/3.0_rp-t)*s-s*t)
     deriv(2, 8)=-3.0_rp*c*(-(1.0_rp/3.0_rp-t)*t&
          - a1*t+a1*(1.0_rp/3.0_rp-t))
     deriv(2, 9)= 3.0_rp*c*(-a1*t-a2*t+a1*a2)
     deriv(2,10)= 6.0_rp*c*(a1*s-s*t)  
  else if(nnode==16) then
     a =81.0_rp/256.0_rp
     c =1.0_rp/3.0_rp
     s1=1.0_rp+s
     s2=c+s
     s3=c-s
     s4=1.0_rp-s
     t1=1.0_rp+t
     t2=c+t
     t3=c-t
     t4=1.0_rp-t
     shapf( 1) =   a*s2*s3*s4*t2*t3*t4                   ! 4    10    9    3
     shapf( 2) =   a*s1*s2*s3*t2*t3*t4                   ! 
     shapf( 3) =   a*s1*s2*s3*t1*t2*t3                   ! 
     shapf( 4) =   a*s2*s3*s4*t1*t2*t3                   ! 11   16   15    8
     shapf( 5) =-3.0_rp*a*s1*s3*s4*t2*t3*t4              !
     shapf( 6) =-3.0_rp*a*s1*s2*s4*t2*t3*t4              !
     shapf( 7) =-3.0_rp*a*s1*s2*s3*t1*t3*t4              ! 12   13   14    7
     shapf( 8) =-3.0_rp*a*s1*s2*s3*t1*t2*t4              !
     shapf( 9) =-3.0_rp*a*s1*s2*s4*t1*t2*t3              !
     shapf(10) =-3.0_rp*a*s1*s3*s4*t1*t2*t3              ! 1     5    6    2
     shapf(11) =-3.0_rp*a*s2*s3*s4*t1*t2*t4                 
     shapf(12) =-3.0_rp*a*s2*s3*s4*t1*t3*t4
     shapf(13) = 9.0_rp*a*s1*s3*s4*t1*t3*t4
     shapf(14) = 9.0_rp*a*s1*s2*s4*t1*t3*t4
     shapf(15) = 9.0_rp*a*s1*s2*s4*t1*t2*t4
     shapf(16) = 9.0_rp*a*s1*s3*s4*t1*t2*t4
     deriv(1, 1)=  a *t2*t3*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1, 2)=  a *t2*t3*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 3)=  a *t1*t2*t3*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 4)=  a *t1*t2*t3*(-s2*s3-s2*s4+s3*s4)
     deriv(1, 5)=-3.0_rp*a*t2*t3*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(1, 6)=-3.0_rp*a*t2*t3*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1, 7)=-3.0_rp*a*t1*t3*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 8)=-3.0_rp*a*t1*t2*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 9)=-3.0_rp*a*t1*t2*t3*(-s1*s2+s1*s4+s2*s4)
     deriv(1,10)=-3.0_rp*a*t1*t2*t3*(-s1*s3-s1*s4+s3*s4)
     deriv(1,11)=-3.0_rp*a*t1*t2*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1,12)=-3.0_rp*a*t1*t3*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1,13)= 9.0_rp*a*t1*t3*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(1,14)= 9.0_rp*a*t1*t3*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1,15)= 9.0_rp*a*t1*t2*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1,16)= 9.0_rp*a*t1*t2*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(2, 1)=  a   *s2*s3*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 2)=  a   *s1*s2*s3*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 3)=  a   *s1*s2*s3*(-t1*t2+t1*t3+t2*t3)
     deriv(2, 4)=  a   *s2*s3*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2, 5)= -3.0_rp*a *s1*s3*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 6)= -3.0_rp*a *s1*s2*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 7)= -3.0_rp*a *s1*s2*s3*(-t1*t3-t1*t4+t3*t4)
     deriv(2, 8)= -3.0_rp*a *s1*s2*s3*(-t1*t2+t1*t4+t2*t4)
     deriv(2, 9)= -3.0_rp*a *s1*s2*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2,10)= -3.0_rp*a *s1*s3*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2,11)= -3.0_rp*a *s2*s3*s4*(-t1*t2+t1*t4+t2*t4)
     deriv(2,12)= -3.0_rp*a *s2*s3*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,13)=  9.0_rp*a *s1*s3*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,14)=  9.0_rp*a *s1*s2*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,15)=  9.0_rp*a *s1*s2*s4*(-t1*t2+t1*t4+t2*t4)
     deriv(2,16)=  9.0_rp*a *s1*s3*s4*(-t1*t2+t1*t4+t2*t4)
    
  else
     ierro=1
  end if
  
end subroutine shape2
!
!
!
subroutine shape3(s,t,z,nnode,shapf,deriv,ierro)

  !-----------------------------------------------------------------------
  !
  ! This routine evaluates shape functions and their first and
  ! second derivatives 3-d standar continuous interpolation
  ! elements.
  ! 
  ! TETRAHEDRA:  4  10  &  20  nodes
  ! HEXAHEDRA:   8  27  &  64  nodes
  ! PRISM:       6             nodes
  !
  !-----------------------------------------------------------------------
  !
  implicit none
  integer(ip), intent(in)  :: nnode
  integer(ip), intent(out) :: ierro
  real(rp),    intent(in)  :: s,t,z
  real(rp),    intent(out) :: deriv(3,nnode),shapf(nnode)
  integer(ip)              :: i,ii,jj
  real(rp)                 :: a1,a2,a3,a4,a,p1,p2,p3,z1,z2,z3,z4,s1,s2,s3,s4
  real(rp)                 :: t1,t2,t3,t4,sm,tm,zm,sq,tp,zp,s11,s21,s31,s41
  real(rp)                 :: t11,t21,t31,t41,z11,z21,z31,s12,s22,s32,s42
  real(rp)                 :: t12,t22,t32,t42,z41,z12,z22,z32,z42,sl,tl,zl
  real(rp)                 :: one8

  ierro =0 
  if(nnode==4) then
     !
     ! Linear tetrahedron 
     !
     shapf(   1) = 1.0_rp-s-t-z
     shapf(   2) = s
     shapf(   3) = t
     shapf(   4) = z
     deriv(1, 1) =-1.0_rp
     deriv(2, 1) =-1.0_rp
     deriv(3, 1) =-1.0_rp
     deriv(1, 2) = 1.0_rp
     deriv(2, 3) = 1.0_rp
     deriv(3, 4) = 1.0_rp

  else if(nnode==5) then
     !
     ! Linear Pyramid
     !
     one8        =  0.125_rp
     shapf(   1) =  one8*(1.0_rp-s)*(1.0_rp-t)*(1.0_rp-z)
     shapf(   2) =  one8*(1.0_rp+s)*(1.0_rp-t)*(1.0_rp-z)
     shapf(   3) =  one8*(1.0_rp+s)*(1.0_rp+t)*(1.0_rp-z) 
     shapf(   4) =  one8*(1.0_rp-s)*(1.0_rp+t)*(1.0_rp-z)
     shapf(   5) =  0.500_rp*(1.0_rp+z)       
     deriv(1, 1) = -one8*(1.0_rp-t)*(1.0_rp-z)
     deriv(2, 1) = -one8*(1.0_rp-s)*(1.0_rp-z)
     deriv(3, 1) = -one8*(1.0_rp-s)*(1.0_rp-t)
     deriv(1, 2) =  one8*(1.0_rp-t)*(1.0_rp-z)
     deriv(2, 2) = -one8*(1.0_rp+s)*(1.0_rp-z)
     deriv(3, 2) = -one8*(1.0_rp+s)*(1.0_rp-t)
     deriv(1, 3) =  one8*(1.0_rp+t)*(1.0_rp-z) 
     deriv(2, 3) =  one8*(1.0_rp+s)*(1.0_rp-z) 
     deriv(3, 3) = -one8*(1.0_rp+s)*(1.0_rp+t)
     deriv(1, 4) = -one8*(1.0_rp+t)*(1.0_rp-z)
     deriv(2, 4) =  one8*(1.0_rp-s)*(1.0_rp-z)
     deriv(3, 4) = -one8*(1.0_rp-s)*(1.0_rp+t)
     deriv(1, 5) =  0.0_rp
     deriv(2, 5) =  0.0_rp
     deriv(3, 5) =  0.5_rp 
   

  else if(nnode==6) then
     !
     ! Linear Prism
     !
     shapf(   1) = (1.0_rp-s-t)*(1.0_rp-z)
     deriv(1, 1) = z-1.0_rp
     deriv(2, 1) = z-1.0_rp
     deriv(3, 1) = s+t-1.0_rp
     shapf(   2) = s*(1.0_rp-z)
     deriv(1, 2) = 1.0_rp-z
     deriv(3, 2) = -s
     shapf(   3) = t*(1.0_rp-z)
     deriv(2, 3) = 1.0_rp-z
     deriv(3, 3) = -t
     shapf(   4) = (1.0_rp-s-t)*z
     deriv(1, 4) = -z
     deriv(2, 4) = -z
     deriv(3, 4) = 1.0_rp-s-t
     shapf(   5) = s*z
     deriv(1, 5) = z
     deriv(3, 5) = s
     shapf(   6) = t*z
     deriv(2, 6) = z
     deriv(3, 6) = t 
 
 

  else if(nnode==8) then
     !
     ! Trilinear brick 
     !   
     sm = 0.5_rp*(1.0_rp-s)
     tm = 0.5_rp*(1.0_rp-t)
     zm = 0.5_rp*(1.0_rp-z)
     sq = 0.5_rp*(1.0_rp+s)
     tp = 0.5_rp*(1.0_rp+t)
     zp = 0.5_rp*(1.0_rp+z)
     shapf(   1) = sm*tm*zm
     deriv(1, 1) =-0.5_rp*tm*zm
     deriv(2, 1) =-0.5_rp*sm*zm
     deriv(3, 1) =-0.5_rp*sm*tm

     shapf(   2) = sq*tm*zm
     deriv(1, 2) = 0.5_rp*tm*zm
     deriv(2, 2) =-0.5_rp*sq*zm
     deriv(3, 2) =-0.5_rp*sq*tm

     shapf(   3) = sq*tp*zm
     deriv(1, 3) = 0.5_rp*tp*zm
     deriv(2, 3) = 0.5_rp*sq*zm
     deriv(3, 3) =-0.5_rp*sq*tp

     shapf(   4) = sm*tp*zm
     deriv(1, 4) =-0.5_rp*tp*zm
     deriv(2, 4) = 0.5_rp*sm*zm
     deriv(3, 4) =-0.5_rp*sm*tp

     shapf(   5) = sm*tm*zp
     deriv(1, 5) =-0.5_rp*tm*zp
     deriv(2, 5) =-0.5_rp*sm*zp
     deriv(3, 5) = 0.5_rp*sm*tm

     shapf(   6) = sq*tm*zp 
     deriv(1, 6) = 0.5_rp*tm*zp
     deriv(2, 6) =-0.5_rp*sq*zp
     deriv(3, 6) = 0.5_rp*sq*tm

     shapf(   7) = sq*tp*zp
     deriv(1, 7) = 0.5_rp*tp*zp
     deriv(2, 7) = 0.5_rp*sq*zp
     deriv(3, 7) = 0.5_rp*sq*tp

     shapf(   8) = sm*tp*zp
     deriv(1, 8) =-0.5_rp*tp*zp
     deriv(2, 8) = 0.5_rp*sm*zp
     deriv(3, 8) = 0.5_rp*sm*tp


  else
     ierro = 1
  end if
  !
  ! Errors
  !      
  if(ierro==1) then
     write (*,*) nnode
     call runend('SHAPE3: NOT AVAILABLE ELEMENT INERPOLATION')
  end if


end subroutine shape3
!
!
!
subroutine elsest_invmtx(a,b,deter,nsize)
  !----------------------------------------------------------------------
  !
  ! This routine inverts a square matrix A -> Mat(nsize,nsize). The
  ! inverse is stored in B. Its determinant is DETER
  !
  !----------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)  :: nsize
  real(rp),    intent(in)  :: a(nsize,*)
  real(rp),    intent(out) :: deter,b(nsize,*)
  integer(ip)              :: jsize,isize
  real(rp)                 :: t1,t2,t3,t4,denom

  if(nsize==1) then
     !
     ! Inverse of a 1*1 matrix
     !
     deter=a(1,1)
     if(deter/=0.0_rp) return
     b(1,1) = 1.0_rp/a(1,1)

  else if(nsize==2) then
     !
     ! Inverse of a 2*2 matrix
     !
     deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
     if(deter==0.0_rp) return
     denom  = 1.0_rp/deter
     b(1,1) = a(2,2)*denom
     b(2,2) = a(1,1)*denom
     b(2,1) =-a(2,1)*denom
     b(1,2) =-a(1,2)*denom

  else if(nsize==3) then
     !
     ! Inverse of a 3*3 matrix
     !
     t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
     t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
     t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
     deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
     if(deter==0.0_rp) return
     denom  = 1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
     b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
     b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
     b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
     b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
     b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

  else if(nsize==4) then
     !
     ! Inverse of a 4*4 matrix
     !
     t1=   a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
          +a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
          -a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
     t2=  -a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
          -a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
          +a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
     t3=   a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
          +a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
          -a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
     t4=  -a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
          -a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
          +a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
     deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
     if(deter==0.0_rp) return
     denom = 1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(4,1) = t4*denom
     b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
          &   - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
          &   + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
     b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
          &   + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
          &   - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
     b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
          &   - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
          &   + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
     b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
          &   + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
          &   - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
     b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
          &   + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
          &   - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
     b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
          &   - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
          &   + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
     b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
          &   + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
          &   - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
     b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
          &   - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
          &   + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
     b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
          &   - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
          &   + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
     b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
          &   + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
          &   - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
     b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
          &   - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
          &   + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
     b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
          &   + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
          &   - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom
  else
     !
     ! Inverse of a nsize*nsize matrix
     !
     do isize=1,nsize
        do jsize=1,nsize
           b(isize,jsize)=a(isize,jsize)
        enddo
     enddo
     call elsest_invert(b,nsize,nsize)
  end if

end subroutine elsest_invmtx
!
!
!
subroutine elsest_invert(a,nmax,ndm)
  !----------------------------------------------------------------------
  !
  ! This routine performs the inversion of a ndm*ndm square matrix
  ! or just part of it (nmax*nmax)
  !
  !----------------------------------------------------------------------
  implicit none
  integer(ip), intent(in)    :: ndm,nmax
  real(rp),    intent(inout) :: a(ndm,ndm)
  integer(ip)                :: n,i,j
  real(rp)                   :: d

  do n = 1,nmax
     d = a(n,n)
     do j = 1,nmax
        a(n,j) = -a(n,j)/d
     end do
     do i = 1,nmax
        if(n.ne.i) then
           do j = 1,nmax
              if(n/=j) a(i,j) = a(i,j) +a(i,n)*a(n,j)
           end do
        end if
        a(i,n) = a(i,n)/d
     end do
     a(n,n) = 1.0_rp/d
  end do

end subroutine elsest_invert




!
!
!
!-----------------------------------------------------------------------
!
! Memory: - allocation check
!         - deallocation check
!         - reallocation
!
!-----------------------------------------------------------------------

  


  subroutine memrp1(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if 
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp1

  subroutine memrp2(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp2

  subroutine memrp3(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp3

  subroutine memrp4(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:,:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp4

  subroutine memip1(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(4)                  :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memip1

  subroutine memip2(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(4)                  :: varia(:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memip2

  subroutine memip3(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Integer(4)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(4)                  :: varia(:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memip3

  subroutine mei1p1(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Type(i1p)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    type(i1p)                   :: varia(:)
    integer(ip)                 :: lbyts,isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbyts=nsize*ip
          do isize=1,nsize
             nullify(varia(isize)%l)
          end do
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*ip
    end if
    call memsum(cumem,ithre,lbyts)

  end subroutine mei1p1

  subroutine mei1p2(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Type(i1p)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    type(i1p)                   :: varia(:,:)
    integer(ip)                 :: lbyts,nsize,isiz1,nsiz1,isiz2,nsiz2
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          nsiz1=size(varia,1)
          nsiz2=size(varia,2)
          lbyts=nsize*ip
          do isiz2=1,nsiz2
             do isiz1=1,nsiz1
                nullify(varia(isiz1,isiz2)%l)
             end do
          end do
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*ip
    end if
    call memsum(cumem,ithre,lbyts)

  end subroutine mei1p2

  subroutine memi81(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(8)                  :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memi81

  subroutine memi82(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(8)                  :: varia(:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memi82

  subroutine memi83(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Integer(8)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(8)                  :: varia(:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memi83

  subroutine memlg1(itask,ithre,istat,cumem,vanam,vacal,varia)
    use def_elsest, only : lg,i1p,memax,memor,nthre
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    logical(lg)                 :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=.false.
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memlg1
  
  subroutine memsum(cumem,ithre,lbyts)

    use def_elsest, only : lg,i1p,memax,memor,nthre
    implicit none


    integer(ip), intent(in)  :: ithre,lbyts
    integer(ip), intent(out) :: cumem
    integer(ip)              :: tomem,i
    
    cumem = cumem+lbyts
    tomem = 0
    do i=1,nthre
       tomem = tomem + memor(1,i)+memor(2,i)+memor(3,i)+memor(4,i)&
            &  +memor(5,i)+memor(6,i)+memor(7,i)+memor(8,i)&
            &  +memor(9,i)+memor(10,i)  
    end do
  !$OMP CRITICAL(chkmemax)
    memax = max(memax,tomem)
  !$OMP END CRITICAL(chkmemax)

  end subroutine memsum

END MODULE Elsest_MOD
