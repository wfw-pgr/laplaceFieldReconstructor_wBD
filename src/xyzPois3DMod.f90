module xyzPois3DMod
  use PBiCGStabMod
  ! --------------------   USAGE   ------------------------- !
  ! source :: include boundary condition of x.           --- !
  ! --- e.g.) x=phi_b at (i=1) --> source(1,j) = phi_b   --- !
  ! -------------------------------------------------------- !
contains
  
  ! =================================================================== !
  ! ===  BiCGSTABCrt3D_MPI  :: BiCGSTAB solver for Cartesian 3D     === !
  ! =================================================================== !
  subroutine BiCGSTABCrt3D_MPI( xglobal, source, bdflag, dx, dy, dz, LI, LJ, LK )
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: LI, LJ, LK
    double precision, intent(in)    :: dx, dy, dz
    integer         , intent(in)    :: bdflag (LI,LJ,LK)
    double precision, intent(in)    :: source (LI,LJ,LK)
    double precision, intent(inout) :: xglobal(LI,LJ,LK)
    integer                         :: ierr, iPE, LIbuf, LJbuf, LKbuf, Nitem_buf
    integer                         :: LIloc, LJloc, LKloc, LJLKloc, LJxLK, LIxLJxLK
    logical, save                   :: Flag__MakeMatrix = .true.
    double precision, allocatable   :: xbuf(:,:,:)
    integer         , parameter     :: Fr_=1, To_=2
    
    ! ------------------------------------------------------ !
    ! --- [1] Make Matrix if No Matrix                   --- !
    ! ------------------------------------------------------ !

    if ( Flag__MakeMatrix ) then
       ! -- [1-1] variables for comm.   --  !
       LJxLK = LJ * LK
       call DefineCommunication( LI, LJxLK )
       LIloc    = Nloc
       LJloc    = LJ
       LKloc    = LK
       LJLKloc  = Mloc
       LIxLJxLK = NxM
       ! -- [1-2] Make Matrix           --  !
       allocate( xloc(LIloc,LJloc,LKloc), sloc(LIloc,LJloc,LKloc) )
       allocate( Amat (LIxLJxLK,Ncomm), pt  (LIxLJxLK,Ncomm),  &
            &    Nelem(LIxLJxLK)      , Minv(LIxLJxLK)         )
       xloc(:,:,:) = 0.d0
       sloc(:,:,:) = 0.d0
       Amat(:,:)   = 0.d0
       pt  (:,:)   = 0.d0
       Nelem(:)    = 0.d0
       Minv (:)    = 0.d0
       call XYZLaplacian( dx, dy, dz, LIloc, LJloc, LKloc )
       call set__boundary( bdflag, LIloc, LJloc, LKloc )
       Flag__MakeMatrix = .false.
    endif
    
    ! ------------------------------------- !
    ! --- [2] call PBiCGSTAB_Engine     --- !
    ! ------------------------------------- !
    sloc(:,:,:) = source( FromTo(myRank,Fr_):FromTo(myRank,To_), 1:LJloc, 1:LKloc )
    call PBiCGSTAB_Engine_MPI( sloc, xloc )
    ! call MPI_Finalize( ierr )
    ! stop
    ! ------------------------------------- !
    ! --- [3] BroadCast Solution        --- !
    ! ------------------------------------- !
    do iPE=0, PEtot-1
       ! -- [3-1] Prepare iPE's Solution -- !
       LJbuf      = LJloc
       LKbuf      = LKloc
       LIbuf      = FromTo(iPE,2) - FromTo(iPE,1) + 1
       Nitem_buf  = LIbuf * LJbuf * LKbuf
       allocate( xbuf(LIbuf,LJbuf,LKbuf) )
       if ( myRank.eq.iPE ) xbuf = xloc
       ! -- [3-2] BroadCast              -- !
       call MPI_Bcast( xbuf, Nitem_buf     , MPI_DOUBLE_PRECISION, &
            &          iPE , MPI_COMM_WORLD, ierr )
       ! -- [3-3] Save in x(From:To)     -- !
       xglobal( FromTo(iPE,1):FromTo(iPE,2), 1:LJloc, 1:LKloc ) = xbuf(1:LIbuf,1:LJbuf,1:LKbuf)
       deallocate( xbuf )
    enddo
    
    return
  end subroutine BiCGSTABCrt3D_MPI


  ! =================================================================== !
  ! ===  XYLaplacian  :: Calculate Laplacian ( cartesian )          === !
  ! =================================================================== !
  subroutine XYZLaplacian( dx, dy, dz, LI, LJ, LK )
    implicit none
    include 'mpif.h'
    integer         , intent(in)  :: LI, LJ, LK
    double precision, intent(in)  :: dx, dy, dz
    integer                       :: i, j, k, m, mL, mR, LIxLJ, LJxLK, LIxLJxLK
    logical                       :: Flag_SizeERROR
    double precision              :: coef(Ncomm)

    ! ------------------------------------------------------ !
    ! --- [1] Preparation                                --- !
    ! ------------------------------------------------------ !
    !  -- [1-1] Display                 --  !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'       ) '[ XYZLaplacian  @ PBiCGStabMod.f90 ]'
       write(6,'(5x,2(a,i12))') '* A-Matrix (local) :: ', LI*LJ*LK, ' x ', LI*LJ*LK
    endif
    !  -- [1-2] Coefficient             --  !
    LIxLJ    = LI * LJ
    LJxLK    = LJ * LJ
    LIxLJxLK = LI * LJ * LK
    coef(1)  = 1.d0 / dz**2
    coef(2)  = 1.d0 / dy**2
    coef(3)  = 1.d0 / dx**2
    coef(4)  = - 2.d0 * ( 1.d0 / dx**2 + 1.d0 / dy**2 + 1.d0 / dz**2 )
    coef(5)  = 1.d0 / dx**2
    coef(6)  = 1.d0 / dy**2
    coef(7)  = 1.d0 / dz**2
    !  -- [1-3] Initialization          --  !
    Amat(:,:) = 0.d0
    pt(:,:)   = 0
    Minv(:)   = 1.d0
    
    ! ------------------------------------------------------ !
    ! --- [2] Make Matrix                                --- !
    ! ------------------------------------------------------ !
    !  -- [2-1] Main Region                              --  !
    m = 1
    do k=1, LK
       do j=1, LJ
          do i=1, LI
             ! - Laplacian - !
             Amat(m,1)   = coef(1)
             Amat(m,2)   = coef(2)
             Amat(m,3)   = coef(3)
             Amat(m,4)   = coef(4)
             Amat(m,5)   = coef(5)
             Amat(m,6)   = coef(6)
             Amat(m,7)   = coef(7)
             pt(m,1)     = m-LIxLJ
             pt(m,2)     = m-LI
             pt(m,3)     = m-1
             pt(m,4)     = m
             pt(m,5)     = m+1
             pt(m,6)     = m+LI
             pt(m,7)     = m+LIxLJ
             Nelem(m)    = 7
             m           = m+1
          enddo
       enddo
    enddo
    write(6,*) m

    ! ------------------------------------------------------ !
    ! --- [3] Dirichlet Boundary condition               --- !
    ! ------------------------------------------------------ !
    m = 1
    do k=1, LK
       do j=1, LJ
          do i=1, LI
             if ( (i.eq.1).or.(i.eq.LI).or.(j.eq.1).or.(j.eq.LJ) &
                  &                    .or.(k.eq.1).or.(k.eq.LK) ) then
                Amat(m,1) = 1.d0
                pt(m,1)   = m
                Nelem(m)  = 1
                Minv(m)   = 1.d0 
             endif
             m = m+1
          enddo
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [4] Boundary ( x=xMin, x=xMax ) ( parallel )   --- !
    ! ------------------------------------------------------ !
    m = LIxLJxLK+1
    if ( PEtot.gt.1 ) then
       ! if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
       if ( myRank.gt.0 ) then
          !  -- (Left) Neibour PE       --  !
          mL = 1
          do k=1, LK
             do j=1, LJ
                pt(mL,3)     = m
                mL           = mL + LI
                m            = m  + 1
             enddo
          enddo
       endif
       if ( myRank.lt.PEtot-1 ) then
          !  -- (Right) Neibour PE      --  !
          mR = LI
          do k=1, LK
             do j=1, LJ
                pt(mR,5)     = m
                mR           = mR + LI
                m            = m  + 1
             enddo
          enddo
       endif
    endif
    
    ! ------------------------------------------------------ !
    ! --- [5] Size Check                                 --- !
    ! ------------------------------------------------------ !
    Flag_SizeERROR = .false.
    if ( PEtot.gt.1 ) then
       if ( ( myRank.eq.0 ).or.( myRank.eq.PEtot-1 ) ) then
          if ( (m-1).ne.(LIxLJxLK+  LJxLK) ) Flag_SizeERROR = .true.
       else
          if ( (m-1).ne.(LIxLJxLK+2*LJxLK) ) Flag_SizeERROR = .true.
       endif
    endif
    if ( PEtot.eq.1 ) then
       if    ( (m-1).ne.(LIxLJxLK        ) ) Flag_SizeERROR = .true.
    endif
    if ( Flag_SizeERROR ) then
       write(6,*) ' [WARNING] Size and m_count are Incompatable.'
       write(6,*) "      myRank / PEtot :: ", myRank, " / ", PEtot
       write(6,*) "    LIxLJxLK         :: ", LIxLJxLK
       write(6,*) "    LIxLJxLK+2*LJxLK :: ", LIxLJxLK+2*LJxLK
       write(6,*) "    LIxLJxLK+  LJxLK :: ", LIxLJxLK+  LJxLK
       write(6,*) "               LIxLJ :: ", LIxLJ
       write(6,*) "                 m-1 :: ", m-1
       stop
    endif

    return
  end subroutine XYZLaplacian


  ! ====================================================== !
  ! === set boundary for flag =1 element               === !
  ! ====================================================== !
  subroutine set__boundary( bdflag, LI, LJ, LK )
    implicit none
    integer, intent(in) :: LI, LJ, LK
    integer, intent(in) :: bdflag(LI,LJ,LK)
    integer             :: i, j, k, m

    m = 1
    do k=1, LK
       do j=1, LJ
          do i=1, LI
             if ( bdflag(i,j,k).eq.1 ) then
                Amat(m,:)   = 0.d0
                Amat(m,1)   = 1.d0
                pt(m,:)     = m
                Minv(m)     = 1.d0 
                Nelem(m)    = 1
             endif
             m           = m+1
          enddo
       enddo
    enddo

    return
  end subroutine set__boundary
  
  
end module xyzPois3DMod
