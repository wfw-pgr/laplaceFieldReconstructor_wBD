module PBiCGStabMod
  implicit none
  ! -- Constants        -- !
  integer         , parameter   :: Ncomm                = 7
  double precision, parameter   :: error_rate           = 1.d-12
  integer         , parameter   :: NeibPEtot_Max        = 2
  double precision, parameter   :: iteration_max_factor = 1.5d0
  ! -- Global Variables -- !
  integer                       :: myRank, PEtot, CommLen, NeibPEtot
  integer                       :: NxM, NxMw, iterMax, Nloc, Mloc
  integer                       :: NeibPE(NeibPEtot_Max)
  ! -- Allocatable Var. -- !
  integer         , allocatable ::       FromTo(:,:)
  integer         , allocatable ::           pt(:,:),     Nelem(:)
  double precision, allocatable ::         Amat(:,:),         Minv(:)
  integer         , allocatable :: export_Index(:)  , import_Index(:)
  integer         , allocatable :: export_Addrs(:)  , import_Addrs(:)
  integer         , allocatable :: request_Send(:)  , request_Recv(:)
  integer         , allocatable ::  status_Send(:,:),  status_Recv(:,:)
  double precision, allocatable ::       xloc(:,:,:),       sloc(:,:,:)
  ! -- suffix           -- !
  integer         , parameter   :: fr_=1, to_=2

contains
  ! --------------------   USAGE   ------------------------- !
  ! source :: include boundary condition of x.           --- !
  ! --- e.g.) x=phi_b at (i=1) --> source(1,j) = phi_b   --- !
  ! -------------------------------------------------------- !

  
  ! =================================================================== !
  ! ===  DefineCommunication  :: Define Comm.Table ( 1D Domain )    === !
  ! =================================================================== !
  subroutine DefineCommunication( N, M )
    implicit none
    include 'mpif.h'
    integer, intent(in) :: N, M
    integer             :: j, jM, ierr, surplus, iPE

    ! ------------------------------------- !
    ! --- [1] Define Grid & Partition   --- !
    ! ------------------------------------- !
    !  -- [1-1] PE infomation           --  !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )

    !  -- [1-2] FromTo :: Partition     --  !
    allocate( FromTo( 0:PEtot-1, 2 ) )
    Nloc    = N / PEtot
    surplus = N - Nloc * PEtot
    do iPE=0, PEtot-1
       FromTo(iPE,to_) = Nloc
    enddo
    do iPE=0, surplus-1
       FromTo(iPE,to_) = FromTo(iPE,to_) + 1
    enddo
    FromTo(0,fr_) = 1
    do iPE=1, PEtot-1
       FromTo(iPE,fr_) = FromTo(iPE-1,to_) + 1
       FromTo(iPE,to_) = FromTo(iPE-1,to_) + FromTo(iPE,to_)
    enddo
    
    !  -- [1-3] Grid Definition         --  !
    Nloc    = FromTo(myRank,to_) - FromTo(myRank,fr_) + 1
    Mloc    = M
    NxM     = Nloc * Mloc
    NxMw    = Nloc * Mloc + 2*Mloc
    iterMax = nint( NxM * iteration_max_factor )
    CommLen = Mloc * NeibPEtot_Max

    ! ------------------------------------- !
    ! --- [2] Allocate MPI variables    --- !
    ! ------------------------------------- !
    allocate( export_Index(0:NeibPEtot_Max), import_Index(0:NeibPEtot_Max) )
    allocate( export_Addrs(CommLen),         import_Addrs(CommLen)         )
    allocate( request_Send(NeibPEtot_Max),   request_Recv(NeibPEtot_Max)   )
    allocate(  status_Send(MPI_STATUS_SIZE,NeibPEtot_Max), &
         &     status_Recv(MPI_STATUS_SIZE,NeibPEtot_Max)  )

    ! ------------------------------------- !
    ! --- [3] Generate MPI Comm. Table  --- !
    ! ------------------------------------- !
    !  -- [3-0] Exception for PEtot=1 -- !
    if ( PEtot.eq.1 ) then
       NeibPEtot = 0
       return
    endif
    !  -- [3-1] Define Index/Address for Comm. -- !
    export_Index = 0
    import_Index = 0
    !  -- (1) Edge.1 myRank=0        -- !
    if ( myRank.eq.0       ) then
       NeibPEtot           = 1
       NeibPE(1)           = 1
       export_Index(1)     = Mloc
       import_Index(1)     = Mloc
       do j=1, Mloc
          export_Addrs(j)  = Nloc* j
          import_Addrs(j)  = NxM + j
       enddo
    endif
    !  -- (2) Edge.2 myRank=PEtot-1  -- !
    if ( myRank.eq.PEtot-1 ) then
       NeibPEtot           = 1
       NeibPE(1)           = PEtot-2
       export_Index(1)     = Mloc
       import_Index(1)     = Mloc
       do j=1, Mloc
          export_Addrs(j)  = Nloc*(j-1) + 1
          import_Addrs(j)  = NxM + j
       enddo
    endif
    !  -- (3) Inner Region           -- !
    if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
       NeibPEtot           = 2
       NeibPE(1)           = myRank-1
       NeibPE(2)           = myRank+1
       export_Index(1)     =   Mloc
       export_Index(2)     = 2*Mloc
       import_Index(1)     =   Mloc
       import_Index(2)     = 2*Mloc
       do j=1, Mloc
          jM               = Mloc+ j
          export_Addrs(j)  = Nloc*(j-1) + 1
          export_Addrs(jM) = Nloc*(j  )
          import_Addrs(j)  = NxM + j
          import_Addrs(jM) = NxM + jM
       enddo
    endif

    return
  end subroutine DefineCommunication


  ! =================================================================== !
  ! ===  MatrixMultiply  :: return Ux                               === !
  ! =================================================================== !
  ! -- it seems to be wrong, however, it works. -- !
  ! -- & Prof.Nakajima recommended this coding. -- !
  function MatrixMultiply( Umat, x, Upt, UNelem, N, M, Nitem )
    implicit none
    integer         , intent(in) :: N, M, Nitem
    integer         , intent(in) :: UNelem(N), Upt(N,Nitem)
    double precision, intent(in) :: Umat(N,Nitem), x(M)
    integer                      :: i, j 
    double precision             :: MatrixMultiply(N), sumi

    !$omp parallel default(none) &
    !$omp shared(N,Umat,UNelem,Upt,x,MatrixMultiply) private(i,j,sumi)
    !$omp do
    do i=1, N
       sumi = 0.d0
       do j=1, UNelem(i)
           sumi = sumi + Umat(i,j) * x( Upt(i,j) )
       enddo
       MatrixMultiply(i) = sumi
    enddo
    !$omp end do
    !$omp end parallel
    return
  end function MatrixMultiply


  ! =================================================================== !
  ! ===  DiagonalMultiply  ::  return u(i)v(i)                      === !
  ! =================================================================== !
  function DiagonalMultiply( u, v, N )
    implicit none
    integer         , intent(in) :: N
    double precision, intent(in) :: u(N), v(N)
    double precision             :: DiagonalMultiply(N)
    integer                      :: i

    !$omp parallel default(none) &
    !$omp shared(N,DiagonalMultiply,u,v) private(i)
    !$omp do
    do i=1, N
       DiagonalMultiply(i) = u(i) * v(i)
    enddo
    !$omp end do
    !$omp end parallel
    return
  end function DiagonalMultiply


  ! =================================================================== !
  ! ===  myDotProduct  :: primitive OpenMP/MPI dot product Function === !
  ! =================================================================== !
  function myDotProduct( r1, r2, N )
    implicit none
    include 'mpif.h'
    integer                      :: i, ierr
    integer         , intent(in) :: N
    double precision, intent(in) :: r1(N), r2(N)
    double precision             :: sum
    double precision             :: myDotProduct
    
    sum = 0.d0
    !$omp parallel default(none) &
    !$omp shared(N,r1,r2,sum) private(i) 
    !$omp do reduction(+:sum)
    do i=1, N
       sum = sum + r1(i)*r2(i)
    enddo
    !$omp end do
    !$omp end parallel
    call MPI_AllReduce( sum, myDotProduct, 1, MPI_DOUBLE_PRECISION, &
         &              MPI_SUM, MPI_COMM_WORLD, ierr )
    return
  end function myDotProduct


  ! =================================================================== !
  ! ===  daxpy  :: Double precision ax + y calculation              === !
  ! =================================================================== !
  function daxpy( acoef, xvec, yvec, N )
    implicit none
    integer         , intent(in)    :: N
    double precision, intent(in)    :: acoef
    double precision, intent(in)    :: xvec(N), yvec(N)
    integer                         :: i
    double precision                :: daxpy(N)
    
    !$omp parallel default(none) &
    !$omp shared(acoef,xvec,yvec,daxpy,N) private(i)
    !$omp do 
    do i=1, N
       daxpy(i) = acoef * xvec(i) + yvec(i)       
    enddo
    !$omp end do
    !$omp end parallel
    return
  end function daxpy


  ! =================================================================== !
  ! ===  BoundaryCommunicate  ::  Exchange Boundary Info.           === !
  ! =================================================================== !
  subroutine BoundaryCommunicate( xvc )
    implicit none
    include 'mpif.h'
    double precision, intent(inout) :: xvc(NxMw)
    integer                         :: j, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(CommLen)
    double precision                :: recvBuff(CommLen)
    if ( neibPEtot.eq.0 ) return
    
    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    do neib=1, NeibPEtot
       do j=export_Index(neib-1)+1, export_Index(neib)
          sendBuff(j) = xvc(export_Addrs(j))
       enddo
    enddo
    
    ! ------------------------------------- !
    ! --- [2] ISend / IRecv             --- !
    ! ------------------------------------- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot
       is    = export_Index(neib-1) + 1
       len_s = export_Index(neib  ) - export_Index(neib-1)
       call MPI_Isend( sendBuff(is)      , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot
       ir    = import_Index(neib-1) + 1
       len_r = import_Index(neib  ) - import_Index(neib-1)
       call MPI_Irecv( recvBuff(ir)      , len_r, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Recv(neib), ierr )
    enddo
    
    ! ------------------------------------- !
    ! --- [3] WaitAll / Update          --- !
    ! ------------------------------------- !
    !  -- [3-1] WaitAll Receive         --  !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )
    !  -- [3-2] Update  xVector         --  !
    do neib=1, NeibPEtot
       do j=import_Index(neib-1)+1, import_Index(neib)
          xvc(import_Addrs(j)) = recvBuff(j)
       enddo
    enddo
    !  -- [3-3] WaitAll Send            --  !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )
    
    return
  end subroutine BoundaryCommunicate


  ! =================================================================== !
  ! ===  PBiCGSTAB_Engine_MPI  ::  Main Routine                     === !
  ! =================================================================== !
  subroutine PBiCGSTAB_Engine_MPI( b, x )
    implicit none
    double precision, intent(in)  :: b(NxM)
    double precision, intent(out) :: x(NxM)
    integer                       :: iter, i
    double precision              :: alpha, beta, c1, c2, c3, rr, bb, criterion
    double precision              :: w(NxMw), r(NxMw), p(NxMw), r0(NxMw), p0(NxMw)
    double precision              :: y(NxMw), z(NxMw), e(NxMw),  v(NxMw),  u(NxMw)
    logical, parameter            :: flag__display_residual = .true.
    ! --- Data Structure ------------------------------------------!
    ! b(NxM ) :: { |NxM| }        = { |Data| }                     !
    ! w(NxMw) :: { |NxM|, |2*M| } = { |Data|, |Commu. Region.| }   !
    ! -------------------------------------------------------------!
    
    ! ------------------------------------- !
    ! --- [1] 1st step                  --- !
    ! ------------------------------------- !
    ! w (:)     =  0.d0
    w (:)     =  b(:)
    r (1:NxM) =  b(1:NxM) - MatrixMultiply( Amat, w, pt, Nelem, NxM, NxMw, Ncomm )
    r0(1:NxM) =  r(1:NxM)
    p (1:NxM) =  r(1:NxM)
    p0(1:NxM) =  r(1:NxM)
    c1        = myDotProduct( r0(1:NxM), r(1:NxM), NxM )
    bb        = myDotProduct(  b(1:NxM), b(1:NxM), NxM )
    criterion = error_rate * bb

    ! ------------------------------------- !
    ! --- [2] Main Loop ( P-BiCGSTAB )  --- !
    ! ------------------------------------- !
    do iter=1, iterMax
       z(1:NxM)  = DiagonalMultiply( Minv(1:NxM), p(1:NxM), NxM )
       call BoundaryCommunicate( z )
       y(1:NxM)  = MatrixMultiply( Amat, z, pt, Nelem, NxM, NxMw, Ncomm )
       c2        = myDotProduct( r0(1:NxM), y(1:NxM), NxM )
       alpha     = c1 / c2

       e(1:NxM)  = daxpy( -alpha, y(1:NxM), r(1:NxM), NxM )
       u(1:NxM)  = DiagonalMultiply( Minv(1:NxM), e(1:NxM), NxM )
       call BoundaryCommunicate( u )
       v(1:NxM)  = MatrixMultiply( Amat, u, pt, Nelem, NxM, NxMw, Ncomm )
       
       c3        =   myDotProduct( e(1:NxM), v(1:NxM), NxM ) &
            &      / myDotProduct( v(1:NxM), v(1:NxM), NxM )
       
       !$omp parallel default(none) &
       !$omp shared(w,alpha,z,c3,u,NxM) private(i)
       !$omp do
       do i=1, NxM
          w(i)   = w(i) + alpha*z(i) + c3*u(i)
       enddo
       !$omp end do
       !$omp end parallel
       
       r (1:NxM) = daxpy( -c3, v(1:NxM), e(1:NxM), NxM )
       rr        = myDotProduct( r(1:NxM) , r(1:NxM), NxM )
       if ( flag__display_residual ) then
          write(6,"(a,i8,a,i8,a,e15.8,a,e15.8)") "  iter / iterMax ::", iter, " / ", iterMax, &
               &    " residual / criterion ::", rr, " / ", criterion
       endif
       if ( rr.lt.criterion ) exit
       c1        = myDotProduct( r0(1:NxM), r(1:NxM), NxM )
       beta      = c1 / ( c2 * c3 )
       
       !$omp parallel default(none) &
       !$omp shared(p,r,beta,c3,y,NxM) private(i)
       !$omp do
       do i=1, NxM
          p(i)   = r(i) + beta * ( p(i) - c3 * y(i) )
       enddo
       !$omp end do
       !$omp end parallel

    enddo

    ! ------------------------------------- !
    ! --- [3] Post Process              --- !
    ! ------------------------------------- !
    !  -- [3-1] store Answer            --  !
    x(1:NxM) = w(1:NxM)
    !  -- [3-2] Display Loop Summary    --  !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '[ PBiCGSTAB    @ myPBICGSTAB ]'
       write(6,'(8x,2(a,i16  ))') 'Iteration  = ', iter, '      / ', iterMax
       write(6,'(8x,2(a,e16.9))') 'Error      = ', rr,   '      / ', criterion
    endif
    return
  end subroutine PBiCGSTAB_Engine_MPI

  
end module PBiCGStabMod
