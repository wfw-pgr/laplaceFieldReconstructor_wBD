module ioUtilityMod
contains

  ! ====================================================== !
  ! === load__source_and_boundary                      === !
  ! ====================================================== !
  subroutine load__source_and_boundary
    use variablesMod
    implicit none
    integer          :: i, j, k
    double precision :: dsurf, dflag
    character(cLen)  :: cmt
    integer, parameter :: x_=1, y_=2, z_=3

    ! ------------------------------------------------------ !
    ! --- [1] load Amatrix                               --- !
    ! ------------------------------------------------------ !
    open(lun,file=trim(srcFile),status="old")
    read(lun,*)
    read(lun,*)
    read(lun,*) cmt, LK, LJ, LI, nCmp
    allocate( source(6,LI,LJ,LK), BField(6,LI,LJ,LK), svec(LI,LJ,LK), xvec(LI,LJ,LK) )
    allocate( boundary_flag(LI,LJ,LK) )
    source(:,:,:,:)      = 0.d0
    BField(:,:,:,:)      = 0.d0
    svec  (:,:,:)        = 0.d0
    xvec  (:,:,:)        = 0.d0
    boundary_flag(:,:,:) = 0.d0
    do k=1, LK
       do j=1, LJ
          do i=1, LI
             read(lun,*) source(1:6,i,j,k), dsurf, dflag
             boundary_flag(i,j,k) = int( dflag )
          enddo
       enddo
    enddo
    close(lun)
    write(6,*) "[load__source_and_boundary] srcFile is loaded...."
    
    ! ------------------------------------------------------ !
    ! --- [2] define delta :: dx, dy, dz                 --- !
    ! ------------------------------------------------------ !
    dx = source(x_,2,1,1) - source(x_,1,1,1)
    dy = source(y_,1,2,1) - source(y_,1,1,1)
    dz = source(z_,1,1,2) - source(z_,1,1,1)

    ! ------------------------------------------------------ !
    ! --- [3] copy coordinate into BField                --- !
    ! ------------------------------------------------------ !
    do k=1, LK
       do j=1, LJ
          do i=1, LI
             BField(x_,i,j,k) = source(x_,i,j,k)
             BField(y_,i,j,k) = source(y_,i,j,k)
             BField(z_,i,j,k) = source(z_,i,j,k)
          enddo
       enddo
    enddo
    
    return
  end subroutine load__source_and_boundary
  

  ! ====================================================== !
  ! === save results in a file                         === !
  ! ====================================================== !
  subroutine save__results
    use variablesMod
    implicit none
    integer            :: i, j, k

    ! ------------------------------------------------------ !
    ! --- [1] save results                               --- !
    ! ------------------------------------------------------ !
    nCmp = 6
    open(lun,file=trim(rslFile),status="replace")
    write(lun,"(a)") "# xg yg zg bx by bz"
    write(lun,"(a,1x,2(i10))") "#", LK*LJ*LI, nCmp
    write(lun,"(a,1x,4(i10))") "#", LK, LJ, LI, nCmp
    do k=1, LK
       do j=1, LJ
          do i=1, LI
             write(lun,"(6(e15.8,1x))") BField(1:6,i,j,k)
          enddo
       enddo
    enddo
    close(lun)
    write(6,*) "[save__results] rslFile is saved...."
    
    return
  end subroutine save__results

  
end module ioUtilityMod
