program main
  use variablesMod
  use ioUtilityMod
  use xyzPois3DMod, only : BiCGStabCrt3D_MPI
  implicit none
  include "mpif.h"
  integer :: cmp
  integer :: i, j, k

  call MPI_Init     ( ierr )
  call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
  call MPI_Comm_Size( MPI_COMM_WORLD, PEtot , ierr )

  write(6,*)
  write(6,*) "[main.f90] myRank == ", myRank
  write(6,*) "[main.f90] PEtot  == ", PEtot
  write(6,*)
  
  call load__source_and_boundary

  do cmp=4, 6
     svec(:,:,:)       = source(cmp,:,:,:)
     call BiCGSTABCrt3D_MPI( xvec, svec, boundary_flag, dx, dy, dz, LI, LJ, LK )
     BField(cmp,:,:,:) = xvec(:,:,:)
  enddo
  
  call save__results
  call MPI_Finalize ( ierr )
  
end program main
  
