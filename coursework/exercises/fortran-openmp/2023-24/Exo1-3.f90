PROGRAM SEQ

   IMPLICIT NONE

   INTEGER :: N=300, i

   REAL, ALLOCATABLE :: A(:,:,:)

   ALLOCATE(A(N,N,N))
   A = 0.

! Test1: compare file dump to the disk or to the RAM, with file in 1 piece
   OPEN(100, FILE='toto', FORM='unformatted')
!   OPEN(100, FILE='/dev/shm/toto', FORM='unformatted')
   WRITE(100) A
   CLOSE(100)

!! Test2: compare file dump to the disk or to the RAM, with file in many piece
!   OPEN(100, FILE='tata', FORM='formatted')
!!   OPEN(100, FILE='/dev/shm/toto', FORM='unformatted')
!   DO i=1,N
!      WRITE(100,*) A(i,:,:)
!   END DO
!   CLOSE(100)


END PROGRAM SEQ
