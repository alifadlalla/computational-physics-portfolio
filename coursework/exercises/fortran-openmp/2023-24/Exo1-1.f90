PROGRAM SEQ

   IMPLICIT NONE

   INTEGER :: i, j, k, N=128, dx, dy, dz

   REAL, ALLOCATABLE :: A(:,:,:), B(:,:,:)

   ALLOCATE(A(N,N,N), B(N,N,N))

   DO k = 1, N
      DO j = 1, N
         DO i = 1, N
            A(i,j,k) = REAL(i+j+k)
         END DO
      END DO
   END DO

   B = 0.


! Invert the k,j,i loops to see the effect of cache miss
   DO i = 3, N-2
      DO j = 3, N-2
         DO k = 3, N-2
! Invert the k,j,i loops to see the effect of cache miss

            DO dx = i-2, i+2
               DO dy = j-2, j+2
                  DO dz = k-2, k+2
                     B(i,j,k) = B(i,j,k) + A(dx,dy,dz)
                  END DO
               END DO
            END DO

         END DO
      END DO
   END DO

   PRINT*, B(N/2,N/2,N/2)

END PROGRAM SEQ
