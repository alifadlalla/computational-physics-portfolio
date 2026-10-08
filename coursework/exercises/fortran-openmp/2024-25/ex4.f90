Program ex1
Implicit None

integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
INTEGER :: i,j,k,N=256,dx,dy,dz
REAL, ALLOCATABLE :: A(:,:,:), B(:,:,:)
ALLOCATE(A(N,N,N), B(N,N,N))

!$OMP PARALLEL shared(N, A) Private(i,j,k)
!$OMP DO SCHEDULE(Guided)
!SCHEDULE(Static, N) where N = total nb of iterations/ nob of threads
Do i=1, N
  Do j=1, N
    Do k=1, N
      A(i,j,k) = real(i+j+k)
    end do
  end do 
end do
!$OMP END DO
!$OMP END PARALLEL

B = 0.0
!$OMP PARALLEL shared(N, A, B) Private(dx,dy,dz,i,j,k)
!$OMP DO SCHEDULE(guided)
Do k=3, N-2
  Do j=3, N-2
    Do i=3, N-2

      do dz=k-2, k+2
        do dy=j-2, j+2
          do dx=i-2, i+2
            B(i,j,k) = B(i,j,k) + A(dx,dy,dz)
          end do
        end do 
      end do
    end do
  end do 
end do
!$OMP END DO
!$OMP END PARALLEL

PRINT*, B(N/2,N/2,N/2)

end program ex1