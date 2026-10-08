Program ex1
Implicit None

INTEGER :: i,j,k,N=256,dx,dy,dz
REAL, ALLOCATABLE :: A(:,:,:), B(:,:,:)
ALLOCATE(A(N,N,N), B(N,N,N))
Do i=1, N
  Do j=1, N
    Do k=1, N
      A(i,j,k) = real(i+j+k)
    end do
  end do 
end do

B = 0.0
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

PRINT*, B(N/2,N/2,N/2)

end program ex1