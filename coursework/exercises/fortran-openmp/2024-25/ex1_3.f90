program ex1_3

implicit None
INTEGER :: i,j,k, N=200
real, ALLOCATABLE :: A(:,:,:)
ALLOCATE(A(N,N,N))
A = 0.

!open(9, file="A_ascii.dat")
! do k = 1, N
!     do j = 1, N
!        do i = 1, N
!           write(9, *) A(i, j, k)
!        end do
!     end do
!  end do
!close(9)

!open(10, file="A_unformatted", Form="unformatted")
! do i = 1, N
!     do j = 1, N
!        do k = 1, N
!           write(10) A(i, j, k)
!        end do
!     end do
!  end do
!close(10)

open(11, file="A_unformatted_dumped", Form="unformatted", access="stream")
  write(11) A
close(11)



end program ex1_3