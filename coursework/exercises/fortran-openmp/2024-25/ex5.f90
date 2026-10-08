Program ex5
Implicit None

integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
INTEGER(kind=8) :: i,N=10000000
real(kind=8) :: dx, pi, x

dx = 1.0/N
pi = 0.0

!$OMP PARALLEL DO shared(dx) reduction(+ : pi) private(x)
!!$OMP PARALLEL DO shared(pi, dx) private(x) 

Do i=1, N
    x = dx*i
    !!$OMP ATOMIC
    pi = pi + 4.0*dx/(1.0 + x*x)
end do
!$OMP END PARALLEL DO

PRINT*, "pi = ",pi

end program ex5