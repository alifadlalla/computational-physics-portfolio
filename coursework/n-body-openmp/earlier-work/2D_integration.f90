program integration
Implicit None
integer :: i, j, k,t,nt=100, N=2000
real :: h=2, dt
real, dimension(2, 100) :: r,v,a
dt = 0.639/nt
r = 0.
v = 0.
a = 0.
r(1,1) = 0.
r(2,1) = 2.

v(1,1) = 1.
v(2,1) = 0.

a(1,:) = 0.
a(2,:) = -9.81

do t = 1, nt-1
    !a(1, t+1) = 0.
    !a(2, t+1) = -9.81
    do i = 1,2
        r(i, t+1) = r(i,t) + v(i,t)*dt + 0.5*a(i,t)*dt**2
        !print*, r(2, t+1)
        v(i, t+1) = v(i, t) + 0.5*(a(i,t)+a(i,t+1))*dt
    end do

end do
!print*, r(2, 3)
!print*, r(1, nt)

 open(unit = 10, file ='2D_evolution_r.dat')
    do t=1, nt
        write(10,*) r(1, t), r(2, t) 
    end do
  close(10)

  open(unit = 11, file ='2D_evolution_v.dat')
    do t=1, nt
        write(11,*) v(1, t), v(2, t) 
    end do
  close(11)

end program integration