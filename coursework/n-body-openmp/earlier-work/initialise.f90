program initialise
Implicit None
integer :: i, j, N=1
real, dimension(3, 2000) :: r_vec, v_vec
real :: r, x, y, z

  do while (N < 2000)
    call random_number(x)
    call random_number(y)
    call random_number(z)
    r = (2*x-1)**2 + (2*y-1)**2 + (2*z-1)**2
    !print*, r
    if (r <=1) then 
        ! positions vector
        r_vec(1, N) = 2*x-1
        r_vec(2, N) = 2*y-1
        r_vec(3, N) = 2*z-1

        ! velocity vector
        v_vec(1, N) = -(2*y-1)
        v_vec(2, N) = 2*x-1
        v_vec(3, N) = 0

        N = N+1
    end if
  end do 
  open(unit = 10, file ='initial_position.dat')
    do i=1, N
        write(10,*) r_vec(1, i), r_vec(2, i), r_vec(3, i) 
    end do
  close(10)

  open(unit = 11, file ='initial_velocity.dat')
    do i=1, N
        write(11,*) v_vec(1, i), v_vec(2, i), v_vec(3, i) 
    end do
  close(11)

  



end program initialise