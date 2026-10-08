program integration
  implicit none
  integer, parameter :: nt = 1000, N = 1000
  integer :: i, j, k, t
  real, dimension(3, N) :: r, v, a, a_new
  real, dimension(3) :: pos_diff

  real :: x, y, z, KE, PE, r_module, dt, epsilon = 0.005, m = 1.0, G = 1.0
  integer :: nseed
  integer, allocatable :: seed(:)

  ! File units for binary output
  integer, parameter :: pos_unit = 20, vel_unit = 21

  ! Seeding for random number generation
  call random_seed(size=nseed)
  allocate(seed(nseed))
  seed = 123456789    ! Arbitrary seed
  call random_seed(put=seed)
  deallocate(seed)

  ! Open binary files for positions and velocities
  open(unit=pos_unit, file="positions.bin", form="unformatted", action="write")
  open(unit=vel_unit, file="velocities.bin", form="unformatted", action="write")
  ! Save energies to files
  open(unit=12, file="KE.dat")
  open(unit=13, file="PE.dat")

  ! Initialize particle positions and velocities
  j = 1
  do while (j <= N)
     call random_number(x)
     call random_number(y)
     call random_number(z)
     if ((2*x-1)**2 + (2*y-1)**2 + (2*z-1)**2 <= 1) then
         r(1, j) = 2*x - 1
         r(2, j) = 2*y - 1
         r(3, j) = 2*z - 1

         ! Initialize velocity (solid body rotation)
         v(1, j) = -r(2, j)
         v(2, j) =  r(1, j)
         v(3, j) = 0.0

         j = j + 1
     end if
  end do

  dt = 5.0/nt

  ! Initialize accelerations
a = 0.0
do i = 1, N
    do j = 1, N
        if (i == j) cycle
        pos_diff = r(:, j) - r(:, i)
        r_module = sqrt(sum(pos_diff**2) + epsilon**2)
        a(:, i) = a(:, i) + G * m * pos_diff / r_module**3
    end do
end do

! Main loop
do t = 1, nt
    ! Update positions
    r = r + v * dt + 0.5 * a * dt**2

    ! Compute new accelerations
    a_new = 0.0
    do i = 1, N
        do j = 1, N
            if (i == j) cycle
            pos_diff = r(:, j) - r(:, i)
            r_module = sqrt(sum(pos_diff**2) + epsilon**2)
            a_new(:, i) = a_new(:, i) + G * m * pos_diff / r_module**3
        end do
    end do

    ! Update velocities
    v = v + 0.5 * (a + a_new) * dt

    ! Save positions and velocities to binary files
    write(pos_unit) r
    write(vel_unit) v

    ! Update accelerations for the next step
    a = a_new

    ! Compute energies
    KE = 0.0
    PE = 0.0
    do i = 1, N
        KE = KE + 0.5 * m * sum(v(:, i)**2)
        do j = 1, N
            if (i == j) cycle
            r_module = sqrt(sum((r(:, j) - r(:, i))**2) + epsilon**2)
            PE = PE - G * m**2 / r_module
        end do
    end do
    write(12,*) KE
    write(13,*) PE
end do


  ! Close binary files
  close(pos_unit)
  close(vel_unit)
  close(12)
  close(13)

end program integration
