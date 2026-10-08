program main
    use settings
    use iso_c_binding, only : c_backspace
    implicit NONE

    integer :: i, j, t
    real(kind=dp), allocatable :: pos(:,:), vel(:,:), acc(:,:), acc_new(:,:), posdiff(:)
    real(kind=dp) :: ec, ep

    ! Read some settings in 'input.ini'
    call readsettings()

    ! Files to save the positions and the energies
    open(unit = 10, file = 'position.dat', action = 'write')
    open(unit = 100, file = 'energy.dat', action="write")

    ! Allocate the main arrays
    allocate(pos(n, 3), vel(n, 3), acc(n, 3), acc_new(n, 3), posdiff(3))

    do i=1, n
        ! Initialize positions
        pos(i, :) = 2

        do while (sqrt(sum(pos(i, :)**2)) > 1.0)
            call random_number(pos(i, :))
            pos(i, :) = 1-2*pos(i, :)
        end do

        ! Initialize velocities based on solid rotation
        vel(i, 1) = -pos(i, 2)
        vel(i, 2) = pos(i, 1)
        vel(i, 3) = 0

    end do

    acc = 0_dp

    write(10, *) n, nt, dp
    write(10, *) pos

    ! Initialize accelerations
    do i=1, n
        do j=1, n

            if (i == j) cycle

            posdiff = pos(j, :) - pos(i, :)
            acc(i, :) = acc(i, :) + (gg*m / sqrt(sum(posdiff**2) + epsilon**2)**3) * posdiff
        end do
    end do
    
    
    ! Main loop
    !$omp parallel private(i, j, t, posdiff, r_module, force_factor) 
    !$OMP DO SCHEDULE(Guided)
    do t=1, nt
        ec = 0_dp
        ep = 0_dp
        acc_new = 0_dp

        ! Update positions
        pos = pos + vel*dt + 0.5*acc*dt**2

        ! Compute new accelerations
        do i=1, n
            do j=1, n

                if (i == j) cycle
 
                posdiff = pos(j, :) - pos(i, :)
                acc_new(i, :) = acc_new(i, :) + (gg*m / sqrt(sum(posdiff**2) + epsilon**2)**3) * posdiff

                ! Compute the potential energy
                ep = ep - 1_dp/sqrt(sum(posdiff**2)+epsilon**2)

            end do
        end do

        ! Update velocities
        vel = vel + 0.5*(acc + acc_new)*dt

        ! Update current accelerations and compute kinetic energy 
        acc = acc_new
        ec = ec + sum(vel**2)

        ! Save kinetic and potantial energies
        write(100, *) 0.5*ep*gg*m**2, ec*0.5*m
        write(10, *) pos
        write(*, "(A,A,I7)", advance="no") repeat(c_backspace, 17), "Iteration ", t
    end do
    !$OMP END DO SCHEDULE(Guided)
    !$omp end parallel

    close(10)
    close(100)

    ! Deallocate the main arrays
    deallocate(pos, vel, acc, acc_new, posdiff)
    
end program main