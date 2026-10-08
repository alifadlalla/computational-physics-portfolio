program main
    use settings
    use iso_c_binding, only : c_backspace
    implicit none

    integer :: i, j, t
    real(kind=dp), allocatable :: pos(:,:), vel(:,:), acc(:,:), acc_new(:,:), r_ij(:)
    real(kind=dp) :: ec, ep, r_module

    ! random seed variables
    integer, parameter :: seed_value = 12345678
    integer, allocatable :: seed(:)
    integer :: seed_size

    ! Parallelization variables
    integer :: num_threads
    
    ! Read some settings in 'input.ini'
    call readsettings()

    !'''''''''Random Seed Logic'''''''''''''''''
    ! Get seed size and allocate
    call random_seed(size = seed_size)
    allocate(seed(seed_size))

    ! Initialize seed array
    do i = 1, seed_size
        seed(i) = i*seed_value
    end do
    call random_seed(put = seed)
    deallocate(seed)
    !'''''''''''''''''''''''''''''''''''''''''

    ! Allocate the main arrays
    allocate(pos(3, n), vel(3, n), acc(3, n), acc_new(3, n), r_ij(3))

    ! Initialize positions and velocities
    !$omp parallel do private(i)
    do i=1, n
        ! Initialize positions
        do
            call random_number(pos(:, i))
            pos(:, i) = 2.0_dp * pos(:, i) - 1.0_dp
            if (sum(pos(:, i)**2) <= 1.0_dp) exit     
        end do

        ! Initialize velocities based on solid rotation
        vel(1, i) = -pos(2, i)
        vel(2, i) = pos(1, i)
        vel(3, i) = 0.0_dp
    end do
    !$omp end parallel do

    

    ! Files to save the positions and the energies
    open(unit = 10, file = 'position.bin', access='direct', recl=3*n*8, status='replace') ! Binary for positions
    open(unit = 100, file = 'energy.bin', access='direct', recl=2*8, status='replace')   ! Binary for energies

    write(10,rec=1) pos

    ! Initialize accelerations
    acc = 0.0_dp
    !$omp parallel do private(i, j, r_ij, r_module)
    do i = 1, n
        do j = 1, n
            if (i == j) cycle
            r_ij = pos(:, j) - pos(:, i)
            r_module =  sqrt(sum(r_ij**2) + epsilon**2)
            acc(:, i) = acc(:, i) + gg * m * r_ij / r_module**3
            end do
    end do
    !$omp end parallel do

    ! Main loop
    do t = 1, nt
        ep = 0.0_dp
        ec = 0.0_dp
        acc_new = 0.0_dp

        ! Position update (parallelized)
        !$omp parallel do private(i)
        do i = 1, n
            pos(:, i) = pos(:, i) + vel(:, i) * dt + 0.5_dp * acc(:, i) * dt**2
        end do
        !$omp end parallel do

        ! position update with Fortran intrinic function
        !pos = pos + vel * dt + 0.5_dp * acc * dt**2

        ! Compute new accelerations
        !$omp parallel do private(i, j, r_ij, r_module) reduction(+:ep)
        do i = 1, n
            do j = 1, n
                if (i == j) cycle
                r_ij = pos(:, j) - pos(:, i)
                r_module =  sqrt(sum(r_ij**2) + epsilon**2)
                acc_new(:, i) = acc_new(:, i) + gg * m * r_ij / r_module**3
                ep = ep - 1.0_dp /r_module
            end do
        end do
        !$omp end parallel do

        ! Update velocities, acc, ec using Fortran intrinic functions
        !vel = vel + 0.5_dp * (acc + acc_new) * dt
        !acc = acc_new
        !ec = 0.5_dp * m * sum(vel**2)

        !$omp parallel do private(i)
        do i = 1, n
            vel(:, i) = vel(:, i) + 0.5_dp * (acc(:, i) + acc_new(:, i)) * dt
        end do
        !$omp end parallel do

        acc = acc_new

        ! Kinetic energy calculation (parallelized with reduction)
        !$omp parallel do private(i) reduction(+:ec)
        do i = 1, n
            ec = ec + 0.5_dp * m * sum(vel(:, i)**2)
        end do
        !$omp end parallel do

        ! Save kinetic and potential energies
        write(100, rec=t) 0.5_dp * ep * gg * m**2, ec
        write(10, rec=t+1) pos !writing data from the second record. the first one contains the initial position
    end do

    close(10)
    close(100)

    ! Deallocate the main arrays
    deallocate(pos, vel, acc, acc_new, r_ij)

end program main