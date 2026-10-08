

program main
	implicit none

	integer :: omp_get_num_threads, omp_get_thread_num
	integer, parameter :: m=20000, n=20000
	real, allocatable :: mat(:,:)
	real :: trace
	integer :: i, j
	integer :: nb_th, rang, t1, t2, rate
	integer, allocatable :: sync(:)
	
	!$OMP PARALLEL
	nb_th = omp_get_num_threads()
	!$OMP END PARALLEL


	allocate(mat(m,n))
	allocate(sync(nb_th))
	sync(:) = 1

	do i=1, n
		do j=1, m
			mat(j,i) = real(i+j)/(m+n)	
		end do
	end do

	call system_clock(t1, rate)

	!$OMP PARALLEL PRIVATE(i,j, rang) shared(sync)
	rang = omp_get_thread_num() + 1
	
	do j=2, n
		if (rang /= 1) then
			do
				!$OMP FLUSH(sync)
				if (sync(rang-1) > sync(rang)) exit
			end do
		!$OMP FLUSH(sync, mat)
		end if

		!$OMP DO SCHEDULE(STATIC)
!		do i=2, m
		do i=m, 2, -1
			mat(i, j) = (mat(i,j) + mat(i-1, j) + mat(i, j-1))/3.0
		end do
		!$OMP END DO NOWAIT
		sync(rang) = sync(rang) + 1
		!$OMP FLUSH(sync, mat)
	end do
	!$OMP END PARALLEL

	call system_clock(t2, rate)

	write(*,*) mat(m/2, n/2)
	write(*,*) "Compute time:", real(t2-t1)/rate, "sec"


end program main





