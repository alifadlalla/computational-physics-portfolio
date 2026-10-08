program main
    implicit None

	integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
	
	!$OMP PARALLEL

	write(*,*) "Hello world thread", omp_get_thread_num(), "out of", omp_get_num_threads()

	!$OMP END PARALLEL

	write(*,*) "end parallel"


end program main