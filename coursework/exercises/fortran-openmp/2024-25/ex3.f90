program main
    implicit None

	integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
    integer :: i, thread_rank, limit
	real :: val

    limit = 3
    i = 9
    val = 1.1

	!$OMP PARALLEL shared(limit, val) Private(i, thread_rank, val )
	thread_rank = omp_get_thread_num()
    do i=0, limit
        write(*,*) "loop index",i, " thread", thread_rank
    end do

    val = val +1
	!$OMP END PARALLEL

    write(*,*) "Outside: ", "i=",i,"val=",val,"thread_rank=", thread_rank

end program main