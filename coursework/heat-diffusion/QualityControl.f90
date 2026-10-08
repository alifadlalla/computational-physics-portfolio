program alphaHeat1D
    
    use functions
    implicit none
    
    real(kind=8), parameter :: PI = 3.14159265358979D0

    real(kind=8) :: bar_length, bar_width, bar_height, alphaCu, lambdaCu,&
    L0,T0, time0, dt, dx, r, du_0, du_1, du_N, cdu, h
    integer(kind=8) :: M, N, i, i_max, j, t, tau, flag, BC_x0, BC_xL
    real(kind=8), dimension(:, :), allocatable :: T_fwd_euler, T_bwd_euler,T_analt,f
    real(kind=8), dimension(1000) :: eps_1, eps_2, dt_array  
    
    data bar_length/0.195D0/, bar_width/0.02D0/, bar_height/0.01D0/, M/1200D0/, N/14D0/,&
    alphaCu/1.11D-04/, lambdaCu/401D0/, L0/0.195/, T0/295/, h/50/

    time0 = (L0**2)/alphaCu     ! t0 = 342.5s
    !dt = 2000.0/(time0*M)
    dx = 0.195/(L0*N)
    cdu = -h*1.0/(lambdaCu)
 
    !r = dt/(dx**2)


    du_0 = 0
    du_1 = 0
    du_N = 0    

    !print*, bar_length*100

    !-------------Block 1: Quality Control---------------!
    N = 14
    !N= 40
    dx = 0.195/(L0*N)
    i = 0
    do M = 2500, 73000, 1000
    !do M = 2500, 20100, 100 
    if (M==20000) then
        print*, "excuse you ?"
       ! N = 40
    else 
        !N = 14
    endif
    !print*, i
    flag = 0
    dt = 2000.0/(time0*M)
    r = dt/(dx**2)
    if (r > 0.5 ) then 
        print*, "Error r = ", r, " > 0.5"
        STOP
    endif

    allocate(T_analt(M+1, N+1))
    allocate(T_fwd_euler(M+1, N+1))
    allocate(T_bwd_euler(M+1, N+1))
    ! M changes during this sweep: the source array must use the current size.
    allocate(f(M+1, N+1))
    f=0
    
    do j= 1, N+1
        T_analt(1,j)     = 1 + sin(PI*(j-1)*dx) 
        T_fwd_euler(1,j) = 1 + sin(PI*(j-1)*dx)
        T_bwd_euler(1,j) = 1 + sin(PI*(j-1)*dx)
    end do
    
    do t = 1, M
        !----Backward Euler----!
        BC_x0 =  0
        BC_xL =  0

        call Backward_Euler(r,N,dx,du_0,du_1,du_N,T_bwd_euler,t,cdu,BC_x0, BC_xL)
        !----------------------!
        
        !-----------Forwad Euler & Analytical--------------!
        !Dirichlet at x=0
        call Dirichlet_0(T_analt, t)
        call Dirichlet_0(T_fwd_euler, t)

        !Dirichlet at x=L
        call Dirichlet_L(T_analt, t, N)
        call Dirichlet_L(T_fwd_euler, t, N)
        
        ! Euler: x= deltaX ----> N-1
        do j = 2, N
            T_analt(t+1, j) = 1 + sin(PI*(j-1)*dx)*exp((-1)*PI*PI*(t)*dt)
            call euler(T_fwd_euler, r, dt, f, t, j)
        end do
        ! error 
        if (ABS(T_analt(t, N/2) - 1.5D0) < 1e-3 .and. flag==0) then
            i = i+1
            tau = t
            flag = 1
            eps_1(i) = ABS(T_fwd_euler(tau, N/2) - T_analt(tau, N/2))*100.0/ (T_analt(tau, N/2)*1.0)
            eps_2(i) = ABS(T_bwd_euler(tau, N/2) - T_analt(tau, N/2))*100.0/ (T_analt(tau, N/2)*1.0)
            dt_array(i) = dt
        end if
    end do
    print*, "calculating ... This is nessecary for the code to run"
    i_max = i
    deallocate(T_analt)
    deallocate(T_fwd_euler)
    deallocate(T_bwd_euler)
    deallocate(f)

    
    end do
    
    
    print*, "Finished"
    ! write T(x,t) to a file
    open(unit=10, file = "error.txt")
    write(10,*) "delta_t,eps_fwd_euler,eps_bwd_euler"
    do i=1, i_max    
        write(10, '(E15.7, ",", E15.7, ",", E15.7)') dt_array(i), eps_1(i), eps_2(i)
    end do
    close (10)
    
end program alphaHeat1D
