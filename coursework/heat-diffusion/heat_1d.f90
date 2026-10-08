program alphaHeat1D
    
    use functions
    implicit none
    
    real(kind=8), parameter :: PI = 3.14159265358979D0

    real(kind=8) :: bar_length, bar_width, bar_height, alphaCu, lambdaCu,&
    L0,T0, time0, dt, dx, r, du_0, du_1, du_N, cdu, h, area, sigma
    integer(kind=8) :: M, N, j, t, BC_x0, BC_xL, method
    real(kind=8), dimension(:, :), allocatable :: u, f, f_EM_radiation  
    
    
    data bar_length/0.195D0/, bar_width/0.02D0/, bar_height/0.01D0/, M/20000D0/, N/40D0/,&
    alphaCu/1.11D-04/, lambdaCu/401D0/, L0/0.195/, T0/295/, h/200/, sigma/5.67037442D-8/

    allocate(u(M+1, N+1))
    allocate(f(M+1, N+1))
    allocate(f_EM_radiation(M+1, N+1))
    
    area = 2*(bar_length*bar_height) + 2*(bar_length*bar_width) + 2*(bar_width*bar_height)

    !--------------------Print relevant info--------------------------!
    print*, "--> System:"
    print  "(' Bar length: ', F5.3, ' m')", bar_length
    !print ('F5.3'), 'Bar length:', bar_length*100, 'cm'
    print"(' Bar width (length along y axis): ',F5.2,' m')", bar_width
    print"(' Bar height (length along z axis): ',F5.2, ' m')", bar_height
    print*, "Material: Copper"
    print"(' Themal Diffusivity: alpha= ',F8.6,' m^2/s')", alphaCu
    print"(' Initial temperature profile T(x,0): [uniform] T = 'F8.2 ,' K')", T0
    print*, "Left boundary condition: [Neumann] phi(x = 0, all t) = 60000 W/m^2"
    print"(' Right boundary condition: [Dirichlet] T(x = L, all t) = 'F8.2 ,' K')", T0
    print*, "--------------------------------------------------------------------"
    
    time0 = (L0**2)/alphaCu     ! t0 = 342.5s
    dt = 2000.0/(time0*M)
    dx = 0.195/(L0*N)
    cdu = -h*1.0/(lambdaCu)  ! h: convection coefficient
    r = dt/(dx**2)
    
    if (r > 0.5 ) then 
        print*, "Error r = ", r, " > 0.5"
        STOP
    endif

    print*, "--> Simulation Parameters:"
    print*, "Total Simulation Time: 2000 s"
    print"(' Nb of points for the discretisation in time: ',I5)", M
    print"(' Nb of points for the discretisation in space: ',I2)", N
    print"(' Characteristic length L0 = ',F8.2 ,' m')", L0
    print"(' Characteristic time t0 = ',F8.2 ,' s')", time0
    print"(' Characteristic temperature T0 = ',F8.2 ,' K')", T0
    print"(' time step: delta_t = ',F8.6, ' in units of t0')", dt 
    print"(' spatial grid spacing: delta_x = ',F8.3, ' in units of L0')", dx
    
    print*, "r=", r    

    !-------------Simulation---------------!
    method = 0  !can be switched between 0 to use the Forward Scheme or 1 for the Backward scheme
    
    do j= 1, N+1
        u(1,j) = 1
        !u(1,j)= 1 + sin(PI*(j-1)*dx)      !initial condition for the anayltical solution
    enddo
    
    do t = 1, M
        
        if (method==1) then
        !----Backward Euler----!
        BC_x0 =  1    ! 0 for Dirchlet, 1 for Neuman
        BC_xL =  2    ! 0 for Dirchlet, 1 for Neuman, 2 for Robin
        
        !parameter for the boundary conditions
        du_0 = 0 
        du_1 = -0.098905 

        ! Simulate the Temperature
        f=0

        call Backward_Euler(r,N,dx,du_0,du_1,du_N,u,t,cdu,BC_x0, BC_xL)
        
        !----------------------!
        
        else if (method == 0) then
        BC_x0 =  1     ! 0 for Dirchlet, 1 for Neuman
        BC_xL =  0     ! 0 for Dirchlet, 1 for Neuman
        
        !parameter for the boundary conditions
        du_0 = -0.098905 
        du_N = 0

        ! Source/Sink Term
        do j = 1, N+1
            f_EM_radiation(t,j) = (time0*alphaCu/(T0*lambdaCu)) * area*sigma*(u(t,j)**4 - T0**4)/(bar_height*bar_length*bar_width)
        enddo
        
        ! choose between No sink Term or With EM radiation
        ! the default is without
        !f = 0
        f = f_EM_radiation

        !-----------Forwad Euler--------------!
        if (BC_x0 == 0) then ! Dirichlet at x=0
        call Dirichlet_0(u, t)

        else if(BC_x0 == 1) then ! Neuman at x=0
        call Neuman_0(u, du_0, r, dx, dt, f, t)
        endif

        if (BC_xL == 0) then ! Dirichlet at x=L
        call Dirichlet_L(u, t, N)

        else if(BC_xL == 1) then ! Neuman at x=L
        call Neuman_L(u, N, du_N, r, dx, dt, f, t)
        endif
        
        ! Euler: x= deltaX ---- N-1
        do j = 2, N
            call euler(u, r, dt, f, t, j)
            !u(t+1, j) = 1 + sin(PI*(j-1)*dx)*exp((-1)*PI*PI*(t)*dt)
        end do
        
        end if

    end do
    print*, "Here"

    ! write T(x,t) to a file
    if (method==0) then
        open(unit=10, file = "Temp_FE.txt")
    else if (method == 1) then
        open(unit=10, file = "Temp_BE.txt")
    endif
    do t=1, M+1
        do j=1, N+1    
            write(10,'(F8.5, " ")', advance='no') u(t,j)
        end do
        write(10,*)""
    end do
    close (10)
    
    ! write T(x,t) to a file
    if (method==0) then
        open(unit=10, file = "Tsensors_FE.txt")
    else if (method == 1) then
        open(unit=10, file = "Tsensors_BE.txt")
    endif
    do t=1, M+1
        do j=8, 26, 3 ! this works only if N=40 otherwise the indexes for the 8 points must be adjusted   
            write(10,'(F8.5, " ")', advance='no') u(t,j)
        end do
        write(10,*)""
    end do
    close(10)

end program alphaHeat1D