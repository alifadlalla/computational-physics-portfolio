module functions
    implicit none
    
contains
    
    !--------------Forward Euler---------------------!
    subroutine Dirichlet_0 (u, i)
        real(kind=8), dimension(:,:) :: u
        integer(kind=8) :: i
        u(i+1, 1) = u(i, 1)
    end subroutine

    subroutine Dirichlet_L (u, i, N)
        real(kind=8), dimension(:,:) :: u
        integer(kind=8) :: i, N
        u(i+1, N+1) = u(i, N+1)
    end subroutine

    subroutine Neuman_0 (u, du_0, r, delta_x, delta_t, f, i)
        real(kind=8), dimension(:,:) :: u, f
        real(kind=8) :: r, du_0, delta_x, delta_t
        integer(kind=8) :: i
        u(i+1, 1) = u(i, 1) + r*(2*u(i, 2) - 2*du_0*delta_x - 2*u(i, 1)) + f(i,1)*delta_t
    end subroutine

    subroutine Neuman_L (u, N, du_N, r, delta_x, delta_t, f, i)
        real(kind=8), dimension(:,:) :: u, f
        real(kind=8) :: r, du_N, delta_x, delta_t
        integer(kind=8) :: i, N
        u(i+1, N+1) = u(i, N+1) + r*(2*u(i, N) + 2*du_N*delta_x - 2*u(i, N+1)) + f(i,N+1)*delta_t
    end subroutine

    subroutine euler (u, r, delta_t, f, i, j)
        real(kind=8), dimension(:,:) :: u, f
        real(kind=8) :: r, delta_t
        integer(kind=8) :: i, j
        u(i+1, j) = u(i, j) + r*(u(i, j-1) - 2*u(i, j) + u(i, j+1)) + f(i,j)*delta_t
    end subroutine
    !--------------End Forward Euler---------------------!

    !--------------Backward Euler------------------------!
    
    subroutine thomas(a, b, c, d, x, n)
        implicit none
        ! a,c - lower and upper diagonal of the tridiagonal matrix
        ! b - main diagonal of the tridiagonal matrix
        ! d - right-hand side of the system of equations
        ! x - solution to the system of equations
        ! n - dimension of the tridiagonal matrix
        real(kind=8), dimension(n-1) :: a, c
        real(kind=8), dimension(n) :: b, d, x
        integer(kind=8) :: i, n
        real(kind=8) :: w
        ! forward elimination step
        do i = 2, n
            ! calculate the ratio between the lower and main diagonals
            w = a(i-1)*1.0 / b(i-1)
            ! update the main diagonal
            b(i) = b(i) - w * c(i-1)
            ! update the right-hand side
            d(i) = d(i) - w * d(i-1)
        end do
        ! back substitution step
        ! calculate the solution for the last equation
        x(n) = d(n)*1.0 / b(n)
        do i = n-1, 1, -1
            ! calculate the solution for the remaining equations
            x(i) = (d(i) - c(i) * x(i+1))*1.0 / b(i)
        end do
    end subroutine thomas
    
    subroutine Backward_Euler(r,N,dx,du_0,du_1,du_N,u,t,cdu,BC_x0, BC_xL)
        
        real(kind=8) :: r, dx,du_0,du_1,du_N, cdu
        real(kind=8), dimension(:,:) :: u
        real(kind=8), dimension(N-1) :: b, d, x
        real(kind=8), dimension(N-2) :: a, c
        integer(kind=8) :: N, i, t, BC_x0, BC_xL

        ! Tridiagonal matrix
        a = -r*1.0 ! left diagonal 
        
        ! main diagonal
        if (BC_x0 == 0) then  
            b(1) = 1+ 2*r*1.0
        else if(BC_x0 == 1) then
            b(1) = 1+ r*1.0    
        end if

        do i = 2, N-2
            b(i) = 1 + 2*r*1.0
        end do
        if (BC_xL == 0) then  
            b(N-1) = 1+ 2*r*1.0
        else
            b(N-1) = 1+ r*1.0    
        end if
        ! right diagonal
        c = -r*1.0

        ! right hand side
        if (BC_x0 == 0) then  
            d(1) = u(t,2) + r*u(t,1)
        else if (BC_x0 == 1) then
            d(1) = u(t,2) - r*du_1*dx
        end if
        
        do i = 2, N-2
            d(i) = u(t, i+1)
        end do

        if (BC_xL == 0) then  
            d(N-1) = u(t,N) + r*u(t,N+1) 
        else if (BC_xL == 1) then 
            du_N = 0
            d(N-1) = u(t,N) + r*du_N*dx
        else if (BC_xL == 2) then
            du_N = cdu*(u(t,N+1) - 1)
            d(N-1) = u(t,N) + r*du_N*dx
        end if

        ! solving the system
        call thomas(a,b,c,d,x,N-1)
        do i=2,N
            u(t+1, i) = x(i-1)
        enddo
        
        if (BC_x0 == 0) then  
            u(t+1, 1) = u(t, 1) ! Dirichlet at 0
        else if (BC_x0 == 1) then
            u(t+1, 1)   = (u(t, 1)   +2*r*(u(t+1,2) - du_0*dx))*1.0/(1+2*r)   ! Neuman at 0    
            !u(t+1, 1)   = u(t, 1) - du_1*dx     ! Neuman at 0
        end if
        
        if (BC_xL == 0) then  
            u(t+1, N+1) = u(t, N+1)  ! Dirchlet at L
        else if (BC_xL == 1) then
            !u(t+1, N+1) = (u(t, N+1) +2*r*(u(t+1,N) + du_N*dx))*1.0/(1+2*r) ! Neuman at L
            u(t+1, N+1) = u(t+1, N) + du_N*dx ! Neuman at L
        else if (BC_xL == 2) then
            u(t+1, N+1) = (u(t, N+1) + 2*r*cdu*dx +2*r*u(t+1,N))*1.0/(1+ 2*r*(1+cdu*dx))  ! Robin at L
        end if

    end subroutine Backward_Euler


end module functions