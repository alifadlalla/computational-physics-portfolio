program main
    use Settings                                      ! Module with the settings of the simulation                  
    use iso_c_binding, only: c_backspace              ! For the backspace character
    implicit none

    real(pr), allocatable, dimension(:, :) :: u, v, u0, v0   ! Current and previous states of the system                                                         
    real(pr) :: du0, dv0, delta_u0, delta_v0
    integer :: i, j, si, sj, t    
    real(pr), dimension(3,3) :: slice_u0, slice_v0                                   ! Loop variables                                            

    ! Read settings from the settings file
    call read_settings()

    !!TO BE FILED BY THE STUDENT!!
    ! Initialize the reaction-diffusion system
    ALLOCATE(u0(nx,ny), v0(nx,ny), u(nx,ny), v(nx,ny))

    u0 = 1
    v0 = 0
    print*, nx, ny
     do i= (nx/2)-15, (nx/2)+15 
         do j= (nx/2)-15, (nx/2)+15 
         u0(i,j) = 0.5
         v0(i,j) = 0.25
         end do
     end do
    !!TO BE FILED BY THE STUDENT!!

    ! Open files to save the system states
    open(unit = 10, file = output_file_u, action = 'write')
    open(unit = 11, file = output_file_v, action = 'write')

    ! Main loop (over time)
     do t = 1, nsteps

         !!TO BE FILED BY THE STUDENT!!
         ! Advance the reaction-diffusion system of one time step
         ! Loop over all elements
        do j = 1, ny
            do i = 1, nx
            ! Extract the 3x3 neighborhood with zero padding
            slice_u0 = 0.0
            slice_v0 = 0.0
            do si = -1, 1
                do sj = -1, 1
                    if (i+si >= 1 .and. i+si <= nx .and. j+sj >= 1 .and. j+sj <= ny) then
                        slice_u0(si+2, sj+2) = u0(i+si, j+sj)
                        slice_v0(si+2, sj+2) = v0(i+si, j+sj)
                    else
                        slice_u0(si+2, sj+2) = 0.0
                        slice_v0(si+2, sj+2) = 0.0
                    end if
                end do
            end do
            delta_u0 = sum(slice_u0*stencil)
            !print*, delta_u0
            delta_v0 = sum(slice_v0*stencil)
            du0 = Du*delta_u0 - u0(i,j)*v0(i,j)**2 + fr*(1-u0(i,j))
            dv0 = Dv*delta_v0 - u0(i,j)*v0(i,j)**2 - (fr+kr)*v0(i,j)
            u(i,j) = u0(i,j) + du0*dt
            v(i,j) = v0(i,j) + dv0*dt
            end do
        end do
         !!TO BE FILED BY THE STUDENT!!

         ! Write the data to the output files every sr time steps
         if (modulo(t, sr) == 0) then
             write(*, '(A30, A11 ,I10)', advance='no') repeat(c_backspace, 30), 'Time step: ', t
             ! Save the data
             write(10, *) u
             write(11, *) v
         end if

         !!TO BE FILED BY THE STUDENT!!
         ! Update the previous state of the system
         u0 = u
         v0 = v
         !!TO BE FILED BY THE STUDENT!!

     end do

    !TO BE FILED BY THE STUDENT!!
    ! Something must be done here...
    DEALLOCATE(u0, v0, u, v)
    !TO BE FILED BY THE STUDENT!!

    close(10)
    close(11)

    write(*, *) 'Done!'

end program main