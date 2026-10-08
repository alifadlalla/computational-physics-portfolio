!################################################################
! K-means exercise for the M2-CompuPhys Machine Learning lectures
! Author: David Cornu - david.cornu@obspm.fr
!################################################################


module utils
	implicit none

contains
	!This Function return a distance in ndim dimension
	!between two points (arguments are table with all the dimensions) 
	function dist(dat, cent, ndim)
		real :: dist
		real, intent(in) :: dat(:), cent(:)
		integer, intent(in) :: ndim
		integer :: i

		dist = 0.0
		do i = 1, ndim
			dist = dist + (dat(i)-cent(i))**2
		end do
	end function dist

end module utils



program kmeans

	use utils
	implicit none

	!Usefull variables, adding some will be necessary
	!for your own implementation of the algorithm
	real, allocatable :: input_data(:,:), centers(:,:)
	integer, allocatable :: nb_points_per_center(:)
	integer :: nb_dim, nb_data
	real :: rand
	integer :: nb_k
	integer :: i, j

	nb_k = 4

	!This entry file must be edited to change to other
	!number of dimension. The code must be re-compiled !
	open(unit = 10, file="kmeans_input_file_2d.dat")

	!Read the dimensions of the data in the file
	read(10, *) nb_dim, nb_data
	allocate(input_data(nb_data, nb_dim))

	write(*,*) nb_dim, nb_data

	!Load all the data
	do i = 1, nb_dim
		read(10, *) input_data(:, i)
	end do

	close(10)

	!allocate the tables according to the dimension
	!gave in the input file
	allocate(centers(nb_k, nb_dim))
	allocate(nb_points_per_center(nb_k))

	!the origin of the centers are selected randomly
	!to the position of some points in the dataset
	do i=1, nb_k
		call random_number(rand)
		centers(i, :) = input_data(int(rand*nb_data),:)
	end do

	!################################################################
	!     Main loop, until the new centers do not move anymore
	!################################################################
	logical :: converged
	integer :: iter
	iter = 0
	converged = .false.

	do while (.not. converged)
		iter = iter + 1
		print *, 'Iteration:', iter


		!################################################################
		!         Association phase, loop on the data points
		!################################################################
		
!################################################################
    !         Association phase, loop on the data points
    !################################################################

    do i = 1, nb_data
        real :: min_dist
        integer :: nearest_center

        min_dist = 1.0e20  ! Set an initial large value
        nearest_center = 0

        ! Loop over centers to find the nearest one to the current data point
        do j = 1, nb_k
            real :: distance
            distance = dist(input_data(i, :), centers(j, :), nb_dim)

            if (distance < min_dist) then
                min_dist = distance
                nearest_center = j
            end if
        end do

        ! Assign the data point to the nearest center
        ! (You may want to keep track of the number of points per center for the update phase)
        nb_points_per_center(nearest_center) = nb_points_per_center(nearest_center) + 1
    end do

		!################################################################
		!           Update phase, calculate the new centers
		!################################################################




		

	!################################################################
	!      Save the ending centroid position for visualisation
	!################################################################	
	open(unit = 10, file="kmeans_output_2d.dat")
	
	write(10,*) nb_dim, nb_k
	
	do i = 1, nb_k
		write (10,*) centers(i, :)
	end do
	
	

end program kmeans




