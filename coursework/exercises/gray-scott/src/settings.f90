module Settings
    ! This module contains all the settings of the simulation

    implicit none

    ! Floating point precision
    integer, public, parameter :: pr = SELECTED_REAL_KIND(P=8)

    ! Grid size nx
    integer :: nx = 1024

    ! Grid size ny
    integer :: ny = 1024

    ! Number of time steps
    integer :: nsteps = 2

    ! Time step
    real(pr) :: dt = 1_pr

    ! Save rate
    integer :: sr = 20

    ! Diffusion rate of u
    real(pr) :: Du = 0.1_pr

    ! Diffusion rate of v
    real(pr) :: Dv = 0.05_pr

    ! Feed rate
    real(pr) :: fr = 0.014_pr

    ! Kill rate
    real(pr) :: kr = 0.054_pr

    ! Moore neighbourhood 
    real(pr), dimension(3, 3), parameter :: stencil = reshape([real(pr) :: &
        1 , 1, 1, &
        1,  -8, 1, &
        1,  1, 1], [3, 3])

    ! output files
    character(len=*), parameter :: output_file_u = 'output_u.txt'
    character(len=*), parameter :: output_file_v = 'output_v.txt'

    ! Gray-Scott parameters
    namelist /GRAYSCOTT/ nx, ny, nsteps, sr, dt, Du, Dv, fr, kr

    contains

    ! Read settings from file
    subroutine read_settings()
        character(len=*), parameter :: filename = '../settings.nml'
        integer :: rc, fu

        ! Open settings file
        open(file=filename, newunit=fu, action='read')

        ! Read settings from file
        read(nml=GRAYSCOTT, iostat=rc, unit=fu)

        ! Close settings file
        close(unit=fu)

    end subroutine read_settings

end module Settings