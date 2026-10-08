module settings
    implicit none

    ! Precision selector
    integer, public, parameter :: dp = SELECTED_REAL_KIND(P=8)
    
    ! Number of particles
    integer, public :: N

    ! Number of time step
    integer, public :: nt

    ! Time step
    real(kind=dp), public :: dt

    ! Mass of each particle
    real(kind=dp), public :: m

    ! Minimal distance between two particles (smoothing lenght)
    real(kind=dp), public :: epsilon

    ! Pi
    real(kind=dp), public, parameter :: pi = acos(-1.0)

    ! Gravitational constant
    real(kind=dp), public, parameter :: gg = 1.0

    contains

    subroutine readsettings()

        integer :: fend
        character(len=20) :: settingname
        real(kind=dp) :: value
    
        fend = 0
        open(unit = 100, file = '../input.yaml', action = 'read')

        do while (fend >= 0)
            read(100, *, IOSTAT=fend) settingname, value
            
            select case (trim(settingname))
                case ("npart:")
                    N = int(value)
                    m = 1/value

                case ("nstep:")
                    nt = int(value)

                case ("dt:")
                    dt = value

                case ("epsilon:")
                    epsilon = value
            end select
        end do
        
        close(100)
    end subroutine readsettings

end module settings