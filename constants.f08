module physical_constants
    implicit none
    
    integer, parameter  :: dp = kind(1.0d0)
    integer, parameter  :: sp = kind(1)
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: pi = 3.141592653589793238460_dp
    real(dp), parameter :: pi_squared_over_six = pi*pi/6.0_dp
    complex(dp), parameter :: i_unit = (0.0_dp, 1.0_dp)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)
end module physical_constants

module grid_parameters
    use physical_constants, only: dp
    implicit none

    integer, parameter :: nr = 401
    real(dp), parameter :: r_min = -40.0_dp
    real(dp), parameter :: r_max = 40.0_dp
    real(dp), parameter :: dr = (r_max - r_min) / real(nr - 1, dp)

end module grid_parameters
