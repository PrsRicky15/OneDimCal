! Created by prash on 9/12/2025.
module const
    integer, parameter  :: dp = kind(1.0d0)
    integer, parameter  :: sp = kind(1)
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: pi_squared_over_six = pi*pi/6.0_dp
    complex(dp), parameter :: i_unit = (0.0_dp, 1.0_dp)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)
end module const