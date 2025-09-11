module potential
    implicit none

    integer, parameter  :: dp = kind(1.0d0)
    integer, parameter  :: sp = kind(1)
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: pi_squared_over_six = pi*pi/6.0_dp
    complex(dp), parameter :: i_unit = (0.0_dp, 1.0_dp)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)

    ! Base potential Definition
    type, abstract :: potential_base
        character(len=40) :: potential_type
    contains
        procedure(evaluate_real), deferred :: evaluate_real
        procedure(evaluate_complex), deferred :: evaluate_complex
        generic :: evaluate => evaluate_real, evaluate_complex
    end type potential_base

    abstract interface
        function evaluate_real(this, x) result(v)
            import :: potential_base, dp
            class(potential_base), intent(in) :: this
            real(dp), intent(in) :: x
            real(dp) :: v
        end function evaluate_real

        function evaluate_complex(this, x) result(v)
            import :: potential_base, dp
            class(potential_base), intent(in) :: this
            complex(dp), intent(in) :: x
            complex(dp) :: v
        end function evaluate_complex
    end interface

    ! Harmonic potential
    type, extends(potential_base) :: harmonic
        real(dp) :: k          ! force constant
        real(dp) :: x0         ! equilibrium position
    contains
        procedure :: evaluate_real => evalHarmonic_real
        procedure :: evaluate_complex => evalHarmonic_complex
    end type harmonic

    ! Polynomial potential
    type, extends(potential_base) :: polynomial
        integer :: order
        real(dp) :: x0         ! reference position
        real(dp), allocatable :: coeffs(:)  ! coefficients array (0:order)
    contains
        procedure :: evaluate_real => evalPolynomial_real
        procedure :: evaluate_complex => evalPolynomial_complex
    end type polynomial

    ! Morse potential
    type, extends(potential_base) :: morse
        real(dp) :: D          ! dissociation energy
        real(dp) :: a          ! width parameter
        real(dp) :: re         ! equilibrium bond length
    contains
        procedure :: evaluate_real => evalMorse_real
        procedure :: evaluate_complex => evalMorse_complex
    end type morse

    ! Gaussian potential - single barrier
    type, extends(potential_base) :: gaussian
        real(dp) :: amplitude  ! peak height
        real(dp) :: sigma      ! width parameter
        real(dp) :: x0         ! center position
    contains
        procedure :: evaluate_real => evalGauss_real
        procedure :: evaluate_complex => evalGauss_complex
    end type gaussian

    ! Multi-Gaussian potential
    type, extends(gaussian) :: multi_gaussian
        integer :: nbarrier    ! number of barriers
        real(dp) :: spacing    ! spacing between barriers
    contains
        procedure :: evaluate_real => evalMultiGauss_real
        procedure :: evaluate_complex => evalMultiGauss_complex
    end type multi_gaussian

    ! Super-Gaussian potential
    type, extends(gaussian) :: SuperGauss
        integer :: norder
    contains
        procedure :: evaluate_real => evalSuperGauss_real
        procedure :: evaluate_complex => evalSuperGauss_complex
    end type SuperGauss

    ! Super Multi-Gaussian potential
    type, extends(SuperGauss) :: superMultiGauss
        integer :: nbarrier    ! number of barriers
        real(dp) :: spacing    ! spacing between barriers
    contains
        procedure :: evaluate_real => evalSuperMultiGauss_real
        procedure :: evaluate_complex => evalSuperMultiGauss_complex
    end type superMultiGauss

    type, extends(potential_base) :: ResonancePoten
        real(dp) :: width = 0.1_dp
        real(dp) :: shift = -0.8_dp
        real(dp) :: k0 = 1.0_dp
    contains
        procedure :: evaluate_real => evalResonancePoten_real
        procedure :: evaluate_complex => evalResonancePoten_complex
    end type ResonancePoten

contains

    ! Harmonic oscillator evaluation
    function evalHarmonic_real(this, x) result(v)
        class(harmonic), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v

        v = 0.5_dp * this%k * (x - this%x0)**2
    end function evalHarmonic_real

    function evalHarmonic_complex(this, x) result(v)
        class(harmonic), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v

        v = 0.5_dp * this%k * (x - this%x0)**2
    end function evalHarmonic_complex

    ! Polynomial evaluation
    function evalPolynomial_real(this, x) result(v)
        class(polynomial), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        integer :: i
        real(dp) :: dx

        dx = x - this%x0
        v = 0.0_dp

        ! Use Horner's method for efficiency
        if (this%order >= 0) then
            v = this%coeffs(this%order)
            do i = this%order - 1, 0, -1
                v = v * dx + this%coeffs(i)
            end do
        end if
    end function evalPolynomial_real

    function evalPolynomial_complex(this, x) result(v)
        class(polynomial), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v
        integer :: i
        complex(dp) :: dx

        dx = x - this%x0
        v = 0.0_dp

        ! Use Horner's method for efficiency
        if (this%order >= 0) then
            v = this%coeffs(this%order)
            do i = this%order - 1, 0, -1
                v = v * dx + this%coeffs(i)
            end do
        end if
    end function evalPolynomial_complex

    ! Morse potential evaluation
    function evalMorse_real(this, x) result(v)
        class(morse), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        real(dp) :: exp_term

        exp_term = exp(-this%a * (x - this%re))
        v = this%D * (1.0_dp - exp_term)**2
    end function evalMorse_real

    function evalMorse_complex(this, x) result(v)
        class(morse), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v
        complex(dp) :: exp_term

        exp_term = exp(-this%a * (x - this%re))
        v = this%D * (1.0_dp - exp_term)**2
    end function evalMorse_complex

    ! Single Gaussian potential evaluation
    function evalGauss_real(this, x) result(v)
        class(gaussian), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v

        v = this%amplitude * exp(-0.5_dp * ((x - this%x0) / this%sigma)**2)
    end function evalGauss_real

    function evalGauss_complex(this, x) result(v)
        class(gaussian), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v

        v = this%amplitude * exp(-0.5_dp * ((x - this%x0) / this%sigma)**2)
    end function evalGauss_complex

    ! Multi-Gaussian potential evaluation
    function evalMultiGauss_real(this, x) result(v)
        class(multi_gaussian), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        integer :: i
        real(dp) :: center

        v = 0.0_dp
        do i = 1, this%nbarrier
            center = this%x0 + real(i-1, dp) * this%spacing
            v = v + this%amplitude * exp(-0.5_dp * ((x - center) / this%sigma)**2)
        end do
    end function evalMultiGauss_real

    function evalMultiGauss_complex(this, x) result(v)
        class(multi_gaussian), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v
        integer :: i
        complex(dp) :: center

        v = 0.0_dp
        do i = 1, this%nbarrier
            center = this%x0 + real(i-1, dp) * this%spacing
            v = v + this%amplitude * exp(-0.5_dp * ((x - center) / this%sigma)**2)
        end do
    end function evalMultiGauss_complex

    ! Super-Gaussian potential evaluation
    function evalSuperGauss_real(this, x) result(v)
        class(SuperGauss), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v

        v = this%amplitude * exp(-(abs(x - this%x0) / this%sigma)**this%norder)
    end function evalSuperGauss_real

    function evalSuperGauss_complex(this, x) result(v)
        class(SuperGauss), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v

        v = this%amplitude * exp(-(abs(x - this%x0) / this%sigma)**this%norder)
    end function evalSuperGauss_complex

    ! Super Multi-Gaussian potential evaluation
    function evalSuperMultiGauss_real(this, x) result(v)
        class(superMultiGauss), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        integer :: i
        real(dp) :: center

        v = 0.0_dp
        do i = 1, this%nbarrier
            center = this%x0 + real(i-1, dp) * this%spacing
            v = v + this%amplitude * exp(-(abs(x - center) / this%sigma)**this%norder)
        end do
    end function evalSuperMultiGauss_real

    function evalSuperMultiGauss_complex(this, x) result(v)
        class(superMultiGauss), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v
        integer :: i
        complex(dp) :: center

        v = 0.0_dp
        do i = 1, this%nbarrier
            center = this%x0 + real(i-1, dp) * this%spacing
            v = v + this%amplitude * exp(-(abs(x - center) / this%sigma)**this%norder)
        end do
    end function evalSuperMultiGauss_complex

    function evalResonancePoten_real(this, x) result(v)
        class(ResonancePoten), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v

        v = dexp(-this%width * x * x) * (x * x * 0.5_dp * this%k0 + this%shift)
    end function evalResonancePoten_real

    function evalResonancePoten_complex(this, x) result(v)
        class(ResonancePoten), intent(in) :: this
        complex(dp), intent(in) :: x
        complex(dp) :: v

        v = cdexp(-this%width * x * x) * (x * x * 0.5_dp * this%k0 + this%shift)
    end function evalResonancePoten_complex

    function create_resonancePoten() result(pot)
        type(ResonancePoten) :: pot
        pot%potential_type = 'One-Bound_Resonances'
    end function create_resonancePoten

    ! Constructor functions
    function create_harmonic(k, x0) result(pot)
        real(dp), intent(in) :: k, x0
        type(harmonic) :: pot

        pot%potential_type = 'harmonic'
        pot%k = k
        pot%x0 = x0
    end function create_harmonic

    function create_polynomial(order, x0, coeffs) result(pot)
        integer, intent(in) :: order
        real(dp), intent(in) :: x0
        real(dp), intent(in) :: coeffs(0:order)
        type(polynomial) :: pot

        pot%potential_type = 'polynomial'
        pot%order = order
        pot%x0 = x0
        allocate(pot%coeffs(0:order))
        pot%coeffs = coeffs
    end function create_polynomial

    function create_morse(D, a, re) result(pot)
        real(dp), intent(in) :: D, a, re
        type(morse) :: pot

        pot%potential_type = 'morse'
        pot%D = D
        pot%a = a
        pot%re = re
    end function create_morse

    function create_gaussian(amplitude, sigma, x0) result(pot)
        real(dp), intent(in) :: amplitude, sigma, x0
        type(gaussian) :: pot

        pot%potential_type = 'gaussian'
        pot%amplitude = amplitude
        pot%sigma = sigma
        pot%x0 = x0
    end function create_gaussian

    function create_multi_gaussian(amplitude, sigma, x0, nbarrier, spacing) result(pot)
        real(dp), intent(in) :: amplitude, sigma, x0, spacing
        integer, intent(in) :: nbarrier
        type(multi_gaussian) :: pot

        pot%potential_type = 'multi_gaussian'
        pot%amplitude = amplitude
        pot%sigma = sigma
        pot%x0 = x0
        pot%nbarrier = nbarrier
        pot%spacing = spacing
    end function create_multi_gaussian

    function create_SuperGauss(amplitude, sigma, x0, norder) result(pot)
        real(dp), intent(in) :: amplitude, sigma, x0
        integer, intent(in) :: norder
        type(SuperGauss) :: pot

        pot%potential_type = 'supergaussian'
        pot%amplitude = amplitude
        pot%sigma = sigma
        pot%x0 = x0
        pot%norder = norder
    end function create_SuperGauss

    function create_superMultiGauss(amplitude, sigma, x0, norder, nbarrier, spacing) result(pot)
        real(dp), intent(in) :: amplitude, sigma, x0, spacing
        integer, intent(in) :: norder, nbarrier
        type(superMultiGauss) :: pot

        pot%potential_type = 'super_multi_gauss'
        pot%amplitude = amplitude
        pot%sigma = sigma
        pot%x0 = x0
        pot%norder = norder
        pot%nbarrier = nbarrier
        pot%spacing = spacing
    end function create_superMultiGauss

    function create_quartic(x0, c0, c2, c4) result(pot)
        real(dp), intent(in) :: x0, c0, c2, c4
        type(polynomial) :: pot
        real(dp) :: coeffs(0:4)

        coeffs = [c0, 0.0_dp, c2, 0.0_dp, c4]  ! V = c0 + c2*(x-x0)^2 + c4*(x-x0)^4
        pot = create_polynomial(4, x0, coeffs)
    end function create_quartic
end module potential