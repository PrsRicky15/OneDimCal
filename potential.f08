module potential
    use physical_constants
    use grid_parameters
    implicit none

    ! Base potential Definition
    type :: potential_base
        character(len=20) :: potential_type
    contains
        procedure :: evaluate => evaluate_base
    end type potential_base

    ! Potential Definition
    type, extends(potential_base) :: harmonic_potential
        real(dp) :: k          ! force constant
        real(dp) :: x0         ! equilibrium position
    contains
        procedure :: evaluate => evaluate_harmonic
    end type harmonic_potential

    type, extends(potential_base) :: polynomial_potential
        integer :: order
        real(dp) :: x0         ! reference position
        real(dp), allocatable :: coeffs(:)  ! coefficients array (0:order)
    contains
        procedure :: evaluate => evaluate_polynomial
    end type polynomial_potential

    type, extends(potential_base) :: morse_potential
        real(dp) :: D          ! dissociation energy
        real(dp) :: a          ! width parameter
        real(dp) :: re         ! equilibrium bond length
    contains
        procedure :: evaluate => evaluate_morse
    end type morse_potential

    type, extends(potential_base) :: gaussian_potential
        real(dp) :: amplitude  ! peak height
        real(dp) :: sigma      ! width parameter
        real(dp) :: x0         ! center position
    contains
        procedure :: evaluate => evaluate_gaussian
    end type gaussian_potential

    type, extends(gaussian_potential) :: multi_gaussian
        integer :: nbarrier    ! number of barriers
        real(dp) :: spacing    ! spacing between barriers
    contains
        procedure :: evaluate => evaluate_multi_gaussian
    end type multi_gaussian

    type, extends(gaussian_potential) :: SuperGauss
        integer :: norder
    contains
        procedure :: evaluate => evaluate_SuperGauss
    end type SuperGauss

    type, extends(SuperGauss) :: superMultiGauss_poten
        integer :: nbarrier    ! number of barriers
        real(dp) :: spacing    ! spacing between barriers
    contains
        procedure :: evaluate => evaluate_superMultiGauss
    end type superMultiGauss_poten

contains

    ! Base evaluation
    function evaluate_base(this, x) result(v)
        class(potential_base), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v

        v = 0.0_dp
        print *, "Warning: Base potential evaluate called"
    end function evaluate_base

    ! Harmonic oscillator evaluation
    function evaluate_harmonic(this, x) result(v)
        class(harmonic_potential), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v

        v = 0.5_dp * this%k * (x - this%x0)**2
    end function evaluate_harmonic

    ! Polynomial evaluation
    function evaluate_polynomial(this, x) result(v)
        class(polynomial_potential), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        integer :: i
        real(dp) :: dx

        dx = x - this%x0
        v = 0.0_dp

        do i = 0, this%order
            v = v + this%coeffs(i) * dx**i
        end do
    end function evaluate_polynomial

    ! Morse potential evaluation
    function evaluate_morse(this, x) result(v)
        class(morse_potential), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        real(dp) :: exp_term

        exp_term = exp(-this%a * (x - this%re))
        v = this%D * (1.0_dp - exp_term)**2
    end function evaluate_morse

    ! Gaussian potential evaluation
    function evaluate_gaussian(this, x) result(v)
        class(gaussian_potential), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        integer :: i
        real(dp) :: center

        v = 0.0_dp
        do i = 1, this%nbarrier
            center = this%x0 + real(i-1, dp) * this%spacing
            v = v + this%amplitude * exp(-0.5_dp * ((x - center) / this%sigma)**2)
        end do
    end function evaluate_gaussian

    ! Super-Gaussian potential evaluation
    function evaluate_supergaussian(this, x) result(v)
        class(supergaussian_potential), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: v
        integer :: i
        real(dp) :: center

        v = 0.0_dp
        do i = 1, this%nbarrier
            center = this%x0 + real(i-1, dp) * this%spacing
            v = v + this%amplitude * exp(-(abs(x - center) / this%sigma)**this%order)
        end do
    end function evaluate_supergaussian

    ! Unified potential energy subroutine
    subroutine potential_onGrid(potential, v_vector)
        class(potential_base), intent(in) :: potential
        real(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: x

        do i = 1, nr
            x = r_min + real(i-1, dp) * dr
            v_vector(i) = potential%evaluate(x)
        end do

    end subroutine potential_onGrid

    ! Constructor functions
    function create_harmonic_potential(k, x0) result(pot)
        real(dp), intent(in) :: k, x0
        type(harmonic_potential) :: pot

        pot%potential_type = 'harmonic'
        pot%k = k
        pot%x0 = x0
    end function create_harmonic_potential

    function create_polynomial_potential(order, x0, coeffs) result(pot)
        integer, intent(in) :: order
        real(dp), intent(in) :: x0
        real(dp), intent(in) :: coeffs(0:order)
        type(polynomial_potential) :: pot

        pot%potential_type = 'polynomial'
        pot%order = order
        pot%x0 = x0
        allocate(pot%coeffs(0:order))
        pot%coeffs = coeffs
    end function create_polynomial_potential

    function create_morse_potential(D, a, re) result(pot)
        real(dp), intent(in) :: D, a, re
        type(morse_potential) :: pot

        pot%potential_type = 'morse'
        pot%D = D
        pot%a = a
        pot%re = re
    end function create_morse_potential

    function create_gaussian_potential(amplitude, sigma, x0, nbarrier, spacing) result(pot)
        real(dp), intent(in) :: amplitude, sigma, x0, spacing
        integer, intent(in) :: nbarrier
        type(gaussian_potential) :: pot

        pot%potential_type = 'gaussian'
        pot%amplitude = amplitude
        pot%sigma = sigma
        pot%x0 = x0
        pot%nbarrier = nbarrier
        pot%spacing = spacing
    end function create_gaussian_potential

    function create_supergaussian_potential(amplitude, sigma, x0, order, nbarrier, spacing) result(pot)
        real(dp), intent(in) :: amplitude, sigma, x0, spacing
        integer, intent(in) :: order, nbarrier
        type(supergaussian_potential) :: pot

        pot%potential_type = 'supergaussian'
        pot%amplitude = amplitude
        pot%sigma = sigma
        pot%x0 = x0
        pot%order = order
        pot%nbarrier = nbarrier
        pot%spacing = spacing
    end function create_supergaussian_potential

    ! Utility function to create quartic potential (common polynomial case)
    function create_quartic_potential(x0, c0, c2, c4) result(pot)
        real(dp), intent(in) :: x0, c0, c2, c4
        type(polynomial_potential) :: pot
        real(dp) :: coeffs(0:4)

        coeffs = [c0, 0.0_dp, c2, 0.0_dp, c4]  ! V = c0 + c2*(x-x0)^2 + c4*(x-x0)^4
        pot = create_polynomial_potential(4, x0, coeffs)
    end function create_quartic_potential

    subroutine test_potentials
        implicit none

        ! Declare potential objects
        type(harmonic_potential) :: harm_pot
        type(morse_potential) :: morse_pot
        type(gaussian_potential) :: gauss_pot
        type(supergaussian_potential) :: super_pot
        type(polynomial_potential) :: poly_pot

        real(dp), allocatable :: v_array(:)
        real(dp) :: coeffs(0:4)

        allocate(v_array(nr))

        ! Harmonic Oscillator: V = 0.5 * k * (x - x0)^2
        harm_pot = harmonic_poten(k=1.0_dp, x0=0.0_dp)
        call potential_onGrid(harm_pot, v_array)

        ! Morse Potential: V = D * (1 - exp(-a*(x-re)))^2
        morse_pot = morse_poten(D=5.0_dp, a=1.5_dp, re=2.0_dp)
        call potential_onGrid(morse_pot, v_array)

        ! Gaussian Barriers
        gauss_pot = gaussian_poten(amplitude=2.0_dp, sigma=0.5_dp, &
                x0=0.0_dp, nbarrier=3, spacing=2.0_dp)
        call potential_onGrid(gauss_pot, v_array)

        ! Super-Gaussian
        super_pot = supergauss_poten(amplitude=1.5_dp, sigma=0.8_dp, &
                x0=0.0_dp, order=4, nbarrier=2, spacing=3.0_dp)
        call potential_onGrid(super_pot, v_array)

        ! Quartic Polynomial: V = c0 + c2*(x-x0)^2 + c4*(x-x0)^4
        poly_pot = quartic_poten(x0=0.0_dp, c0=0.0_dp, c2=1.0_dp, c4=0.1_dp)
        call potential_onGrid(poly_pot, v_array)

        ! General Polynomial: V = sum(coeffs(i) * (x-x0)^i)
        coeffs = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.2_dp]  ! quadratic + quartic
        poly_pot = polynomial_poten(4, 0.0_dp, coeffs)
        call potential_onGrid(poly_pot, v_array)

    end subroutine test_potentials
end module potential