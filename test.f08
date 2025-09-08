module test
    use grid_module
    implicit none
contains    ! Test subroutine
    subroutine test_potentials
    implicit none

    type(rgrid) :: grid

    ! Declare potential objects
    type(harmonic) :: harm_pot
    type(morse) :: morse_pot
    type(multi_gaussian) :: multi_gauss_pot
    type(SuperGauss) :: super_pot
    type(superMultiGauss) :: super_multi_pot
    type(polynomial) :: poly_pot

    real(dp), allocatable :: v_array(:)
    complex(dp), allocatable :: vz_array(:)
    real(dp) :: coeffs(0:4)

    call grid%init_rgrid(-5.0_dp, 5.0_dp, 100)
    allocate(v_array(grid%nr), vz_array(grid%nr))

    ! Test Harmonic Oscillator: V = 0.5 * k * (x - x0)^2
    harm_pot = create_harmonic(k=1.0_dp, x0=0.0_dp)
    call grid%potential_onGrid_real(harm_pot, v_array)
    call grid%printPotToFile_real(harm_pot, "harmonic.dat")

    ! Test Morse Potential: V = D * (1 - exp(-a*(x-re)))^2
    morse_pot = create_morse(D=5.0_dp, a=1.5_dp, re=2.0_dp)
    call grid%potential_onGrid_real(morse_pot, v_array)
    call grid%printPotToFile_real(morse_pot, "morse_potential.dat")

    ! Test Multi-Gaussian Barriers
    multi_gauss_pot = create_multi_gaussian(amplitude=2.0_dp, sigma=0.5_dp, &
            x0=0.0_dp, nbarrier=3, spacing=2.0_dp)
    call grid%potential_onGrid_real(multi_gauss_pot, v_array)
    call grid%printPotToFile_real(multi_gauss_pot, "multi_gaussian.dat")

    ! Test Super-Gaussian
    super_pot = create_SuperGauss(amplitude=1.5_dp, sigma=0.8_dp, &
            x0=0.0_dp, norder=4)
    call grid%potential_onGrid_real(super_pot, v_array)
    call grid%printPotToFile_real(super_pot, "super_gaussian.dat")

    ! Test Super Multi-Gaussian
    super_multi_pot = create_superMultiGauss(amplitude=1.5_dp, sigma=0.8_dp, &
            x0=0.0_dp, norder=4, nbarrier=2, spacing=3.0_dp)
    call grid%potential_onGrid_real(super_multi_pot, v_array)
    call grid%printPotToFile_real(super_multi_pot, "super_multiGaussian.dat")

    ! Test Quartic Polynomial: V = c0 + c2*(x-x0)^2 + c4*(x-x0)^4
    poly_pot = create_quartic(x0=0.0_dp, c0=0.0_dp, c2=1.0_dp, c4=0.1_dp)
    call grid%potential_onGrid_real(poly_pot, v_array)
    call grid%printPotToFile_real(poly_pot, "quadratic_poten.dat")

    ! Test General Polynomial: V = sum(coeffs(i) * (x-x0)^i)
    coeffs = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.2_dp]
    poly_pot = create_polynomial(4, 0.0_dp, coeffs)
    call grid%potential_onGrid_real(poly_pot, v_array)
    call grid%printPotToFile_real(poly_pot, "polynomial_poten.dat")

    deallocate(v_array, vz_array)
end subroutine test_potentials
end module test