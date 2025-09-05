program modern_quantum_solver
!    use physical_constants, only: dp
    use quantum_solver
    use potential
    implicit none
!    real(dp), parameter :: sigma = 10._dp, coef = 5
!    real(dp), parameter :: displace = 10._dp
!    integer, parameter :: nbarrier =  2, order = 12
!    call Gauss_real(sigma, coef, displace, nbarrier)
!    call superGauss_(sigma, coef, order, displace, nbarrier)
    call test_potentials
end program modern_quantum_solver