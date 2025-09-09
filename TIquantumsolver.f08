module quantum_solver
    use grid_module
    use matrix_algebra
    implicit none
    private

    public :: kinetic_energy_matrix, potential_energy_vector
    public :: real_Hamiltonian, print_vector, printPoten_grid, printEigvecs

    type, abstract :: kinetic_base
        class(rgrid) :: grid
        real(dp), allocatable :: real_mat(:,:)
        complex(dp), allocatable :: complex_mat(:,:)
    contains
        generic :: init => initial_real_ke, initial_complex_ke
        procedure :: initial_real_ke, initial_complex_ke
    end type kinetic_base

    abstract interface
        subroutine getKe(this)
            import kinetic_base, dp
            class(kinetic_base), intent(in) :: this
            integer :: i, j
        end subroutine getKe

        subroutine getKe(this, theta)
            import kinetic_base, dp
            class(kinetic_base), intent(in) :: this
            real(dp), intent(in) :: theta
            integer :: i, j
        end subroutine getKe
    end interface

    interface print_vector
        procedure printVec_real, printVec_complex
    end interface print_vector

    interface printPoten_grid
        procedure printPoten_real_grid, printPoten_complex_grid
    end interface printPoten_grid

    interface printEigvecs
        procedure printEigvecs_real, printEigvecs_complex
    end interface printEigvecs

contains
    
    subroutine kinetic_energy_matrix(t_matrix)
        real(dp), intent(out) :: t_matrix(:,:)
        integer :: i, j
        
        do i = 1, nr
            t_matrix(i, i) = pi_squared_over_six
            do j = 1, i-1
                t_matrix(i, j) = ((-one)**(i-j)) / real((i-j)**2, dp)
                t_matrix(j, i) = t_matrix(i, j)
            end do
        end do
        
        t_matrix = t_matrix / (dr * dr)
    end subroutine kinetic_energy_matrix

    subroutine printPoten_real_grid(filename, vVec)
        character(*), intent(in) :: filename
        real(dp), intent(in) :: vVec(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, size(vVec)
                write(11,*) r_min + (ir-1) * dr, vVec(ir)
            enddo
        close(11)
    end subroutine printPoten_real_grid

    subroutine printPoten_complex_grid(filename, vVec)
        character(*), intent(in) :: filename
        complex(dp), intent(in) :: vVec(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, size(vVec)
                write(11,*) r_min + (ir-1) * dr, real(vVec(ir), dp), aimag(vVec(ir))
            enddo
        close(11)
    end subroutine printPoten_complex_grid

    subroutine printVec_real(filename, eigvals)
        character(*), intent(in) :: filename
        real(dp), intent(in) :: eigvals(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, nr
                write(11,*) eigvals(ir)
            enddo
        close(11)
    end subroutine printVec_real

    subroutine printVec_complex(filename, eigvals)
        character(*), intent(in) :: filename
        complex(dp), intent(in) :: eigvals(:)
        integer :: ir

        open(11, file=filename, status='replace')
        do ir = 1, nr
            write(11,*) real(eigvals(ir), dp), aimag(eigvals(ir))
        enddo
        close(11)
    end subroutine printVec_complex

    subroutine printEigvecs_real(filename, eigvecs, evals, nstate)
        character(*), intent(in) :: filename
        integer, intent(in) :: nstate
        real(dp), intent(in) :: eigvecs(:,:), evals(:)
        integer :: ir, istate

        open(11, file=filename, status='replace')
            do ir = 1, nr

                write(11,*) r_min + (ir - 1) * dr, &
                    (eigvecs(ir, istate) + evals(istate), istate = 1, nstate)
            enddo
        close(11)
    end subroutine printEigvecs_real

    subroutine printEigvecs_complex(filename, eigvecs, evals, nstate)
        character(*), intent(in) :: filename
        integer, intent(in) :: nstate
        complex(dp), intent(in) :: eigvecs(:,:)
        real(dp), intent(in) :: evals(:)
        integer :: ir, istate

        open(11, file=filename, status='replace')
        do ir = 1, nr
            write(11,*) r_min + (ir - 1) * dr, &
                    (real(eigvecs(ir, istate), dp) + evals(istate), &
            aimag(eigvecs(ir, istate)) + evals(istate), istate = 1, nstate)
        enddo
        close(11)
    end subroutine printEigvecs_complex

    subroutine real_Hamiltonian(potential)
        class(potential_base), intent(in) :: potential
        real(dp), allocatable :: hmat(:, :), vVec(:)
        real(dp) :: eigval(nr)
        integer :: ir

        allocate(hmat(nr, nr), vVec(nr))
        call kinetic_energy_matrix(hmat)

        do ir = 1, nr
            hmat(ir,ir) = hmat(ir,ir) + vVec(ir)
        enddo

        call solve_Hamil(hmat, eigval)
        call print_vector("gaussEigval.dat", eigval)
        call printEigvecs("gauss_eigvecs.dat", hmat, eigval, 20)
        deallocate(hmat, vVec)
    end subroutine real_Hamiltonian

    subroutine complex_Hamiltonian(potential, theta)
        class(potential_base), intent(in) :: potential
        real(dp), intent(in) :: theta
        real(dp), allocatable :: tmat(:,:)
        complex(dp), allocatable :: hmat(:, :), vVec(:)
        complex(dp), allocatable :: vecsL(:, :), vecsR(:)
        complex(dp) :: eigval(nr)
        integer :: ir

        allocate(tmat)
        call kinetic_energy_matrix(tmat)

        allocate(hmat(nr, nr), vVec(nr))

        do ir = 1, nr
            hmat(ir,ir) = hmat(ir,ir) + vVec(ir)
        enddo

        call solve_Hamil(hmat, eigval)
        call print_vector("gaussEigval.dat", eigval)
        call printEigvecs("gauss_eigvecs.dat", hmat, eigval, 20)
        deallocate(hmat, vVec)
    end subroutine complex_Hamiltonian

    subroutine cap_trajectory(lmda_max, nlmda, Vcap)
        use physical_constants
        use grid_parameters
        use quantum_solver
        use matrix_algebra
        implicit none

        real(dp), intent(in) :: lmda_max, Vcap(:)
        integer(sp), intent(in) :: nlmda

        real(dp), allocatable :: Hmatr(:,:)
        complex(dp), allocatable :: Vx(:), Hmatc(:,:)
        complex(dp), allocatable :: eigenvalues(:)
        complex(dp), allocatable :: left_vecs(:,:), right_vecs(:,:)
        complex(dp), allocatable :: prev_right_vecs(:,:)

        real(dp) :: lambda, dlambda
        complex(dp) :: phase_factor
        integer :: i, j, ilmda, ndim

        dlambda = lmda_max/(nlmda - 1._dp)
        ndim = size(Vcap)

        allocate(Hmatr(ndim, ndim), Hmatc(ndim, ndim))
        allocate(eigenvalues(ndim))
        allocate(left_vecs(ndim, ndim), right_vecs(ndim, ndim))
        allocate(prev_right_vecs(ndim, ndim))

        ! Calculate kinetic energy matrix (constant for all theta)
        call kinetic_energy_matrix(Hmatr)

        write(*,'(a)') 'Starting eigenvalue calculation...'

        ! Loop over theta values
        open(22, file='eigenvalues_vs_theta.dat', status='replace')
        write(22, '(a)') '# Lambda-dependent eigenvalues: Real and Imaginary parts'

        do ilmda = 1, nlmda
            lambda = real(ilmda - 1, dp) * dlambda

            write(*,'(a,i0,a,f8.4)') 'Theta step ', ilmda, ', theta = ', lambda

            call solve_hamil(Hmatc, eigenvalues, left_vecs, right_vecs)

            if (ilmda == 1) then
                call sort_eigenvalues_real(eigenvalues, left_vecs, right_vecs)
                call print_vector('real_eigenvals.dat', eigenvalues)
            else
                call sort_eigenvalues_overlap(eigenvalues, left_vecs, right_vecs, &
                        prev_right_vecs)
            end if

            prev_right_vecs = right_vecs
            write(22, '(*(es16.8,1x))') (real(eigenvalues(i), dp), &
                    aimag(eigenvalues(i)), i = 1, nr)
        end do
        close(22)

        write(*,'(a)') 'Calculation completed successfully!'
        write(*,'(a)') 'Output files:'
        write(*,'(a)') '  - initial_eigenvalues.dat: First theta eigenvalues'
        write(*,'(a)') '  - eigenvalues_vs_theta.dat: All eigenvalues vs theta'

        deallocate(Vx, Hmatc, Hmatr)
        deallocate(eigenvalues)
        deallocate(left_vecs, right_vecs)
        deallocate(prev_right_vecs)

    end subroutine cap_trajectory

    subroutine thetatraj_complex_scale(theta_max, ntheta)
        use physical_constants
        use grid_parameters
        use quantum_solver
        use matrix_algebra
        implicit none

        real(dp), intent(in) :: theta_max
        integer(sp), intent(in) :: ntheta
        real(dp), allocatable :: t_matrix(:,:)
        complex(dp), allocatable :: v_vector(:), h_matrix(:,:)
        complex(dp), allocatable :: eigenvalues(:)
        complex(dp), allocatable :: left_vecs(:,:), right_vecs(:,:)
        complex(dp), allocatable :: prev_right_vecs(:,:)

        real(dp) :: theta, dtheta
        complex(dp) :: phase_factor
        integer :: i, j, ith

        dtheta = theta_max / real(ntheta - 1, dp)

        allocate(t_matrix(nr, nr))
        allocate(v_vector(nr), h_matrix(nr, nr))
        allocate(eigenvalues(nr))
        allocate(left_vecs(nr, nr), right_vecs(nr, nr))
        allocate(prev_right_vecs(nr, nr))

        ! Calculate kinetic energy matrix (constant for all theta)
        call kinetic_energy_matrix(t_matrix)

        write(*,'(a)') 'Starting eigenvalue calculation...'

        ! Loop over theta values
        open(22, file='eigenvalues_vs_theta.dat', status='replace')
        write(22, '(a)') '# Theta-dependent eigenvalues: Real and Imaginary parts'

        do ith = 1, ntheta
            theta = real(ith - 1, dp) * dtheta

            write(*,'(a,i0,a,f8.4)') 'Theta step ', ith, ', theta = ', theta

            call potential_energy_vector(theta, v_vector)

            phase_factor = exp(-2.0_dp * i_unit * theta)
            do i = 1, nr
                h_matrix(i, i) = t_matrix(i, i) * phase_factor + v_vector(i)
                do j = 1, i-1
                    h_matrix(i, j) = t_matrix(i, j) * phase_factor
                    h_matrix(j, i) = h_matrix(i, j)
                end do
            end do

            call solve_hamil(h_matrix, eigenvalues, left_vecs, right_vecs)

            if (ith == 1) then
                call sort_eigenvalues_real(eigenvalues, left_vecs, right_vecs)
                call print_vector('real_eigenvals.dat', eigenvalues)
            else
                call sort_eigenvalues_overlap(eigenvalues, left_vecs, right_vecs, &
                        prev_right_vecs)
            end if

            prev_right_vecs = right_vecs
            write(22, '(*(es16.8,1x))') (real(eigenvalues(i), dp), &
                    aimag(eigenvalues(i)), i = 1, nr)
        end do
        close(22)

        write(*,'(a)') 'Calculation completed successfully!'
        write(*,'(a)') 'Output files:'
        write(*,'(a)') '  - initial_eigenvalues.dat: First theta eigenvalues'
        write(*,'(a)') '  - eigenvalues_vs_theta.dat: All eigenvalues vs theta'

        deallocate(t_matrix, v_vector, h_matrix)
        deallocate(eigenvalues)
        deallocate(left_vecs, right_vecs)
        deallocate(prev_right_vecs)

    end subroutine thetatraj_complex_scale
end module quantum_solver
