module quantum_solver
    use derivative_operators
    implicit none

    type, extends(kinetic_base) :: Hamiltonian
        real(dp), allocatable :: hmatr(:,:), eigr(:), eigvecs(:,:)
        complex(dp), allocatable :: hmatz(:,:), eigz(:)
        complex(dp), allocatable :: eigvecsL(:,:), eigvecsR(:,:)
    contains

        ! Generic variable name for subroutine call
        generic :: construct_hamil => construct_rhamil, construct_zhamil
        generic :: expHOp => expHOp_real, expHOp_real_mut, expHOp_complex, expHOp_complex_mut
        generic :: get_hamil => getH_real, getH_complex
        generic :: getH_eig => getH_eigr, getH_eigz

        ! Procedures
        procedure :: construct_realKE, construct_complexKE
        procedure :: expOp_real, expOp_real_mut, expOp_complex, expOp_complex_mut
        procedure :: get_eigr, get_eigz
        procedure :: delete
        procedure :: dimension
        procedure :: printVec_real, printVec_complex
        procedure :: printPoten_real_grid, printPoten_complex_grid
        procedure :: printEigvecs_real, printEigvecs_complex
        procedure :: thetatraj_complex_scale, cap_trajectory

        ! Deferred procedures for inheritance for different basis
        procedure(hamil_interface_real), deferred :: getH_real
        procedure(hamil_interface_complex), deferred :: getH_complex
    end type Hamiltonian

contains

    subroutine printPoten_real_grid(this, filename, vVec)
        class(Hamiltonian), intent(inout) :: this
        character(*), intent(in) :: filename
        real(dp), intent(in) :: vVec(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, size(vVec)
                write(11,*) this%grid%rmin + (ir-1) * this%grid%dr, vVec(ir)
            enddo
        close(11)
    end subroutine printPoten_real_grid

    subroutine printPoten_complex_grid(this, filename, vVec)
        class(Hamiltonian), intent(inout) :: this
        character(*), intent(in) :: filename
        complex(dp), intent(in) :: vVec(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, size(vVec)
                write(11,*) this%grid%rmin + (ir-1) * this%grid%dr, real(vVec(ir), dp), aimag(vVec(ir))
            enddo
        close(11)
    end subroutine printPoten_complex_grid

    subroutine printVec_real(this, filename, eigvals)
        class(Hamiltonian), intent(inout) :: this
        character(*), intent(in) :: filename
        real(dp), intent(in) :: eigvals(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, this%grid%nr
                write(11,*) eigvals(ir)
            enddo
        close(11)
    end subroutine printVec_real

    subroutine printVec_complex(this, filename, eigvals)
        class(Hamiltonian), intent(inout) :: this
        character(*), intent(in) :: filename
        complex(dp), intent(in) :: eigvals(:)
        integer :: ir

        open(11, file=filename, status='replace')
        do ir = 1, this%grid%nr
            write(11,*) real(eigvals(ir), dp), aimag(eigvals(ir))
        enddo
        close(11)
    end subroutine printVec_complex

    subroutine printEigvecs_real(this, filename, eigvecs, evals, nstate)
        class(Hamiltonian), intent(inout) :: this
        character(*), intent(in) :: filename
        integer, intent(in) :: nstate
        real(dp), intent(in) :: eigvecs(:,:), evals(:)
        integer :: ir, istate

        open(11, file=filename, status='replace')
            do ir = 1, this%grid%nr

                write(11,*) this%grid%rmin + (ir - 1) * this%grid%dr, &
                    (eigvecs(ir, istate) + evals(istate), istate = 1, nstate)
            enddo
        close(11)
    end subroutine printEigvecs_real

    subroutine printEigvecs_complex(this, filename, eigvecs, evals, nstate)
        class(Hamiltonian), intent(inout) :: this
        character(*), intent(in) :: filename
        integer, intent(in) :: nstate
        complex(dp), intent(in) :: eigvecs(:,:)
        real(dp), intent(in) :: evals(:)
        integer :: ir, istate

        open(11, file=filename, status='replace')
        do ir = 1, this%grid%nr
            write(11,*) this%grid%rmin + (ir - 1) * this%grid%dr, &
                    (real(eigvecs(ir, istate), dp) + evals(istate), &
            aimag(eigvecs(ir, istate)) + evals(istate), istate = 1, nstate)
        enddo
        close(11)
    end subroutine printEigvecs_complex

    subroutine getH_real(this, potential)
        class(Hamiltonian), intent(inout) :: this
        class(potential_base), intent(in) :: potential
        real(dp), allocatable :: hmat(:, :), vVec(:)
        real(dp) :: eigval(this%grid%nr)
        integer :: ir

        allocate(hmat(this%grid%nr, this%grid%nr), vVec(this%grid%nr))
        call kinetic_energy_matrix(hmat)

        do ir = 1, this%grid%nr
            hmat(ir,ir) = hmat(ir,ir) + vVec(ir)
        enddo

        call solve_Hamil(hmat, eigval)
        call print_vector("gaussEigval.dat", eigval)
        call printEigvecs("gauss_eigvecs.dat", hmat, eigval, 20)
        deallocate(hmat, vVec)
    end subroutine getH_real

    subroutine getH_complex(self, potential, theta)
        class(Hamiltonian), intent(inout) :: this
        class(potential_base), intent(in) :: potential
        real(dp), intent(in) :: theta
        real(dp), allocatable :: tmat(:,:)
        complex(dp), allocatable :: hmat(:, :), vVec(:)
        complex(dp), allocatable :: vecsL(:, :), vecsR(:)
        complex(dp), allocatable :: eigval(:)
        integer :: ir

        allocate(tmat(this%grid%nr, this%grid%nr))
        call kinetic_energy_matrix(tmat)

        allocate(hmat(this%grid%nr, this%grid%nr), vVec(this%grid%nr))

        do ir = 1, this%grid%nr
            hmat(ir,ir) = hmat(ir,ir) + vVec(ir)
        enddo

        allocate(eigval(this%grid%nr))
        call solve_Hamil(hmat, eigval)
        call print_vector("gaussEigval.dat", eigval)
        call printEigvecs("gauss_eigvecs.dat", hmat, eigval, 20)
        deallocate(hmat, vVec)
    end subroutine getH_complex

    subroutine cap_trajectory(lmda_max, nlmda, Vcap)
        class(Hamiltonian), intent(inout) :: this
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
                    aimag(eigenvalues(i)), i = 1, this%grid%nr)
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

    subroutine thetatraj_complex_scale(this, theta_max, ntheta)
        class(Hamiltonian), intent(inout) :: this
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

        allocate(t_matrix(this%grid%nr, this%grid%nr))
        allocate(v_vector(this%grid%nr), h_matrix(this%grid%nr, this%grid%nr))
        allocate(eigenvalues(this%grid%nr))
        allocate(left_vecs(this%grid%nr, this%grid%nr), right_vecs(this%grid%nr, this%grid%nr))
        allocate(prev_right_vecs(this%grid%nr, this%grid%nr))

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
            do i = 1, this%grid%nr
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
                    aimag(eigenvalues(i)), i = 1, this%grid%nr)
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
