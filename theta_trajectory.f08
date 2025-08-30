subroutine thetatraj_complex_scale
    use physical_constants
    use grid_parameters
    use quantum_solver
    use eigenvalue_sorting
    implicit none

    real(dp), allocatable :: t_matrix(:,:)
    complex(dp), allocatable :: v_vector(:), h_matrix(:,:)
    complex(dp), allocatable :: eigenvalues(:)
    complex(dp), allocatable :: left_vecs(:,:), right_vecs(:,:)
    complex(dp), allocatable :: prev_right_vecs(:,:)

    real(dp) :: theta
    complex(dp) :: phase_factor
    integer :: i, j, ith

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

    do ith = 1, nth
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

        call solve_hamil_complex(h_matrix, eigenvalues, left_vecs, right_vecs)

        if (ith == 1) then

            call sort_eigenvalues_real(eigenvalues, left_vecs, right_vecs)

            open(11, file='initial_eigenvalues.dat', status='replace')
            do i = 1, nr
                write(11, '(2es16.8)') real(eigenvalues(i), dp), aimag(eigenvalues(i))
            end do
            close(11)
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