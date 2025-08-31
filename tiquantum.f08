module quantum_solver
    use physical_constants
    use grid_parameters
    implicit none
    private

    public :: solve_hamil_real, solve_hamil_complex
    public :: kinetic_energy_matrix, potential_energy_vector
    public :: superGauss_, Gauss_real

    interface print_vector
        procedure printVec_real, printVec_complex
    end interface print_vector

    interface printPoten_grid
        procedure printPoten_real_grid, printPoten_complex_grid
    end interface printPoten_grid

    interface printEigvecs
        procedure printEigvecs_real, printEigvecs_complex
    end interface printEigvecs

    interface Gauss
        procedure gaussian_real, gaussian_complex
    end interface

    interface model_jol
        procedure potential_real, potential_complex
    end interface model_jol

    interface superGauss
        procedure superGauss_real, superGauss_complex
    end interface superGauss

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

    pure elemental real(dp) function gaussian_real(x, sigma, coef, displace, nbarrier) result(res)
        real(dp), intent(in) :: x, sigma, coef, displace
        integer, intent(in) :: nbarrier
        integer :: igauss
        real(dp) :: val, two_sigma2

        two_sigma2 = 2 * sigma * sigma
        res = zero

        if ( mod(nbarrier, 2) == 0 ) then
            do igauss = 1, nbarrier / 2
                val = 0.5 * igauss * displace
                res = res + dexp(-(x - val)**2/two_sigma2) &
                                    + dexp(-(x + val)**2/two_sigma2)
            enddo
        else
            res = dexp(-(x*x)/(2*sigma*sigma))
            do igauss = 1, (nbarrier - 1)/2
                val = igauss * displace
                res = res + dexp(-(x - val)**2/two_sigma2) &
                                    + dexp(-(x + val)**2/two_sigma2)
            enddo
        endif

        res = res * coef
        return
    end function gaussian_real

    pure elemental complex(dp) function gaussian_complex(x, sigma, coef, displace, nbarrier)
        complex(dp), intent(in) :: x
        real(dp), intent(in) :: sigma, coef, displace
        integer, intent(in) :: nbarrier
        integer :: igauss
        real(dp) :: val, two_sigma2 

        two_sigma2 = 2 * sigma * sigma
        gaussian_complex = zero

        if ( mod(nbarrier, 2) == 0 ) then
            do igauss = 1, nbarrier / 2
                val = 0.5 * igauss * displace
                gaussian_complex = gaussian_complex + cdexp(-(x - val)**2/two_sigma2) &
                                    + cdexp(-(x + val)**2/two_sigma2)
            enddo
        else
            gaussian_complex = cdexp(-(x*x)/(2*sigma*sigma))
            do igauss = 1, (nbarrier - 1)/2
                val = igauss * displace
                gaussian_complex = gaussian_complex + cdexp(-(x - val)**2/two_sigma2) &
                                    + cdexp(-(x + val)**2/two_sigma2)
            enddo
        endif

        gaussian_complex = gaussian_complex * coef
        return
    end function gaussian_complex

    pure elemental real(dp) function superGauss_real(x, nbarrier, order, sigma, coef, displace)
        real(dp), intent(in) :: x, sigma, coef, displace
        integer, intent(in) :: order, nbarrier
        integer :: igauss
        real(dp) :: val, two_sigma2 

        two_sigma2 = 2 * sigma * sigma
        superGauss_real = zero

        if ( mod(nbarrier, 2) == 0 ) then
            do igauss = 1, nbarrier / 2
                val = 0.5 * igauss * displace
                superGauss_real = superGauss_real + dexp(-(x - val)**order/two_sigma2) &
                        + dexp(-(x + val)**order/two_sigma2)
            enddo
        else
            superGauss_real = dexp(-(x**order)/(2*sigma*sigma))
            do igauss = 1, (nbarrier - 1) / 2
                val = igauss * displace
                superGauss_real = superGauss_real + dexp(-(x - val)**order/two_sigma2) &
                        + dexp(-(x + val)**order/two_sigma2)
            enddo
        endif

            superGauss_real = superGauss_real * coef
        return
    end function superGauss_real

    pure elemental complex(dp) function superGauss_complex(x, nbarrier, order, sigma, coef, displace)
        complex(dp), intent(in) :: x
        real(dp), intent(in) :: sigma, coef, displace
        integer, intent(in) :: order, nbarrier
        integer :: igauss
        real(dp) :: val, two_sigma2

        two_sigma2 = 2 * sigma * sigma
        superGauss_complex = zero

        if ( mod(nbarrier, 2) == 0 ) then
            do igauss = 1, nbarrier / 2
                val = 0.5 * igauss * displace
                superGauss_complex = superGauss_complex + cdexp(-(x - val)**order/two_sigma2) &
                        + cdexp(-(x + val)**order/two_sigma2)
            enddo
        else
            superGauss_complex = cdexp(-(x**order)/(2*sigma*sigma))
            do igauss = 1, (nbarrier - 1) / 2
                val = igauss * displace
                superGauss_complex = superGauss_complex + cdexp(-(x - val)**order/two_sigma2) &
                        + cdexp(-(x + val)**order/two_sigma2)
            enddo
        endif

        superGauss_complex = superGauss_complex * coef
        return
    end function superGauss_complex
    
    pure elemental real(dp) function potential_real(x)
        real(dp), intent(in) :: x
        potential_real = dexp(-0.1_dp * x * x) * (x * x * 0.5_dp - 0.8_dp)
        return 
    end function potential_real

    pure elemental complex(dp) function potential_complex(cx)
        complex(dp), intent(in) :: cx
        potential_complex = cdexp(-0.1_dp * cx * cx) * (cx * cx * 0.5_dp - 0.8_dp)
       return
    end function potential_complex

    subroutine potential_energy_real(sigma, coef, displace, nbarrier, v_vector)
        integer, intent(in) :: nbarrier
        real(dp), intent(in) :: sigma, coef, displace
        real(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: x
        
        do i = 1, nr
            x = r_min + real(i-1, dp) * dr
!            v_vector(i) = model_jol(x)
            v_vector(i) = gauss(x, sigma, coef, displace, nbarrier)
        end do

    end subroutine potential_energy_real

    subroutine potenSuperGauss_real(sigma, coef, order, displace, nbarrier, v_vector)
        integer, intent(in) :: nbarrier, order
        real(dp), intent(in) :: sigma, coef, displace
        real(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: x

        do i = 1, nr
            x = r_min + real(i-1, dp) * dr
            v_vector(i) = superGauss(x, nbarrier, order, sigma, coef, displace)
        end do
    end subroutine potenSuperGauss_real
    
    subroutine potential_energy_vector(theta, v_vector)
        real(dp), intent(in) :: theta
        complex(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: x
        complex(dp) :: cx
        
        do i = 1, nr
            x = r_min + real(i-1, dp) * dr
            cx = x * cdexp(i_unit * theta)
            v_vector(i) = model_jol(cx)
        end do
    end subroutine potential_energy_vector

    subroutine printPoten_real_grid(filename, vVec)
        character, intent(in) :: filename
        real(dp), intent(in) :: vVec(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, size(vVec)
                write(11,*) r_min + (ir-1) * dr, vVec(ir)
            enddo
        close(11)
    end subroutine printPoten_real_grid

    subroutine printPoten_complex_grid(filename, vVec)
        character, intent(in) :: filename
        complex(dp), intent(in) :: vVec(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, size(vVec)
                write(11,*) r_min + (ir-1) * dr, real(vVec(ir), dp), aimag(vVec(ir))
            enddo
        close(11)
    end subroutine printPoten_complex_grid

    subroutine printVec_real(filename, eigvals)
        character, intent(in) :: filename
        real(dp), intent(in) :: eigvals(:)
        integer :: ir

        open(11, file=filename, status='replace')
            do ir = 1, nr
                write(11,*) eigvals(ir)
            enddo
        close(11)
    end subroutine printVec_real

    subroutine printVec_complex(filename, eigvals)
        character, intent(in) :: filename
        complex(dp), intent(in) :: eigvals(:)
        integer :: ir

        open(11, file=filename, status='replace')
        do ir = 1, nr
            write(11,*) real(eigvals(ir), dp), aimag(eigvals(ir))
        enddo
        close(11)
    end subroutine printVec_complex

    subroutine printEigvecs_real(filename, eigvecs, evals, nstate)
        character, intent(in) :: filename
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
        character, intent(in) :: filename
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

    subroutine Gauss_real(sigma, coef, displace, nbarrier)
        integer, intent(in) :: nbarrier
        real(dp), intent(in) :: sigma, coef, displace
        real(dp), allocatable :: hmat(:, :), vVec(:)
        real(dp) :: eigval(nr)
        integer :: ir

        allocate(hmat(nr, nr), vVec(nr))
        call kinetic_energy_matrix(hmat)
        call potential_energy_real(sigma, coef, displace, nbarrier, vVec)
        call printPoten_grid("gaussPoten.dat", vVec)

        do ir = 1, nr
            hmat(ir,ir) = hmat(ir,ir) + vVec(ir)
        enddo

        call solve_hamil_real(hmat, eigval)
        call print_vector("gaussEigval.dat", eigval)
        call printEigvecs("gauss_eigvecs.dat", hmat, eigval, 20)
        deallocate(hmat, vVec)
    end subroutine Gauss_real

    subroutine superGauss_(sigma, coef, order, displace, nbarrier)
        integer, intent(in) :: nbarrier, order
        real(dp), intent(in) :: sigma, coef, displace
        real(dp), allocatable :: hmat(:, :), vVec(:)
        real(dp) :: eigval(nr)
        integer :: ir

        allocate(hmat(nr, nr), vVec(nr))
        call kinetic_energy_matrix(hmat)
        call potenSuperGauss_real(sigma, coef, order, displace, nbarrier, vVec)
        call printPoten_grid("SuperGauss_Poten.dat", vVec)

        do ir = 1, nr
            hmat(ir,ir) = hmat(ir,ir) + vVec(ir)
        enddo

        call solve_hamil_real(hmat, eigval)
        call print_vector("SuperGauss_Eigval.dat", eigval)
        call printEigvecs("SuperGauss_eigvecs.dat", hmat, eigval, 20)

        deallocate(hmat, vVec)
    end subroutine superGauss_
    
    SUBROUTINE solve_hamil_real(H, eigvals)
      real(dp), intent(inout) :: H(:, :)
      real(dp), intent(out) :: eigvals(:)
      real(dp), allocatable :: WORK(:)
      integer :: ndim, LDA, ldwork, info

      external :: zgeev

        ndim = size(h, 1)
        LDA = ndim
        LDWORK = 3 * ndim - 1
        allocate(WORK(LDWORK))

        CALL DSYEV('V', 'L', ndim, H, LDA, EIGVALS, WORK, LDWORK, info)

        if (info /= 0) then
            write(*,'(a,i0)') 'Error in ZGEEV: info = ', info
            stop
        end if
        
        deallocate(WORK)
    end subroutine solve_hamil_real

    subroutine solve_hamil_complex(h_matrix, eigenvalues, left_vectors, right_vectors)
        complex(dp), intent(inout) :: h_matrix(:,:)
        complex(dp), intent(out) :: eigenvalues(:)
        complex(dp), intent(out) :: left_vectors(:,:)
        complex(dp), intent(out) :: right_vectors(:,:)
        
        external :: zgeev
        
        integer :: n, lda, ldvl, ldvr, lwork, info
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)
        
        n = size(h_matrix, 1)
        lda = n
        ldvl = n
        ldvr = n
        lwork = 2*n
        
        allocate(work(lwork), rwork(2*n))
        
        call zgeev('V', 'V', n, h_matrix, lda, eigenvalues, &
                   left_vectors, ldvl, right_vectors, ldvr, &
                   work, lwork, rwork, info)
        
        if (info /= 0) then
            write(*,'(a,i0)') 'Error in ZGEEV: info = ', info
            stop
        end if
        
        deallocate(work, rwork)
    end subroutine solve_hamil_complex
    
end module quantum_solver