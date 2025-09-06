module matrix_algebra
    use grid_module
    implicit none

    public :: sort_eigenvalues_real, sort_eigenvalues_overlap
    public :: solve_Hamil, matmul

    interface matmul
        procedure complex_matmul, real_matmul
    end interface matmul

    interface solve_Hamil
        procedure solve_hamil_complex, solve_hamil_real
    end interface solve_Hamil
    
contains
    
    subroutine sort_eigenvalues_real(eigenvalues, left_vecs, right_vecs)
        complex(dp), intent(inout) :: eigenvalues(:)
        complex(dp), intent(inout) :: left_vecs(:,:), right_vecs(:,:)
        complex(dp), allocatable :: sorted_eigenvals(:)
        complex(dp), allocatable :: sorted_left(:,:), sorted_right(:,:)
        
        real(dp), allocatable :: real_parts(:)
        integer, allocatable :: indices(:)
        integer :: i, ndim
        
        ndim = size(eigenvalues)

        allocate(real_parts(ndim), indices(ndim))
        allocate(sorted_eigenvals(ndim), sorted_left(ndim, ndim))
        allocate(sorted_right(ndim,ndim))
        
        real_parts = real(eigenvalues, dp)

        do i = 1, ndim
            indices(i) = minloc(real_parts, 1)
            sorted_eigenvals(i) = eigenvalues(indices(i))
            real_parts(indices(i)) = 1.0D+7
        end do
        
        sorted_left = left_vecs(:, indices)
        sorted_right = right_vecs(:, indices)
        
        eigenvalues = sorted_eigenvals
        left_vecs = sorted_left
        right_vecs = sorted_right

        deallocate(real_parts, indices, sorted_eigenvals, sorted_left, sorted_right)
    end subroutine sort_eigenvalues_real
    
    subroutine sort_eigenvalues_overlap(eigenvalues, left_vecs, right_vecs, &
                                      prev_right_vecs)
        complex(dp), intent(inout) :: eigenvalues(:)
        complex(dp), intent(inout) :: left_vecs(:,:), right_vecs(:,:)
        complex(dp), intent(inout) :: prev_right_vecs(:,:)
        complex(dp), allocatable :: sorted_eigenvals(:)
        complex(dp), allocatable :: sorted_left(:,:), sorted_right(:,:)
        
        complex(dp), allocatable :: overlap_matrix(:,:)
        integer, allocatable :: indices(:)
        integer :: ndim
        
        ndim = size(eigenvalues)
        allocate(overlap_matrix(ndim, ndim), indices(ndim))
        allocate(sorted_eigenvals(ndim), sorted_left(ndim,ndim), sorted_right(ndim,ndim))
        
        call complex_matmul(left_vecs, 'c', prev_right_vecs, 'N', overlap_matrix)
        indices = maxloc(cdabs(overlap_matrix), 1)
        
        sorted_eigenvals = eigenvalues(indices)
        sorted_left = left_vecs(:, indices)
        sorted_right = right_vecs(:, indices)

        eigenvalues = sorted_eigenvals
        left_vecs =  sorted_left
        right_vecs = sorted_right
        
        deallocate(overlap_matrix, indices, sorted_eigenvals)
        deallocate(sorted_left, sorted_right)
    end subroutine sort_eigenvalues_overlap

    subroutine complex_matmul(a, chara, b, charb, c)
        character(len=1) :: chara, charb
        complex(dp), intent(in) :: a(:,:), b(:,:)
        complex(dp), intent(out) :: c(:,:)
        integer :: ndim
        real(dp) :: aval = cone, bval = czero
        external :: zgemm
        ndim = size(a, 1)
        call zgemm (chara, charb, ndim, ndim, ndim, AVAL, a, ndim, &
                      b, ndim, BVAL, c, ndim)
    end subroutine complex_matmul

    subroutine real_matmul(a, chara, b, charb, c)
        character(len=1) :: chara, charb
        real(dp), intent(in) :: a(:,:), b(:,:)
        real(dp), intent(out) :: c(:,:)
        integer :: ndim
        real(dp) :: aval = cone, bval = czero
        external :: zgemm
        ndim = size(a, 1)
        call dgemm (chara, charb, ndim, ndim, ndim, AVAL, a, ndim, &
                b, ndim, BVAL, c, ndim)
    end subroutine real_matmul

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

end module matrix_algebra

