module eigenvalue_sorting
    use physical_constants
    use grid_parameters
    implicit none
    
    private :: complex_matmul
    public :: sort_eigenvalues_real, sort_eigenvalues_overlap
    
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
    
end module eigenvalue_sorting

