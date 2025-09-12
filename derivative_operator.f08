! Created by prash on 9/10/2025.
module derivative_operators
    use grid_module
    use matrix_algebra
    implicit none

    type, abstract :: kinetic_base
        class(rgrid), allocatable :: grid
        real(dp), allocatable :: RKEmat(:,:), RKEvec(:)
        complex(dp), allocatable :: ZKEmat(:,:), ZKEvec(:)
        integer(sp) :: size = 0, error_unit = 201
        integer(sp) :: iComplex = 0
        real(dp) :: mass = 1.
        logical :: is_constructed = .false.
        logical :: complex = .false.
    contains
        ! Generic variable name for subroutine call
        generic :: construct_ke => construct_realKE, construct_complexKE
        generic :: expOp => expOp_real, expOp_real_mut, expOp_complex, expOp_complex_mut
        generic :: get_ke => getKE_real, getKE_complex
        generic :: get_eig => get_eigr, get_eigz

        ! Procedures
        procedure :: construct_realKE, construct_complexKE
        procedure :: expOp_real, expOp_real_mut, expOp_complex, expOp_complex_mut
        procedure :: get_eigr, get_eigz
        procedure :: delete
        procedure :: dimension

        ! Deferred procedures for inheritance for different basis
        procedure(ke_interface_real), deferred :: getKE_real
        procedure(ke_interface_complex), deferred :: getKE_complex
    end type kinetic_base

    type, extends(kinetic_base) :: sinc_dvr
    contains
        procedure :: getKE_real => sinc_get_real
        procedure :: getKE_complex => sinc_get_complex
        procedure, private :: compute_sinc_matrix
        ! Automatic deallocate the memory
        final :: sinc_destruct
    end type sinc_dvr

    abstract interface
        subroutine ke_interface_real(this)
            import :: kinetic_base, dp
            class(kinetic_base), intent(inout) :: this
        end subroutine ke_interface_real

        subroutine ke_interface_complex(this, theta)
            import :: kinetic_base, dp
            class(kinetic_base), intent(inout) :: this
            real(dp), intent(in) :: theta
        end subroutine ke_interface_complex
    end interface

contains

    subroutine delete(this)
        class(kinetic_base), intent(inout) :: this

        if (allocated(this%RKEvec)) deallocate(this%RKEvec)
        if (allocated(this%ZKEvec)) deallocate(this%ZKEvec)
        if (allocated(this%RKEmat)) deallocate(this%RKEmat)
        if (allocated(this%ZKEmat)) deallocate(this%ZKEmat)
        if (allocated(this%grid)) deallocate(this%grid)

        this%size = 0
        this%is_constructed = .false.
    end subroutine delete

    subroutine sinc_destruct(this)
        type(sinc_dvr), intent(inout) :: this
        call this%delete
    end subroutine sinc_destruct

    function dimension(this) result(dim)
        class(kinetic_base), intent(in) :: this
        integer :: dim
        dim = this%size
    end function dimension

    subroutine construct_realKE(this, igrid, imass)
        class(kinetic_base), intent(inout) :: this
        class(rgrid), intent(in) :: igrid
        real(dp), intent(in) :: imass
        integer :: alloc_stat

        call this%delete()

        ! Allocate grid (more efficient than pointer)
        allocate(this%grid, source=igrid, stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(this%error_unit, *) 'Error: Failed to allocate grid'
            return
        end if

        this%mass = imass
        this%complex = .false.
        this%size = this%grid%nr

        allocate(this%RKEmat(this%size, this%size), stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(this%error_unit, *) 'Error: Failed to allocate kinetic energy matrix'
            return
        end if

        this%is_constructed = .true.

        ! Calculate the kinetic energy matrix
        call this%getKE_real
    end subroutine construct_realKE

    subroutine construct_complexKE(this, igrid, imass, iComplex)
        class(kinetic_base), intent(inout) :: this
        class(rgrid), intent(in) :: igrid
        real(dp), intent(in) :: imass
        integer(sp), intent(in) :: iComplex
        integer :: alloc_stat

        if (iComplex == 0) then
            write(this%error_unit, *) 'Error: 0 is for setting up real Ke-Matrix'
            return
        end if

        call this%delete

        ! Allocate grid (more efficient than pointer)
        allocate(this%grid, source=igrid, stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(this%error_unit, *) 'Error: Failed to allocate grid'
            return
        end if

        this%mass = imass
        this%complex = .true.
        this%size = this%grid%nr

        allocate(this%ZKEmat(this%size, this%size), stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(this%error_unit, *) 'Error: Failed to allocate kinetic energy matrix'
            return
        end if

        this%is_constructed = .true.
    end subroutine construct_complexKE

    subroutine get_eigr(this, evals, eigvecs)
        class(kinetic_base), intent(in) :: this
        real(dp), intent(out) :: evals(:), eigvecs(:,:)

        if (.not. this%is_constructed) then
            write(this%error_unit, *) 'Kinetic energy is not constructed'
            return
        end if

        eigvecs = this%RKEmat
        call solve_Hamil(eigvecs, evals)
    end subroutine get_eigr

    subroutine get_eigz(this, evals, eigvecsL, eigvecsR)
        class(kinetic_base), intent(in) :: this
        complex(dp), intent(out) :: evals(:)
        complex(dp), intent(out) :: eigvecsL(:,:), eigvecsR(:,:)
        complex(dp), allocatable :: kemat(:,:)

        if (.not. this%is_constructed) then
            write(this%error_unit, *) 'Kinetic energy is not constructed'
            return
        end if

        allocate(kemat(this%size, this%size))
        call solve_Hamil(kemat, evals, eigvecsL, eigvecsR)
        deallocate(kemat)
    end subroutine get_eigz

    ! Concrete implementations for sinc_dvr type
    subroutine sinc_get_real(this)
        class(sinc_dvr), intent(inout) :: this

        if (.not. this%is_constructed) then
            write(this%error_unit, *) 'Error: Kinetic energy object not constructed'
            return
        end if

        call this%compute_sinc_matrix()
        this%RKEmat = this%RKEmat
    end subroutine sinc_get_real

    subroutine sinc_get_complex(this, theta)
        class(sinc_dvr), intent(inout) :: this
        real(dp), intent(in) :: theta

        if (.not. this%is_constructed) then
            write(this%error_unit, *) 'Error: Kinetic energy object not constructed'
            return
        end if

        call this%compute_sinc_matrix()
        this%ZKEmat = this%RKEmat * cdexp(-2*i_unit * theta)
    end subroutine sinc_get_complex

    subroutine compute_sinc_matrix(this)
        class(sinc_dvr), intent(inout) :: this
        integer :: i, j

        if (.not. allocated(this%grid)) then
            write(this%error_unit, *) 'Error: Grid not allocated'
            return
        end if

        do i = 1, this%size
            this%RKEmat(i, i) = pi_squared_over_six / (this%mass * this%grid%dr ** 2)
            do j = 1, i - 1
                this%RKEmat(i, j) = ((-one)**(i-j)) / (this%mass * ((i-j) * this%grid%dr)**2)
                this%RKEmat(j, i) = this%RKEmat(i, j)
            end do
        end do
    end subroutine compute_sinc_matrix

    ! Placeholder implementations for expKE procedures
    subroutine expOp_real(this, dt, psi_in, psi_out)
        class(kinetic_base), intent(in) :: this
        real(dp), intent(in) :: dt
        complex(dp), intent(in) :: psi_in(:)
        complex(dp), intent(out) :: psi_out(:)
        psi_out = psi_in  ! Placeholder
    end subroutine expOp_real

    subroutine expOp_real_mut(this, dt, psi_inout)
        class(kinetic_base), intent(in) :: this
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: psi_inout(:)
        psi_inout = 0.  ! Placeholder
    end subroutine expOp_real_mut

    subroutine expOp_complex(this, dt, psi_in, psi_out, dummy)
        class(kinetic_base), intent(in) :: this
        complex(dp), intent(in) :: dt
        complex(dp), intent(in) :: psi_in(:)
        complex(dp), intent(out) :: psi_out(:)
        integer, intent(in) :: dummy
        psi_out = psi_in  ! Placeholder
    end subroutine expOp_complex

    subroutine expOp_complex_mut(this, dt, psi_inout, dummy)
        class(kinetic_base), intent(in) :: this
        real(dp), intent(in) :: dt
        complex(dp), intent(inout) :: psi_inout(:)
        integer, intent(in) :: dummy
        psi_inout = 0.  ! Placeholder
    end subroutine expOp_complex_mut

end module derivative_operators