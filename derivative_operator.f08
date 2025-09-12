! Created by prash on 9/10/2025.

module derivative_operators
    use grid_module
    use matrix_algebra
    implicit none

    type, abstract :: momentum_base
        class(rgrid), allocatable :: grid
        real(dp), allocatable :: rmom_mat(:,:), rpvec(:)
        complex(dp), allocatable :: zmom_mat(:,:), zpvec(:)
        integer(sp) :: size = 0, error_unit = 201
        logical :: is_constructed = .false.
        logical :: complex = .false.
    contains
        ! Generic variable name for subroutine call
        generic :: construct => construct_rKE, construct_zKE
        generic :: exp => exprKE, expzKE
        generic :: get_ke => getrKE, getzKE

        ! Procedures
        procedure :: construct_rmom, construct_zmom
        procedure :: exp_rp, exp_zp
        procedure :: delete, get_eig
        procedure :: dimension

        ! Deferred procedures for inheritance for different basis
        procedure(ke_interface_real), deferred :: getrKE
        procedure(ke_interface_complex), deferred :: getzKE
    end type momentum_base

    type, abstract :: kinetic_base
        class(rgrid), allocatable :: grid
        real(dp), allocatable :: RKEmat(:,:), RKEvec(:)
        complex(dp), allocatable :: ZKEmat(:,:), ZKEvec(:)
        integer(sp) :: size = 0, error_unit = 201
        real(dp) :: mass
        logical :: is_constructed = .false.
        logical :: complex = .false.
    contains
        ! Generic variable name for subroutine call
        generic :: construct => construct_realKE, construct_complexKE
        generic :: exp => expKE_real, expKE_complex
        generic :: get_ke => getKE_real, getKE_complex

        ! Procedures
        procedure :: construct_realKE, construct_complexKE
        procedure :: expKE_real, expKE_complex
        procedure :: delete, get_eig
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
        subroutine ke_interface_real(this, mass)
            import :: kinetic_base, dp
            class(kinetic_base), intent(inout) :: this
            real(dp), intent(in) :: mass
        end subroutine ke_interface_real

        subroutine ke_interface_complex(this, mass, theta)
            import :: kinetic_base, dp
            class(kinetic_base), intent(inout) :: this
            real(dp), intent(in) :: mass, theta
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
        call this%getKE_real(this%mass)
    end subroutine construct_realKE

    subroutine construct_complexKE(this, igrid, theta, imass)
        class(kinetic_base), intent(inout) :: this
        class(rgrid), intent(in) :: igrid
        real(dp), intent(in) :: theta
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
        this%complex = .true.
        this%size = this%grid%nr

        allocate(this%ZKEmat(this%size, this%size), stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(this%error_unit, *) 'Error: Failed to allocate kinetic energy matrix'
            return
        end if

        this%is_constructed = .true.
    end subroutine construct_complexKE

    subroutine get_eig_real(this)
        class(kinetic_base), intent(in) :: this
        call solve_Hamil(this)
    end subroutine get_eig_real

    subroutine get_eig_complex(this)
        class(kinetic_base), intent(in) :: this
        call solve_Hamil(this)
    end subroutine get_eig_complex

    ! Placeholder implementations for expKE procedures
    subroutine expKE_real(this, dt, psi_in, psi_out)
        class(kinetic_base), intent(in) :: this
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: psi_in(:)
        real(dp), intent(out) :: psi_out(:)

        ! Implementation would go here
        ! This is a placeholder - actual implementation depends on your needs
        psi_out = psi_in  ! Placeholder
    end subroutine expKE_real

    subroutine expKE_complex(this, dt, psi_in, psi_out)
        class(kinetic_base), intent(in) :: this
        complex(dp), intent(in) :: dt
        complex(dp), intent(in) :: psi_in(:)
        complex(dp), intent(out) :: psi_out(:)

        ! Implementation would go here
        ! This is a placeholder - actual implementation depends on your needs
        psi_out = psi_in  ! Placeholder
    end subroutine expKE_complex

    ! Concrete implementations for sinc_dvr type
    subroutine sinc_get_real(this, mass)
        class(sinc_dvr), intent(inout) :: this
        real(dp), intent(in) :: mass

        if (.not. this%is_constructed) then
            write(this%error_unit, *) 'Error: Kinetic energy object not constructed'
            return
        end if

        call this%compute_sinc_matrix()

        ! Scale by mass factor
        this%RKEmat = this%RKEmat / (2.0_dp * mass)
    end subroutine sinc_get_real

    subroutine sinc_get_complex(this, mass, theta)
        class(sinc_dvr), intent(inout) :: this
        real(dp), intent(in) :: mass, theta

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
            this%RKEmat(i, i) = pi_squared_over_six / (this%grid%dr ** 2)
            do j = 1, i - 1
                this%RKEmat(i, j) = ((-one)**(i-j)) / ((i-j) * this%grid%dr)**2
                this%RKEmat(j, i) = this%RKEmat(i, j)
            end do
        end do
    end subroutine compute_sinc_matrix

end module derivative_operators