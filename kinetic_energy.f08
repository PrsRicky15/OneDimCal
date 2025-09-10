! Created by prash on 9/10/2025.

module kinetic_energy
    use grid_module

    type, abstract :: kinetic_base
        class(rgrid), allocatable :: grid
        real(dp), allocatable :: real_mat(:,:), ke_real(:)
        complex(dp), allocatable :: complex_mat(:,:), ke_complex(:)
        real(dp) :: mass = 1.0_dp
        integer(sp) :: size = 0
    contains
        generic :: construct => construct_realKE, construct_complexKE
        generic :: exp => expKE_real, expKE_complex
        generic :: get_ke => getKE_real, getKE_complex
        procedure :: initiate_real_ke, initiate_complex_ke
        procedure :: expKE_real, expKE_complex
        procedure :: delete
        procedure :: dimension

        ! Deferred procedures (must be implemented by child types)
        procedure(ke_interface_real), deferred :: getKe_real
        procedure(ke_interface_complex), deferred :: getKe_theta
        final :: destruct
    end type kinetic_base

    abstract interface
        subroutine ke_interface_real(this, mass)
            import kinetic_base, dp
            class(kinetic_base), intent(in) :: this
            real(dp), intent(in) :: mass
            integer :: i, j
        end subroutine ke_interface_real

        subroutine ke_interface_complex(this, mass, theta)
            import kinetic_base, dp
            class(kinetic_base), intent(in) :: this
            real(dp), intent(in) :: mass, theta
            integer :: i, j
        end subroutine ke_interface_complex
    end interface

contains

    subroutine delete(this)
        class(kinetic_base), intent(inout) :: this

        if (allocated(this%ke_real)) deallocate(this%ke_real)
        if (allocated(this%ke_complex)) deallocate(this%ke_complex)
        if (allocated(this%real_mat)) deallocate(this%real_mat)
        if (allocated(this%complex_mat)) deallocate(this%complex_mat)
        if (allocated(this%grid)) deallocate(this%grid)

        this%size = 0
    end subroutine delete

    subroutine destruct(this)
        type(kinetic_base), intent(inout) :: this
        call this%delete()
    end subroutine destruct

end module kinetic_energy