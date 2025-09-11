! Created by prash on 9/10/2025.

module kinetic_energy
    use grid_module

    type, abstract :: kinetic_base
        class(rgrid), allocatable :: grid
        real(dp), allocatable :: RKEmat(:,:), RKEvec(:)
        complex(dp), allocatable :: ZKEmat(:,:), ZKEvec(:)
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

        ! Deferred procedures for inheritance for different basis
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

        if (allocated(this%RKEvec)) deallocate(this%RKEvec)
        if (allocated(this%ZKEvec)) deallocate(this%ZKEvec)
        if (allocated(this%RKEmat)) deallocate(this%RKEmat)
        if (allocated(this%ZKEmat)) deallocate(this%ZKEmat)
        if (allocated(this%grid)) deallocate(this%grid)

        this%size = 0
    end subroutine delete

    subroutine destruct(this)
        type(kinetic_base), intent(inout) :: this
        call this%delete()
    end subroutine destruct

end module kinetic_energy