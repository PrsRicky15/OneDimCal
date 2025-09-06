module grid_module
    implicit none

    integer, parameter  :: dp = kind(1.0d0)
    integer, parameter  :: sp = kind(1)
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: pi_squared_over_six = pi*pi/6.0_dp
    complex(dp), parameter :: i_unit = (0.0_dp, 1.0_dp)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)

    ! Abstract interfaces for potential evaluation (to be implemented by potential module)
    abstract interface
        function evaluate_real_interface(pot, x) result(v)
            import :: dp
            class(*), intent(in) :: pot
            real(dp), intent(in) :: x
            real(dp) :: v
        end function evaluate_real_interface

        function evaluate_complex_interface(pot, x) result(v)
            import :: dp
            class(*), intent(in) :: pot
            complex(dp), intent(in) :: x
            complex(dp) :: v
        end function evaluate_complex_interface
    end interface

    ! RGrid type definition
    type :: rgrid
        real(dp) :: rmin, rmax, dr, len, dk, kmin, kmax, cutoff_e
        integer(sp) :: nr
        real(dp), allocatable :: r_grid(:), k_grid(:)
        logical :: grid_generated = .false.
        logical :: k_grid_generated = .false.
    contains
        procedure :: init_rgrid
        procedure :: potential_onGrid_real
        procedure :: potential_onGrid_complex
        procedure :: printPotToFile_real
        procedure :: printPotToFile_complex
        generic :: genPoten => potential_onGrid_real, potential_onGrid_complex
        generic :: printPotToFile => printPotToFile_real, printPotToFile_complex
        procedure, nopass :: from_length => rgrid_from_length
        procedure, nopass :: from_center_width => rgrid_from_center_width
        procedure :: generate_grid => rgrid_generate_grid
        procedure :: get_rgrid => rgrid_get_rgrid
        procedure :: get_kgrid => rgrid_get_kgrid
        procedure :: get_shifted_momentum_grid => rgrid_get_shifted_momentum_grid
        procedure :: display_grid => rgrid_display_grid
        procedure :: display_momentum_info => rgrid_display_momentum_info
        procedure :: save_grid_to_file => rgrid_save_grid_to_file
        procedure :: cleanup => rgrid_cleanup
    end type rgrid

    ! TGrid type definition
    type :: tgrid
        real(dp) :: tmin, tmax, dt, duration, sample_rate, nyquist_freq, freq_resolution
        integer(sp) :: nt
        real(dp), allocatable :: grid(:), omega_grid(:)
        logical :: grid_generated = .false.
        logical :: omega_grid_generated = .false.
    contains
        procedure :: init_tgrid
        procedure :: potential_onGrid_real_t
        procedure :: potential_onGrid_complex_t
        procedure :: printPotToFile_real_t
        procedure :: printPotToFile_complex_t
        generic :: genPoten => potential_onGrid_real_t, potential_onGrid_complex_t
        generic :: printPotToFile => printPotToFile_real_t, printPotToFile_complex_t
        procedure, nopass :: from_duration => tgrid_from_duration
        procedure, nopass :: from_sample_rate => tgrid_from_sample_rate
        procedure :: generate_grid => tgrid_generate_grid
        procedure :: get_grid => tgrid_get_grid
        procedure :: get_frequency_grid => tgrid_get_frequency_grid
        procedure :: get_shifted_frequency_grid => tgrid_get_shifted_frequency_grid
        procedure :: display_grid => tgrid_display_grid
        procedure :: display_frequency_info => tgrid_display_frequency_info
        procedure :: save_grid_to_file => tgrid_save_grid_to_file
        procedure :: cleanup => tgrid_cleanup
    end type tgrid

contains

    !====================================================================
    ! RGrid procedures
    !====================================================================

    subroutine init_rgrid(this, rmin, rmax, ngrid)
        class(rgrid), intent(inout) :: this
        real(dp), intent(in) :: rmin, rmax
        integer(sp), intent(in) :: ngrid

        this%rmin = rmin
        this%rmax = rmax
        this%nr = ngrid
        this%dr = (rmax - rmin) / real(ngrid, dp)
        this%len = rmax - rmin
        this%dk = 2.0_dp * pi / this%len
        this%kmin = -pi / this%dr
        this%kmax = pi / this%dr
        this%cutoff_e = this%kmax**2 / 2.0_dp
        this%grid_generated = .false.
        this%k_grid_generated = .false.

        if (allocated(this%r_grid)) deallocate(this%r_grid)
        if (allocated(this%k_grid)) deallocate(this%k_grid)
        allocate(this%r_grid(ngrid))
        allocate(this%k_grid(ngrid))
    end subroutine init_rgrid

    function rgrid_from_length(length, ngrid) result(grid)
        real(dp), intent(in) :: length
        integer(sp), intent(in) :: ngrid
        type(rgrid) :: grid

        real(dp) :: grd_min, grd_max

        grd_min = -length / 2.0_dp
        grd_max = length / 2.0_dp
        call grid%init_rgrid(grd_min, grd_max, ngrid)
    end function rgrid_from_length

    function rgrid_from_center_width(center, width, ngrid) result(grid)
        real(dp), intent(in) :: center, width
        integer(sp), intent(in) :: ngrid
        type(rgrid) :: grid

        real(dp) :: grd_min, grd_max

        grd_min = center - width / 2.0_dp
        grd_max = center + width / 2.0_dp
        call grid%init_rgrid(grd_min, grd_max, ngrid)
    end function rgrid_from_center_width

    subroutine rgrid_generate_grid(this)
        class(rgrid), intent(inout) :: this
        integer(sp) :: i

        do i = 1, this%nr
            this%r_grid(i) = this%rmin + real(i-1, dp) * this%dr
        end do
        this%grid_generated = .true.
    end subroutine rgrid_generate_grid

    function rgrid_get_rgrid(this) result(grid)
        class(rgrid), intent(inout) :: this
        real(dp) :: grid(this%nr)

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if
        grid = this%r_grid
    end function rgrid_get_rgrid

    subroutine rgrid_get_kgrid(this, k_grid)
        class(rgrid), intent(inout) :: this
        real(dp), intent(out) :: k_grid(:)
        integer(sp) :: i, n

        if (.not. this%k_grid_generated) then
            n = this%nr
            ! Generate FFT frequencies
            do i = 1, n/2
                this%k_grid(i) = 2.0_dp * pi * real(i-1, dp) / (this%dr * real(n, dp))
            end do
            do i = n/2 + 1, n
                this%k_grid(i) = 2.0_dp * pi * real(i-1-n, dp) / (this%dr * real(n, dp))
            end do
            this%k_grid_generated = .true.
        end if
        k_grid = this%k_grid
    end subroutine rgrid_get_kgrid

    subroutine rgrid_get_shifted_momentum_grid(this, k_shifted)
        class(rgrid), intent(inout) :: this
        real(dp), intent(out) :: k_shifted(:)
        real(dp) :: k_temp(this%nr)

        call this%get_kgrid(k_temp)
        call fftshift_real(k_temp, k_shifted)
    end subroutine rgrid_get_shifted_momentum_grid

    subroutine rgrid_display_grid(this)
        class(rgrid), intent(inout) :: this

        if (.not. this%grid_generated) then
            print *, "Grid not generated yet."
        else
            print '(A,I0,A,F0.6,A,F0.6)', "Spatial Grid: ", this%nr, &
                    " points from ", this%rmin, " to ", this%rmax
            print '(A,F0.6)', "Grid spacing: ", this%dr
            print '(A,5F10.6)', "First few points: ", this%r_grid(1:min(5,this%nr))
            print '(A,5F10.6)', "Last few points: ", &
                    this%r_grid(max(1,this%nr-4):this%nr)
        end if
    end subroutine rgrid_display_grid

    subroutine rgrid_display_momentum_info(this)
        class(rgrid), intent(inout) :: this
        real(dp) :: k_temp(this%nr)

        call this%get_kgrid(k_temp)
        print '(A,I0,A)', "Momentum Grid: ", size(k_temp), " points"
        print '(A,F0.6)', "dk: ", this%dk
        print '(A,F0.6)', "kMax: ", this%kmax
        print '(A,F0.6,A,F0.6,A)', "k range: [", minval(k_temp), ", ", maxval(k_temp), "]"
    end subroutine rgrid_display_momentum_info

    ! Evaluate potential on rgrid (real) - using procedure pointer
    subroutine potential_onGrid_real(this, inPot, eval_proc, v_vector)
        class(rgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_real_interface) :: eval_proc
        real(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: x

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        do i = 1, this%nr
            x = this%rmin + real(i-1, dp) * this%dr
            v_vector(i) = eval_proc(inPot, x)
        end do
    end subroutine potential_onGrid_real

    ! Evaluate potential on rgrid (complex) - using procedure pointer
    subroutine potential_onGrid_complex(this, inPot, eval_proc, theta, v_vector)
        class(rgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_complex_interface) :: eval_proc
        real(dp), intent(in) :: theta
        complex(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: x

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        do i = 1, this%nr
            x = this%rmin + real(i-1, dp) * this%dr
            v_vector(i) = eval_proc(inPot, x * exp(cmplx(0.0_dp, theta, dp)))
        end do
    end subroutine potential_onGrid_complex

    ! Print potential to file (real) for rgrid
    subroutine printPotToFile_real(this, inPot, eval_proc, filename)
        class(rgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_real_interface) :: eval_proc
        character(len=*), intent(in) :: filename
        integer :: i, unit
        real(dp) :: x

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        open(newunit=unit, file=filename, status='replace')
        do i = 1, this%nr
            x = this%rmin + real(i-1, dp) * this%dr
            write(unit,'(2ES15.6)') x, eval_proc(inPot, x)
        end do
        close(unit)
    end subroutine printPotToFile_real

    ! Print potential to file (complex) for rgrid
    subroutine printPotToFile_complex(this, inPot, eval_proc, theta, filename)
        class(rgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_complex_interface) :: eval_proc
        real(dp), intent(in) :: theta
        character(len=*), intent(in) :: filename
        integer :: i, unit
        real(dp) :: x
        complex(dp) :: poten

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        open(newunit=unit, file=filename, status='replace')
        do i = 1, this%nr
            x = this%rmin + real(i-1, dp) * this%dr
            poten = eval_proc(inPot, x * exp(cmplx(0.0_dp, theta, dp)))
            write(unit,'(3ES15.6)') x, real(poten, dp), aimag(poten)
        end do
        close(unit)
    end subroutine printPotToFile_complex

    subroutine rgrid_save_grid_to_file(this, filename)
        class(rgrid), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: unit
        integer(sp) :: i

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        open(newunit=unit, file=filename, status='replace')
        write(unit, '(A,I0,A,F0.6)') "# Spatial grid: ", this%nr, " points, dx=", this%dr
        do i = 1, this%nr
            write(unit, '(F0.12)') this%r_grid(i)
        end do
        close(unit)
        print *, "Grid saved to ", filename
    end subroutine rgrid_save_grid_to_file

    subroutine rgrid_cleanup(this)
        class(rgrid), intent(inout) :: this
        if (allocated(this%r_grid)) deallocate(this%r_grid)
        if (allocated(this%k_grid)) deallocate(this%k_grid)
    end subroutine rgrid_cleanup

    !====================================================================
    ! TGrid procedures
    !====================================================================

    subroutine init_tgrid(this, tmin, tmax, nt)
        class(tgrid), intent(inout) :: this
        real(dp), intent(in) :: tmin, tmax
        integer(sp), intent(in) :: nt

        this%tmin = tmin
        this%tmax = tmax
        this%nt = nt
        this%dt = (tmax - tmin) / real(nt - 1, dp)
        this%duration = tmax - tmin
        this%sample_rate = 1.0_dp / this%dt
        this%nyquist_freq = pi / this%dt
        this%freq_resolution = 2.0_dp * pi / this%duration
        this%grid_generated = .false.
        this%omega_grid_generated = .false.

        if (allocated(this%grid)) deallocate(this%grid)
        if (allocated(this%omega_grid)) deallocate(this%omega_grid)
        allocate(this%grid(nt))
        allocate(this%omega_grid(nt))
    end subroutine init_tgrid

    function tgrid_from_duration(duration, nt) result(grid)
        real(dp), intent(in) :: duration
        integer(sp), intent(in) :: nt
        type(tgrid) :: grid

        call grid%init_tgrid(0.0_dp, duration, nt)
    end function tgrid_from_duration

    function tgrid_from_sample_rate(sample_rate, duration) result(grid)
        real(dp), intent(in) :: sample_rate, duration
        type(tgrid) :: grid
        integer(sp) :: nt

        nt = int(sample_rate * duration, sp) + 1
        call grid%init_tgrid(0.0_dp, duration, nt)
    end function tgrid_from_sample_rate

    subroutine tgrid_generate_grid(this)
        class(tgrid), intent(inout) :: this
        integer(sp) :: i

        do i = 1, this%nt
            this%grid(i) = this%tmin + real(i-1, dp) * this%dt
        end do
        this%grid_generated = .true.
    end subroutine tgrid_generate_grid

    function tgrid_get_grid(this) result(grid)
        class(tgrid), intent(inout) :: this
        real(dp) :: grid(this%nt)

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if
        grid = this%grid
    end function tgrid_get_grid

    subroutine tgrid_get_frequency_grid(this, omega_grid)
        class(tgrid), intent(inout) :: this
        real(dp), intent(out) :: omega_grid(:)
        integer(sp) :: i, n

        if (.not. this%omega_grid_generated) then
            n = this%nt
            ! Generate FFT frequencies
            do i = 1, n/2
                this%omega_grid(i) = 2.0_dp * pi * real(i-1, dp) / (this%dt * real(n, dp))
            end do
            do i = n/2 + 1, n
                this%omega_grid(i) = 2.0_dp * pi * real(i-1-n, dp) / (this%dt * real(n, dp))
            end do
            this%omega_grid_generated = .true.
        end if
        omega_grid = this%omega_grid
    end subroutine tgrid_get_frequency_grid

    subroutine tgrid_get_shifted_frequency_grid(this, omega_shifted)
        class(tgrid), intent(inout) :: this
        real(dp), intent(out) :: omega_shifted(:)
        real(dp) :: omega_temp(this%nt)

        call this%get_frequency_grid(omega_temp)
        call fftshift_real(omega_temp, omega_shifted)
    end subroutine tgrid_get_shifted_frequency_grid

    ! Evaluate potential on tgrid (real)
    subroutine potential_onGrid_real_t(this, inPot, eval_proc, v_vector)
        class(tgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_real_interface) :: eval_proc
        real(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: t

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        do i = 1, this%nt
            t = this%tmin + real(i-1, dp) * this%dt
            v_vector(i) = eval_proc(inPot, t)
        end do
    end subroutine potential_onGrid_real_t

    ! Evaluate potential on tgrid (complex)
    subroutine potential_onGrid_complex_t(this, inPot, eval_proc, theta, v_vector)
        class(tgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_complex_interface) :: eval_proc
        real(dp), intent(in) :: theta
        complex(dp), intent(out) :: v_vector(:)
        integer :: i
        real(dp) :: t

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        do i = 1, this%nt
            t = this%tmin + real(i-1, dp) * this%dt
            v_vector(i) = eval_proc(inPot, t * exp(cmplx(0.0_dp, theta, dp)))
        end do
    end subroutine potential_onGrid_complex_t

    ! Print potential to file (real) for tgrid
    subroutine printPotToFile_real_t(this, inPot, eval_proc, filename)
        class(tgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_real_interface) :: eval_proc
        character(len=*), intent(in) :: filename
        integer :: i, unit
        real(dp) :: t

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        open(newunit=unit, file=filename, status='replace')
        do i = 1, this%nt
            t = this%tmin + real(i-1, dp) * this%dt
            write(unit,'(2ES15.6)') t, eval_proc(inPot, t)
        end do
        close(unit)
    end subroutine printPotToFile_real_t

    ! Print potential to file (complex) for tgrid
    subroutine printPotToFile_complex_t(this, inPot, eval_proc, theta, filename)
        class(tgrid), intent(inout) :: this
        class(*), intent(in) :: inPot
        procedure(evaluate_complex_interface) :: eval_proc
        real(dp), intent(in) :: theta
        character(len=*), intent(in) :: filename
        integer :: i, unit
        real(dp) :: t
        complex(dp) :: poten

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        open(newunit=unit, file=filename, status='replace')
        do i = 1, this%nt
            t = this%tmin + real(i-1, dp) * this%dt
            poten = eval_proc(inPot, t * exp(cmplx(0.0_dp, theta, dp)))
            write(unit,'(3ES15.6)') t, real(poten, dp), aimag(poten)
        end do
        close(unit)
    end subroutine printPotToFile_complex_t

    subroutine tgrid_display_grid(this)
        class(tgrid), intent(inout) :: this

        if (.not. this%grid_generated) then
            print *, "Grid not generated yet."
        else
            print '(A,I0,A,F0.6,A,F0.6)', "Time Grid: ", this%nt, &
                    " points from ", this%tmin, " to ", this%tmax
            print '(A,F0.6)', "Time step: ", this%dt
            print '(A,F0.2,A)', "Sample rate: ", this%sample_rate, " Hz"
            print '(A,5F10.6)', "First few points: ", this%grid(1:min(5,this%nt))
            print '(A,5F10.6)', "Last few points: ", &
                    this%grid(max(1,this%nt-4):this%nt)
        end if
    end subroutine tgrid_display_grid

    subroutine tgrid_display_frequency_info(this)
        class(tgrid), intent(inout) :: this
        real(dp) :: omega_temp(this%nt)

        call this%get_frequency_grid(omega_temp)
        print '(A,I0,A)', "Frequency Grid: ", size(omega_temp), " points"
        print '(A,F0.6)', "dω: ", this%freq_resolution
        print '(A,F0.6)', "ω_max (Nyquist): ", this%nyquist_freq
        print '(A,F0.6,A,F0.6,A)', "ω range: [", minval(omega_temp), ", ", maxval(omega_temp), "]"
    end subroutine tgrid_display_frequency_info

    subroutine tgrid_save_grid_to_file(this, filename)
        class(tgrid), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: unit
        integer(sp) :: i

        if (.not. this%grid_generated) then
            call this%generate_grid()
        end if

        open(newunit=unit, file=filename, status='replace')
        write(unit, '(A,I0,A,F0.6)') "# Time grid: ", this%nt, " points, dt=", this%dt
        do i = 1, this%nt
            write(unit, '(F0.12)') this%grid(i)
        end do
        close(unit)
        print *, "Grid saved to ", filename
    end subroutine tgrid_save_grid_to_file

    subroutine tgrid_cleanup(this)
        class(tgrid), intent(inout) :: this
        if (allocated(this%grid)) deallocate(this%grid)
        if (allocated(this%omega_grid)) deallocate(this%omega_grid)
    end subroutine tgrid_cleanup

    ! Helper subroutine for fftshift
    subroutine fftshift_real(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: n, mid

        n = size(input)
        mid = n/2

        if (mod(n, 2) == 0) then
            output(1:mid) = input(mid+1:n)
            output(mid+1:n) = input(1:mid)
        else
            output(1:mid) = input(mid+2:n)
            output(mid+1:n) = input(1:mid+1)
        end if
    end subroutine fftshift_real

end module grid_module