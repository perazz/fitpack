! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_surfaces_c.f90 (module fitpack_surfaces_c)
!> @brief C interface to fitpack_surface
!
!   @author Federico Perini
!   @date   2026-08-27
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************

module fitpack_surfaces_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_surfaces, only: fitpack_surface
    use fitpack_surfaces_c_types, only: fitpack_surface_c, fitpack_surface_c_null, fitpack_surface_c_typename
    use fitpack_fx_status, only: fx_status, fx_status_ok, FX_ERROR_ALLOCATION, &
                        FX_ERROR_DEALLOCATION, handle_error
    use fitpack_curves, only: fitpack_curve
    use fitpack_curves_c, only: f_pointer
    use fitpack_curves_c_types, only: fitpack_curve_c, fitpack_curve_c_null, f_type_name

    use, intrinsic :: iso_c_binding
    implicit none(type, external)
    private

    ! ===========================================================================================
    ! SECTION 2: Public Exports
    ! ===========================================================================================

    ! Core interfaces (ALWAYS export these)
    public :: fitpack_surface_c            ! Type and constructor
    public :: fitpack_surface_c_null       ! Null constant
    public :: f_pointer              ! Get Fortran pointer (non-allocating)
    public :: f_associated           ! Pointer-identity check (non bind(C))
    public :: fitpack_surface_c_is_same   ! Pointer-identity check [bind(C)]
    public :: fitpack_surface_c_allocate   ! Allocate new object [bind(C)]
    public :: fitpack_surface_c_destroy    ! Deallocate object [bind(C)]
    public :: fitpack_surface_c_copy       ! Deep copy [bind(C)]
    public :: fitpack_surface_c_associate  ! Shallow copy/pointer [bind(C)]
    public :: fitpack_surface_c_move_alloc ! Move ownership [bind(C)]

    ! Method wrappers
    public :: fitpack_surface_c_fit
    public :: fitpack_surface_c_interpolate
    public :: fitpack_surface_c_least_squares
    public :: fitpack_surface_c_surface_eval_one
    public :: fitpack_surface_c_surface_derivatives_one
    public :: fitpack_surface_c_cross_section
    public :: fitpack_surface_c_derivative_spline
    public :: fitpack_surface_c_comm_size
    public :: fitpack_surface_c_mse
    public :: fitpack_surface_c_core_comm_size
    public :: fitpack_surface_c_destroy_base
    public :: fitpack_surface_c_new_points
    public :: fitpack_surface_c_new_fit
    public :: fitpack_surface_c_surface_eval_many
    public :: fitpack_surface_c_surface_eval_gridded
    public :: fitpack_surface_c_surface_derivatives_many
    public :: fitpack_surface_c_surface_derivatives_gridded
    public :: fitpack_surface_c_integral
    public :: fitpack_surface_c_comm_pack
    public :: fitpack_surface_c_comm_expand
    public :: fitpack_surface_c_core_comm_pack
    public :: fitpack_surface_c_core_comm_expand
    public :: fitpack_surface_c_getcomp_x
    public :: fitpack_surface_c_getcomp_y
    public :: fitpack_surface_c_getcomp_z
    public :: fitpack_surface_c_getcomp_w
    public :: fitpack_surface_c_getcomp_wrk2
    public :: fitpack_surface_c_getcomp_t
    public :: fitpack_surface_c_ref_m
    public :: fitpack_surface_c_ref_nmax
    public :: fitpack_surface_c_ref_lwrk2
    public :: fitpack_surface_c_ref_bc

    ! ===========================================================================================
    ! SECTION 4: Fortran-Side Interfaces (NOT bind(C))
    ! ===========================================================================================

    !> Get pointer to embedded Fortran object
    interface f_pointer
        module procedure fitpack_surface_c_get_pointer
    end interface f_pointer

    !> Pointer-identity check (Fortran-side wrapper of the bind(C) function).
    interface f_associated
        module procedure fitpack_surface_c_is_same
    end interface f_associated

    !> Construct from Fortran object
    interface fitpack_surface_c
        module procedure fitpack_surface_c_f2c
        module procedure fitpack_surface_c_fptr2c
    end interface fitpack_surface_c

contains

    ! ===========================================================================================
    ! SECTION 5: Core Memory Management (ALWAYS IMPLEMENT)
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_surface_c_is_same(this, that) &
            bind(C, name='fitpack_surface_c_is_same') result(same)
        type(fitpack_surface_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_surface_c_is_same

    !> Get pointer to embedded Fortran object (non-allocating)
    function fitpack_surface_c_get_pointer(this) result(fptr)
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_surface_c_get_pointer

    !> Get/allocate embedded Fortran object (internal, not bind(C))
    !> Returns .true. on success, .false. on allocation failure
    function fitpack_surface_c_pointer(this, fthis) result(success)
        type(fitpack_surface_c), intent(inout) :: this
        type(fitpack_surface), pointer, intent(inout) :: fthis
        logical :: success
        integer :: ierr

        success = .true.

        ! Try non-allocating get first
        fthis => f_pointer(this)

        ! Already associated: do nothing
        if (associated(fthis)) return

        ! Allocate if necessary
        allocate(fthis, stat=ierr)
        if (ierr /= 0) then
            success = .false.
            return
        end if

        ! Initialize and locate
        this%cptr = c_loc(fthis)
        this%is_pointer = .false._c_bool
        this%name_cptr = c_loc(fitpack_surface_c_typename)
    end function fitpack_surface_c_pointer

    !> [bind(C)] Allocate new object
    !> @param status Optional. If absent, allocation failure triggers error stop.
    subroutine fitpack_surface_c_allocate(this, status) bind(C, name='fitpack_surface_c_allocate')
        type(fitpack_surface_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_surface), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_surface_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_surface_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_surface_c_allocate

    !> [bind(C)] Deallocate object
    !> @param status Optional. If absent, deallocation failure triggers error stop.
    subroutine fitpack_surface_c_destroy(this, status) bind(C, name='fitpack_surface_c_destroy')
        type(fitpack_surface_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_surface), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => f_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_surface_c_destroy: deallocate failed')
        end if
        this = fitpack_surface_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_surface_c_destroy

    !> [bind(C)] Copy. By default, mirrors source ownership: a view source
    !> yields a view (shallow handle copy); an owned source yields a deep
    !> copy via Fortran intrinsic assignment. Pass deep_copy=.true. to force
    !> a deep copy regardless of the source's ownership.
    !> @param deep_copy When .true., always allocate a fresh Fortran object
    !>                  and deep-copy data, even if the source is a view.
    !> @param status    Optional. If absent, errors trigger error stop.
    subroutine fitpack_surface_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_surface_c_copy')
        type(fitpack_surface_c), intent(inout) :: this
        type(fitpack_surface_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_surface), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => f_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_surface_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_surface_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_surface_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_surface_c_copy

    !> [bind(C)] Shallow copy (pointer semantics — Fortran "associate" construct)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_surface_c_associate(this, that, status) &
            bind(C, name='fitpack_surface_c_associate')
        type(fitpack_surface_c), intent(inout) :: this
        type(fitpack_surface_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_surface_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_surface_c_associate

    !> [bind(C)] Move allocation (transfer ownership)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_surface_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_surface_c_move_alloc')
        type(fitpack_surface_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_surface_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_surface_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_surface_c_move_alloc

    ! ===========================================================================================
    ! SECTION 6: Fortran Constructors (NOT bind(C))
    ! ===========================================================================================

    !> Create from Fortran object (owns copy)
    type(fitpack_surface_c) function fitpack_surface_c_f2c(f) result(c)
        type(fitpack_surface), intent(in) :: f
        type(fitpack_surface), pointer :: fptr
        logical :: ok

        ok = fitpack_surface_c_pointer(c, fptr)
        if (ok) fptr = f
    end function fitpack_surface_c_f2c

    !> Create from Fortran object (optionally as pointer)
    type(fitpack_surface_c) function fitpack_surface_c_fptr2c(f, want_pointer) result(c)
        type(fitpack_surface), intent(inout), target :: f
        logical, intent(in) :: want_pointer

        if (want_pointer) then
            c%cptr = c_loc(f)
            c%is_pointer = .true._c_bool
        else
            c = fitpack_surface_c_f2c(f)
        end if
    end function fitpack_surface_c_fptr2c

    ! ===========================================================================================
    ! SECTION 7: Method Wrappers (bind(C))
    ! ===========================================================================================
    !> fit
    integer(c_int32_t) function fitpack_surface_c_fit(this, smoothing, order, keep_knots) bind(C, name='fitpack_surface_c_fit')
        type(fitpack_surface_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_surface), pointer :: fthis
        logical :: f_keep_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_surface_c_fit = fthis%fit(smoothing=smoothing, order=order, keep_knots=f_keep_knots)
            else
                fitpack_surface_c_fit = fthis%fit(smoothing=smoothing, order=order)
            end if
        else
            fitpack_surface_c_fit = 0
        end if
    end function fitpack_surface_c_fit

    !> interpolate
    integer(c_int32_t) function fitpack_surface_c_interpolate(this, reset_knots) bind(C, name='fitpack_surface_c_interpolate')
        type(fitpack_surface_c), intent(in) :: this
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_surface), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_surface_c_interpolate = fthis%interpolate(reset_knots=f_reset_knots)
            else
                fitpack_surface_c_interpolate = fthis%interpolate()
            end if
        else
            fitpack_surface_c_interpolate = 0
        end if
    end function fitpack_surface_c_interpolate

    !> least_squares
    integer(c_int32_t) function fitpack_surface_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_surface_c_least_squares')
        type(fitpack_surface_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_surface), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_surface_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_surface_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_surface_c_least_squares = 0
        end if
    end function fitpack_surface_c_least_squares

    !> surface_eval_one
    real(c_double) function fitpack_surface_c_surface_eval_one(this, x, y, ierr) bind(C, name='fitpack_surface_c_surface_eval_one')
        type(fitpack_surface_c), intent(in) :: this
        real(c_double), intent(in), value :: x, y
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_surface_eval_one = fthis%eval(x, y, ierr)
        else
            fitpack_surface_c_surface_eval_one = 0.0_c_double
        end if
    end function fitpack_surface_c_surface_eval_one

    !> surface_derivatives_one
    real(c_double) function fitpack_surface_c_surface_derivatives_one(this, x, y, dx, dy, ierr) bind(C,  &
        & name='fitpack_surface_c_surface_derivatives_one')
        type(fitpack_surface_c), intent(in) :: this
        real(c_double), intent(in), value :: x, y
        integer(c_int32_t), intent(in), value :: dx, dy
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_surface_derivatives_one = fthis%dfdx(x, y, dx, dy, ierr)
        else
            fitpack_surface_c_surface_derivatives_one = 0.0_c_double
        end if
    end function fitpack_surface_c_surface_derivatives_one

    !> cross_section
    subroutine fitpack_surface_c_cross_section(this, u, along_y, ierr, result) bind(C, name='fitpack_surface_c_cross_section')
        type(fitpack_surface_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        logical(c_bool), intent(in), value :: along_y
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_curve_c), intent(out) :: result
        type(fitpack_curve), pointer :: fresult
        type(fitpack_surface), pointer :: fthis
        logical :: f_along_y

        fthis => f_pointer(this)
        if (associated(fthis)) then
            f_along_y = logical(along_y)
            ! Sourced allocation (F2018 6.7.1.2) bypasses user-defined assignment(=); see CLAUDE.md §8.10.
            allocate(fresult, source=fthis%cross_section(u=u, along_y=f_along_y, ierr=ierr))
            result%cptr = c_loc(fresult)
            result%is_pointer = .false._c_bool
            result%name_cptr = f_type_name(fitpack_curve_c_null)
        end if
    end subroutine fitpack_surface_c_cross_section

    !> derivative_spline
    subroutine fitpack_surface_c_derivative_spline(this, nux, nuy, ierr, result) bind(C, name='fitpack_surface_c_derivative_spline')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: nux, nuy
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_surface_c), intent(out) :: result
        type(fitpack_surface), pointer :: fresult, fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            ! Sourced allocation (F2018 6.7.1.2) bypasses user-defined assignment(=); see CLAUDE.md §8.10.
            allocate(fresult, source=fthis%derivative_spline(nux=nux, nuy=nuy, ierr=ierr))
            result%cptr = c_loc(fresult)
            result%is_pointer = .false._c_bool
            result%name_cptr = c_loc(fitpack_surface_c_typename)
        end if
    end subroutine fitpack_surface_c_derivative_spline

    !> comm_size
    integer(c_int32_t) function fitpack_surface_c_comm_size(this) bind(C, name='fitpack_surface_c_comm_size')
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_comm_size = fthis%comm_size()
        else
            fitpack_surface_c_comm_size = 0
        end if
    end function fitpack_surface_c_comm_size

    !> mse
    real(c_double) function fitpack_surface_c_mse(this) bind(C, name='fitpack_surface_c_mse')
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_mse = fthis%mse()
        else
            fitpack_surface_c_mse = 0.0_c_double
        end if
    end function fitpack_surface_c_mse

    !> core_comm_size
    integer(c_int32_t) function fitpack_surface_c_core_comm_size(this) bind(C, name='fitpack_surface_c_core_comm_size')
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_surface_c_core_comm_size = 0
        end if
    end function fitpack_surface_c_core_comm_size

    !> destroy_base
    subroutine fitpack_surface_c_destroy_base(this) bind(C, name='fitpack_surface_c_destroy_base')
        type(fitpack_surface_c), intent(inout) :: this
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_surface_c_destroy_base

    !> new_points
    subroutine fitpack_surface_c_new_points(this, n, x, y, z, w) bind(C, name='fitpack_surface_c_new_points')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), intent(in) :: z(n)
        real(c_double), optional, intent(in) :: w(n)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(x, y, z, w)
        end if
    end subroutine fitpack_surface_c_new_points

    !> new_fit
    integer(c_int32_t) function fitpack_surface_c_new_fit(this, n, x, y, z, w, smoothing, order) bind(C,  &
        & name='fitpack_surface_c_new_fit')
        type(fitpack_surface_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), intent(in) :: z(n)
        real(c_double), optional, intent(in) :: w(n)
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_new_fit = fthis%new_fit(x, y, z, w, smoothing, order)
        else
            fitpack_surface_c_new_fit = 0
        end if
    end function fitpack_surface_c_new_fit

    !> surface_eval_many
    subroutine fitpack_surface_c_surface_eval_many(this, n, x, y, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_surface_c_surface_eval_many')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(x=x, y=y, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_surface_c_surface_eval_many

    !> surface_eval_gridded
    subroutine fitpack_surface_c_surface_eval_gridded(this, n, x, y, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_surface_c_surface_eval_gridded')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval_ongrid(x=x, y=y, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = reshape(fresult, [min(n_result, max_size)])
        end if
    end subroutine fitpack_surface_c_surface_eval_gridded

    !> surface_derivatives_many
    subroutine fitpack_surface_c_surface_derivatives_many(this, n, x, y, dx, dy, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_surface_c_surface_derivatives_many')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, dx, dy
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(x=x, y=y, dx=dx, dy=dy, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_surface_c_surface_derivatives_many

    !> surface_derivatives_gridded
    subroutine fitpack_surface_c_surface_derivatives_gridded(this, n, x, y, dx, dy, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_surface_c_surface_derivatives_gridded')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, dx, dy
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_ongrid(x=x, y=y, dx=dx, dy=dy, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = reshape(fresult, [min(n_result, max_size)])
        end if
    end subroutine fitpack_surface_c_surface_derivatives_gridded

    !> integral
    real(c_double) function fitpack_surface_c_integral(this, n, lower, upper) bind(C, name='fitpack_surface_c_integral')
        type(fitpack_surface_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: lower(n)
        real(c_double), intent(in) :: upper(n)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_integral = fthis%integral(lower, upper)
        else
            fitpack_surface_c_integral = 0.0_c_double
        end if
    end function fitpack_surface_c_integral

    !> comm_pack
    subroutine fitpack_surface_c_comm_pack(this, n, buffer) bind(C, name='fitpack_surface_c_comm_pack')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_surface_c_comm_pack

    !> comm_expand
    subroutine fitpack_surface_c_comm_expand(this, n, buffer) bind(C, name='fitpack_surface_c_comm_expand')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_surface_c_comm_expand

    !> core_comm_pack
    subroutine fitpack_surface_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_surface_c_core_comm_pack')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_surface_c_core_comm_pack

    !> core_comm_expand
    subroutine fitpack_surface_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_surface_c_core_comm_expand')
        type(fitpack_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_surface_c_core_comm_expand

    ! ===========================================================================================
    ! SECTION 8: Component Array Accessors (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'x': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_surface_c_getcomp_x(this, ptr, extents) &
        bind(C, name='fitpack_surface_c_getcomp_x')
        type(fitpack_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%x)) then
                extents = int(shape(fthis%x), kind=c_int64_t)
                if (size(fthis%x) > 0) ptr = c_loc(fthis%x)
            end if
        end if
    end subroutine fitpack_surface_c_getcomp_x

    !> Raw view of component 'y': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_surface_c_getcomp_y(this, ptr, extents) &
        bind(C, name='fitpack_surface_c_getcomp_y')
        type(fitpack_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%y)) then
                extents = int(shape(fthis%y), kind=c_int64_t)
                if (size(fthis%y) > 0) ptr = c_loc(fthis%y)
            end if
        end if
    end subroutine fitpack_surface_c_getcomp_y

    !> Raw view of component 'z': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_surface_c_getcomp_z(this, ptr, extents) &
        bind(C, name='fitpack_surface_c_getcomp_z')
        type(fitpack_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%z)) then
                extents = int(shape(fthis%z), kind=c_int64_t)
                if (size(fthis%z) > 0) ptr = c_loc(fthis%z)
            end if
        end if
    end subroutine fitpack_surface_c_getcomp_z

    !> Raw view of component 'w': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_surface_c_getcomp_w(this, ptr, extents) &
        bind(C, name='fitpack_surface_c_getcomp_w')
        type(fitpack_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%w)) then
                extents = int(shape(fthis%w), kind=c_int64_t)
                if (size(fthis%w) > 0) ptr = c_loc(fthis%w)
            end if
        end if
    end subroutine fitpack_surface_c_getcomp_w

    !> Raw view of component 'wrk2': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_surface_c_getcomp_wrk2(this, ptr, extents) &
        bind(C, name='fitpack_surface_c_getcomp_wrk2')
        type(fitpack_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%wrk2)) then
                extents = int(shape(fthis%wrk2), kind=c_int64_t)
                if (size(fthis%wrk2) > 0) ptr = c_loc(fthis%wrk2)
            end if
        end if
    end subroutine fitpack_surface_c_getcomp_wrk2

    !> Raw view of component 't': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_surface_c_getcomp_t(this, ptr, extents) &
        bind(C, name='fitpack_surface_c_getcomp_t')
        type(fitpack_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%t)) then
                extents = int(shape(fthis%t), kind=c_int64_t)
                if (size(fthis%t) > 0) ptr = c_loc(fthis%t)
            end if
        end if
    end subroutine fitpack_surface_c_getcomp_t

    ! ===========================================================================================
    ! SECTION 9: Scalar Property Accessors
    ! ===========================================================================================

    !> Get pointer to scalar property 'm'
    type(c_ptr) function fitpack_surface_c_ref_m(this) &
        bind(C, name='fitpack_surface_c_ref_m')
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_ref_m = c_loc(fthis%m)
        else
            fitpack_surface_c_ref_m = c_null_ptr
        end if
    end function fitpack_surface_c_ref_m

    !> Get pointer to scalar property 'nmax'
    type(c_ptr) function fitpack_surface_c_ref_nmax(this) &
        bind(C, name='fitpack_surface_c_ref_nmax')
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_ref_nmax = c_loc(fthis%nmax)
        else
            fitpack_surface_c_ref_nmax = c_null_ptr
        end if
    end function fitpack_surface_c_ref_nmax

    !> Get pointer to scalar property 'lwrk2'
    type(c_ptr) function fitpack_surface_c_ref_lwrk2(this) &
        bind(C, name='fitpack_surface_c_ref_lwrk2')
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_ref_lwrk2 = c_loc(fthis%lwrk2)
        else
            fitpack_surface_c_ref_lwrk2 = c_null_ptr
        end if
    end function fitpack_surface_c_ref_lwrk2

    !> Get pointer to scalar property 'bc'
    type(c_ptr) function fitpack_surface_c_ref_bc(this) &
        bind(C, name='fitpack_surface_c_ref_bc')
        type(fitpack_surface_c), intent(in) :: this
        type(fitpack_surface), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_surface_c_ref_bc = c_loc(fthis%bc)
        else
            fitpack_surface_c_ref_bc = c_null_ptr
        end if
    end function fitpack_surface_c_ref_bc

end module fitpack_surfaces_c
