! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_curves_c.f90 (module fitpack_curves_c)
!> @brief C interface to module fitpack_curves
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

module fitpack_curves_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_curves, only: fitpack_curve, fitpack_periodic_curve
    use fitpack_curves_c_types
    use fitpack_fx_status, only: fx_status, fx_status_ok, FX_ERROR_ALLOCATION, &
                        FX_ERROR_DEALLOCATION, handle_error
    use, intrinsic :: iso_c_binding
    implicit none(type, external)
    private

    ! ===========================================================================================
    ! SECTION 2: Public Exports
    ! ===========================================================================================

    public :: f_pointer
    public :: f_associated

    ! --- fitpack_curve ---
    public :: fitpack_curve_c_is_same
    public :: fitpack_curve_c_allocate
    public :: fitpack_curve_c_destroy
    public :: fitpack_curve_c_copy
    public :: fitpack_curve_c_associate
    public :: fitpack_curve_c_move_alloc
    public :: fitpack_curve_c_fit
    public :: fitpack_curve_c_interpolate
    public :: fitpack_curve_c_least_squares
    public :: fitpack_curve_c_curve_eval_one
    public :: fitpack_curve_c_curve_eval_one_noerr
    public :: fitpack_curve_c_integral
    public :: fitpack_curve_c_curve_derivative
    public :: fitpack_curve_c_curve_insert_knot_one
    public :: fitpack_curve_c_comm_size
    public :: fitpack_curve_c_mse
    public :: fitpack_curve_c_core_comm_size
    public :: fitpack_curve_c_destroy_base
    public :: fitpack_curve_c_new_points
    public :: fitpack_curve_c_new_fit
    public :: fitpack_curve_c_curve_eval_many
    public :: fitpack_curve_c_curve_eval_many_pure
    public :: fitpack_curve_c_fourier_coefficients
    public :: fitpack_curve_c_zeros
    public :: fitpack_curve_c_curve_derivatives
    public :: fitpack_curve_c_curve_all_derivatives
    public :: fitpack_curve_c_curve_all_derivatives_pure
    public :: fitpack_curve_c_curve_insert_knot_many
    public :: fitpack_curve_c_comm_pack
    public :: fitpack_curve_c_comm_expand
    public :: fitpack_curve_c_core_comm_pack
    public :: fitpack_curve_c_core_comm_expand
    public :: fitpack_curve_c_getcomp_x
    public :: fitpack_curve_c_getcomp_y
    public :: fitpack_curve_c_getcomp_sp
    public :: fitpack_curve_c_getcomp_w
    public :: fitpack_curve_c_getcomp_wrk_fou
    public :: fitpack_curve_c_getcomp_t
    public :: fitpack_curve_c_ref_m
    public :: fitpack_curve_c_ref_order
    public :: fitpack_curve_c_ref_xleft
    public :: fitpack_curve_c_ref_xright
    public :: fitpack_curve_c_ref_nest
    public :: fitpack_curve_c_ref_bc
    public :: fitpack_curve_c_ref_knots

    ! --- fitpack_periodic_curve ---
    public :: fitpack_periodic_curve_c_is_same
    public :: fitpack_periodic_curve_c_allocate
    public :: fitpack_periodic_curve_c_destroy
    public :: fitpack_periodic_curve_c_copy
    public :: fitpack_periodic_curve_c_associate
    public :: fitpack_periodic_curve_c_move_alloc
    public :: fitpack_periodic_curve_c_fit
    public :: fitpack_periodic_curve_c_interpolate
    public :: fitpack_periodic_curve_c_least_squares
    public :: fitpack_periodic_curve_c_curve_eval_one
    public :: fitpack_periodic_curve_c_curve_eval_one_noerr
    public :: fitpack_periodic_curve_c_integral
    public :: fitpack_periodic_curve_c_curve_derivative
    public :: fitpack_periodic_curve_c_curve_insert_knot_one
    public :: fitpack_periodic_curve_c_comm_size
    public :: fitpack_periodic_curve_c_mse
    public :: fitpack_periodic_curve_c_core_comm_size
    public :: fitpack_periodic_curve_c_destroy_base
    public :: fitpack_periodic_curve_c_new_points
    public :: fitpack_periodic_curve_c_new_fit
    public :: fitpack_periodic_curve_c_curve_eval_many
    public :: fitpack_periodic_curve_c_curve_eval_many_pure
    public :: fitpack_periodic_curve_c_fourier_coefficients
    public :: fitpack_periodic_curve_c_zeros
    public :: fitpack_periodic_curve_c_curve_derivatives
    public :: fitpack_periodic_curve_c_curve_all_derivatives
    public :: fitpack_periodic_curve_c_curve_all_derivatives_pure
    public :: fitpack_periodic_curve_c_curve_insert_knot_many
    public :: fitpack_periodic_curve_c_comm_pack
    public :: fitpack_periodic_curve_c_comm_expand
    public :: fitpack_periodic_curve_c_core_comm_pack
    public :: fitpack_periodic_curve_c_core_comm_expand

    interface f_pointer
        module procedure fitpack_curve_c_get_pointer
        module procedure fitpack_periodic_curve_c_get_pointer
    end interface f_pointer

    !> Pointer-identity check. Returns true iff `a` and `b` view the same
    !> underlying Fortran object (their internal c_ptrs match and are non-null).
    !> Useful for C-API consumers that want object identity without exposing
    !> raw pointers — e.g. GUI panels checking whether an input view points
    !> at the same Fortran allocation as the project that owns it.
    interface f_associated
        module procedure fitpack_curve_c_is_same
        module procedure fitpack_periodic_curve_c_is_same
    end interface f_associated

contains

    ! ===========================================================================================
    ! FITPACK_CURVE: Core Memory Management
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_curve_c_is_same(this, that) &
            bind(C, name='fitpack_curve_c_is_same') result(same)
        type(fitpack_curve_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_curve_c_is_same

    function fitpack_curve_c_get_pointer(this) result(fptr)
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_curve_c_get_pointer

    function fitpack_curve_c_pointer(this, fthis) result(success)
        type(fitpack_curve_c), intent(inout) :: this
        type(fitpack_curve), pointer, intent(inout) :: fthis
        logical :: success
        integer :: ierr

        success = .true.
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) return

        allocate(fthis, stat=ierr)
        if (ierr /= 0) then
            success = .false.
            return
        end if

        this%cptr = c_loc(fthis)
        this%is_pointer = .false._c_bool
        this%name_cptr = c_loc(fitpack_curve_c_typename)
    end function fitpack_curve_c_pointer

    subroutine fitpack_curve_c_allocate(this, status) &
            bind(C, name='fitpack_curve_c_allocate')
        type(fitpack_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_curve), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_curve_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_curve_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_curve_c_allocate

    subroutine fitpack_curve_c_destroy(this, status) &
            bind(C, name='fitpack_curve_c_destroy')
        type(fitpack_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_curve), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_curve_c_destroy: deallocate failed')
        end if
        this = fitpack_curve_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_curve_c_destroy

    subroutine fitpack_curve_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_curve_c_copy')
        type(fitpack_curve_c), intent(inout) :: this
        type(fitpack_curve_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_curve), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => fitpack_curve_c_get_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_curve_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_curve_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_curve_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_curve_c_copy

    subroutine fitpack_curve_c_associate(this, that, status) &
            bind(C, name='fitpack_curve_c_associate')
        type(fitpack_curve_c), intent(inout) :: this
        type(fitpack_curve_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_curve_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_curve_c_associate

    subroutine fitpack_curve_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_curve_c_move_alloc')
        type(fitpack_curve_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_curve_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_curve_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_curve_c_move_alloc

    ! ===========================================================================================
    ! FITPACK_CURVE: Method Wrappers
    ! ===========================================================================================
    integer(c_int32_t) function fitpack_curve_c_fit(this, smoothing, order, keep_knots) bind(C, name='fitpack_curve_c_fit')
        type(fitpack_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_curve), pointer :: fthis
        logical :: f_keep_knots

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_curve_c_fit = fthis%fit(smoothing=smoothing, order=order, keep_knots=f_keep_knots)
            else
                fitpack_curve_c_fit = fthis%fit(smoothing=smoothing, order=order)
            end if
        else
            fitpack_curve_c_fit = 0
        end if
    end function fitpack_curve_c_fit

    integer(c_int32_t) function fitpack_curve_c_interpolate(this, order, reset_knots) bind(C, name='fitpack_curve_c_interpolate')
        type(fitpack_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_curve_c_interpolate = fthis%interpolate(order=order, reset_knots=f_reset_knots)
            else
                fitpack_curve_c_interpolate = fthis%interpolate(order=order)
            end if
        else
            fitpack_curve_c_interpolate = 0
        end if
    end function fitpack_curve_c_interpolate

    integer(c_int32_t) function fitpack_curve_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_curve_c_least_squares')
        type(fitpack_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_curve_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_curve_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_curve_c_least_squares = 0
        end if
    end function fitpack_curve_c_least_squares

    real(c_double) function fitpack_curve_c_curve_eval_one(this, x, ierr) bind(C, name='fitpack_curve_c_curve_eval_one')
        type(fitpack_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(out) :: ierr
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_curve_eval_one = fthis%eval(x, ierr)
        else
            fitpack_curve_c_curve_eval_one = 0.0_c_double
        end if
    end function fitpack_curve_c_curve_eval_one

    real(c_double) function fitpack_curve_c_curve_eval_one_noerr(this, x) bind(C, name='fitpack_curve_c_curve_eval_one_noerr')
        type(fitpack_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: x
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_curve_eval_one_noerr = fthis%eval(x)
        else
            fitpack_curve_c_curve_eval_one_noerr = 0.0_c_double
        end if
    end function fitpack_curve_c_curve_eval_one_noerr

    real(c_double) function fitpack_curve_c_integral(this, from, to) bind(C, name='fitpack_curve_c_integral')
        type(fitpack_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: from, to
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_integral = fthis%integral(from, to)
        else
            fitpack_curve_c_integral = 0.0_c_double
        end if
    end function fitpack_curve_c_integral

    real(c_double) function fitpack_curve_c_curve_derivative(this, x, order, ierr) bind(C, name='fitpack_curve_c_curve_derivative')
        type(fitpack_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(in), value :: order
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_curve_derivative = fthis%dfdx(x, order, ierr)
        else
            fitpack_curve_c_curve_derivative = 0.0_c_double
        end if
    end function fitpack_curve_c_curve_derivative

    subroutine fitpack_curve_c_curve_insert_knot_one(this, x, ierr) bind(C, name='fitpack_curve_c_curve_insert_knot_one')
        type(fitpack_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%insert_knot(x, ierr)
        end if
    end subroutine fitpack_curve_c_curve_insert_knot_one

    integer(c_int32_t) function fitpack_curve_c_comm_size(this) bind(C, name='fitpack_curve_c_comm_size')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_comm_size = fthis%comm_size()
        else
            fitpack_curve_c_comm_size = 0
        end if
    end function fitpack_curve_c_comm_size

    real(c_double) function fitpack_curve_c_mse(this) bind(C, name='fitpack_curve_c_mse')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_mse = fthis%mse()
        else
            fitpack_curve_c_mse = 0.0_c_double
        end if
    end function fitpack_curve_c_mse

    integer(c_int32_t) function fitpack_curve_c_core_comm_size(this) bind(C, name='fitpack_curve_c_core_comm_size')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_curve_c_core_comm_size = 0
        end if
    end function fitpack_curve_c_core_comm_size

    subroutine fitpack_curve_c_destroy_base(this) bind(C, name='fitpack_curve_c_destroy_base')
        type(fitpack_curve_c), intent(inout) :: this
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_curve_c_destroy_base

    subroutine fitpack_curve_c_new_points(this, n, x, y, w) bind(C, name='fitpack_curve_c_new_points')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), optional, intent(in) :: w(n)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(x, y, w)
        end if
    end subroutine fitpack_curve_c_new_points

    integer(c_int32_t) function fitpack_curve_c_new_fit(this, n, x, y, w, smoothing, order) bind(C, name='fitpack_curve_c_new_fit')
        type(fitpack_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), optional, intent(in) :: w(n)
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_new_fit = fthis%new_fit(x, y, w, smoothing, order)
        else
            fitpack_curve_c_new_fit = 0
        end if
    end function fitpack_curve_c_new_fit

    subroutine fitpack_curve_c_curve_eval_many(this, n, x, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_curve_c_curve_eval_many')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        integer(c_int32_t), intent(out) :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(x, ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_curve_c_curve_eval_many

    subroutine fitpack_curve_c_curve_eval_many_pure(this, n, x, result, n_result, max_size) bind(C,  &
        & name='fitpack_curve_c_curve_eval_many_pure')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(x)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_curve_c_curve_eval_many_pure

    subroutine fitpack_curve_c_fourier_coefficients(this, n, alpha, a, b, ierr) bind(C, name='fitpack_curve_c_fourier_coefficients')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: alpha(n)
        real(c_double), intent(inout) :: a(n)
        real(c_double), intent(inout) :: b(n)
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%fourier_coefficients(alpha, a, b, ierr)
        end if
    end subroutine fitpack_curve_c_fourier_coefficients

    subroutine fitpack_curve_c_zeros(this, ierr, result, n_result, max_size) bind(C, name='fitpack_curve_c_zeros')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%zeros(ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_curve_c_zeros

    subroutine fitpack_curve_c_curve_derivatives(this, n, x, order, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_curve_c_curve_derivatives')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, order
        real(c_double), intent(in) :: x(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(x, order, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_curve_c_curve_derivatives

    subroutine fitpack_curve_c_curve_all_derivatives(this, x, ierr, result, n_result) bind(C,  &
        & name='fitpack_curve_c_curve_all_derivatives')
        type(fitpack_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(out) :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_all(x, ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_curve_c_curve_all_derivatives

    subroutine fitpack_curve_c_curve_all_derivatives_pure(this, x, result, n_result) bind(C,  &
        & name='fitpack_curve_c_curve_all_derivatives_pure')
        type(fitpack_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: x
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_all(x)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_curve_c_curve_all_derivatives_pure

    subroutine fitpack_curve_c_curve_insert_knot_many(this, n, x, ierr) bind(C, name='fitpack_curve_c_curve_insert_knot_many')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%insert_knot(x, ierr)
        end if
    end subroutine fitpack_curve_c_curve_insert_knot_many

    subroutine fitpack_curve_c_comm_pack(this, n, buffer) bind(C, name='fitpack_curve_c_comm_pack')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_curve_c_comm_pack

    subroutine fitpack_curve_c_comm_expand(this, n, buffer) bind(C, name='fitpack_curve_c_comm_expand')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_curve_c_comm_expand

    subroutine fitpack_curve_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_curve_c_core_comm_pack')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_curve_c_core_comm_pack

    subroutine fitpack_curve_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_curve_c_core_comm_expand')
        type(fitpack_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_curve), pointer :: fthis

        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_curve_c_core_comm_expand

    ! ===========================================================================================
    ! Component Array Accessors for fitpack_curve (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'x': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_curve_c_getcomp_x(this, ptr, extents) &
        bind(C, name='fitpack_curve_c_getcomp_x')
        type(fitpack_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%x)) then
                extents = int(shape(fthis%x), kind=c_int64_t)
                if (size(fthis%x) > 0) ptr = c_loc(fthis%x)
            end if
        end if
    end subroutine fitpack_curve_c_getcomp_x

    !> Raw view of component 'y': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_curve_c_getcomp_y(this, ptr, extents) &
        bind(C, name='fitpack_curve_c_getcomp_y')
        type(fitpack_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%y)) then
                extents = int(shape(fthis%y), kind=c_int64_t)
                if (size(fthis%y) > 0) ptr = c_loc(fthis%y)
            end if
        end if
    end subroutine fitpack_curve_c_getcomp_y

    !> Raw view of component 'sp': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_curve_c_getcomp_sp(this, ptr, extents) &
        bind(C, name='fitpack_curve_c_getcomp_sp')
        type(fitpack_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%sp)) then
                extents = int(shape(fthis%sp), kind=c_int64_t)
                if (size(fthis%sp) > 0) ptr = c_loc(fthis%sp)
            end if
        end if
    end subroutine fitpack_curve_c_getcomp_sp

    !> Raw view of component 'w': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_curve_c_getcomp_w(this, ptr, extents) &
        bind(C, name='fitpack_curve_c_getcomp_w')
        type(fitpack_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%w)) then
                extents = int(shape(fthis%w), kind=c_int64_t)
                if (size(fthis%w) > 0) ptr = c_loc(fthis%w)
            end if
        end if
    end subroutine fitpack_curve_c_getcomp_w

    !> Raw view of component 'wrk_fou': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_curve_c_getcomp_wrk_fou(this, ptr, extents) &
        bind(C, name='fitpack_curve_c_getcomp_wrk_fou')
        type(fitpack_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%wrk_fou)) then
                extents = int(shape(fthis%wrk_fou), kind=c_int64_t)
                if (size(fthis%wrk_fou) > 0) ptr = c_loc(fthis%wrk_fou)
            end if
        end if
    end subroutine fitpack_curve_c_getcomp_wrk_fou

    !> Raw view of component 't': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_curve_c_getcomp_t(this, ptr, extents) &
        bind(C, name='fitpack_curve_c_getcomp_t')
        type(fitpack_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%t)) then
                extents = int(shape(fthis%t), kind=c_int64_t)
                if (size(fthis%t) > 0) ptr = c_loc(fthis%t)
            end if
        end if
    end subroutine fitpack_curve_c_getcomp_t

    ! ===========================================================================================
    ! Scalar Property Accessors for fitpack_curve
    ! ===========================================================================================

    !> Get pointer to scalar property 'm'
    type(c_ptr) function fitpack_curve_c_ref_m(this) &
        bind(C, name='fitpack_curve_c_ref_m')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_ref_m = c_loc(fthis%m)
        else
            fitpack_curve_c_ref_m = c_null_ptr
        end if
    end function fitpack_curve_c_ref_m

    !> Get pointer to scalar property 'order'
    type(c_ptr) function fitpack_curve_c_ref_order(this) &
        bind(C, name='fitpack_curve_c_ref_order')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_ref_order = c_loc(fthis%order)
        else
            fitpack_curve_c_ref_order = c_null_ptr
        end if
    end function fitpack_curve_c_ref_order

    !> Get pointer to scalar property 'xleft'
    type(c_ptr) function fitpack_curve_c_ref_xleft(this) &
        bind(C, name='fitpack_curve_c_ref_xleft')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_ref_xleft = c_loc(fthis%xleft)
        else
            fitpack_curve_c_ref_xleft = c_null_ptr
        end if
    end function fitpack_curve_c_ref_xleft

    !> Get pointer to scalar property 'xright'
    type(c_ptr) function fitpack_curve_c_ref_xright(this) &
        bind(C, name='fitpack_curve_c_ref_xright')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_ref_xright = c_loc(fthis%xright)
        else
            fitpack_curve_c_ref_xright = c_null_ptr
        end if
    end function fitpack_curve_c_ref_xright

    !> Get pointer to scalar property 'nest'
    type(c_ptr) function fitpack_curve_c_ref_nest(this) &
        bind(C, name='fitpack_curve_c_ref_nest')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_ref_nest = c_loc(fthis%nest)
        else
            fitpack_curve_c_ref_nest = c_null_ptr
        end if
    end function fitpack_curve_c_ref_nest

    !> Get pointer to scalar property 'bc'
    type(c_ptr) function fitpack_curve_c_ref_bc(this) &
        bind(C, name='fitpack_curve_c_ref_bc')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_ref_bc = c_loc(fthis%bc)
        else
            fitpack_curve_c_ref_bc = c_null_ptr
        end if
    end function fitpack_curve_c_ref_bc

    !> Get pointer to scalar property 'knots'
    type(c_ptr) function fitpack_curve_c_ref_knots(this) &
        bind(C, name='fitpack_curve_c_ref_knots')
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_curve_c_ref_knots = c_loc(fthis%knots)
        else
            fitpack_curve_c_ref_knots = c_null_ptr
        end if
    end function fitpack_curve_c_ref_knots

    ! ===========================================================================================
    ! FITPACK_PERIODIC_CURVE: Core Memory Management
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_periodic_curve_c_is_same(this, that) &
            bind(C, name='fitpack_periodic_curve_c_is_same') result(same)
        type(fitpack_periodic_curve_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_periodic_curve_c_is_same

    function fitpack_periodic_curve_c_get_pointer(this) result(fptr)
        type(fitpack_periodic_curve_c), intent(in) :: this
        type(fitpack_periodic_curve), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_periodic_curve_c_get_pointer

    function fitpack_periodic_curve_c_pointer(this, fthis) result(success)
        type(fitpack_periodic_curve_c), intent(inout) :: this
        type(fitpack_periodic_curve), pointer, intent(inout) :: fthis
        logical :: success
        integer :: ierr

        success = .true.
        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) return

        allocate(fthis, stat=ierr)
        if (ierr /= 0) then
            success = .false.
            return
        end if

        this%cptr = c_loc(fthis)
        this%is_pointer = .false._c_bool
        this%name_cptr = c_loc(fitpack_periodic_curve_c_typename)
    end function fitpack_periodic_curve_c_pointer

    subroutine fitpack_periodic_curve_c_allocate(this, status) &
            bind(C, name='fitpack_periodic_curve_c_allocate')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_periodic_curve), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_periodic_curve_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_periodic_curve_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_periodic_curve_c_allocate

    subroutine fitpack_periodic_curve_c_destroy(this, status) &
            bind(C, name='fitpack_periodic_curve_c_destroy')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_periodic_curve), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_periodic_curve_c_destroy: deallocate failed')
        end if
        this = fitpack_periodic_curve_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_periodic_curve_c_destroy

    subroutine fitpack_periodic_curve_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_periodic_curve_c_copy')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        type(fitpack_periodic_curve_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_periodic_curve), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => fitpack_periodic_curve_c_get_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_periodic_curve_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_periodic_curve_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_periodic_curve_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_periodic_curve_c_copy

    subroutine fitpack_periodic_curve_c_associate(this, that, status) &
            bind(C, name='fitpack_periodic_curve_c_associate')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        type(fitpack_periodic_curve_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_periodic_curve_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_periodic_curve_c_associate

    subroutine fitpack_periodic_curve_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_periodic_curve_c_move_alloc')
        type(fitpack_periodic_curve_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_periodic_curve_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_periodic_curve_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_periodic_curve_c_move_alloc

    ! ===========================================================================================
    ! FITPACK_PERIODIC_CURVE: Method Wrappers
    ! ===========================================================================================
    integer(c_int32_t) function fitpack_periodic_curve_c_fit(this, smoothing, order, keep_knots) bind(C,  &
        & name='fitpack_periodic_curve_c_fit')
        type(fitpack_periodic_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_periodic_curve), pointer :: fthis
        logical :: f_keep_knots

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_periodic_curve_c_fit = fthis%fit(smoothing=smoothing, order=order, keep_knots=f_keep_knots)
            else
                fitpack_periodic_curve_c_fit = fthis%fit(smoothing=smoothing, order=order)
            end if
        else
            fitpack_periodic_curve_c_fit = 0
        end if
    end function fitpack_periodic_curve_c_fit

    integer(c_int32_t) function fitpack_periodic_curve_c_interpolate(this, order, reset_knots) bind(C,  &
        & name='fitpack_periodic_curve_c_interpolate')
        type(fitpack_periodic_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_periodic_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_periodic_curve_c_interpolate = fthis%interpolate(order=order, reset_knots=f_reset_knots)
            else
                fitpack_periodic_curve_c_interpolate = fthis%interpolate(order=order)
            end if
        else
            fitpack_periodic_curve_c_interpolate = 0
        end if
    end function fitpack_periodic_curve_c_interpolate

    integer(c_int32_t) function fitpack_periodic_curve_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_periodic_curve_c_least_squares')
        type(fitpack_periodic_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_periodic_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_periodic_curve_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_periodic_curve_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_periodic_curve_c_least_squares = 0
        end if
    end function fitpack_periodic_curve_c_least_squares

    real(c_double) function fitpack_periodic_curve_c_curve_eval_one(this, x, ierr) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_eval_one')
        type(fitpack_periodic_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(out) :: ierr
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_curve_eval_one = fthis%eval(x, ierr)
        else
            fitpack_periodic_curve_c_curve_eval_one = 0.0_c_double
        end if
    end function fitpack_periodic_curve_c_curve_eval_one

    real(c_double) function fitpack_periodic_curve_c_curve_eval_one_noerr(this, x) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_eval_one_noerr')
        type(fitpack_periodic_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: x
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_curve_eval_one_noerr = fthis%eval(x)
        else
            fitpack_periodic_curve_c_curve_eval_one_noerr = 0.0_c_double
        end if
    end function fitpack_periodic_curve_c_curve_eval_one_noerr

    real(c_double) function fitpack_periodic_curve_c_integral(this, from, to) bind(C, name='fitpack_periodic_curve_c_integral')
        type(fitpack_periodic_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: from, to
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_integral = fthis%integral(from, to)
        else
            fitpack_periodic_curve_c_integral = 0.0_c_double
        end if
    end function fitpack_periodic_curve_c_integral

    real(c_double) function fitpack_periodic_curve_c_curve_derivative(this, x, order, ierr) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_derivative')
        type(fitpack_periodic_curve_c), intent(in) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(in), value :: order
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_curve_derivative = fthis%dfdx(x, order, ierr)
        else
            fitpack_periodic_curve_c_curve_derivative = 0.0_c_double
        end if
    end function fitpack_periodic_curve_c_curve_derivative

    subroutine fitpack_periodic_curve_c_curve_insert_knot_one(this, x, ierr) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_insert_knot_one')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%insert_knot(x, ierr)
        end if
    end subroutine fitpack_periodic_curve_c_curve_insert_knot_one

    integer(c_int32_t) function fitpack_periodic_curve_c_comm_size(this) bind(C, name='fitpack_periodic_curve_c_comm_size')
        type(fitpack_periodic_curve_c), intent(in) :: this
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_comm_size = fthis%comm_size()
        else
            fitpack_periodic_curve_c_comm_size = 0
        end if
    end function fitpack_periodic_curve_c_comm_size

    real(c_double) function fitpack_periodic_curve_c_mse(this) bind(C, name='fitpack_periodic_curve_c_mse')
        type(fitpack_periodic_curve_c), intent(in) :: this
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_mse = fthis%mse()
        else
            fitpack_periodic_curve_c_mse = 0.0_c_double
        end if
    end function fitpack_periodic_curve_c_mse

    integer(c_int32_t) function fitpack_periodic_curve_c_core_comm_size(this) bind(C,  &
        & name='fitpack_periodic_curve_c_core_comm_size')
        type(fitpack_periodic_curve_c), intent(in) :: this
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_periodic_curve_c_core_comm_size = 0
        end if
    end function fitpack_periodic_curve_c_core_comm_size

    subroutine fitpack_periodic_curve_c_destroy_base(this) bind(C, name='fitpack_periodic_curve_c_destroy_base')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_periodic_curve_c_destroy_base

    subroutine fitpack_periodic_curve_c_new_points(this, n, x, y, w) bind(C, name='fitpack_periodic_curve_c_new_points')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), optional, intent(in) :: w(n)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(x, y, w)
        end if
    end subroutine fitpack_periodic_curve_c_new_points

    integer(c_int32_t) function fitpack_periodic_curve_c_new_fit(this, n, x, y, w, smoothing, order) bind(C,  &
        & name='fitpack_periodic_curve_c_new_fit')
        type(fitpack_periodic_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), optional, intent(in) :: w(n)
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_periodic_curve_c_new_fit = fthis%new_fit(x, y, w, smoothing, order)
        else
            fitpack_periodic_curve_c_new_fit = 0
        end if
    end function fitpack_periodic_curve_c_new_fit

    subroutine fitpack_periodic_curve_c_curve_eval_many(this, n, x, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_eval_many')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        integer(c_int32_t), intent(out) :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(x, ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_periodic_curve_c_curve_eval_many

    subroutine fitpack_periodic_curve_c_curve_eval_many_pure(this, n, x, result, n_result, max_size) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_eval_many_pure')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(x)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_periodic_curve_c_curve_eval_many_pure

    subroutine fitpack_periodic_curve_c_fourier_coefficients(this, n, alpha, a, b, ierr) bind(C,  &
        & name='fitpack_periodic_curve_c_fourier_coefficients')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: alpha(n)
        real(c_double), intent(inout) :: a(n)
        real(c_double), intent(inout) :: b(n)
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%fourier_coefficients(alpha, a, b, ierr)
        end if
    end subroutine fitpack_periodic_curve_c_fourier_coefficients

    subroutine fitpack_periodic_curve_c_zeros(this, ierr, result, n_result, max_size) bind(C, name='fitpack_periodic_curve_c_zeros')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%zeros(ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_periodic_curve_c_zeros

    subroutine fitpack_periodic_curve_c_curve_derivatives(this, n, x, order, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_derivatives')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, order
        real(c_double), intent(in) :: x(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(x, order, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = fresult(1:min(n_result, max_size))
        end if
    end subroutine fitpack_periodic_curve_c_curve_derivatives

    subroutine fitpack_periodic_curve_c_curve_all_derivatives(this, x, ierr, result, n_result) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_all_derivatives')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: x
        integer(c_int32_t), intent(out) :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_all(x, ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_periodic_curve_c_curve_all_derivatives

    subroutine fitpack_periodic_curve_c_curve_all_derivatives_pure(this, x, result, n_result) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_all_derivatives_pure')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: x
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_all(x)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_periodic_curve_c_curve_all_derivatives_pure

    subroutine fitpack_periodic_curve_c_curve_insert_knot_many(this, n, x, ierr) bind(C,  &
        & name='fitpack_periodic_curve_c_curve_insert_knot_many')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%insert_knot(x, ierr)
        end if
    end subroutine fitpack_periodic_curve_c_curve_insert_knot_many

    subroutine fitpack_periodic_curve_c_comm_pack(this, n, buffer) bind(C, name='fitpack_periodic_curve_c_comm_pack')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_periodic_curve_c_comm_pack

    subroutine fitpack_periodic_curve_c_comm_expand(this, n, buffer) bind(C, name='fitpack_periodic_curve_c_comm_expand')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_periodic_curve_c_comm_expand

    subroutine fitpack_periodic_curve_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_periodic_curve_c_core_comm_pack')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_periodic_curve_c_core_comm_pack

    subroutine fitpack_periodic_curve_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_periodic_curve_c_core_comm_expand')
        type(fitpack_periodic_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_periodic_curve), pointer :: fthis

        fthis => fitpack_periodic_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_periodic_curve_c_core_comm_expand

end module fitpack_curves_c
