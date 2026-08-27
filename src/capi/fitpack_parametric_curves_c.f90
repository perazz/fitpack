! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_parametric_curves_c.f90 (module fitpack_parametric_curves_c)
!> @brief C interface to module fitpack_parametric_curves
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

module fitpack_parametric_curves_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_parametric_curves, only: fitpack_parametric_curve, fitpack_closed_curve, fitpack_constrained_curve
    use fitpack_parametric_curves_c_types
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

    ! --- fitpack_parametric_curve ---
    public :: fitpack_parametric_curve_c_is_same
    public :: fitpack_parametric_curve_c_allocate
    public :: fitpack_parametric_curve_c_destroy
    public :: fitpack_parametric_curve_c_copy
    public :: fitpack_parametric_curve_c_associate
    public :: fitpack_parametric_curve_c_move_alloc
    public :: fitpack_parametric_curve_c_set_default_parameters
    public :: fitpack_parametric_curve_c_fit
    public :: fitpack_parametric_curve_c_interpolate
    public :: fitpack_parametric_curve_c_least_squares
    public :: fitpack_parametric_curve_c_comm_size
    public :: fitpack_parametric_curve_c_mse
    public :: fitpack_parametric_curve_c_core_comm_size
    public :: fitpack_parametric_curve_c_destroy_base
    public :: fitpack_parametric_curve_c_new_points
    public :: fitpack_parametric_curve_c_new_fit
    public :: fitpack_parametric_curve_c_eval_one
    public :: fitpack_parametric_curve_c_eval_many
    public :: fitpack_parametric_curve_c_curve_derivative
    public :: fitpack_parametric_curve_c_curve_derivatives
    public :: fitpack_parametric_curve_c_curve_all_derivatives
    public :: fitpack_parametric_curve_c_comm_pack
    public :: fitpack_parametric_curve_c_comm_expand
    public :: fitpack_parametric_curve_c_core_comm_pack
    public :: fitpack_parametric_curve_c_core_comm_expand
    public :: fitpack_parametric_curve_c_getcomp_x
    public :: fitpack_parametric_curve_c_getcomp_u
    public :: fitpack_parametric_curve_c_getcomp_sp
    public :: fitpack_parametric_curve_c_getcomp_w
    public :: fitpack_parametric_curve_c_getcomp_dd
    public :: fitpack_parametric_curve_c_getcomp_t
    public :: fitpack_parametric_curve_c_ref_m
    public :: fitpack_parametric_curve_c_ref_idim
    public :: fitpack_parametric_curve_c_get_has_params
    public :: fitpack_parametric_curve_c_set_has_params
    public :: fitpack_parametric_curve_c_ref_order
    public :: fitpack_parametric_curve_c_ref_ubegin
    public :: fitpack_parametric_curve_c_ref_uend
    public :: fitpack_parametric_curve_c_ref_nest
    public :: fitpack_parametric_curve_c_ref_knots

    ! --- fitpack_closed_curve ---
    public :: fitpack_closed_curve_c_is_same
    public :: fitpack_closed_curve_c_allocate
    public :: fitpack_closed_curve_c_destroy
    public :: fitpack_closed_curve_c_copy
    public :: fitpack_closed_curve_c_associate
    public :: fitpack_closed_curve_c_move_alloc
    public :: fitpack_closed_curve_c_set_default_parameters
    public :: fitpack_closed_curve_c_fit
    public :: fitpack_closed_curve_c_interpolate
    public :: fitpack_closed_curve_c_least_squares
    public :: fitpack_closed_curve_c_comm_size
    public :: fitpack_closed_curve_c_mse
    public :: fitpack_closed_curve_c_core_comm_size
    public :: fitpack_closed_curve_c_destroy_base
    public :: fitpack_closed_curve_c_new_points
    public :: fitpack_closed_curve_c_new_fit
    public :: fitpack_closed_curve_c_eval_one
    public :: fitpack_closed_curve_c_eval_many
    public :: fitpack_closed_curve_c_curve_derivative
    public :: fitpack_closed_curve_c_curve_derivatives
    public :: fitpack_closed_curve_c_curve_all_derivatives
    public :: fitpack_closed_curve_c_comm_pack
    public :: fitpack_closed_curve_c_comm_expand
    public :: fitpack_closed_curve_c_core_comm_pack
    public :: fitpack_closed_curve_c_core_comm_expand

    ! --- fitpack_constrained_curve ---
    public :: fitpack_constrained_curve_c_is_same
    public :: fitpack_constrained_curve_c_allocate
    public :: fitpack_constrained_curve_c_destroy
    public :: fitpack_constrained_curve_c_copy
    public :: fitpack_constrained_curve_c_associate
    public :: fitpack_constrained_curve_c_move_alloc
    public :: fitpack_constrained_curve_c_clean_constraints
    public :: fitpack_constrained_curve_c_comm_size
    public :: fitpack_constrained_curve_c_set_default_parameters
    public :: fitpack_constrained_curve_c_fit
    public :: fitpack_constrained_curve_c_interpolate
    public :: fitpack_constrained_curve_c_least_squares
    public :: fitpack_constrained_curve_c_mse
    public :: fitpack_constrained_curve_c_core_comm_size
    public :: fitpack_constrained_curve_c_destroy_base
    public :: fitpack_constrained_curve_c_set_constraints
    public :: fitpack_constrained_curve_c_comm_pack
    public :: fitpack_constrained_curve_c_comm_expand
    public :: fitpack_constrained_curve_c_new_points
    public :: fitpack_constrained_curve_c_new_fit
    public :: fitpack_constrained_curve_c_eval_one
    public :: fitpack_constrained_curve_c_eval_many
    public :: fitpack_constrained_curve_c_curve_derivative
    public :: fitpack_constrained_curve_c_curve_derivatives
    public :: fitpack_constrained_curve_c_curve_all_derivatives
    public :: fitpack_constrained_curve_c_core_comm_pack
    public :: fitpack_constrained_curve_c_core_comm_expand
    public :: fitpack_constrained_curve_c_getcomp_deriv_begin
    public :: fitpack_constrained_curve_c_getcomp_deriv_end
    public :: fitpack_constrained_curve_c_getcomp_xx
    public :: fitpack_constrained_curve_c_getcomp_cp
    public :: fitpack_constrained_curve_c_ref_ib
    public :: fitpack_constrained_curve_c_ref_ie

    interface f_pointer
        module procedure fitpack_parametric_curve_c_get_pointer
        module procedure fitpack_closed_curve_c_get_pointer
        module procedure fitpack_constrained_curve_c_get_pointer
    end interface f_pointer

    !> Pointer-identity check. Returns true iff `a` and `b` view the same
    !> underlying Fortran object (their internal c_ptrs match and are non-null).
    !> Useful for C-API consumers that want object identity without exposing
    !> raw pointers — e.g. GUI panels checking whether an input view points
    !> at the same Fortran allocation as the project that owns it.
    interface f_associated
        module procedure fitpack_parametric_curve_c_is_same
        module procedure fitpack_closed_curve_c_is_same
        module procedure fitpack_constrained_curve_c_is_same
    end interface f_associated

contains

    ! ===========================================================================================
    ! FITPACK_PARAMETRIC_CURVE: Core Memory Management
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_parametric_curve_c_is_same(this, that) &
            bind(C, name='fitpack_parametric_curve_c_is_same') result(same)
        type(fitpack_parametric_curve_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_parametric_curve_c_is_same

    function fitpack_parametric_curve_c_get_pointer(this) result(fptr)
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_parametric_curve_c_get_pointer

    function fitpack_parametric_curve_c_pointer(this, fthis) result(success)
        type(fitpack_parametric_curve_c), intent(inout) :: this
        type(fitpack_parametric_curve), pointer, intent(inout) :: fthis
        logical :: success
        integer :: ierr

        success = .true.
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) return

        allocate(fthis, stat=ierr)
        if (ierr /= 0) then
            success = .false.
            return
        end if

        this%cptr = c_loc(fthis)
        this%is_pointer = .false._c_bool
        this%name_cptr = c_loc(fitpack_parametric_curve_c_typename)
    end function fitpack_parametric_curve_c_pointer

    subroutine fitpack_parametric_curve_c_allocate(this, status) &
            bind(C, name='fitpack_parametric_curve_c_allocate')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_parametric_curve), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_parametric_curve_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_parametric_curve_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_parametric_curve_c_allocate

    subroutine fitpack_parametric_curve_c_destroy(this, status) &
            bind(C, name='fitpack_parametric_curve_c_destroy')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_parametric_curve), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_parametric_curve_c_destroy: deallocate failed')
        end if
        this = fitpack_parametric_curve_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_parametric_curve_c_destroy

    subroutine fitpack_parametric_curve_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_parametric_curve_c_copy')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        type(fitpack_parametric_curve_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_parametric_curve), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => fitpack_parametric_curve_c_get_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_parametric_curve_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_parametric_curve_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_parametric_curve_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_parametric_curve_c_copy

    subroutine fitpack_parametric_curve_c_associate(this, that, status) &
            bind(C, name='fitpack_parametric_curve_c_associate')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        type(fitpack_parametric_curve_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_parametric_curve_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_parametric_curve_c_associate

    subroutine fitpack_parametric_curve_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_parametric_curve_c_move_alloc')
        type(fitpack_parametric_curve_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_parametric_curve_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_parametric_curve_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_parametric_curve_c_move_alloc

    ! ===========================================================================================
    ! FITPACK_PARAMETRIC_CURVE: Method Wrappers
    ! ===========================================================================================
    subroutine fitpack_parametric_curve_c_set_default_parameters(this) bind(C,  &
        & name='fitpack_parametric_curve_c_set_default_parameters')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%set_default_parameters()
        end if
    end subroutine fitpack_parametric_curve_c_set_default_parameters

    integer(c_int32_t) function fitpack_parametric_curve_c_fit(this, smoothing, order, keep_knots) bind(C,  &
        & name='fitpack_parametric_curve_c_fit')
        type(fitpack_parametric_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_parametric_curve), pointer :: fthis
        logical :: f_keep_knots

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_parametric_curve_c_fit = fthis%fit(smoothing=smoothing, order=order, keep_knots=f_keep_knots)
            else
                fitpack_parametric_curve_c_fit = fthis%fit(smoothing=smoothing, order=order)
            end if
        else
            fitpack_parametric_curve_c_fit = 0
        end if
    end function fitpack_parametric_curve_c_fit

    integer(c_int32_t) function fitpack_parametric_curve_c_interpolate(this, order, reset_knots) bind(C,  &
        & name='fitpack_parametric_curve_c_interpolate')
        type(fitpack_parametric_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_parametric_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_parametric_curve_c_interpolate = fthis%interpolate(order=order, reset_knots=f_reset_knots)
            else
                fitpack_parametric_curve_c_interpolate = fthis%interpolate(order=order)
            end if
        else
            fitpack_parametric_curve_c_interpolate = 0
        end if
    end function fitpack_parametric_curve_c_interpolate

    integer(c_int32_t) function fitpack_parametric_curve_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_parametric_curve_c_least_squares')
        type(fitpack_parametric_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_parametric_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_parametric_curve_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_parametric_curve_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_parametric_curve_c_least_squares = 0
        end if
    end function fitpack_parametric_curve_c_least_squares

    integer(c_int32_t) function fitpack_parametric_curve_c_comm_size(this) bind(C, name='fitpack_parametric_curve_c_comm_size')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_comm_size = fthis%comm_size()
        else
            fitpack_parametric_curve_c_comm_size = 0
        end if
    end function fitpack_parametric_curve_c_comm_size

    real(c_double) function fitpack_parametric_curve_c_mse(this) bind(C, name='fitpack_parametric_curve_c_mse')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_mse = fthis%mse()
        else
            fitpack_parametric_curve_c_mse = 0.0_c_double
        end if
    end function fitpack_parametric_curve_c_mse

    integer(c_int32_t) function fitpack_parametric_curve_c_core_comm_size(this) bind(C,  &
        & name='fitpack_parametric_curve_c_core_comm_size')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_parametric_curve_c_core_comm_size = 0
        end if
    end function fitpack_parametric_curve_c_core_comm_size

    subroutine fitpack_parametric_curve_c_destroy_base(this) bind(C, name='fitpack_parametric_curve_c_destroy_base')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_parametric_curve_c_destroy_base

    subroutine fitpack_parametric_curve_c_new_points(this, x_n1, x_n2, x, u, w) bind(C,  &
        & name='fitpack_parametric_curve_c_new_points')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: x_n1, x_n2
        real(c_double), intent(in) :: x(x_n1, x_n2)
        real(c_double), optional, intent(in) :: u(x_n2)
        real(c_double), optional, intent(in) :: w(x_n2)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(x, u, w)
        end if
    end subroutine fitpack_parametric_curve_c_new_points

    integer(c_int32_t) function fitpack_parametric_curve_c_new_fit(this, x_n1, x_n2, x, u, w, smoothing, order) bind(C,  &
        & name='fitpack_parametric_curve_c_new_fit')
        type(fitpack_parametric_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: x_n1, x_n2
        real(c_double), intent(in) :: x(x_n1, x_n2)
        real(c_double), optional, intent(in) :: u(x_n2)
        real(c_double), optional, intent(in) :: w(x_n2)
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_new_fit = fthis%new_fit(x, u, w, smoothing, order)
        else
            fitpack_parametric_curve_c_new_fit = 0
        end if
    end function fitpack_parametric_curve_c_new_fit

    subroutine fitpack_parametric_curve_c_eval_one(this, u, ierr, result, n_result) bind(C,  &
        & name='fitpack_parametric_curve_c_eval_one')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval_one(u, ierr=ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_parametric_curve_c_eval_one

    subroutine fitpack_parametric_curve_c_eval_many(this, n, u, ierr, result, n_result) bind(C,  &
        & name='fitpack_parametric_curve_c_eval_many')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: u(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval_many(u, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_parametric_curve_c_eval_many

    subroutine fitpack_parametric_curve_c_curve_derivative(this, u, order, ierr, result, n_result) bind(C,  &
        & name='fitpack_parametric_curve_c_curve_derivative')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(in), value :: order
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(u, order, ierr=ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_parametric_curve_c_curve_derivative

    subroutine fitpack_parametric_curve_c_curve_derivatives(this, n, u, order, ierr, result, n_result) bind(C,  &
        & name='fitpack_parametric_curve_c_curve_derivatives')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, order
        real(c_double), intent(in) :: u(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(u, order, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_parametric_curve_c_curve_derivatives

    subroutine fitpack_parametric_curve_c_curve_all_derivatives(this, u, ierr, result, n_result) bind(C,  &
        & name='fitpack_parametric_curve_c_curve_all_derivatives')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_all(u, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_parametric_curve_c_curve_all_derivatives

    subroutine fitpack_parametric_curve_c_comm_pack(this, n, buffer) bind(C, name='fitpack_parametric_curve_c_comm_pack')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_parametric_curve_c_comm_pack

    subroutine fitpack_parametric_curve_c_comm_expand(this, n, buffer) bind(C, name='fitpack_parametric_curve_c_comm_expand')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_parametric_curve_c_comm_expand

    subroutine fitpack_parametric_curve_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_parametric_curve_c_core_comm_pack')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_parametric_curve_c_core_comm_pack

    subroutine fitpack_parametric_curve_c_core_comm_expand(this, n, buffer) bind(C,  &
        & name='fitpack_parametric_curve_c_core_comm_expand')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_parametric_curve), pointer :: fthis

        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_parametric_curve_c_core_comm_expand

    ! ===========================================================================================
    ! Component Array Accessors for fitpack_parametric_curve (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'x': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_curve_c_getcomp_x(this, ptr, extents) &
        bind(C, name='fitpack_parametric_curve_c_getcomp_x')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_parametric_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%x)) then
                extents = int(shape(fthis%x), kind=c_int64_t)
                if (size(fthis%x) > 0) ptr = c_loc(fthis%x)
            end if
        end if
    end subroutine fitpack_parametric_curve_c_getcomp_x

    !> Raw view of component 'u': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_curve_c_getcomp_u(this, ptr, extents) &
        bind(C, name='fitpack_parametric_curve_c_getcomp_u')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_parametric_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%u)) then
                extents = int(shape(fthis%u), kind=c_int64_t)
                if (size(fthis%u) > 0) ptr = c_loc(fthis%u)
            end if
        end if
    end subroutine fitpack_parametric_curve_c_getcomp_u

    !> Raw view of component 'sp': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_curve_c_getcomp_sp(this, ptr, extents) &
        bind(C, name='fitpack_parametric_curve_c_getcomp_sp')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_parametric_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%sp)) then
                extents = int(shape(fthis%sp), kind=c_int64_t)
                if (size(fthis%sp) > 0) ptr = c_loc(fthis%sp)
            end if
        end if
    end subroutine fitpack_parametric_curve_c_getcomp_sp

    !> Raw view of component 'w': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_curve_c_getcomp_w(this, ptr, extents) &
        bind(C, name='fitpack_parametric_curve_c_getcomp_w')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_parametric_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%w)) then
                extents = int(shape(fthis%w), kind=c_int64_t)
                if (size(fthis%w) > 0) ptr = c_loc(fthis%w)
            end if
        end if
    end subroutine fitpack_parametric_curve_c_getcomp_w

    !> Raw view of component 'dd': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_curve_c_getcomp_dd(this, ptr, extents) &
        bind(C, name='fitpack_parametric_curve_c_getcomp_dd')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_parametric_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%dd)) then
                extents = int(shape(fthis%dd), kind=c_int64_t)
                if (size(fthis%dd) > 0) ptr = c_loc(fthis%dd)
            end if
        end if
    end subroutine fitpack_parametric_curve_c_getcomp_dd

    !> Raw view of component 't': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_curve_c_getcomp_t(this, ptr, extents) &
        bind(C, name='fitpack_parametric_curve_c_getcomp_t')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_parametric_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%t)) then
                extents = int(shape(fthis%t), kind=c_int64_t)
                if (size(fthis%t) > 0) ptr = c_loc(fthis%t)
            end if
        end if
    end subroutine fitpack_parametric_curve_c_getcomp_t

    ! ===========================================================================================
    ! Scalar Property Accessors for fitpack_parametric_curve
    ! ===========================================================================================

    !> Get pointer to scalar property 'm'
    type(c_ptr) function fitpack_parametric_curve_c_ref_m(this) &
        bind(C, name='fitpack_parametric_curve_c_ref_m')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_ref_m = c_loc(fthis%m)
        else
            fitpack_parametric_curve_c_ref_m = c_null_ptr
        end if
    end function fitpack_parametric_curve_c_ref_m

    !> Get pointer to scalar property 'idim'
    type(c_ptr) function fitpack_parametric_curve_c_ref_idim(this) &
        bind(C, name='fitpack_parametric_curve_c_ref_idim')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_ref_idim = c_loc(fthis%idim)
        else
            fitpack_parametric_curve_c_ref_idim = c_null_ptr
        end if
    end function fitpack_parametric_curve_c_ref_idim

    !> Get logical property 'has_params'
    logical(c_bool) function fitpack_parametric_curve_c_get_has_params(this) &
        bind(C, name='fitpack_parametric_curve_c_get_has_params')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_get_has_params = logical(fthis%has_params, c_bool)
        else
            fitpack_parametric_curve_c_get_has_params = .false._c_bool
        end if
    end function fitpack_parametric_curve_c_get_has_params

    !> Set logical property 'has_params'
    subroutine fitpack_parametric_curve_c_set_has_params(this, value) &
        bind(C, name='fitpack_parametric_curve_c_set_has_params')
        type(fitpack_parametric_curve_c), intent(inout) :: this
        logical(c_bool), intent(in), value :: value
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) fthis%has_params = logical(value)
    end subroutine fitpack_parametric_curve_c_set_has_params

    !> Get pointer to scalar property 'order'
    type(c_ptr) function fitpack_parametric_curve_c_ref_order(this) &
        bind(C, name='fitpack_parametric_curve_c_ref_order')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_ref_order = c_loc(fthis%order)
        else
            fitpack_parametric_curve_c_ref_order = c_null_ptr
        end if
    end function fitpack_parametric_curve_c_ref_order

    !> Get pointer to scalar property 'ubegin'
    type(c_ptr) function fitpack_parametric_curve_c_ref_ubegin(this) &
        bind(C, name='fitpack_parametric_curve_c_ref_ubegin')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_ref_ubegin = c_loc(fthis%ubegin)
        else
            fitpack_parametric_curve_c_ref_ubegin = c_null_ptr
        end if
    end function fitpack_parametric_curve_c_ref_ubegin

    !> Get pointer to scalar property 'uend'
    type(c_ptr) function fitpack_parametric_curve_c_ref_uend(this) &
        bind(C, name='fitpack_parametric_curve_c_ref_uend')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_ref_uend = c_loc(fthis%uend)
        else
            fitpack_parametric_curve_c_ref_uend = c_null_ptr
        end if
    end function fitpack_parametric_curve_c_ref_uend

    !> Get pointer to scalar property 'nest'
    type(c_ptr) function fitpack_parametric_curve_c_ref_nest(this) &
        bind(C, name='fitpack_parametric_curve_c_ref_nest')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_ref_nest = c_loc(fthis%nest)
        else
            fitpack_parametric_curve_c_ref_nest = c_null_ptr
        end if
    end function fitpack_parametric_curve_c_ref_nest

    !> Get pointer to scalar property 'knots'
    type(c_ptr) function fitpack_parametric_curve_c_ref_knots(this) &
        bind(C, name='fitpack_parametric_curve_c_ref_knots')
        type(fitpack_parametric_curve_c), intent(in) :: this
        type(fitpack_parametric_curve), pointer :: fthis
        fthis => fitpack_parametric_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_curve_c_ref_knots = c_loc(fthis%knots)
        else
            fitpack_parametric_curve_c_ref_knots = c_null_ptr
        end if
    end function fitpack_parametric_curve_c_ref_knots

    ! ===========================================================================================
    ! FITPACK_CLOSED_CURVE: Core Memory Management
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_closed_curve_c_is_same(this, that) &
            bind(C, name='fitpack_closed_curve_c_is_same') result(same)
        type(fitpack_closed_curve_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_closed_curve_c_is_same

    function fitpack_closed_curve_c_get_pointer(this) result(fptr)
        type(fitpack_closed_curve_c), intent(in) :: this
        type(fitpack_closed_curve), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_closed_curve_c_get_pointer

    function fitpack_closed_curve_c_pointer(this, fthis) result(success)
        type(fitpack_closed_curve_c), intent(inout) :: this
        type(fitpack_closed_curve), pointer, intent(inout) :: fthis
        logical :: success
        integer :: ierr

        success = .true.
        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) return

        allocate(fthis, stat=ierr)
        if (ierr /= 0) then
            success = .false.
            return
        end if

        this%cptr = c_loc(fthis)
        this%is_pointer = .false._c_bool
        this%name_cptr = c_loc(fitpack_closed_curve_c_typename)
    end function fitpack_closed_curve_c_pointer

    subroutine fitpack_closed_curve_c_allocate(this, status) &
            bind(C, name='fitpack_closed_curve_c_allocate')
        type(fitpack_closed_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_closed_curve), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_closed_curve_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_closed_curve_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_closed_curve_c_allocate

    subroutine fitpack_closed_curve_c_destroy(this, status) &
            bind(C, name='fitpack_closed_curve_c_destroy')
        type(fitpack_closed_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_closed_curve), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_closed_curve_c_destroy: deallocate failed')
        end if
        this = fitpack_closed_curve_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_closed_curve_c_destroy

    subroutine fitpack_closed_curve_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_closed_curve_c_copy')
        type(fitpack_closed_curve_c), intent(inout) :: this
        type(fitpack_closed_curve_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_closed_curve), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => fitpack_closed_curve_c_get_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_closed_curve_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_closed_curve_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_closed_curve_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_closed_curve_c_copy

    subroutine fitpack_closed_curve_c_associate(this, that, status) &
            bind(C, name='fitpack_closed_curve_c_associate')
        type(fitpack_closed_curve_c), intent(inout) :: this
        type(fitpack_closed_curve_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_closed_curve_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_closed_curve_c_associate

    subroutine fitpack_closed_curve_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_closed_curve_c_move_alloc')
        type(fitpack_closed_curve_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_closed_curve_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_closed_curve_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_closed_curve_c_move_alloc

    ! ===========================================================================================
    ! FITPACK_CLOSED_CURVE: Method Wrappers
    ! ===========================================================================================
    subroutine fitpack_closed_curve_c_set_default_parameters(this) bind(C, name='fitpack_closed_curve_c_set_default_parameters')
        type(fitpack_closed_curve_c), intent(inout) :: this
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%set_default_parameters()
        end if
    end subroutine fitpack_closed_curve_c_set_default_parameters

    integer(c_int32_t) function fitpack_closed_curve_c_fit(this, smoothing, order, keep_knots) bind(C,  &
        & name='fitpack_closed_curve_c_fit')
        type(fitpack_closed_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_closed_curve), pointer :: fthis
        logical :: f_keep_knots

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_closed_curve_c_fit = fthis%fit(smoothing=smoothing, order=order, keep_knots=f_keep_knots)
            else
                fitpack_closed_curve_c_fit = fthis%fit(smoothing=smoothing, order=order)
            end if
        else
            fitpack_closed_curve_c_fit = 0
        end if
    end function fitpack_closed_curve_c_fit

    integer(c_int32_t) function fitpack_closed_curve_c_interpolate(this, order, reset_knots) bind(C,  &
        & name='fitpack_closed_curve_c_interpolate')
        type(fitpack_closed_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_closed_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_closed_curve_c_interpolate = fthis%interpolate(order=order, reset_knots=f_reset_knots)
            else
                fitpack_closed_curve_c_interpolate = fthis%interpolate(order=order)
            end if
        else
            fitpack_closed_curve_c_interpolate = 0
        end if
    end function fitpack_closed_curve_c_interpolate

    integer(c_int32_t) function fitpack_closed_curve_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_closed_curve_c_least_squares')
        type(fitpack_closed_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_closed_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_closed_curve_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_closed_curve_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_closed_curve_c_least_squares = 0
        end if
    end function fitpack_closed_curve_c_least_squares

    integer(c_int32_t) function fitpack_closed_curve_c_comm_size(this) bind(C, name='fitpack_closed_curve_c_comm_size')
        type(fitpack_closed_curve_c), intent(in) :: this
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_closed_curve_c_comm_size = fthis%comm_size()
        else
            fitpack_closed_curve_c_comm_size = 0
        end if
    end function fitpack_closed_curve_c_comm_size

    real(c_double) function fitpack_closed_curve_c_mse(this) bind(C, name='fitpack_closed_curve_c_mse')
        type(fitpack_closed_curve_c), intent(in) :: this
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_closed_curve_c_mse = fthis%mse()
        else
            fitpack_closed_curve_c_mse = 0.0_c_double
        end if
    end function fitpack_closed_curve_c_mse

    integer(c_int32_t) function fitpack_closed_curve_c_core_comm_size(this) bind(C, name='fitpack_closed_curve_c_core_comm_size')
        type(fitpack_closed_curve_c), intent(in) :: this
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_closed_curve_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_closed_curve_c_core_comm_size = 0
        end if
    end function fitpack_closed_curve_c_core_comm_size

    subroutine fitpack_closed_curve_c_destroy_base(this) bind(C, name='fitpack_closed_curve_c_destroy_base')
        type(fitpack_closed_curve_c), intent(inout) :: this
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_closed_curve_c_destroy_base

    subroutine fitpack_closed_curve_c_new_points(this, x_n1, x_n2, x, u, w) bind(C, name='fitpack_closed_curve_c_new_points')
        type(fitpack_closed_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: x_n1, x_n2
        real(c_double), intent(in) :: x(x_n1, x_n2)
        real(c_double), optional, intent(in) :: u(x_n2)
        real(c_double), optional, intent(in) :: w(x_n2)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(x, u, w)
        end if
    end subroutine fitpack_closed_curve_c_new_points

    integer(c_int32_t) function fitpack_closed_curve_c_new_fit(this, x_n1, x_n2, x, u, w, smoothing, order) bind(C,  &
        & name='fitpack_closed_curve_c_new_fit')
        type(fitpack_closed_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: x_n1, x_n2
        real(c_double), intent(in) :: x(x_n1, x_n2)
        real(c_double), optional, intent(in) :: u(x_n2)
        real(c_double), optional, intent(in) :: w(x_n2)
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_closed_curve_c_new_fit = fthis%new_fit(x, u, w, smoothing, order)
        else
            fitpack_closed_curve_c_new_fit = 0
        end if
    end function fitpack_closed_curve_c_new_fit

    subroutine fitpack_closed_curve_c_eval_one(this, u, ierr, result, n_result) bind(C, name='fitpack_closed_curve_c_eval_one')
        type(fitpack_closed_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval_one(u, ierr=ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_closed_curve_c_eval_one

    subroutine fitpack_closed_curve_c_eval_many(this, n, u, ierr, result, n_result) bind(C, name='fitpack_closed_curve_c_eval_many')
        type(fitpack_closed_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: u(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval_many(u, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_closed_curve_c_eval_many

    subroutine fitpack_closed_curve_c_curve_derivative(this, u, order, ierr, result, n_result) bind(C,  &
        & name='fitpack_closed_curve_c_curve_derivative')
        type(fitpack_closed_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(in), value :: order
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(u, order, ierr=ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_closed_curve_c_curve_derivative

    subroutine fitpack_closed_curve_c_curve_derivatives(this, n, u, order, ierr, result, n_result) bind(C,  &
        & name='fitpack_closed_curve_c_curve_derivatives')
        type(fitpack_closed_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, order
        real(c_double), intent(in) :: u(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(u, order, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_closed_curve_c_curve_derivatives

    subroutine fitpack_closed_curve_c_curve_all_derivatives(this, u, ierr, result, n_result) bind(C,  &
        & name='fitpack_closed_curve_c_curve_all_derivatives')
        type(fitpack_closed_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_all(u, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_closed_curve_c_curve_all_derivatives

    subroutine fitpack_closed_curve_c_comm_pack(this, n, buffer) bind(C, name='fitpack_closed_curve_c_comm_pack')
        type(fitpack_closed_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_closed_curve_c_comm_pack

    subroutine fitpack_closed_curve_c_comm_expand(this, n, buffer) bind(C, name='fitpack_closed_curve_c_comm_expand')
        type(fitpack_closed_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_closed_curve_c_comm_expand

    subroutine fitpack_closed_curve_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_closed_curve_c_core_comm_pack')
        type(fitpack_closed_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_closed_curve_c_core_comm_pack

    subroutine fitpack_closed_curve_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_closed_curve_c_core_comm_expand')
        type(fitpack_closed_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_closed_curve), pointer :: fthis

        fthis => fitpack_closed_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_closed_curve_c_core_comm_expand

    ! ===========================================================================================
    ! FITPACK_CONSTRAINED_CURVE: Core Memory Management
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_constrained_curve_c_is_same(this, that) &
            bind(C, name='fitpack_constrained_curve_c_is_same') result(same)
        type(fitpack_constrained_curve_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_constrained_curve_c_is_same

    function fitpack_constrained_curve_c_get_pointer(this) result(fptr)
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(fitpack_constrained_curve), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_constrained_curve_c_get_pointer

    function fitpack_constrained_curve_c_pointer(this, fthis) result(success)
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve), pointer, intent(inout) :: fthis
        logical :: success
        integer :: ierr

        success = .true.
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) return

        allocate(fthis, stat=ierr)
        if (ierr /= 0) then
            success = .false.
            return
        end if

        this%cptr = c_loc(fthis)
        this%is_pointer = .false._c_bool
        this%name_cptr = c_loc(fitpack_constrained_curve_c_typename)
    end function fitpack_constrained_curve_c_pointer

    subroutine fitpack_constrained_curve_c_allocate(this, status) &
            bind(C, name='fitpack_constrained_curve_c_allocate')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_constrained_curve), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_constrained_curve_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_constrained_curve_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_constrained_curve_c_allocate

    subroutine fitpack_constrained_curve_c_destroy(this, status) &
            bind(C, name='fitpack_constrained_curve_c_destroy')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_constrained_curve), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_constrained_curve_c_destroy: deallocate failed')
        end if
        this = fitpack_constrained_curve_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_constrained_curve_c_destroy

    subroutine fitpack_constrained_curve_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_constrained_curve_c_copy')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_constrained_curve), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => fitpack_constrained_curve_c_get_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_constrained_curve_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_constrained_curve_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_constrained_curve_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_constrained_curve_c_copy

    subroutine fitpack_constrained_curve_c_associate(this, that, status) &
            bind(C, name='fitpack_constrained_curve_c_associate')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_constrained_curve_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_constrained_curve_c_associate

    subroutine fitpack_constrained_curve_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_constrained_curve_c_move_alloc')
        type(fitpack_constrained_curve_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_constrained_curve_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_constrained_curve_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_constrained_curve_c_move_alloc

    ! ===========================================================================================
    ! FITPACK_CONSTRAINED_CURVE: Method Wrappers
    ! ===========================================================================================
    subroutine fitpack_constrained_curve_c_clean_constraints(this) bind(C, name='fitpack_constrained_curve_c_clean_constraints')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%clean_constraints()
        end if
    end subroutine fitpack_constrained_curve_c_clean_constraints

    integer(c_int32_t) function fitpack_constrained_curve_c_comm_size(this) bind(C, name='fitpack_constrained_curve_c_comm_size')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_constrained_curve_c_comm_size = fthis%comm_size()
        else
            fitpack_constrained_curve_c_comm_size = 0
        end if
    end function fitpack_constrained_curve_c_comm_size

    subroutine fitpack_constrained_curve_c_set_default_parameters(this) bind(C,  &
        & name='fitpack_constrained_curve_c_set_default_parameters')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%set_default_parameters()
        end if
    end subroutine fitpack_constrained_curve_c_set_default_parameters

    integer(c_int32_t) function fitpack_constrained_curve_c_fit(this, smoothing, order, keep_knots) bind(C,  &
        & name='fitpack_constrained_curve_c_fit')
        type(fitpack_constrained_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_constrained_curve), pointer :: fthis
        logical :: f_keep_knots

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_constrained_curve_c_fit = fthis%fit(smoothing=smoothing, order=order, keep_knots=f_keep_knots)
            else
                fitpack_constrained_curve_c_fit = fthis%fit(smoothing=smoothing, order=order)
            end if
        else
            fitpack_constrained_curve_c_fit = 0
        end if
    end function fitpack_constrained_curve_c_fit

    integer(c_int32_t) function fitpack_constrained_curve_c_interpolate(this, order, reset_knots) bind(C,  &
        & name='fitpack_constrained_curve_c_interpolate')
        type(fitpack_constrained_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_constrained_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_constrained_curve_c_interpolate = fthis%interpolate(order=order, reset_knots=f_reset_knots)
            else
                fitpack_constrained_curve_c_interpolate = fthis%interpolate(order=order)
            end if
        else
            fitpack_constrained_curve_c_interpolate = 0
        end if
    end function fitpack_constrained_curve_c_interpolate

    integer(c_int32_t) function fitpack_constrained_curve_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_constrained_curve_c_least_squares')
        type(fitpack_constrained_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_constrained_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_constrained_curve_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_constrained_curve_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_constrained_curve_c_least_squares = 0
        end if
    end function fitpack_constrained_curve_c_least_squares

    real(c_double) function fitpack_constrained_curve_c_mse(this) bind(C, name='fitpack_constrained_curve_c_mse')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_constrained_curve_c_mse = fthis%mse()
        else
            fitpack_constrained_curve_c_mse = 0.0_c_double
        end if
    end function fitpack_constrained_curve_c_mse

    integer(c_int32_t) function fitpack_constrained_curve_c_core_comm_size(this) bind(C,  &
        & name='fitpack_constrained_curve_c_core_comm_size')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_constrained_curve_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_constrained_curve_c_core_comm_size = 0
        end if
    end function fitpack_constrained_curve_c_core_comm_size

    subroutine fitpack_constrained_curve_c_destroy_base(this) bind(C, name='fitpack_constrained_curve_c_destroy_base')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_constrained_curve_c_destroy_base

    subroutine fitpack_constrained_curve_c_set_constraints(this, ddx_begin_n1, ddx_begin_n2, ddx_begin, ddx_end_n1, ddx_end_n2,  &
        & ddx_end, ierr) bind(C, name='fitpack_constrained_curve_c_set_constraints')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: ddx_begin_n1, ddx_begin_n2, ddx_end_n1, ddx_end_n2
        real(c_double), optional, intent(in) :: ddx_begin(ddx_begin_n1, ddx_begin_n2)
        real(c_double), optional, intent(in) :: ddx_end(ddx_end_n1, ddx_end_n2)
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%set_constraints(ddx_begin, ddx_end, ierr)
        end if
    end subroutine fitpack_constrained_curve_c_set_constraints

    subroutine fitpack_constrained_curve_c_comm_pack(this, n, buffer) bind(C, name='fitpack_constrained_curve_c_comm_pack')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_constrained_curve_c_comm_pack

    subroutine fitpack_constrained_curve_c_comm_expand(this, n, buffer) bind(C, name='fitpack_constrained_curve_c_comm_expand')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_constrained_curve_c_comm_expand

    subroutine fitpack_constrained_curve_c_new_points(this, x_n1, x_n2, x, u, w) bind(C,  &
        & name='fitpack_constrained_curve_c_new_points')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: x_n1, x_n2
        real(c_double), intent(in) :: x(x_n1, x_n2)
        real(c_double), optional, intent(in) :: u(x_n2)
        real(c_double), optional, intent(in) :: w(x_n2)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(x, u, w)
        end if
    end subroutine fitpack_constrained_curve_c_new_points

    integer(c_int32_t) function fitpack_constrained_curve_c_new_fit(this, x_n1, x_n2, x, u, w, smoothing, order) bind(C,  &
        & name='fitpack_constrained_curve_c_new_fit')
        type(fitpack_constrained_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: x_n1, x_n2
        real(c_double), intent(in) :: x(x_n1, x_n2)
        real(c_double), optional, intent(in) :: u(x_n2)
        real(c_double), optional, intent(in) :: w(x_n2)
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_constrained_curve_c_new_fit = fthis%new_fit(x, u, w, smoothing, order)
        else
            fitpack_constrained_curve_c_new_fit = 0
        end if
    end function fitpack_constrained_curve_c_new_fit

    subroutine fitpack_constrained_curve_c_eval_one(this, u, ierr, result, n_result) bind(C,  &
        & name='fitpack_constrained_curve_c_eval_one')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval_one(u, ierr=ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_constrained_curve_c_eval_one

    subroutine fitpack_constrained_curve_c_eval_many(this, n, u, ierr, result, n_result) bind(C,  &
        & name='fitpack_constrained_curve_c_eval_many')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: u(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval_many(u, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_constrained_curve_c_eval_many

    subroutine fitpack_constrained_curve_c_curve_derivative(this, u, order, ierr, result, n_result) bind(C,  &
        & name='fitpack_constrained_curve_c_curve_derivative')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(in), value :: order
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(u, order, ierr=ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_constrained_curve_c_curve_derivative

    subroutine fitpack_constrained_curve_c_curve_derivatives(this, n, u, order, ierr, result, n_result) bind(C,  &
        & name='fitpack_constrained_curve_c_curve_derivatives')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, order
        real(c_double), intent(in) :: u(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx(u, order, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_constrained_curve_c_curve_derivatives

    subroutine fitpack_constrained_curve_c_curve_all_derivatives(this, u, ierr, result, n_result) bind(C,  &
        & name='fitpack_constrained_curve_c_curve_all_derivatives')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        real(c_double), intent(in), value :: u
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%dfdx_all(u, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_constrained_curve_c_curve_all_derivatives

    subroutine fitpack_constrained_curve_c_core_comm_pack(this, n, buffer) bind(C,  &
        & name='fitpack_constrained_curve_c_core_comm_pack')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_constrained_curve_c_core_comm_pack

    subroutine fitpack_constrained_curve_c_core_comm_expand(this, n, buffer) bind(C,  &
        & name='fitpack_constrained_curve_c_core_comm_expand')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_constrained_curve), pointer :: fthis

        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_constrained_curve_c_core_comm_expand

    ! ===========================================================================================
    ! Component Array Accessors for fitpack_constrained_curve (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'deriv_begin': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_constrained_curve_c_getcomp_deriv_begin(this, ptr, extents) &
        bind(C, name='fitpack_constrained_curve_c_getcomp_deriv_begin')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_constrained_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%deriv_begin)) then
                extents = int(shape(fthis%deriv_begin), kind=c_int64_t)
                if (size(fthis%deriv_begin) > 0) ptr = c_loc(fthis%deriv_begin)
            end if
        end if
    end subroutine fitpack_constrained_curve_c_getcomp_deriv_begin

    !> Raw view of component 'deriv_end': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_constrained_curve_c_getcomp_deriv_end(this, ptr, extents) &
        bind(C, name='fitpack_constrained_curve_c_getcomp_deriv_end')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_constrained_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%deriv_end)) then
                extents = int(shape(fthis%deriv_end), kind=c_int64_t)
                if (size(fthis%deriv_end) > 0) ptr = c_loc(fthis%deriv_end)
            end if
        end if
    end subroutine fitpack_constrained_curve_c_getcomp_deriv_end

    !> Raw view of component 'xx': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_constrained_curve_c_getcomp_xx(this, ptr, extents) &
        bind(C, name='fitpack_constrained_curve_c_getcomp_xx')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_constrained_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%xx)) then
                extents = int(shape(fthis%xx), kind=c_int64_t)
                if (size(fthis%xx) > 0) ptr = c_loc(fthis%xx)
            end if
        end if
    end subroutine fitpack_constrained_curve_c_getcomp_xx

    !> Raw view of component 'cp': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_constrained_curve_c_getcomp_cp(this, ptr, extents) &
        bind(C, name='fitpack_constrained_curve_c_getcomp_cp')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_constrained_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%cp)) then
                extents = int(shape(fthis%cp), kind=c_int64_t)
                if (size(fthis%cp) > 0) ptr = c_loc(fthis%cp)
            end if
        end if
    end subroutine fitpack_constrained_curve_c_getcomp_cp

    ! ===========================================================================================
    ! Scalar Property Accessors for fitpack_constrained_curve
    ! ===========================================================================================

    !> Get pointer to scalar property 'ib'
    type(c_ptr) function fitpack_constrained_curve_c_ref_ib(this) &
        bind(C, name='fitpack_constrained_curve_c_ref_ib')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(fitpack_constrained_curve), pointer :: fthis
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_constrained_curve_c_ref_ib = c_loc(fthis%ib)
        else
            fitpack_constrained_curve_c_ref_ib = c_null_ptr
        end if
    end function fitpack_constrained_curve_c_ref_ib

    !> Get pointer to scalar property 'ie'
    type(c_ptr) function fitpack_constrained_curve_c_ref_ie(this) &
        bind(C, name='fitpack_constrained_curve_c_ref_ie')
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(fitpack_constrained_curve), pointer :: fthis
        fthis => fitpack_constrained_curve_c_get_pointer(this)
        if (associated(fthis)) then
            fitpack_constrained_curve_c_ref_ie = c_loc(fthis%ie)
        else
            fitpack_constrained_curve_c_ref_ie = c_null_ptr
        end if
    end function fitpack_constrained_curve_c_ref_ie

end module fitpack_parametric_curves_c
