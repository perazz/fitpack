!   ***********************************************************************************************
!   **                                        fxArray                                          **
!   **                                  Fortran Arrays for C++                                 **
!   ***********************************************************************************************
!   **    fitpack_convex_curves_c                                                                       **
!> @brief C interface to fitpack_convex_curve
!   ***********************************************************************************************
!> @author Binding Generator
!   ***********************************************************************************************

module fitpack_convex_curves_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_convex_curves, only: fitpack_convex_curve
    use fitpack_convex_curves_c_types, only: fitpack_convex_curve_c, fitpack_convex_curve_c_null, fitpack_convex_curve_c_typename
    use fitpack_fx_status, only: fx_status, fx_status_ok, FX_ERROR_ALLOCATION, &
                        FX_ERROR_DEALLOCATION, handle_error
    use, intrinsic :: iso_c_binding
    implicit none(type, external)
    private

    ! ===========================================================================================
    ! SECTION 2: Public Exports
    ! ===========================================================================================

    ! Core interfaces (ALWAYS export these)
    public :: fitpack_convex_curve_c            ! Type and constructor
    public :: fitpack_convex_curve_c_null       ! Null constant
    public :: f_pointer              ! Get Fortran pointer (non-allocating)
    public :: f_associated           ! Pointer-identity check (non bind(C))
    public :: fitpack_convex_curve_c_is_same   ! Pointer-identity check [bind(C)]
    public :: fitpack_convex_curve_c_allocate   ! Allocate new object [bind(C)]
    public :: fitpack_convex_curve_c_destroy    ! Deallocate object [bind(C)]
    public :: fitpack_convex_curve_c_copy       ! Deep copy [bind(C)]
    public :: fitpack_convex_curve_c_associate  ! Shallow copy/pointer [bind(C)]
    public :: fitpack_convex_curve_c_move_alloc ! Move ownership [bind(C)]

    ! Method wrappers
    public :: fitpack_convex_curve_c_fit
    public :: fitpack_convex_curve_c_least_squares
    public :: fitpack_convex_curve_c_comm_size
    public :: fitpack_convex_curve_c_mse
    public :: fitpack_convex_curve_c_core_comm_size
    public :: fitpack_convex_curve_c_destroy_base
    public :: fitpack_convex_curve_c_new_points
    public :: fitpack_convex_curve_c_set_convexity
    public :: fitpack_convex_curve_c_comm_pack
    public :: fitpack_convex_curve_c_comm_expand
    public :: fitpack_convex_curve_c_core_comm_pack
    public :: fitpack_convex_curve_c_core_comm_expand
    public :: fitpack_convex_curve_c_getcomp_v
    public :: fitpack_convex_curve_c_getcomp_sx
    public :: fitpack_convex_curve_c_ref_maxtr
    public :: fitpack_convex_curve_c_ref_maxbin

    ! ===========================================================================================
    ! SECTION 4: Fortran-Side Interfaces (NOT bind(C))
    ! ===========================================================================================

    !> Get pointer to embedded Fortran object
    interface f_pointer
        module procedure fitpack_convex_curve_c_get_pointer
    end interface f_pointer

    !> Pointer-identity check (Fortran-side wrapper of the bind(C) function).
    interface f_associated
        module procedure fitpack_convex_curve_c_is_same
    end interface f_associated

    !> Construct from Fortran object
    interface fitpack_convex_curve_c
        module procedure fitpack_convex_curve_c_f2c
        module procedure fitpack_convex_curve_c_fptr2c
    end interface fitpack_convex_curve_c

contains

    ! ===========================================================================================
    ! SECTION 5: Core Memory Management (ALWAYS IMPLEMENT)
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_convex_curve_c_is_same(this, that) &
            bind(C, name='fitpack_convex_curve_c_is_same') result(same)
        type(fitpack_convex_curve_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_convex_curve_c_is_same

    !> Get pointer to embedded Fortran object (non-allocating)
    function fitpack_convex_curve_c_get_pointer(this) result(fptr)
        type(fitpack_convex_curve_c), intent(in) :: this
        type(fitpack_convex_curve), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_convex_curve_c_get_pointer

    !> Get/allocate embedded Fortran object (internal, not bind(C))
    !> Returns .true. on success, .false. on allocation failure
    function fitpack_convex_curve_c_pointer(this, fthis) result(success)
        type(fitpack_convex_curve_c), intent(inout) :: this
        type(fitpack_convex_curve), pointer, intent(inout) :: fthis
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
        this%name_cptr = c_loc(fitpack_convex_curve_c_typename)
    end function fitpack_convex_curve_c_pointer

    !> [bind(C)] Allocate new object
    !> @param status Optional. If absent, allocation failure triggers error stop.
    subroutine fitpack_convex_curve_c_allocate(this, status) bind(C, name='fitpack_convex_curve_c_allocate')
        type(fitpack_convex_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_convex_curve), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_convex_curve_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_convex_curve_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_convex_curve_c_allocate

    !> [bind(C)] Deallocate object
    !> @param status Optional. If absent, deallocation failure triggers error stop.
    subroutine fitpack_convex_curve_c_destroy(this, status) bind(C, name='fitpack_convex_curve_c_destroy')
        type(fitpack_convex_curve_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_convex_curve), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => f_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_convex_curve_c_destroy: deallocate failed')
        end if
        this = fitpack_convex_curve_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_convex_curve_c_destroy

    !> [bind(C)] Copy. By default, mirrors source ownership: a view source
    !> yields a view (shallow handle copy); an owned source yields a deep
    !> copy via Fortran intrinsic assignment. Pass deep_copy=.true. to force
    !> a deep copy regardless of the source's ownership.
    !> @param deep_copy When .true., always allocate a fresh Fortran object
    !>                  and deep-copy data, even if the source is a view.
    !> @param status    Optional. If absent, errors trigger error stop.
    subroutine fitpack_convex_curve_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_convex_curve_c_copy')
        type(fitpack_convex_curve_c), intent(inout) :: this
        type(fitpack_convex_curve_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_convex_curve), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => f_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_convex_curve_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_convex_curve_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_convex_curve_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_convex_curve_c_copy

    !> [bind(C)] Shallow copy (pointer semantics — Fortran "associate" construct)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_convex_curve_c_associate(this, that, status) &
            bind(C, name='fitpack_convex_curve_c_associate')
        type(fitpack_convex_curve_c), intent(inout) :: this
        type(fitpack_convex_curve_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_convex_curve_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_convex_curve_c_associate

    !> [bind(C)] Move allocation (transfer ownership)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_convex_curve_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_convex_curve_c_move_alloc')
        type(fitpack_convex_curve_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_convex_curve_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_convex_curve_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_convex_curve_c_move_alloc

    ! ===========================================================================================
    ! SECTION 6: Fortran Constructors (NOT bind(C))
    ! ===========================================================================================

    !> Create from Fortran object (owns copy)
    type(fitpack_convex_curve_c) function fitpack_convex_curve_c_f2c(f) result(c)
        type(fitpack_convex_curve), intent(in) :: f
        type(fitpack_convex_curve), pointer :: fptr
        logical :: ok

        ok = fitpack_convex_curve_c_pointer(c, fptr)
        if (ok) fptr = f
    end function fitpack_convex_curve_c_f2c

    !> Create from Fortran object (optionally as pointer)
    type(fitpack_convex_curve_c) function fitpack_convex_curve_c_fptr2c(f, want_pointer) result(c)
        type(fitpack_convex_curve), intent(inout), target :: f
        logical, intent(in) :: want_pointer

        if (want_pointer) then
            c%cptr = c_loc(f)
            c%is_pointer = .true._c_bool
        else
            c = fitpack_convex_curve_c_f2c(f)
        end if
    end function fitpack_convex_curve_c_fptr2c

    ! ===========================================================================================
    ! SECTION 7: Method Wrappers (bind(C))
    ! ===========================================================================================
    !> fit
    integer(c_int32_t) function fitpack_convex_curve_c_fit(this, smoothing, order, keep_knots) bind(C,  &
        & name='fitpack_convex_curve_c_fit')
        type(fitpack_convex_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        integer(c_int32_t), intent(in), optional :: order
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_convex_curve), pointer :: fthis
        logical :: f_keep_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_convex_curve_c_fit = fthis%fit(smoothing=smoothing, order=order, keep_knots=f_keep_knots)
            else
                fitpack_convex_curve_c_fit = fthis%fit(smoothing=smoothing, order=order)
            end if
        else
            fitpack_convex_curve_c_fit = 0
        end if
    end function fitpack_convex_curve_c_fit

    !> least_squares
    integer(c_int32_t) function fitpack_convex_curve_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_convex_curve_c_least_squares')
        type(fitpack_convex_curve_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_convex_curve), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_convex_curve_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_convex_curve_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_convex_curve_c_least_squares = 0
        end if
    end function fitpack_convex_curve_c_least_squares

    !> comm_size
    integer(c_int32_t) function fitpack_convex_curve_c_comm_size(this) bind(C, name='fitpack_convex_curve_c_comm_size')
        type(fitpack_convex_curve_c), intent(in) :: this
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_convex_curve_c_comm_size = fthis%comm_size()
        else
            fitpack_convex_curve_c_comm_size = 0
        end if
    end function fitpack_convex_curve_c_comm_size

    !> mse
    real(c_double) function fitpack_convex_curve_c_mse(this) bind(C, name='fitpack_convex_curve_c_mse')
        type(fitpack_convex_curve_c), intent(in) :: this
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_convex_curve_c_mse = fthis%mse()
        else
            fitpack_convex_curve_c_mse = 0.0_c_double
        end if
    end function fitpack_convex_curve_c_mse

    !> core_comm_size
    integer(c_int32_t) function fitpack_convex_curve_c_core_comm_size(this) bind(C, name='fitpack_convex_curve_c_core_comm_size')
        type(fitpack_convex_curve_c), intent(in) :: this
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_convex_curve_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_convex_curve_c_core_comm_size = 0
        end if
    end function fitpack_convex_curve_c_core_comm_size

    !> destroy_base
    subroutine fitpack_convex_curve_c_destroy_base(this) bind(C, name='fitpack_convex_curve_c_destroy_base')
        type(fitpack_convex_curve_c), intent(inout) :: this
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_convex_curve_c_destroy_base

    !> new_points
    subroutine fitpack_convex_curve_c_new_points(this, n, x, y, w) bind(C, name='fitpack_convex_curve_c_new_points')
        type(fitpack_convex_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: y(n)
        real(c_double), optional, intent(in) :: w(n)
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(x, y, w)
        end if
    end subroutine fitpack_convex_curve_c_new_points

    !> set_convexity
    integer(c_int32_t) function fitpack_convex_curve_c_set_convexity(this, n, v) bind(C,  &
        & name='fitpack_convex_curve_c_set_convexity')
        type(fitpack_convex_curve_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: v(n)
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_convex_curve_c_set_convexity = fthis%set_convexity(v)
        else
            fitpack_convex_curve_c_set_convexity = 0
        end if
    end function fitpack_convex_curve_c_set_convexity

    !> comm_pack
    subroutine fitpack_convex_curve_c_comm_pack(this, n, buffer) bind(C, name='fitpack_convex_curve_c_comm_pack')
        type(fitpack_convex_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_convex_curve_c_comm_pack

    !> comm_expand
    subroutine fitpack_convex_curve_c_comm_expand(this, n, buffer) bind(C, name='fitpack_convex_curve_c_comm_expand')
        type(fitpack_convex_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_convex_curve_c_comm_expand

    !> core_comm_pack
    subroutine fitpack_convex_curve_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_convex_curve_c_core_comm_pack')
        type(fitpack_convex_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_convex_curve_c_core_comm_pack

    !> core_comm_expand
    subroutine fitpack_convex_curve_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_convex_curve_c_core_comm_expand')
        type(fitpack_convex_curve_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_convex_curve), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_convex_curve_c_core_comm_expand

    ! ===========================================================================================
    ! SECTION 8: Component Array Accessors (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'v': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_convex_curve_c_getcomp_v(this, ptr, extents) &
        bind(C, name='fitpack_convex_curve_c_getcomp_v')
        type(fitpack_convex_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_convex_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%v)) then
                extents = int(shape(fthis%v), kind=c_int64_t)
                if (size(fthis%v) > 0) ptr = c_loc(fthis%v)
            end if
        end if
    end subroutine fitpack_convex_curve_c_getcomp_v

    !> Raw view of component 'sx': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_convex_curve_c_getcomp_sx(this, ptr, extents) &
        bind(C, name='fitpack_convex_curve_c_getcomp_sx')
        type(fitpack_convex_curve_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_convex_curve), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%sx)) then
                extents = int(shape(fthis%sx), kind=c_int64_t)
                if (size(fthis%sx) > 0) ptr = c_loc(fthis%sx)
            end if
        end if
    end subroutine fitpack_convex_curve_c_getcomp_sx

    ! ===========================================================================================
    ! SECTION 9: Scalar Property Accessors
    ! ===========================================================================================

    !> Get pointer to scalar property 'maxtr'
    type(c_ptr) function fitpack_convex_curve_c_ref_maxtr(this) &
        bind(C, name='fitpack_convex_curve_c_ref_maxtr')
        type(fitpack_convex_curve_c), intent(in) :: this
        type(fitpack_convex_curve), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_convex_curve_c_ref_maxtr = c_loc(fthis%maxtr)
        else
            fitpack_convex_curve_c_ref_maxtr = c_null_ptr
        end if
    end function fitpack_convex_curve_c_ref_maxtr

    !> Get pointer to scalar property 'maxbin'
    type(c_ptr) function fitpack_convex_curve_c_ref_maxbin(this) &
        bind(C, name='fitpack_convex_curve_c_ref_maxbin')
        type(fitpack_convex_curve_c), intent(in) :: this
        type(fitpack_convex_curve), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_convex_curve_c_ref_maxbin = c_loc(fthis%maxbin)
        else
            fitpack_convex_curve_c_ref_maxbin = c_null_ptr
        end if
    end function fitpack_convex_curve_c_ref_maxbin

end module fitpack_convex_curves_c
