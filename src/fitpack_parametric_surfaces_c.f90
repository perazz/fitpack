!   ***********************************************************************************************
!   **                                        fxArray                                          **
!   **                                  Fortran Arrays for C++                                 **
!   ***********************************************************************************************
!   **    fitpack_parametric_surfaces_c                                                                       **
!> @brief C interface to fitpack_parametric_surface
!   ***********************************************************************************************
!> @author Binding Generator
!   ***********************************************************************************************

module fitpack_parametric_surfaces_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_parametric_surfaces, only: fitpack_parametric_surface
    use fitpack_parametric_surfaces_c_types, only: fitpack_parametric_surface_c, fitpack_parametric_surface_c_null,  &
        & fitpack_parametric_surface_c_typename
    use fitpack_fx_status, only: fx_status, fx_status_ok, FX_ERROR_ALLOCATION, &
                        FX_ERROR_DEALLOCATION, handle_error
    use, intrinsic :: iso_c_binding
    implicit none(type, external)
    private

    ! ===========================================================================================
    ! SECTION 2: Public Exports
    ! ===========================================================================================

    ! Core interfaces (ALWAYS export these)
    public :: fitpack_parametric_surface_c            ! Type and constructor
    public :: fitpack_parametric_surface_c_null       ! Null constant
    public :: f_pointer              ! Get Fortran pointer (non-allocating)
    public :: f_associated           ! Pointer-identity check (non bind(C))
    public :: fitpack_parametric_surface_c_is_same   ! Pointer-identity check [bind(C)]
    public :: fitpack_parametric_surface_c_allocate   ! Allocate new object [bind(C)]
    public :: fitpack_parametric_surface_c_destroy    ! Deallocate object [bind(C)]
    public :: fitpack_parametric_surface_c_copy       ! Deep copy [bind(C)]
    public :: fitpack_parametric_surface_c_associate  ! Shallow copy/pointer [bind(C)]
    public :: fitpack_parametric_surface_c_move_alloc ! Move ownership [bind(C)]

    ! Method wrappers
    public :: fitpack_parametric_surface_c_interpolate
    public :: fitpack_parametric_surface_c_comm_size
    public :: fitpack_parametric_surface_c_mse
    public :: fitpack_parametric_surface_c_core_comm_size
    public :: fitpack_parametric_surface_c_destroy_base
    public :: fitpack_parametric_surface_c_fit
    public :: fitpack_parametric_surface_c_least_squares
    public :: fitpack_parametric_surface_c_surf_eval_one
    public :: fitpack_parametric_surface_c_surf_eval_grid
    public :: fitpack_parametric_surface_c_comm_pack
    public :: fitpack_parametric_surface_c_comm_expand
    public :: fitpack_parametric_surface_c_core_comm_pack
    public :: fitpack_parametric_surface_c_core_comm_expand
    public :: fitpack_parametric_surface_c_getcomp_u
    public :: fitpack_parametric_surface_c_getcomp_v
    public :: fitpack_parametric_surface_c_getcomp_z
    public :: fitpack_parametric_surface_c_getcomp_t
    public :: fitpack_parametric_surface_c_ref_idim
    public :: fitpack_parametric_surface_c_ref_nmax

    ! ===========================================================================================
    ! SECTION 4: Fortran-Side Interfaces (NOT bind(C))
    ! ===========================================================================================

    !> Get pointer to embedded Fortran object
    interface f_pointer
        module procedure fitpack_parametric_surface_c_get_pointer
    end interface f_pointer

    !> Pointer-identity check (Fortran-side wrapper of the bind(C) function).
    interface f_associated
        module procedure fitpack_parametric_surface_c_is_same
    end interface f_associated

    !> Construct from Fortran object
    interface fitpack_parametric_surface_c
        module procedure fitpack_parametric_surface_c_f2c
        module procedure fitpack_parametric_surface_c_fptr2c
    end interface fitpack_parametric_surface_c

contains

    ! ===========================================================================================
    ! SECTION 5: Core Memory Management (ALWAYS IMPLEMENT)
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_parametric_surface_c_is_same(this, that) &
            bind(C, name='fitpack_parametric_surface_c_is_same') result(same)
        type(fitpack_parametric_surface_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_parametric_surface_c_is_same

    !> Get pointer to embedded Fortran object (non-allocating)
    function fitpack_parametric_surface_c_get_pointer(this) result(fptr)
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(fitpack_parametric_surface), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_parametric_surface_c_get_pointer

    !> Get/allocate embedded Fortran object (internal, not bind(C))
    !> Returns .true. on success, .false. on allocation failure
    function fitpack_parametric_surface_c_pointer(this, fthis) result(success)
        type(fitpack_parametric_surface_c), intent(inout) :: this
        type(fitpack_parametric_surface), pointer, intent(inout) :: fthis
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
        this%name_cptr = c_loc(fitpack_parametric_surface_c_typename)
    end function fitpack_parametric_surface_c_pointer

    !> [bind(C)] Allocate new object
    !> @param status Optional. If absent, allocation failure triggers error stop.
    subroutine fitpack_parametric_surface_c_allocate(this, status) bind(C, name='fitpack_parametric_surface_c_allocate')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_parametric_surface), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_parametric_surface_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_parametric_surface_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_parametric_surface_c_allocate

    !> [bind(C)] Deallocate object
    !> @param status Optional. If absent, deallocation failure triggers error stop.
    subroutine fitpack_parametric_surface_c_destroy(this, status) bind(C, name='fitpack_parametric_surface_c_destroy')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_parametric_surface), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => f_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_parametric_surface_c_destroy: deallocate failed')
        end if
        this = fitpack_parametric_surface_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_parametric_surface_c_destroy

    !> [bind(C)] Copy. By default, mirrors source ownership: a view source
    !> yields a view (shallow handle copy); an owned source yields a deep
    !> copy via Fortran intrinsic assignment. Pass deep_copy=.true. to force
    !> a deep copy regardless of the source's ownership.
    !> @param deep_copy When .true., always allocate a fresh Fortran object
    !>                  and deep-copy data, even if the source is a view.
    !> @param status    Optional. If absent, errors trigger error stop.
    subroutine fitpack_parametric_surface_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_parametric_surface_c_copy')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        type(fitpack_parametric_surface_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_parametric_surface), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => f_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_parametric_surface_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_parametric_surface_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_parametric_surface_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_parametric_surface_c_copy

    !> [bind(C)] Shallow copy (pointer semantics — Fortran "associate" construct)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_parametric_surface_c_associate(this, that, status) &
            bind(C, name='fitpack_parametric_surface_c_associate')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        type(fitpack_parametric_surface_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_parametric_surface_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_parametric_surface_c_associate

    !> [bind(C)] Move allocation (transfer ownership)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_parametric_surface_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_parametric_surface_c_move_alloc')
        type(fitpack_parametric_surface_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_parametric_surface_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_parametric_surface_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_parametric_surface_c_move_alloc

    ! ===========================================================================================
    ! SECTION 6: Fortran Constructors (NOT bind(C))
    ! ===========================================================================================

    !> Create from Fortran object (owns copy)
    type(fitpack_parametric_surface_c) function fitpack_parametric_surface_c_f2c(f) result(c)
        type(fitpack_parametric_surface), intent(in) :: f
        type(fitpack_parametric_surface), pointer :: fptr
        logical :: ok

        ok = fitpack_parametric_surface_c_pointer(c, fptr)
        if (ok) fptr = f
    end function fitpack_parametric_surface_c_f2c

    !> Create from Fortran object (optionally as pointer)
    type(fitpack_parametric_surface_c) function fitpack_parametric_surface_c_fptr2c(f, want_pointer) result(c)
        type(fitpack_parametric_surface), intent(inout), target :: f
        logical, intent(in) :: want_pointer

        if (want_pointer) then
            c%cptr = c_loc(f)
            c%is_pointer = .true._c_bool
        else
            c = fitpack_parametric_surface_c_f2c(f)
        end if
    end function fitpack_parametric_surface_c_fptr2c

    ! ===========================================================================================
    ! SECTION 7: Method Wrappers (bind(C))
    ! ===========================================================================================
    !> interpolate
    integer(c_int32_t) function fitpack_parametric_surface_c_interpolate(this, reset_knots) bind(C,  &
        & name='fitpack_parametric_surface_c_interpolate')
        type(fitpack_parametric_surface_c), intent(in) :: this
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_parametric_surface), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_parametric_surface_c_interpolate = fthis%interpolate(reset_knots=f_reset_knots)
            else
                fitpack_parametric_surface_c_interpolate = fthis%interpolate()
            end if
        else
            fitpack_parametric_surface_c_interpolate = 0
        end if
    end function fitpack_parametric_surface_c_interpolate

    !> comm_size
    integer(c_int32_t) function fitpack_parametric_surface_c_comm_size(this) bind(C, name='fitpack_parametric_surface_c_comm_size')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_surface_c_comm_size = fthis%comm_size()
        else
            fitpack_parametric_surface_c_comm_size = 0
        end if
    end function fitpack_parametric_surface_c_comm_size

    !> mse
    real(c_double) function fitpack_parametric_surface_c_mse(this) bind(C, name='fitpack_parametric_surface_c_mse')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_surface_c_mse = fthis%mse()
        else
            fitpack_parametric_surface_c_mse = 0.0_c_double
        end if
    end function fitpack_parametric_surface_c_mse

    !> core_comm_size
    integer(c_int32_t) function fitpack_parametric_surface_c_core_comm_size(this) bind(C,  &
        & name='fitpack_parametric_surface_c_core_comm_size')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_surface_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_parametric_surface_c_core_comm_size = 0
        end if
    end function fitpack_parametric_surface_c_core_comm_size

    !> destroy_base
    subroutine fitpack_parametric_surface_c_destroy_base(this) bind(C, name='fitpack_parametric_surface_c_destroy_base')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_parametric_surface_c_destroy_base

    !> fit
    integer(c_int32_t) function fitpack_parametric_surface_c_fit(this, smoothing, n, periodic, keep_knots) bind(C,  &
        & name='fitpack_parametric_surface_c_fit')
        type(fitpack_parametric_surface_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), optional, intent(in) :: periodic(n)
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_parametric_surface), pointer :: fthis
        logical, allocatable :: fperiodic(:)
        logical :: f_keep_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fperiodic = logical(periodic)
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_parametric_surface_c_fit = fthis%fit(smoothing=smoothing, periodic=fperiodic, keep_knots=f_keep_knots)
            else
                fitpack_parametric_surface_c_fit = fthis%fit(smoothing=smoothing, periodic=fperiodic)
            end if
        else
            fitpack_parametric_surface_c_fit = 0
        end if
    end function fitpack_parametric_surface_c_fit

    !> least_squares
    integer(c_int32_t) function fitpack_parametric_surface_c_least_squares(this, n, u_knots, v_knots, smoothing, reset_knots)  &
        & bind(C, name='fitpack_parametric_surface_c_least_squares')
        type(fitpack_parametric_surface_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), optional, intent(in) :: u_knots(n)
        real(c_double), optional, intent(in) :: v_knots(n)
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_parametric_surface), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_parametric_surface_c_least_squares = fthis%least_squares(u_knots=u_knots, v_knots=v_knots,  &
                    & smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_parametric_surface_c_least_squares = fthis%least_squares(u_knots=u_knots, v_knots=v_knots,  &
                    & smoothing=smoothing)
            end if
        else
            fitpack_parametric_surface_c_least_squares = 0
        end if
    end function fitpack_parametric_surface_c_least_squares

    !> surf_eval_one
    subroutine fitpack_parametric_surface_c_surf_eval_one(this, u, v, ierr, result, n_result) bind(C,  &
        & name='fitpack_parametric_surface_c_surf_eval_one')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        real(c_double), intent(in), value :: u, v
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double) :: fresult(n_result)
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(u=u, v=v, ierr=ierr)
            result(1:n_result) = fresult
        end if
    end subroutine fitpack_parametric_surface_c_surf_eval_one

    !> surf_eval_grid
    subroutine fitpack_parametric_surface_c_surf_eval_grid(this, n, u, v, ierr, result, n_result) bind(C,  &
        & name='fitpack_parametric_surface_c_surf_eval_grid')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: u(n)
        real(c_double), intent(in) :: v(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(in), value :: n_result
        real(c_double), allocatable :: fresult(:,:,:)
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(u=u, v=v, ierr=ierr)
            result(1:n_result) = reshape(fresult, [n_result])
        end if
    end subroutine fitpack_parametric_surface_c_surf_eval_grid

    !> comm_pack
    subroutine fitpack_parametric_surface_c_comm_pack(this, n, buffer) bind(C, name='fitpack_parametric_surface_c_comm_pack')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_parametric_surface_c_comm_pack

    !> comm_expand
    subroutine fitpack_parametric_surface_c_comm_expand(this, n, buffer) bind(C, name='fitpack_parametric_surface_c_comm_expand')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_parametric_surface_c_comm_expand

    !> core_comm_pack
    subroutine fitpack_parametric_surface_c_core_comm_pack(this, n, buffer) bind(C,  &
        & name='fitpack_parametric_surface_c_core_comm_pack')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_parametric_surface_c_core_comm_pack

    !> core_comm_expand
    subroutine fitpack_parametric_surface_c_core_comm_expand(this, n, buffer) bind(C,  &
        & name='fitpack_parametric_surface_c_core_comm_expand')
        type(fitpack_parametric_surface_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_parametric_surface), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_parametric_surface_c_core_comm_expand

    ! ===========================================================================================
    ! SECTION 8: Component Array Accessors (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'u': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_surface_c_getcomp_u(this, ptr, extents) &
        bind(C, name='fitpack_parametric_surface_c_getcomp_u')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_parametric_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%u)) then
                extents = int(shape(fthis%u), kind=c_int64_t)
                if (size(fthis%u) > 0) ptr = c_loc(fthis%u)
            end if
        end if
    end subroutine fitpack_parametric_surface_c_getcomp_u

    !> Raw view of component 'v': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_surface_c_getcomp_v(this, ptr, extents) &
        bind(C, name='fitpack_parametric_surface_c_getcomp_v')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_parametric_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%v)) then
                extents = int(shape(fthis%v), kind=c_int64_t)
                if (size(fthis%v) > 0) ptr = c_loc(fthis%v)
            end if
        end if
    end subroutine fitpack_parametric_surface_c_getcomp_v

    !> Raw view of component 'z': address of the first element, plus
    !> the component's 3 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_surface_c_getcomp_z(this, ptr, extents) &
        bind(C, name='fitpack_parametric_surface_c_getcomp_z')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(3)
        type(fitpack_parametric_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%z)) then
                extents = int(shape(fthis%z), kind=c_int64_t)
                if (size(fthis%z) > 0) ptr = c_loc(fthis%z)
            end if
        end if
    end subroutine fitpack_parametric_surface_c_getcomp_z

    !> Raw view of component 't': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_parametric_surface_c_getcomp_t(this, ptr, extents) &
        bind(C, name='fitpack_parametric_surface_c_getcomp_t')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_parametric_surface), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%t)) then
                extents = int(shape(fthis%t), kind=c_int64_t)
                if (size(fthis%t) > 0) ptr = c_loc(fthis%t)
            end if
        end if
    end subroutine fitpack_parametric_surface_c_getcomp_t

    ! ===========================================================================================
    ! SECTION 9: Scalar Property Accessors
    ! ===========================================================================================

    !> Get pointer to scalar property 'idim'
    type(c_ptr) function fitpack_parametric_surface_c_ref_idim(this) &
        bind(C, name='fitpack_parametric_surface_c_ref_idim')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(fitpack_parametric_surface), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_surface_c_ref_idim = c_loc(fthis%idim)
        else
            fitpack_parametric_surface_c_ref_idim = c_null_ptr
        end if
    end function fitpack_parametric_surface_c_ref_idim

    !> Get pointer to scalar property 'nmax'
    type(c_ptr) function fitpack_parametric_surface_c_ref_nmax(this) &
        bind(C, name='fitpack_parametric_surface_c_ref_nmax')
        type(fitpack_parametric_surface_c), intent(in) :: this
        type(fitpack_parametric_surface), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_parametric_surface_c_ref_nmax = c_loc(fthis%nmax)
        else
            fitpack_parametric_surface_c_ref_nmax = c_null_ptr
        end if
    end function fitpack_parametric_surface_c_ref_nmax

end module fitpack_parametric_surfaces_c
