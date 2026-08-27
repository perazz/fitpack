! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_sphere_domains_c.f90 (module fitpack_sphere_domains_c)
!> @brief C interface to fitpack_sphere
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

module fitpack_sphere_domains_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_sphere_domains, only: fitpack_sphere
    use fitpack_sphere_domains_c_types, only: fitpack_sphere_c, fitpack_sphere_c_null, fitpack_sphere_c_typename
    use fitpack_fx_status, only: fx_status, fx_status_ok, FX_ERROR_ALLOCATION, &
                        FX_ERROR_DEALLOCATION, handle_error
    use, intrinsic :: iso_c_binding
    implicit none(type, external)
    private

    ! ===========================================================================================
    ! SECTION 2: Public Exports
    ! ===========================================================================================

    ! Core interfaces (ALWAYS export these)
    public :: fitpack_sphere_c            ! Type and constructor
    public :: fitpack_sphere_c_null       ! Null constant
    public :: f_pointer              ! Get Fortran pointer (non-allocating)
    public :: f_associated           ! Pointer-identity check (non bind(C))
    public :: fitpack_sphere_c_is_same   ! Pointer-identity check [bind(C)]
    public :: fitpack_sphere_c_allocate   ! Allocate new object [bind(C)]
    public :: fitpack_sphere_c_destroy    ! Deallocate object [bind(C)]
    public :: fitpack_sphere_c_copy       ! Deep copy [bind(C)]
    public :: fitpack_sphere_c_associate  ! Shallow copy/pointer [bind(C)]
    public :: fitpack_sphere_c_move_alloc ! Move ownership [bind(C)]

    ! Method wrappers
    public :: fitpack_sphere_c_fit
    public :: fitpack_sphere_c_least_squares
    public :: fitpack_sphere_c_interpolate
    public :: fitpack_sphere_c_sphere_eval_one
    public :: fitpack_sphere_c_comm_size
    public :: fitpack_sphere_c_mse
    public :: fitpack_sphere_c_core_comm_size
    public :: fitpack_sphere_c_destroy_base
    public :: fitpack_sphere_c_new_points
    public :: fitpack_sphere_c_new_fit
    public :: fitpack_sphere_c_sphere_eval_many
    public :: fitpack_sphere_c_comm_pack
    public :: fitpack_sphere_c_comm_expand
    public :: fitpack_sphere_c_core_comm_pack
    public :: fitpack_sphere_c_core_comm_expand
    public :: fitpack_sphere_c_getcomp_theta
    public :: fitpack_sphere_c_getcomp_phi
    public :: fitpack_sphere_c_getcomp_r
    public :: fitpack_sphere_c_getcomp_w
    public :: fitpack_sphere_c_getcomp_wrk2
    public :: fitpack_sphere_c_getcomp_t
    public :: fitpack_sphere_c_ref_m
    public :: fitpack_sphere_c_ref_lwrk2
    public :: fitpack_sphere_c_ref_nmax

    ! ===========================================================================================
    ! SECTION 4: Fortran-Side Interfaces (NOT bind(C))
    ! ===========================================================================================

    !> Get pointer to embedded Fortran object
    interface f_pointer
        module procedure fitpack_sphere_c_get_pointer
    end interface f_pointer

    !> Pointer-identity check (Fortran-side wrapper of the bind(C) function).
    interface f_associated
        module procedure fitpack_sphere_c_is_same
    end interface f_associated

    !> Construct from Fortran object
    interface fitpack_sphere_c
        module procedure fitpack_sphere_c_f2c
        module procedure fitpack_sphere_c_fptr2c
    end interface fitpack_sphere_c

contains

    ! ===========================================================================================
    ! SECTION 5: Core Memory Management (ALWAYS IMPLEMENT)
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_sphere_c_is_same(this, that) &
            bind(C, name='fitpack_sphere_c_is_same') result(same)
        type(fitpack_sphere_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_sphere_c_is_same

    !> Get pointer to embedded Fortran object (non-allocating)
    function fitpack_sphere_c_get_pointer(this) result(fptr)
        type(fitpack_sphere_c), intent(in) :: this
        type(fitpack_sphere), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_sphere_c_get_pointer

    !> Get/allocate embedded Fortran object (internal, not bind(C))
    !> Returns .true. on success, .false. on allocation failure
    function fitpack_sphere_c_pointer(this, fthis) result(success)
        type(fitpack_sphere_c), intent(inout) :: this
        type(fitpack_sphere), pointer, intent(inout) :: fthis
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
        this%name_cptr = c_loc(fitpack_sphere_c_typename)
    end function fitpack_sphere_c_pointer

    !> [bind(C)] Allocate new object
    !> @param status Optional. If absent, allocation failure triggers error stop.
    subroutine fitpack_sphere_c_allocate(this, status) bind(C, name='fitpack_sphere_c_allocate')
        type(fitpack_sphere_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_sphere), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_sphere_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_sphere_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_sphere_c_allocate

    !> [bind(C)] Deallocate object
    !> @param status Optional. If absent, deallocation failure triggers error stop.
    subroutine fitpack_sphere_c_destroy(this, status) bind(C, name='fitpack_sphere_c_destroy')
        type(fitpack_sphere_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_sphere), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => f_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_sphere_c_destroy: deallocate failed')
        end if
        this = fitpack_sphere_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_sphere_c_destroy

    !> [bind(C)] Copy. By default, mirrors source ownership: a view source
    !> yields a view (shallow handle copy); an owned source yields a deep
    !> copy via Fortran intrinsic assignment. Pass deep_copy=.true. to force
    !> a deep copy regardless of the source's ownership.
    !> @param deep_copy When .true., always allocate a fresh Fortran object
    !>                  and deep-copy data, even if the source is a view.
    !> @param status    Optional. If absent, errors trigger error stop.
    subroutine fitpack_sphere_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_sphere_c_copy')
        type(fitpack_sphere_c), intent(inout) :: this
        type(fitpack_sphere_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_sphere), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => f_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_sphere_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_sphere_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_sphere_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_sphere_c_copy

    !> [bind(C)] Shallow copy (pointer semantics — Fortran "associate" construct)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_sphere_c_associate(this, that, status) &
            bind(C, name='fitpack_sphere_c_associate')
        type(fitpack_sphere_c), intent(inout) :: this
        type(fitpack_sphere_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_sphere_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_sphere_c_associate

    !> [bind(C)] Move allocation (transfer ownership)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_sphere_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_sphere_c_move_alloc')
        type(fitpack_sphere_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_sphere_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_sphere_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_sphere_c_move_alloc

    ! ===========================================================================================
    ! SECTION 6: Fortran Constructors (NOT bind(C))
    ! ===========================================================================================

    !> Create from Fortran object (owns copy)
    type(fitpack_sphere_c) function fitpack_sphere_c_f2c(f) result(c)
        type(fitpack_sphere), intent(in) :: f
        type(fitpack_sphere), pointer :: fptr
        logical :: ok

        ok = fitpack_sphere_c_pointer(c, fptr)
        if (ok) fptr = f
    end function fitpack_sphere_c_f2c

    !> Create from Fortran object (optionally as pointer)
    type(fitpack_sphere_c) function fitpack_sphere_c_fptr2c(f, want_pointer) result(c)
        type(fitpack_sphere), intent(inout), target :: f
        logical, intent(in) :: want_pointer

        if (want_pointer) then
            c%cptr = c_loc(f)
            c%is_pointer = .true._c_bool
        else
            c = fitpack_sphere_c_f2c(f)
        end if
    end function fitpack_sphere_c_fptr2c

    ! ===========================================================================================
    ! SECTION 7: Method Wrappers (bind(C))
    ! ===========================================================================================
    !> fit
    integer(c_int32_t) function fitpack_sphere_c_fit(this, smoothing, keep_knots) bind(C, name='fitpack_sphere_c_fit')
        type(fitpack_sphere_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_sphere), pointer :: fthis
        logical :: f_keep_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_sphere_c_fit = fthis%fit(smoothing=smoothing, keep_knots=f_keep_knots)
            else
                fitpack_sphere_c_fit = fthis%fit(smoothing=smoothing)
            end if
        else
            fitpack_sphere_c_fit = 0
        end if
    end function fitpack_sphere_c_fit

    !> least_squares
    integer(c_int32_t) function fitpack_sphere_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_sphere_c_least_squares')
        type(fitpack_sphere_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_sphere), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_sphere_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_sphere_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_sphere_c_least_squares = 0
        end if
    end function fitpack_sphere_c_least_squares

    !> interpolate
    integer(c_int32_t) function fitpack_sphere_c_interpolate(this, reset_knots) bind(C, name='fitpack_sphere_c_interpolate')
        type(fitpack_sphere_c), intent(in) :: this
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_sphere), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_sphere_c_interpolate = fthis%interpolate(reset_knots=f_reset_knots)
            else
                fitpack_sphere_c_interpolate = fthis%interpolate()
            end if
        else
            fitpack_sphere_c_interpolate = 0
        end if
    end function fitpack_sphere_c_interpolate

    !> sphere_eval_one
    real(c_double) function fitpack_sphere_c_sphere_eval_one(this, theta, phi, ierr) bind(C,  &
        & name='fitpack_sphere_c_sphere_eval_one')
        type(fitpack_sphere_c), intent(in) :: this
        real(c_double), intent(in), value :: theta, phi
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_sphere_eval_one = fthis%eval(theta, phi, ierr)
        else
            fitpack_sphere_c_sphere_eval_one = 0.0_c_double
        end if
    end function fitpack_sphere_c_sphere_eval_one

    !> comm_size
    integer(c_int32_t) function fitpack_sphere_c_comm_size(this) bind(C, name='fitpack_sphere_c_comm_size')
        type(fitpack_sphere_c), intent(in) :: this
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_comm_size = fthis%comm_size()
        else
            fitpack_sphere_c_comm_size = 0
        end if
    end function fitpack_sphere_c_comm_size

    !> mse
    real(c_double) function fitpack_sphere_c_mse(this) bind(C, name='fitpack_sphere_c_mse')
        type(fitpack_sphere_c), intent(in) :: this
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_mse = fthis%mse()
        else
            fitpack_sphere_c_mse = 0.0_c_double
        end if
    end function fitpack_sphere_c_mse

    !> core_comm_size
    integer(c_int32_t) function fitpack_sphere_c_core_comm_size(this) bind(C, name='fitpack_sphere_c_core_comm_size')
        type(fitpack_sphere_c), intent(in) :: this
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_sphere_c_core_comm_size = 0
        end if
    end function fitpack_sphere_c_core_comm_size

    !> destroy_base
    subroutine fitpack_sphere_c_destroy_base(this) bind(C, name='fitpack_sphere_c_destroy_base')
        type(fitpack_sphere_c), intent(inout) :: this
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_sphere_c_destroy_base

    !> new_points
    subroutine fitpack_sphere_c_new_points(this, n, theta, phi, r, w) bind(C, name='fitpack_sphere_c_new_points')
        type(fitpack_sphere_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: theta(n)
        real(c_double), intent(in) :: phi(n)
        real(c_double), intent(in) :: r(n)
        real(c_double), optional, intent(in) :: w(n)
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(theta, phi, r, w)
        end if
    end subroutine fitpack_sphere_c_new_points

    !> new_fit
    integer(c_int32_t) function fitpack_sphere_c_new_fit(this, n, theta, phi, r, w, smoothing) bind(C,  &
        & name='fitpack_sphere_c_new_fit')
        type(fitpack_sphere_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: theta(n)
        real(c_double), intent(in) :: phi(n)
        real(c_double), intent(in) :: r(n)
        real(c_double), optional, intent(in) :: w(n)
        real(c_double), intent(in), optional :: smoothing
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_new_fit = fthis%new_fit(theta, phi, r, w, smoothing)
        else
            fitpack_sphere_c_new_fit = 0
        end if
    end function fitpack_sphere_c_new_fit

    !> sphere_eval_many
    subroutine fitpack_sphere_c_sphere_eval_many(this, n, theta, phi, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_sphere_c_sphere_eval_many')
        type(fitpack_sphere_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: theta(n)
        real(c_double), intent(in) :: phi(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(theta=theta, phi=phi, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = reshape(fresult, [min(n_result, max_size)])
        end if
    end subroutine fitpack_sphere_c_sphere_eval_many

    !> comm_pack
    subroutine fitpack_sphere_c_comm_pack(this, n, buffer) bind(C, name='fitpack_sphere_c_comm_pack')
        type(fitpack_sphere_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_sphere_c_comm_pack

    !> comm_expand
    subroutine fitpack_sphere_c_comm_expand(this, n, buffer) bind(C, name='fitpack_sphere_c_comm_expand')
        type(fitpack_sphere_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_sphere_c_comm_expand

    !> core_comm_pack
    subroutine fitpack_sphere_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_sphere_c_core_comm_pack')
        type(fitpack_sphere_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_sphere_c_core_comm_pack

    !> core_comm_expand
    subroutine fitpack_sphere_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_sphere_c_core_comm_expand')
        type(fitpack_sphere_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_sphere), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_sphere_c_core_comm_expand

    ! ===========================================================================================
    ! SECTION 8: Component Array Accessors (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'theta': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_sphere_c_getcomp_theta(this, ptr, extents) &
        bind(C, name='fitpack_sphere_c_getcomp_theta')
        type(fitpack_sphere_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_sphere), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%theta)) then
                extents = int(shape(fthis%theta), kind=c_int64_t)
                if (size(fthis%theta) > 0) ptr = c_loc(fthis%theta)
            end if
        end if
    end subroutine fitpack_sphere_c_getcomp_theta

    !> Raw view of component 'phi': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_sphere_c_getcomp_phi(this, ptr, extents) &
        bind(C, name='fitpack_sphere_c_getcomp_phi')
        type(fitpack_sphere_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_sphere), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%phi)) then
                extents = int(shape(fthis%phi), kind=c_int64_t)
                if (size(fthis%phi) > 0) ptr = c_loc(fthis%phi)
            end if
        end if
    end subroutine fitpack_sphere_c_getcomp_phi

    !> Raw view of component 'r': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_sphere_c_getcomp_r(this, ptr, extents) &
        bind(C, name='fitpack_sphere_c_getcomp_r')
        type(fitpack_sphere_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_sphere), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%r)) then
                extents = int(shape(fthis%r), kind=c_int64_t)
                if (size(fthis%r) > 0) ptr = c_loc(fthis%r)
            end if
        end if
    end subroutine fitpack_sphere_c_getcomp_r

    !> Raw view of component 'w': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_sphere_c_getcomp_w(this, ptr, extents) &
        bind(C, name='fitpack_sphere_c_getcomp_w')
        type(fitpack_sphere_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_sphere), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%w)) then
                extents = int(shape(fthis%w), kind=c_int64_t)
                if (size(fthis%w) > 0) ptr = c_loc(fthis%w)
            end if
        end if
    end subroutine fitpack_sphere_c_getcomp_w

    !> Raw view of component 'wrk2': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_sphere_c_getcomp_wrk2(this, ptr, extents) &
        bind(C, name='fitpack_sphere_c_getcomp_wrk2')
        type(fitpack_sphere_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_sphere), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%wrk2)) then
                extents = int(shape(fthis%wrk2), kind=c_int64_t)
                if (size(fthis%wrk2) > 0) ptr = c_loc(fthis%wrk2)
            end if
        end if
    end subroutine fitpack_sphere_c_getcomp_wrk2

    !> Raw view of component 't': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_sphere_c_getcomp_t(this, ptr, extents) &
        bind(C, name='fitpack_sphere_c_getcomp_t')
        type(fitpack_sphere_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_sphere), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%t)) then
                extents = int(shape(fthis%t), kind=c_int64_t)
                if (size(fthis%t) > 0) ptr = c_loc(fthis%t)
            end if
        end if
    end subroutine fitpack_sphere_c_getcomp_t

    ! ===========================================================================================
    ! SECTION 9: Scalar Property Accessors
    ! ===========================================================================================

    !> Get pointer to scalar property 'm'
    type(c_ptr) function fitpack_sphere_c_ref_m(this) &
        bind(C, name='fitpack_sphere_c_ref_m')
        type(fitpack_sphere_c), intent(in) :: this
        type(fitpack_sphere), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_ref_m = c_loc(fthis%m)
        else
            fitpack_sphere_c_ref_m = c_null_ptr
        end if
    end function fitpack_sphere_c_ref_m

    !> Get pointer to scalar property 'lwrk2'
    type(c_ptr) function fitpack_sphere_c_ref_lwrk2(this) &
        bind(C, name='fitpack_sphere_c_ref_lwrk2')
        type(fitpack_sphere_c), intent(in) :: this
        type(fitpack_sphere), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_ref_lwrk2 = c_loc(fthis%lwrk2)
        else
            fitpack_sphere_c_ref_lwrk2 = c_null_ptr
        end if
    end function fitpack_sphere_c_ref_lwrk2

    !> Get pointer to scalar property 'nmax'
    type(c_ptr) function fitpack_sphere_c_ref_nmax(this) &
        bind(C, name='fitpack_sphere_c_ref_nmax')
        type(fitpack_sphere_c), intent(in) :: this
        type(fitpack_sphere), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_sphere_c_ref_nmax = c_loc(fthis%nmax)
        else
            fitpack_sphere_c_ref_nmax = c_null_ptr
        end if
    end function fitpack_sphere_c_ref_nmax

end module fitpack_sphere_domains_c
