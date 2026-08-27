! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_gridded_polar_c.f90 (module fitpack_gridded_polar_c)
!> @brief C interface to fitpack_grid_polar
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

module fitpack_gridded_polar_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_gridded_polar, only: fitpack_grid_polar
    use fitpack_gridded_polar_c_types, only: fitpack_grid_polar_c, fitpack_grid_polar_c_null, fitpack_grid_polar_c_typename
    use fitpack_fx_status, only: fx_status, fx_status_ok, FX_ERROR_ALLOCATION, &
                        FX_ERROR_DEALLOCATION, handle_error
    use, intrinsic :: iso_c_binding
    implicit none(type, external)
    private

    ! ===========================================================================================
    ! SECTION 2: Public Exports
    ! ===========================================================================================

    ! Core interfaces (ALWAYS export these)
    public :: fitpack_grid_polar_c            ! Type and constructor
    public :: fitpack_grid_polar_c_null       ! Null constant
    public :: f_pointer              ! Get Fortran pointer (non-allocating)
    public :: f_associated           ! Pointer-identity check (non bind(C))
    public :: fitpack_grid_polar_c_is_same   ! Pointer-identity check [bind(C)]
    public :: fitpack_grid_polar_c_allocate   ! Allocate new object [bind(C)]
    public :: fitpack_grid_polar_c_destroy    ! Deallocate object [bind(C)]
    public :: fitpack_grid_polar_c_copy       ! Deep copy [bind(C)]
    public :: fitpack_grid_polar_c_associate  ! Shallow copy/pointer [bind(C)]
    public :: fitpack_grid_polar_c_move_alloc ! Move ownership [bind(C)]

    ! Method wrappers
    public :: fitpack_grid_polar_c_set_origin_BC
    public :: fitpack_grid_polar_c_fit
    public :: fitpack_grid_polar_c_least_squares
    public :: fitpack_grid_polar_c_interpolate
    public :: fitpack_grid_polar_c_gridded_eval_one
    public :: fitpack_grid_polar_c_comm_size
    public :: fitpack_grid_polar_c_write
    public :: fitpack_grid_polar_c_mse
    public :: fitpack_grid_polar_c_core_comm_size
    public :: fitpack_grid_polar_c_destroy_base
    public :: fitpack_grid_polar_c_new_points
    public :: fitpack_grid_polar_c_new_fit
    public :: fitpack_grid_polar_c_gridded_eval_many
    public :: fitpack_grid_polar_c_comm_pack
    public :: fitpack_grid_polar_c_comm_expand
    public :: fitpack_grid_polar_c_core_comm_pack
    public :: fitpack_grid_polar_c_core_comm_expand
    public :: fitpack_grid_polar_c_getcomp_u
    public :: fitpack_grid_polar_c_getcomp_v
    public :: fitpack_grid_polar_c_getcomp_z
    public :: fitpack_grid_polar_c_getcomp_t
    public :: fitpack_grid_polar_c_ref_r
    public :: fitpack_grid_polar_c_ref_z0
    public :: fitpack_grid_polar_c_get_z0_present
    public :: fitpack_grid_polar_c_set_z0_present
    public :: fitpack_grid_polar_c_get_z0_exact
    public :: fitpack_grid_polar_c_set_z0_exact
    public :: fitpack_grid_polar_c_get_z0_zero_gradient
    public :: fitpack_grid_polar_c_set_z0_zero_gradient
    public :: fitpack_grid_polar_c_ref_nmax
    public :: fitpack_grid_polar_c_ref_bc_continuity_origin
    public :: fitpack_grid_polar_c_ref_bc_boundary

    ! ===========================================================================================
    ! SECTION 4: Fortran-Side Interfaces (NOT bind(C))
    ! ===========================================================================================

    !> Get pointer to embedded Fortran object
    interface f_pointer
        module procedure fitpack_grid_polar_c_get_pointer
    end interface f_pointer

    !> Pointer-identity check (Fortran-side wrapper of the bind(C) function).
    interface f_associated
        module procedure fitpack_grid_polar_c_is_same
    end interface f_associated

    !> Construct from Fortran object
    interface fitpack_grid_polar_c
        module procedure fitpack_grid_polar_c_f2c
        module procedure fitpack_grid_polar_c_fptr2c
    end interface fitpack_grid_polar_c

    !> Pointer view of a null-terminated C string: a `character(:,c_char),
    !> pointer` result aliasing the C buffer — no copy, no allocation.
    !> Always associated: a null, zero-length or ABSENT input yields a
    !> 0-length module sentinel, safe to pass to required `character(*)`
    !> dummies and to `trim()`. Two overloads:
    !>   * `f_char(c_str)` with `c_str` a `type(c_ptr)`
    !>   * `f_char(c_str)` with `c_str` a `character(*)` buffer (optional)
    !> Where the callee dummy is `optional`, use `f_char_opt` instead: it
    !> forwards an absent argument as a disassociated pointer, which an
    !> optional `character(*)` dummy reads as "not present" (F2018 15.5.2.12).
    interface f_char
        module procedure from_c_string_ptr_cptr
        module procedure from_c_string_ptr_arr
    end interface f_char

contains

    ! ===========================================================================================
    ! SECTION 5: Core Memory Management (ALWAYS IMPLEMENT)
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_grid_polar_c_is_same(this, that) &
            bind(C, name='fitpack_grid_polar_c_is_same') result(same)
        type(fitpack_grid_polar_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_grid_polar_c_is_same

    !> Get pointer to embedded Fortran object (non-allocating)
    function fitpack_grid_polar_c_get_pointer(this) result(fptr)
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr, fptr)
        else
            nullify(fptr)
        end if
    end function fitpack_grid_polar_c_get_pointer

    !> Get/allocate embedded Fortran object (internal, not bind(C))
    !> Returns .true. on success, .false. on allocation failure
    function fitpack_grid_polar_c_pointer(this, fthis) result(success)
        type(fitpack_grid_polar_c), intent(inout) :: this
        type(fitpack_grid_polar), pointer, intent(inout) :: fthis
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
        this%name_cptr = c_loc(fitpack_grid_polar_c_typename)
    end function fitpack_grid_polar_c_pointer

    !> [bind(C)] Allocate new object
    !> @param status Optional. If absent, allocation failure triggers error stop.
    subroutine fitpack_grid_polar_c_allocate(this, status) bind(C, name='fitpack_grid_polar_c_allocate')
        type(fitpack_grid_polar_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_grid_polar), pointer :: fthis
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        ok = fitpack_grid_polar_c_pointer(this, fthis)
        if (.not. ok) stat0 = fx_status(FX_ERROR_ALLOCATION, &
            'fitpack_grid_polar_c_allocate: allocation failed')
        call handle_error(stat0, status)
    end subroutine fitpack_grid_polar_c_allocate

    !> [bind(C)] Deallocate object
    !> @param status Optional. If absent, deallocation failure triggers error stop.
    subroutine fitpack_grid_polar_c_destroy(this, status) bind(C, name='fitpack_grid_polar_c_destroy')
        type(fitpack_grid_polar_c), intent(inout) :: this
        type(fx_status), intent(out), optional :: status
        type(fitpack_grid_polar), pointer :: fthis
        type(fx_status) :: stat0
        integer :: ierr

        stat0 = fx_status_ok
        fthis => f_pointer(this)
        if (associated(fthis) .and. .not. this%is_pointer) then
            deallocate(fthis, stat=ierr)
            if (ierr /= 0) stat0 = fx_status(FX_ERROR_DEALLOCATION, &
                'fitpack_grid_polar_c_destroy: deallocate failed')
        end if
        this = fitpack_grid_polar_c_null
        call handle_error(stat0, status)
    end subroutine fitpack_grid_polar_c_destroy

    !> [bind(C)] Copy. By default, mirrors source ownership: a view source
    !> yields a view (shallow handle copy); an owned source yields a deep
    !> copy via Fortran intrinsic assignment. Pass deep_copy=.true. to force
    !> a deep copy regardless of the source's ownership.
    !> @param deep_copy When .true., always allocate a fresh Fortran object
    !>                  and deep-copy data, even if the source is a view.
    !> @param status    Optional. If absent, errors trigger error stop.
    subroutine fitpack_grid_polar_c_copy(this, that, deep_copy, status) &
            bind(C, name='fitpack_grid_polar_c_copy')
        type(fitpack_grid_polar_c), intent(inout) :: this
        type(fitpack_grid_polar_c), intent(in)    :: that
        logical(c_bool),    intent(in), value :: deep_copy
        type(fx_status), intent(out), optional :: status
        type(fitpack_grid_polar), pointer :: fthis, fthat
        type(fx_status) :: stat0
        logical :: ok

        stat0 = fx_status_ok
        fthat => f_pointer(that)
        if (that%is_pointer .and. .not. deep_copy) then
            this = that
        elseif (associated(fthat)) then
            ok = fitpack_grid_polar_c_pointer(this, fthis)
            if (.not. ok) then
                stat0 = fx_status(FX_ERROR_ALLOCATION, &
                    'fitpack_grid_polar_c_copy: allocation failed')
                goto 100
            end if
            fthis = fthat
            this%name_cptr = that%name_cptr
        else
            call fitpack_grid_polar_c_destroy(this, stat0)
        end if
100     call handle_error(stat0, status)
    end subroutine fitpack_grid_polar_c_copy

    !> [bind(C)] Shallow copy (pointer semantics — Fortran "associate" construct)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_grid_polar_c_associate(this, that, status) &
            bind(C, name='fitpack_grid_polar_c_associate')
        type(fitpack_grid_polar_c), intent(inout) :: this
        type(fitpack_grid_polar_c), intent(in)    :: that
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_grid_polar_c_destroy(this, stat0)
        if (.not. stat0%ok) goto 100
        this%cptr = that%cptr
        this%is_pointer = .true._c_bool
        this%name_cptr = that%name_cptr
100     call handle_error(stat0, status)
    end subroutine fitpack_grid_polar_c_associate

    !> [bind(C)] Move allocation (transfer ownership)
    !> @param status Optional. If absent, errors trigger error stop.
    subroutine fitpack_grid_polar_c_move_alloc(to, from, status) &
            bind(C, name='fitpack_grid_polar_c_move_alloc')
        type(fitpack_grid_polar_c), intent(inout) :: to, from
        type(fx_status), intent(out), optional :: status
        type(fx_status) :: stat0

        stat0 = fx_status_ok
        call fitpack_grid_polar_c_destroy(to, stat0)
        if (.not. stat0%ok) goto 100
        to = from
        from = fitpack_grid_polar_c_null
100     call handle_error(stat0, status)
    end subroutine fitpack_grid_polar_c_move_alloc

    ! ===========================================================================================
    ! SECTION 6: Fortran Constructors (NOT bind(C))
    ! ===========================================================================================

    !> Create from Fortran object (owns copy)
    type(fitpack_grid_polar_c) function fitpack_grid_polar_c_f2c(f) result(c)
        type(fitpack_grid_polar), intent(in) :: f
        type(fitpack_grid_polar), pointer :: fptr
        logical :: ok

        ok = fitpack_grid_polar_c_pointer(c, fptr)
        if (ok) fptr = f
    end function fitpack_grid_polar_c_f2c

    !> Create from Fortran object (optionally as pointer)
    type(fitpack_grid_polar_c) function fitpack_grid_polar_c_fptr2c(f, want_pointer) result(c)
        type(fitpack_grid_polar), intent(inout), target :: f
        logical, intent(in) :: want_pointer

        if (want_pointer) then
            c%cptr = c_loc(f)
            c%is_pointer = .true._c_bool
        else
            c = fitpack_grid_polar_c_f2c(f)
        end if
    end function fitpack_grid_polar_c_fptr2c

    ! ===========================================================================================
    ! SECTION 6b: C String Utilities (emitted on demand — only the helpers
    ! this module's wrappers actually call)
    ! ===========================================================================================

    !> Pointer view from a `type(c_ptr)` C string. Aliases the C buffer;
    !> always associated (0-length sentinel for null/empty input). The
    !> sentinel is a NUL byte viewed as zero-length, so it is simultaneously
    !> an empty Fortran string and a valid empty C string: a callee that
    !> re-derives the length with `strlen` (any `character(len=1) :: s(*)`
    !> dummy) stops at the NUL instead of reading on into adjacent statics.
    function from_c_string_ptr_cptr(c_str) result(fstring)
        character(len=:, kind=c_char), pointer :: fstring
        type(c_ptr), intent(in) :: c_str

        interface
            integer(c_size_t) function strlen(s) bind(C, name='strlen')
                import :: c_ptr, c_size_t
                type(c_ptr), intent(in), value :: s
            end function strlen
        end interface

        character(len=1, kind=c_char), target, save :: NO_STRING = c_null_char

        character(kind=c_char), pointer :: arr(:)
        integer(c_size_t) :: ltot

        fstring => NO_STRING(1:0)
        if (.not. c_associated(c_str)) return

        ltot = strlen(c_str)
        if (ltot == 0) return

        call c_f_pointer(c_str, arr, [ltot])
        if (.not. associated(arr)) return

        call c_array_to_string(int(ltot), arr, fstring)
    end function from_c_string_ptr_cptr

    !> Turn array of characters to a character string
    subroutine c_array_to_string(scalar_len, scalar, ptr)
        integer, intent(in) :: scalar_len
        character(kind=c_char, len=scalar_len), intent(in), target :: scalar(1)
        character(:, kind=c_char), intent(out), pointer :: ptr
        ptr => scalar(1)
    end subroutine c_array_to_string

    !> Pointer view from a `character(*)` assumed-size C buffer. The dummy is
    !> optional: an absent argument (NULL passed to the optional bind(C)
    !> dummy) maps to the 0-length sentinel — always associated, safe for
    !> required `character(*)` callee dummies.
    function from_c_string_ptr_arr(c_str) result(fstring)
        character(len=:, kind=c_char), pointer :: fstring
        character(len=1, kind=c_char), dimension(*), intent(in), optional, target :: c_str

        if (present(c_str)) then
            fstring => from_c_string_ptr_cptr(c_loc(c_str))
        else
            fstring => from_c_string_ptr_cptr(c_null_ptr)
        end if
    end function from_c_string_ptr_arr

    ! ===========================================================================================
    ! SECTION 7: Method Wrappers (bind(C))
    ! ===========================================================================================
    !> set_origin_BC
    subroutine fitpack_grid_polar_c_set_origin_BC(this, z0, exact, differentiable) bind(C,  &
        & name='fitpack_grid_polar_c_set_origin_BC')
        type(fitpack_grid_polar_c), intent(inout) :: this
        real(c_double), intent(in), optional :: z0
        logical(c_bool), intent(in), optional :: exact, differentiable
        type(fitpack_grid_polar), pointer :: fthis
        logical :: f_exact, f_differentiable

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(exact)) f_exact = logical(exact)
            if (present(differentiable)) f_differentiable = logical(differentiable)
            if (present(exact) .and. present(differentiable)) then
                call fthis%set_origin_BC(z0=z0, exact=f_exact, differentiable=f_differentiable)
            else
                call fthis%set_origin_BC(z0=z0)
            end if
        end if
    end subroutine fitpack_grid_polar_c_set_origin_BC

    !> fit
    integer(c_int32_t) function fitpack_grid_polar_c_fit(this, smoothing, keep_knots) bind(C, name='fitpack_grid_polar_c_fit')
        type(fitpack_grid_polar_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: keep_knots
        type(fitpack_grid_polar), pointer :: fthis
        logical :: f_keep_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(keep_knots)) f_keep_knots = logical(keep_knots)
            if (present(keep_knots)) then
                fitpack_grid_polar_c_fit = fthis%fit(smoothing=smoothing, keep_knots=f_keep_knots)
            else
                fitpack_grid_polar_c_fit = fthis%fit(smoothing=smoothing)
            end if
        else
            fitpack_grid_polar_c_fit = 0
        end if
    end function fitpack_grid_polar_c_fit

    !> least_squares
    integer(c_int32_t) function fitpack_grid_polar_c_least_squares(this, smoothing, reset_knots) bind(C,  &
        & name='fitpack_grid_polar_c_least_squares')
        type(fitpack_grid_polar_c), intent(in) :: this
        real(c_double), intent(in), optional :: smoothing
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_grid_polar), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_grid_polar_c_least_squares = fthis%least_squares(smoothing=smoothing, reset_knots=f_reset_knots)
            else
                fitpack_grid_polar_c_least_squares = fthis%least_squares(smoothing=smoothing)
            end if
        else
            fitpack_grid_polar_c_least_squares = 0
        end if
    end function fitpack_grid_polar_c_least_squares

    !> interpolate
    integer(c_int32_t) function fitpack_grid_polar_c_interpolate(this, reset_knots) bind(C, name='fitpack_grid_polar_c_interpolate')
        type(fitpack_grid_polar_c), intent(in) :: this
        logical(c_bool), intent(in), optional :: reset_knots
        type(fitpack_grid_polar), pointer :: fthis
        logical :: f_reset_knots

        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (present(reset_knots)) f_reset_knots = logical(reset_knots)
            if (present(reset_knots)) then
                fitpack_grid_polar_c_interpolate = fthis%interpolate(reset_knots=f_reset_knots)
            else
                fitpack_grid_polar_c_interpolate = fthis%interpolate()
            end if
        else
            fitpack_grid_polar_c_interpolate = 0
        end if
    end function fitpack_grid_polar_c_interpolate

    !> gridded_eval_one
    real(c_double) function fitpack_grid_polar_c_gridded_eval_one(this, u, v, ierr) bind(C,  &
        & name='fitpack_grid_polar_c_gridded_eval_one')
        type(fitpack_grid_polar_c), intent(in) :: this
        real(c_double), intent(in), value :: u, v
        integer(c_int32_t), intent(out), optional :: ierr
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_gridded_eval_one = fthis%eval(u, v, ierr)
        else
            fitpack_grid_polar_c_gridded_eval_one = 0.0_c_double
        end if
    end function fitpack_grid_polar_c_gridded_eval_one

    !> comm_size
    integer(c_int32_t) function fitpack_grid_polar_c_comm_size(this) bind(C, name='fitpack_grid_polar_c_comm_size')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_comm_size = fthis%comm_size()
        else
            fitpack_grid_polar_c_comm_size = 0
        end if
    end function fitpack_grid_polar_c_comm_size

    !> write
    subroutine fitpack_grid_polar_c_write(this, filename) bind(C, name='fitpack_grid_polar_c_write')
        type(fitpack_grid_polar_c), intent(inout) :: this
        character(kind=c_char), dimension(*), optional, intent(in), target :: filename
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%write(f_char(filename))
        end if
    end subroutine fitpack_grid_polar_c_write

    !> mse
    real(c_double) function fitpack_grid_polar_c_mse(this) bind(C, name='fitpack_grid_polar_c_mse')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_mse = fthis%mse()
        else
            fitpack_grid_polar_c_mse = 0.0_c_double
        end if
    end function fitpack_grid_polar_c_mse

    !> core_comm_size
    integer(c_int32_t) function fitpack_grid_polar_c_core_comm_size(this) bind(C, name='fitpack_grid_polar_c_core_comm_size')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_grid_polar_c_core_comm_size = 0
        end if
    end function fitpack_grid_polar_c_core_comm_size

    !> destroy_base
    subroutine fitpack_grid_polar_c_destroy_base(this) bind(C, name='fitpack_grid_polar_c_destroy_base')
        type(fitpack_grid_polar_c), intent(inout) :: this
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_grid_polar_c_destroy_base

    !> new_points
    subroutine fitpack_grid_polar_c_new_points(this, n, u, v, r, z_n1, z_n2, z, z0) bind(C, name='fitpack_grid_polar_c_new_points')
        type(fitpack_grid_polar_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n, z_n1, z_n2
        real(c_double), intent(in) :: u(n)
        real(c_double), intent(in) :: v(n)
        real(c_double), intent(in), value :: r
        real(c_double), intent(in) :: z(z_n1, z_n2)
        real(c_double), intent(in), optional :: z0
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%new_points(u, v, r, z, z0)
        end if
    end subroutine fitpack_grid_polar_c_new_points

    !> new_fit
    integer(c_int32_t) function fitpack_grid_polar_c_new_fit(this, n, u, v, r, z_n1, z_n2, z, z0, smoothing) bind(C,  &
        & name='fitpack_grid_polar_c_new_fit')
        type(fitpack_grid_polar_c), intent(in) :: this
        integer(c_int32_t), intent(in), value :: n, z_n1, z_n2
        real(c_double), intent(in) :: u(n)
        real(c_double), intent(in) :: v(n)
        real(c_double), intent(in), value :: r
        real(c_double), intent(in) :: z(z_n1, z_n2)
        real(c_double), intent(in), optional :: z0, smoothing
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_new_fit = fthis%new_fit(u, v, r, z, z0, smoothing)
        else
            fitpack_grid_polar_c_new_fit = 0
        end if
    end function fitpack_grid_polar_c_new_fit

    !> gridded_eval_many
    subroutine fitpack_grid_polar_c_gridded_eval_many(this, n, u, v, ierr, result, n_result, max_size) bind(C,  &
        & name='fitpack_grid_polar_c_gridded_eval_many')
        type(fitpack_grid_polar_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: u(n)
        real(c_double), intent(in) :: v(n)
        integer(c_int32_t), intent(out), optional :: ierr
        real(c_double), intent(inout) :: result(*)
        integer(c_int32_t), intent(out) :: n_result
        integer(c_int32_t), intent(in), value :: max_size
        real(c_double), allocatable :: fresult(:,:)
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            fresult = fthis%eval(u=u, v=v, ierr=ierr)
            n_result = size(fresult)
            result(1:min(n_result, max_size)) = reshape(fresult, [min(n_result, max_size)])
        end if
    end subroutine fitpack_grid_polar_c_gridded_eval_many

    !> comm_pack
    subroutine fitpack_grid_polar_c_comm_pack(this, n, buffer) bind(C, name='fitpack_grid_polar_c_comm_pack')
        type(fitpack_grid_polar_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_pack(buffer)
        end if
    end subroutine fitpack_grid_polar_c_comm_pack

    !> comm_expand
    subroutine fitpack_grid_polar_c_comm_expand(this, n, buffer) bind(C, name='fitpack_grid_polar_c_comm_expand')
        type(fitpack_grid_polar_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%comm_expand(buffer)
        end if
    end subroutine fitpack_grid_polar_c_comm_expand

    !> core_comm_pack
    subroutine fitpack_grid_polar_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_grid_polar_c_core_comm_pack')
        type(fitpack_grid_polar_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_grid_polar_c_core_comm_pack

    !> core_comm_expand
    subroutine fitpack_grid_polar_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_grid_polar_c_core_comm_expand')
        type(fitpack_grid_polar_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        type(fitpack_grid_polar), pointer :: fthis

        fthis => f_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_grid_polar_c_core_comm_expand

    ! ===========================================================================================
    ! SECTION 8: Component Array Accessors (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'u': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_grid_polar_c_getcomp_u(this, ptr, extents) &
        bind(C, name='fitpack_grid_polar_c_getcomp_u')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_grid_polar), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%u)) then
                extents = int(shape(fthis%u), kind=c_int64_t)
                if (size(fthis%u) > 0) ptr = c_loc(fthis%u)
            end if
        end if
    end subroutine fitpack_grid_polar_c_getcomp_u

    !> Raw view of component 'v': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_grid_polar_c_getcomp_v(this, ptr, extents) &
        bind(C, name='fitpack_grid_polar_c_getcomp_v')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        type(fitpack_grid_polar), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%v)) then
                extents = int(shape(fthis%v), kind=c_int64_t)
                if (size(fthis%v) > 0) ptr = c_loc(fthis%v)
            end if
        end if
    end subroutine fitpack_grid_polar_c_getcomp_v

    !> Raw view of component 'z': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_grid_polar_c_getcomp_z(this, ptr, extents) &
        bind(C, name='fitpack_grid_polar_c_getcomp_z')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_grid_polar), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%z)) then
                extents = int(shape(fthis%z), kind=c_int64_t)
                if (size(fthis%z) > 0) ptr = c_loc(fthis%z)
            end if
        end if
    end subroutine fitpack_grid_polar_c_getcomp_z

    !> Raw view of component 't': address of the first element, plus
    !> the component's 2 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_grid_polar_c_getcomp_t(this, ptr, extents) &
        bind(C, name='fitpack_grid_polar_c_getcomp_t')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(2)
        type(fitpack_grid_polar), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => f_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%t)) then
                extents = int(shape(fthis%t), kind=c_int64_t)
                if (size(fthis%t) > 0) ptr = c_loc(fthis%t)
            end if
        end if
    end subroutine fitpack_grid_polar_c_getcomp_t

    ! ===========================================================================================
    ! SECTION 9: Scalar Property Accessors
    ! ===========================================================================================

    !> Get pointer to scalar property 'r'
    type(c_ptr) function fitpack_grid_polar_c_ref_r(this) &
        bind(C, name='fitpack_grid_polar_c_ref_r')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_ref_r = c_loc(fthis%r)
        else
            fitpack_grid_polar_c_ref_r = c_null_ptr
        end if
    end function fitpack_grid_polar_c_ref_r

    !> Get pointer to scalar property 'z0'
    type(c_ptr) function fitpack_grid_polar_c_ref_z0(this) &
        bind(C, name='fitpack_grid_polar_c_ref_z0')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_ref_z0 = c_loc(fthis%z0)
        else
            fitpack_grid_polar_c_ref_z0 = c_null_ptr
        end if
    end function fitpack_grid_polar_c_ref_z0

    !> Get logical property 'z0_present'
    logical(c_bool) function fitpack_grid_polar_c_get_z0_present(this) &
        bind(C, name='fitpack_grid_polar_c_get_z0_present')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_get_z0_present = logical(fthis%z0_present, c_bool)
        else
            fitpack_grid_polar_c_get_z0_present = .false._c_bool
        end if
    end function fitpack_grid_polar_c_get_z0_present

    !> Set logical property 'z0_present'
    subroutine fitpack_grid_polar_c_set_z0_present(this, value) &
        bind(C, name='fitpack_grid_polar_c_set_z0_present')
        type(fitpack_grid_polar_c), intent(inout) :: this
        logical(c_bool), intent(in), value :: value
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fthis%z0_present = logical(value)
        end if
    end subroutine fitpack_grid_polar_c_set_z0_present

    !> Get logical property 'z0_exact'
    logical(c_bool) function fitpack_grid_polar_c_get_z0_exact(this) &
        bind(C, name='fitpack_grid_polar_c_get_z0_exact')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_get_z0_exact = logical(fthis%z0_exact, c_bool)
        else
            fitpack_grid_polar_c_get_z0_exact = .false._c_bool
        end if
    end function fitpack_grid_polar_c_get_z0_exact

    !> Set logical property 'z0_exact'
    subroutine fitpack_grid_polar_c_set_z0_exact(this, value) &
        bind(C, name='fitpack_grid_polar_c_set_z0_exact')
        type(fitpack_grid_polar_c), intent(inout) :: this
        logical(c_bool), intent(in), value :: value
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fthis%z0_exact = logical(value)
        end if
    end subroutine fitpack_grid_polar_c_set_z0_exact

    !> Get logical property 'z0_zero_gradient'
    logical(c_bool) function fitpack_grid_polar_c_get_z0_zero_gradient(this) &
        bind(C, name='fitpack_grid_polar_c_get_z0_zero_gradient')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_get_z0_zero_gradient = logical(fthis%z0_zero_gradient, c_bool)
        else
            fitpack_grid_polar_c_get_z0_zero_gradient = .false._c_bool
        end if
    end function fitpack_grid_polar_c_get_z0_zero_gradient

    !> Set logical property 'z0_zero_gradient'
    subroutine fitpack_grid_polar_c_set_z0_zero_gradient(this, value) &
        bind(C, name='fitpack_grid_polar_c_set_z0_zero_gradient')
        type(fitpack_grid_polar_c), intent(inout) :: this
        logical(c_bool), intent(in), value :: value
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fthis%z0_zero_gradient = logical(value)
        end if
    end subroutine fitpack_grid_polar_c_set_z0_zero_gradient

    !> Get pointer to scalar property 'nmax'
    type(c_ptr) function fitpack_grid_polar_c_ref_nmax(this) &
        bind(C, name='fitpack_grid_polar_c_ref_nmax')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_ref_nmax = c_loc(fthis%nmax)
        else
            fitpack_grid_polar_c_ref_nmax = c_null_ptr
        end if
    end function fitpack_grid_polar_c_ref_nmax

    !> Get pointer to scalar property 'bc_continuity_origin'
    type(c_ptr) function fitpack_grid_polar_c_ref_bc_continuity_origin(this) &
        bind(C, name='fitpack_grid_polar_c_ref_bc_continuity_origin')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_ref_bc_continuity_origin = c_loc(fthis%bc_continuity_origin)
        else
            fitpack_grid_polar_c_ref_bc_continuity_origin = c_null_ptr
        end if
    end function fitpack_grid_polar_c_ref_bc_continuity_origin

    !> Get pointer to scalar property 'bc_boundary'
    type(c_ptr) function fitpack_grid_polar_c_ref_bc_boundary(this) &
        bind(C, name='fitpack_grid_polar_c_ref_bc_boundary')
        type(fitpack_grid_polar_c), intent(in) :: this
        type(fitpack_grid_polar), pointer :: fthis
        fthis => f_pointer(this)
        if (associated(fthis)) then
            fitpack_grid_polar_c_ref_bc_boundary = c_loc(fthis%bc_boundary)
        else
            fitpack_grid_polar_c_ref_bc_boundary = c_null_ptr
        end if
    end function fitpack_grid_polar_c_ref_bc_boundary

end module fitpack_gridded_polar_c
