!   ***********************************************************************************************
!   **                                        fxArray                                          **
!   **                                  Fortran Arrays for C++                                 **
!   ***********************************************************************************************
!   **    fitpack_fitters_c                                                                       **
!> @brief C interface to module fitpack_fitters
!   ***********************************************************************************************
!> @author Binding Generator
!   ***********************************************************************************************

module fitpack_fitters_c
    ! ===========================================================================================
    ! SECTION 1: Imports
    ! ===========================================================================================
    use fitpack_fitters, only: fitpack_fitter
    use fitpack_fitters_c_types
    use fitpack_curves, only: fitpack_curve
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

    ! --- fitpack_fitter ---
    public :: fitpack_fitter_c_is_same
    public :: fitpack_fitter_c_mse
    public :: fitpack_fitter_c_core_comm_size
    public :: fitpack_fitter_c_destroy_base
    public :: fitpack_fitter_c_core_comm_pack
    public :: fitpack_fitter_c_core_comm_expand
    public :: fitpack_fitter_c_getcomp_c
    public :: fitpack_fitter_c_getcomp_wrk
    public :: fitpack_fitter_c_getcomp_iwrk
    public :: fitpack_fitter_c_ref_iopt
    public :: fitpack_fitter_c_ref_smoothing
    public :: fitpack_fitter_c_ref_fp
    public :: fitpack_fitter_c_ref_lwrk
    public :: fitpack_fitter_c_ref_liwrk

    interface f_pointer
        module procedure fitpack_fitter_c_get_pointer
    end interface f_pointer

    !> Per-consumer-module resolver(s) for polymorphic procedure args.
    !> One module procedure per abstract base any wrapper in this module
    !> consumes — colocated with the consumer so the resolver sees every
    !> bound concrete subtype known at THIS module's regen time
    !> (cross-package coverage included via the cumulative
    !> `all_type_parents` registry). Local (not exported); consumer
    !> wrappers call it via the local generic below.
    interface f_poly_pointer
        module procedure fitpack_fitter_c_get_poly_pointer
    end interface f_poly_pointer

    !> Pointer-identity check. Returns true iff `a` and `b` view the same
    !> underlying Fortran object (their internal c_ptrs match and are non-null).
    !> Useful for C-API consumers that want object identity without exposing
    !> raw pointers — e.g. GUI panels checking whether an input view points
    !> at the same Fortran allocation as the project that owns it.
    interface f_associated
        module procedure fitpack_fitter_c_is_same
    end interface f_associated

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
    ! C String Utilities (emitted on demand — only the helpers this module's
    ! wrappers actually call)
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
    ! FITPACK_FITTER: Core Memory Management
    ! ===========================================================================================

    !> Pointer-identity check: true iff both wrappers point to the same
    !> underlying Fortran object (and that object is allocated).
    logical(c_bool) function fitpack_fitter_c_is_same(this, that) &
            bind(C, name='fitpack_fitter_c_is_same') result(same)
        type(fitpack_fitter_c), intent(in) :: this, that
        same = logical(c_associated(this%cptr, that%cptr), kind=c_bool)
    end function fitpack_fitter_c_is_same

    !> Get opaque pointer to embedded Fortran object (non-allocating)
    !> Abstract types cannot use c_f_pointer; return the raw c_ptr instead
    function fitpack_fitter_c_get_pointer(this) result(fptr)
        type(fitpack_fitter_c), intent(in) :: this
        type(c_ptr) :: fptr

        if (c_associated(this%cptr)) then
            fptr = this%cptr
        else
            fptr = c_null_ptr
        end if
    end function fitpack_fitter_c_get_pointer

    !> Concrete-typed pointer for reading base-declared components on this
    !> proxy-less abstract base. Casts the opaque handle to the first concrete
    !> child (fitpack_curve); a base component sits at the same
    !> offset in every descendant, so the read is valid for any runtime subtype.
    function fitpack_fitter_c_get_accessor_pointer(this) result(fptr)
        type(fitpack_fitter_c), intent(in) :: this
        type(fitpack_curve), pointer :: fptr

        nullify(fptr)
        if (c_associated(this%cptr)) call c_f_pointer(this%cptr, fptr)
    end function fitpack_fitter_c_get_accessor_pointer

    ! ===========================================================================================
    ! FITPACK_FITTER: Method Wrappers
    ! ===========================================================================================
    real(c_double) function fitpack_fitter_c_mse(this) bind(C, name='fitpack_fitter_c_mse')
        type(fitpack_fitter_c), intent(in) :: this
        class(fitpack_fitter), pointer :: fthis

        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            fitpack_fitter_c_mse = fthis%mse()
        else
            fitpack_fitter_c_mse = 0.0_c_double
        end if
    end function fitpack_fitter_c_mse

    integer(c_int32_t) function fitpack_fitter_c_core_comm_size(this) bind(C, name='fitpack_fitter_c_core_comm_size')
        type(fitpack_fitter_c), intent(in) :: this
        class(fitpack_fitter), pointer :: fthis

        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            fitpack_fitter_c_core_comm_size = fthis%core_comm_size()
        else
            fitpack_fitter_c_core_comm_size = 0
        end if
    end function fitpack_fitter_c_core_comm_size

    subroutine fitpack_fitter_c_destroy_base(this) bind(C, name='fitpack_fitter_c_destroy_base')
        type(fitpack_fitter_c), intent(inout) :: this
        class(fitpack_fitter), pointer :: fthis

        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            call fthis%destroy_base()
        end if
    end subroutine fitpack_fitter_c_destroy_base

    subroutine fitpack_fitter_c_core_comm_pack(this, n, buffer) bind(C, name='fitpack_fitter_c_core_comm_pack')
        type(fitpack_fitter_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(inout) :: buffer(n)
        class(fitpack_fitter), pointer :: fthis

        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_pack(buffer)
        end if
    end subroutine fitpack_fitter_c_core_comm_pack

    subroutine fitpack_fitter_c_core_comm_expand(this, n, buffer) bind(C, name='fitpack_fitter_c_core_comm_expand')
        type(fitpack_fitter_c), intent(inout) :: this
        integer(c_int32_t), intent(in), value :: n
        real(c_double), intent(in) :: buffer(n)
        class(fitpack_fitter), pointer :: fthis

        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            call fthis%core_comm_expand(buffer)
        end if
    end subroutine fitpack_fitter_c_core_comm_expand

    ! ===========================================================================================
    ! Component Array Accessors for fitpack_fitter (raw pointer + extents)
    ! ===========================================================================================

    !> Raw view of component 'c': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_fitter_c_getcomp_c(this, ptr, extents) &
        bind(C, name='fitpack_fitter_c_getcomp_c')
        type(fitpack_fitter_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        class(fitpack_fitter), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%c)) then
                extents = int(shape(fthis%c), kind=c_int64_t)
                if (size(fthis%c) > 0) ptr = c_loc(fthis%c)
            end if
        end if
    end subroutine fitpack_fitter_c_getcomp_c

    !> Raw view of component 'wrk': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_fitter_c_getcomp_wrk(this, ptr, extents) &
        bind(C, name='fitpack_fitter_c_getcomp_wrk')
        type(fitpack_fitter_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        class(fitpack_fitter), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%wrk)) then
                extents = int(shape(fthis%wrk), kind=c_int64_t)
                if (size(fthis%wrk) > 0) ptr = c_loc(fthis%wrk)
            end if
        end if
    end subroutine fitpack_fitter_c_getcomp_wrk

    !> Raw view of component 'iwrk': address of the first element, plus
    !> the component's 1 extent(s). Null pointer and zero extents when the
    !> object or the component is unallocated; an allocated but zero-sized
    !> component still reports its true shape, but has no representable address.
    !> Symbol uses the `getcomp_` infix (not `get_`) so it can never collide
    !> with a type-bound procedure named `get_<component>`.
    subroutine fitpack_fitter_c_getcomp_iwrk(this, ptr, extents) &
        bind(C, name='fitpack_fitter_c_getcomp_iwrk')
        type(fitpack_fitter_c), intent(in) :: this
        type(c_ptr), intent(out) :: ptr
        integer(c_int64_t), intent(out) :: extents(1)
        class(fitpack_fitter), pointer :: fthis

        ptr = c_null_ptr
        extents = 0_c_int64_t
        fthis => fitpack_fitter_c_get_poly_pointer(this)
        if (associated(fthis)) then
            if (allocated(fthis%iwrk)) then
                extents = int(shape(fthis%iwrk), kind=c_int64_t)
                if (size(fthis%iwrk) > 0) ptr = c_loc(fthis%iwrk)
            end if
        end if
    end subroutine fitpack_fitter_c_getcomp_iwrk

    ! ===========================================================================================
    ! Scalar Property Accessors for fitpack_fitter
    ! ===========================================================================================

    !> Get pointer to scalar property 'iopt'
    type(c_ptr) function fitpack_fitter_c_ref_iopt(this) &
        bind(C, name='fitpack_fitter_c_ref_iopt')
        type(fitpack_fitter_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_fitter_c_get_accessor_pointer(this)
        if (associated(fthis)) then
            fitpack_fitter_c_ref_iopt = c_loc(fthis%iopt)
        else
            fitpack_fitter_c_ref_iopt = c_null_ptr
        end if
    end function fitpack_fitter_c_ref_iopt

    !> Get pointer to scalar property 'smoothing'
    type(c_ptr) function fitpack_fitter_c_ref_smoothing(this) &
        bind(C, name='fitpack_fitter_c_ref_smoothing')
        type(fitpack_fitter_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_fitter_c_get_accessor_pointer(this)
        if (associated(fthis)) then
            fitpack_fitter_c_ref_smoothing = c_loc(fthis%smoothing)
        else
            fitpack_fitter_c_ref_smoothing = c_null_ptr
        end if
    end function fitpack_fitter_c_ref_smoothing

    !> Get pointer to scalar property 'fp'
    type(c_ptr) function fitpack_fitter_c_ref_fp(this) &
        bind(C, name='fitpack_fitter_c_ref_fp')
        type(fitpack_fitter_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_fitter_c_get_accessor_pointer(this)
        if (associated(fthis)) then
            fitpack_fitter_c_ref_fp = c_loc(fthis%fp)
        else
            fitpack_fitter_c_ref_fp = c_null_ptr
        end if
    end function fitpack_fitter_c_ref_fp

    !> Get pointer to scalar property 'lwrk'
    type(c_ptr) function fitpack_fitter_c_ref_lwrk(this) &
        bind(C, name='fitpack_fitter_c_ref_lwrk')
        type(fitpack_fitter_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_fitter_c_get_accessor_pointer(this)
        if (associated(fthis)) then
            fitpack_fitter_c_ref_lwrk = c_loc(fthis%lwrk)
        else
            fitpack_fitter_c_ref_lwrk = c_null_ptr
        end if
    end function fitpack_fitter_c_ref_lwrk

    !> Get pointer to scalar property 'liwrk'
    type(c_ptr) function fitpack_fitter_c_ref_liwrk(this) &
        bind(C, name='fitpack_fitter_c_ref_liwrk')
        type(fitpack_fitter_c), intent(in) :: this
        type(fitpack_curve), pointer :: fthis
        fthis => fitpack_fitter_c_get_accessor_pointer(this)
        if (associated(fthis)) then
            fitpack_fitter_c_ref_liwrk = c_loc(fthis%liwrk)
        else
            fitpack_fitter_c_ref_liwrk = c_null_ptr
        end if
    end function fitpack_fitter_c_ref_liwrk

    function fitpack_fitter_c_get_poly_pointer(this) result(fptr)
        type(fitpack_fitter_c), intent(in) :: this
        class(fitpack_fitter), pointer :: fptr

        nullify(fptr)
        if (.not. c_associated(this%cptr)) return
        if (.not. c_associated(this%name_cptr)) return
        select case (f_char(this%name_cptr))
        case ('fitpack_periodic_curve')
            block
                use fitpack_curves, only: fitpack_periodic_curve
                type(fitpack_periodic_curve), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_convex_curve')
            block
                use fitpack_convex_curves, only: fitpack_convex_curve
                type(fitpack_convex_curve), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_curve')
            block
                use fitpack_curves, only: fitpack_curve
                type(fitpack_curve), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_closed_curve')
            block
                use fitpack_parametric_curves, only: fitpack_closed_curve
                type(fitpack_closed_curve), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_constrained_curve')
            block
                use fitpack_parametric_curves, only: fitpack_constrained_curve
                type(fitpack_constrained_curve), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_parametric_curve')
            block
                use fitpack_parametric_curves, only: fitpack_parametric_curve
                type(fitpack_parametric_curve), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_surface')
            block
                use fitpack_surfaces, only: fitpack_surface
                type(fitpack_surface), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_grid_surface')
            block
                use fitpack_gridded_surfaces, only: fitpack_grid_surface
                type(fitpack_grid_surface), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_gridded_spline')
            block
                use fitpack_gridded_splines, only: fitpack_gridded_spline
                type(fitpack_gridded_spline), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_parametric_surface')
            block
                use fitpack_parametric_surfaces, only: fitpack_parametric_surface
                type(fitpack_parametric_surface), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_polar')
            block
                use fitpack_polar_domains, only: fitpack_polar
                type(fitpack_polar), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_grid_polar')
            block
                use fitpack_gridded_polar, only: fitpack_grid_polar
                type(fitpack_grid_polar), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_sphere')
            block
                use fitpack_sphere_domains, only: fitpack_sphere
                type(fitpack_sphere), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case ('fitpack_grid_sphere')
            block
                use fitpack_gridded_sphere, only: fitpack_grid_sphere
                type(fitpack_grid_sphere), pointer :: p
                call c_f_pointer(this%cptr, p)
                fptr => p
            end block
        case default
            nullify(fptr)
        end select
    end function fitpack_fitter_c_get_poly_pointer

end module fitpack_fitters_c
