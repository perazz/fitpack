! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_curves_c_types.f90 (module fitpack_curves_c_types)
!> @brief C wrapper type definitions for fitpack_curves
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

module fitpack_curves_c_types
    use, intrinsic :: iso_c_binding
    implicit none(type, external)

    !> Opaque C wrapper for type(fitpack_curve).
    type, public, bind(C) :: fitpack_curve_c
        type(c_ptr)      :: cptr = c_null_ptr
        logical(c_bool)  :: is_pointer = .false._c_bool
        type(c_ptr)      :: name_cptr = c_null_ptr
    end type fitpack_curve_c

    type(fitpack_curve_c), parameter, public :: fitpack_curve_c_null = &
        fitpack_curve_c(cptr=c_null_ptr, is_pointer=.false._c_bool, name_cptr=c_null_ptr)

    character(len=1, kind=c_char), target, protected, &
            bind(C, name='fitpack_curve_c_typename') :: &
        fitpack_curve_c_typename(14) = &
            transfer("fitpack_curve" // c_null_char, [character(c_char)::], &
                     size=14)

    public :: fitpack_curve_c_typename
    public :: fitpack_curve_c_fortran_type_name
    public :: fitpack_curve_c_c_type_name

    !> Opaque C wrapper for type(fitpack_periodic_curve).
    type, public, bind(C) :: fitpack_periodic_curve_c
        type(c_ptr)      :: cptr = c_null_ptr
        logical(c_bool)  :: is_pointer = .false._c_bool
        type(c_ptr)      :: name_cptr = c_null_ptr
    end type fitpack_periodic_curve_c

    type(fitpack_periodic_curve_c), parameter, public :: fitpack_periodic_curve_c_null = &
        fitpack_periodic_curve_c(cptr=c_null_ptr, is_pointer=.false._c_bool, name_cptr=c_null_ptr)

    character(len=1, kind=c_char), target, protected, &
            bind(C, name='fitpack_periodic_curve_c_typename') :: &
        fitpack_periodic_curve_c_typename(23) = &
            transfer("fitpack_periodic_curve" // c_null_char, [character(c_char)::], &
                     size=23)

    public :: fitpack_periodic_curve_c_typename
    public :: fitpack_periodic_curve_c_fortran_type_name
    public :: fitpack_periodic_curve_c_c_type_name

    interface f_type_name
        module procedure fitpack_curve_c_fortran_type_name
        module procedure fitpack_periodic_curve_c_fortran_type_name
    end interface f_type_name
    public :: f_type_name

contains

    !> [bind(C)] Return a C pointer to the dynamic Fortran type-name string of `this`.
    function fitpack_curve_c_fortran_type_name(this) result(name) &
            bind(C, name='fitpack_curve_c_fortran_type_name')
        type(fitpack_curve_c), intent(in), value :: this
        type(c_ptr) :: name
        if (c_associated(this%name_cptr)) then
            name = this%name_cptr
        else
            name = c_loc(fitpack_curve_c_typename(1))
        end if
    end function fitpack_curve_c_fortran_type_name

    !> [bind(C)] Return a C pointer to the static null-terminated C wrapper type name.
    !> The C wrapper struct lives in the same binding layer as the Fortran code,
    !> so its name is owned here. Higher-level language names (C++ class, Python
    !> class, ...) are owned by their respective generators.
    function fitpack_curve_c_c_type_name(this) result(name) &
            bind(C, name='fitpack_curve_c_c_type_name')
        type(fitpack_curve_c), intent(in), value :: this
        type(c_ptr) :: name
        character(len=16, kind=c_char), target, save :: tname = &
            "fitpack_curve_c" // c_null_char
        name = c_loc(tname)
    end function fitpack_curve_c_c_type_name

    !> [bind(C)] Return a C pointer to the dynamic Fortran type-name string of `this`.
    function fitpack_periodic_curve_c_fortran_type_name(this) result(name) &
            bind(C, name='fitpack_periodic_curve_c_fortran_type_name')
        type(fitpack_periodic_curve_c), intent(in), value :: this
        type(c_ptr) :: name
        if (c_associated(this%name_cptr)) then
            name = this%name_cptr
        else
            name = c_loc(fitpack_periodic_curve_c_typename(1))
        end if
    end function fitpack_periodic_curve_c_fortran_type_name

    !> [bind(C)] Return a C pointer to the static null-terminated C wrapper type name.
    !> The C wrapper struct lives in the same binding layer as the Fortran code,
    !> so its name is owned here. Higher-level language names (C++ class, Python
    !> class, ...) are owned by their respective generators.
    function fitpack_periodic_curve_c_c_type_name(this) result(name) &
            bind(C, name='fitpack_periodic_curve_c_c_type_name')
        type(fitpack_periodic_curve_c), intent(in), value :: this
        type(c_ptr) :: name
        character(len=25, kind=c_char), target, save :: tname = &
            "fitpack_periodic_curve_c" // c_null_char
        name = c_loc(tname)
    end function fitpack_periodic_curve_c_c_type_name

end module fitpack_curves_c_types
