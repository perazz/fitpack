!   ***********************************************************************************************
!   **                                        fxArray                                          **
!   **                                  Fortran Arrays for C++                                 **
!   ***********************************************************************************************
!   **    fitpack_gridded_polar_c_types                                                                   **
!> @brief C wrapper type definitions for fitpack_gridded_polar
!   ***********************************************************************************************
!> @author Binding Generator
!   ***********************************************************************************************

module fitpack_gridded_polar_c_types
    use, intrinsic :: iso_c_binding
    implicit none(type, external)

    !> Opaque C wrapper for type(fitpack_grid_polar).
    type, public, bind(C) :: fitpack_grid_polar_c
        type(c_ptr)      :: cptr = c_null_ptr
        logical(c_bool)  :: is_pointer = .false._c_bool
        type(c_ptr)      :: name_cptr = c_null_ptr
    end type fitpack_grid_polar_c

    type(fitpack_grid_polar_c), parameter, public :: fitpack_grid_polar_c_null = &
        fitpack_grid_polar_c(cptr=c_null_ptr, is_pointer=.false._c_bool, name_cptr=c_null_ptr)

    character(len=1, kind=c_char), target, protected, &
            bind(C, name='fitpack_grid_polar_c_typename') :: &
        fitpack_grid_polar_c_typename(19) = &
            transfer("fitpack_grid_polar" // c_null_char, [character(c_char)::], &
                     size=19)

    public :: fitpack_grid_polar_c_typename
    public :: fitpack_grid_polar_c_fortran_type_name
    public :: fitpack_grid_polar_c_c_type_name

    interface f_type_name
        module procedure fitpack_grid_polar_c_fortran_type_name
    end interface f_type_name
    public :: f_type_name

contains

    !> [bind(C)] Return a C pointer to the dynamic Fortran type-name string of `this`.
    function fitpack_grid_polar_c_fortran_type_name(this) result(name) &
            bind(C, name='fitpack_grid_polar_c_fortran_type_name')
        type(fitpack_grid_polar_c), intent(in), value :: this
        type(c_ptr) :: name
        if (c_associated(this%name_cptr)) then
            name = this%name_cptr
        else
            name = c_loc(fitpack_grid_polar_c_typename(1))
        end if
    end function fitpack_grid_polar_c_fortran_type_name

    !> [bind(C)] Return a C pointer to the static null-terminated C wrapper type name.
    !> The C wrapper struct lives in the same binding layer as the Fortran code,
    !> so its name is owned here. Higher-level language names (C++ class, Python
    !> class, ...) are owned by their respective generators.
    function fitpack_grid_polar_c_c_type_name(this) result(name) &
            bind(C, name='fitpack_grid_polar_c_c_type_name')
        type(fitpack_grid_polar_c), intent(in), value :: this
        type(c_ptr) :: name
        character(len=21, kind=c_char), target, save :: tname = &
            "fitpack_grid_polar_c" // c_null_char
        name = c_loc(tname)
    end function fitpack_grid_polar_c_c_type_name

end module fitpack_gridded_polar_c_types
