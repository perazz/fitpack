! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   Refactored by Federico Perini, 10/6/2022
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************
module fitpack_cpp_tests

    use iso_c_binding
    use fitpack_core

    implicit none
    private

    public :: test_cpp_sine_fit
    public :: test_cpp_periodic_fit
    public :: test_cpp_parametric_fit
    public :: test_cpp_closed_fit

    interface
        logical(FP_BOOL) function test_cpp_sine_fit() bind(C,name="test_cpp_sine_fit")
           import FP_BOOL
        end function test_cpp_sine_fit
        logical(FP_BOOL) function test_cpp_periodic_fit() bind(C,name="test_cpp_periodic_fit")
           import FP_BOOL
        end function test_cpp_periodic_fit
        logical(FP_BOOL) function test_cpp_parametric_fit() bind(C,name="test_cpp_parametric_fit")
           import FP_BOOL
        end function test_cpp_parametric_fit
        logical(FP_BOOL) function test_cpp_closed_fit() bind(C,name="test_cpp_closed_fit")
           import FP_BOOL
        end function test_cpp_closed_fit
    end interface

    contains

end module fitpack_cpp_tests
