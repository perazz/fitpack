/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_config.h
!> @brief Build-feature configuration for the generated C++ wrappers
!
!   Single source of truth for the optional-feature macros. Every generated .hpp includes
!   this header before it uses a gate, so the wrappers and their consumers can never
!   disagree on conditionally-compiled class layout.
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
! **************************************************************************************************/

#ifndef FITPACK_CONFIG_H
#define FITPACK_CONFIG_H

/* fxArray (Fortran no-copy views). Auto-enabled wherever fxArrays.hpp is on the include
 * path — the wrappers then expose a zero-copy fxArray view of every array component
 * alongside the std::vector copy. Off in a build with no fortran-arrays dependency, where
 * the copy faces are the whole surface. Define FITPACK_NO_FXARRAY to opt out.
 * The C ABI is identical either way: the gate is C++-side only. */
#if !defined(FITPACK_NO_FXARRAY) && defined(__has_include)
#  if __has_include("fxArrays.hpp")
#    define HAVE_FXARRAY 1
#  endif
#endif

//! @return true when the fxArray (Fortran no-copy views) surface is available.
inline constexpr bool fitpackHasFxArray()
{
#if HAVE_FXARRAY
    return true;
#else
    return false;
#endif
}

#endif // FITPACK_CONFIG_H
