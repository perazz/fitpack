/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fitpack_config.h                                                                      **
 *   ** @brief Build-feature configuration for the generated C++ wrappers
 *   **
 *   ** Single source of truth for the optional-feature macros. Every generated .hpp includes
 *   ** this header before it uses a gate, so the wrappers and their consumers can never
 *   ** disagree on conditionally-compiled class layout.
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

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
