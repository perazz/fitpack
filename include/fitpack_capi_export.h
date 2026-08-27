/* fitpack_capi_export.h
 *
 * Cross-platform export-visibility macro for the FITPACK C API.
 *
 * - On Windows / Cygwin, FITPACK_CAPI_EXPORT expands to dllexport when the
 *   umbrella shared library (fitpack_capi.dll) is being built (the build
 *   system defines FITPACK_CAPI_BUILDING), and to dllimport for consumers of
 *   the DLL. When FITPACK_CAPI_STATIC is defined it expands to nothing, so a
 *   consumer linking the merged static archive (libfitpack_capi.a) sees plain
 *   symbols with neither dllimport nor dllexport.
 * - On GCC / Clang (including Apple clang), FITPACK_CAPI_EXPORT expands to
 *   __attribute__((visibility("default"))) so that the symbol survives
 *   the umbrella's -fvisibility=hidden build setting.
 * - On other compilers / systems the macro expands to nothing.
 */

#ifndef FITPACK_CAPI_EXPORT_H
#define FITPACK_CAPI_EXPORT_H

#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(FITPACK_CAPI_STATIC)
#    define FITPACK_CAPI_EXPORT
#  elif defined(FITPACK_CAPI_BUILDING)
#    define FITPACK_CAPI_EXPORT __declspec(dllexport)
#  else
#    define FITPACK_CAPI_EXPORT __declspec(dllimport)
#  endif
#elif defined(__GNUC__) && (__GNUC__ >= 4 || defined(__clang__))
#  define FITPACK_CAPI_EXPORT __attribute__((visibility("default")))
#else
#  define FITPACK_CAPI_EXPORT
#endif

#endif /* FITPACK_CAPI_EXPORT_H */
