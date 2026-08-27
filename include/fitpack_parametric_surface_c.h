/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_parametric_surface_c.h (module fitpack_parametric_surfaces)
!> @brief Standalone C interface to fitpack_parametric_surface (no fortran-arrays dependency)
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

#ifndef FITPACK_PARAMETRIC_SURFACE_C_H_INCLUDED
#define FITPACK_PARAMETRIC_SURFACE_C_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

/* Minimal fp_status definition for standalone use.
 * Layout contract: this struct is the C half of type(fx_status) in the
 * generated <project>_fx_status module (templates/fortran_fx_status.f90.jinja2),
 * which in turn mirrors fortran-arrays' arrays_c. Field order, widths and
 * FX_LEN_STATUS_MSG must change in all three places at once.
 * ARRAYS_C_H_INCLUDED: a gated build pulls the real arrays_c.h into the
 * translation unit ahead of this header (via fxArrays.hpp), and it defines the
 * same struct and macro behind no FX_STATUS_DEFINED guard — so where both are
 * present the real definition wins and this copy stands down. */
#ifndef FP_STATUS_TYPEDEF_INCLUDED
#define FP_STATUS_TYPEDEF_INCLUDED
#if defined(FX_STATUS_DEFINED) || defined(ARRAYS_C_H_INCLUDED)
typedef fx_status fp_status;               /* the real fortran-arrays struct */
#else
#define FX_LEN_STATUS_MSG 248
typedef struct fp_status { bool ok; int code; char message[FX_LEN_STATUS_MSG]; } fp_status;
typedef fp_status fx_status;               /* layout-identical interop alias */
#define FX_STATUS_DEFINED
#endif
#endif /* FP_STATUS_TYPEDEF_INCLUDED */

#include "fitpack_capi_export.h"
#include "fitpack_parametric_surfaces_c_types.h"  /* For fitpack_parametric_surface_c, fitpack_parametric_surface_c_null */
#include "fitpack_fitter_c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================================
 * Core Memory Management
 *
 * All functions accept an optional status parameter (pass NULL to use error stop on failure).
 * If status is provided, caller is responsible for checking status->ok after the call.
 * =========================================================================================== */

/**
 * @brief Allocate new fitpack_parametric_surface object
 * @param self Pointer to wrapper (will be initialized)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_allocate(fitpack_parametric_surface_c* self, fp_status* status);

/**
 * @brief Deallocate fitpack_parametric_surface object
 * @param self Pointer to wrapper (will be nullified)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_destroy(fitpack_parametric_surface_c* self, fp_status* status);

/**
 * @brief Copy fitpack_parametric_surface object.
 *
 * Default behavior mirrors the source ownership: a view source yields a
 * view (shallow handle copy), an owned source yields a deep copy via
 * Fortran intrinsic assignment. Pass `deep_copy=true` to force a deep
 * copy regardless of the source's ownership — useful for capturing an
 * independent snapshot of a member-view accessor's return value.
 *
 * @param self      Destination wrapper
 * @param other     Source wrapper
 * @param deep_copy When true, always allocate a fresh Fortran object and
 *                  deep-copy data, even if the source is a view.
 * @param status    Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_copy(fitpack_parametric_surface_c* self, const fitpack_parametric_surface_c* other, bool deep_copy, fp_status* status);

/**
 * @brief Shallow copy (pointer semantics — Fortran "associate" construct)
 * @param self Destination wrapper (will point to same object as 'other')
 * @param other Source wrapper (read-only)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_associate(fitpack_parametric_surface_c* self, const fitpack_parametric_surface_c* other, fp_status* status);

/**
 * @brief Move allocation (transfer ownership)
 * @param to Destination wrapper (receives ownership)
 * @param from Source wrapper (becomes null)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_move_alloc(fitpack_parametric_surface_c* to, fitpack_parametric_surface_c* from, fp_status* status);

/* ===========================================================================================
 * Method Wrappers (standalone-compatible only)
 * =========================================================================================== */

/**
 * @brief interpolate
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_parametric_surface_c_interpolate(fitpack_parametric_surface_c* self, bool* reset_knots);

/**
 * @brief comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_parametric_surface_c_comm_size(const fitpack_parametric_surface_c* self);

/**
 * @brief mse
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_parametric_surface_c_mse(const fitpack_parametric_surface_c* self);

/**
 * @brief core_comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_parametric_surface_c_core_comm_size(const fitpack_parametric_surface_c* self);

/**
 * @brief destroy_base
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_destroy_base(fitpack_parametric_surface_c* self);

/**
 * @brief fit
 * @param smoothing 
 * @param periodic 
 * @param keep_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_parametric_surface_c_fit(fitpack_parametric_surface_c* self, double* smoothing, int32_t n, bool* periodic, bool* keep_knots);

/**
 * @brief least_squares
 * @param u_knots 
 * @param v_knots 
 * @param smoothing 
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_parametric_surface_c_least_squares(fitpack_parametric_surface_c* self, int32_t n, double* u_knots, double* v_knots, double* smoothing, bool* reset_knots);

/**
 * @brief surf_eval_one
 * @param u 
 * @param v 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of (this % idim) — known to the caller from object state / inputs
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_surf_eval_one(fitpack_parametric_surface_c* self, double u, double v, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief surf_eval_grid
 * @param u 
 * @param v 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of (SIZE(v), SIZE(u), this % idim) — known to the caller from object state / inputs
 * @note result is the rank-3 Fortran result flattened column-major
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_surf_eval_grid(fitpack_parametric_surface_c* self, int32_t n, const double* u, const double* v, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_comm_pack(const fitpack_parametric_surface_c* self, int32_t n, double* buffer);

/**
 * @brief comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_comm_expand(fitpack_parametric_surface_c* self, int32_t n, const double* buffer);

/**
 * @brief core_comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_core_comm_pack(const fitpack_parametric_surface_c* self, int32_t n, double* buffer);

/**
 * @brief core_comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_core_comm_expand(fitpack_parametric_surface_c* self, int32_t n, const double* buffer);

/* ===========================================================================================
 * Component Array Accessors (raw pointer + extents)
 * =========================================================================================== */

/**
 * @brief Borrow the storage of allocatable array component 'u'.
 *
 * @param self    Object handle
 * @param ptr     Receives the address of the first element, or NULL when the
 *                object or the component is unallocated (also NULL for an
 *                allocated but zero-sized component).
 * @param extents Receives the component's 1 extent(s); all zero when the
 *                component is unallocated.
 * @note The pointer is borrowed, not owned: it is invalidated by any fitting,
 *       assignment or destroy call on the object. Do not free it, and re-read
 *       it after every call that can reallocate the component.
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_getcomp_u(const fitpack_parametric_surface_c* self, double** ptr, int64_t extents[1]);

/**
 * @brief Borrow the storage of allocatable array component 'v'.
 *
 * @param self    Object handle
 * @param ptr     Receives the address of the first element, or NULL when the
 *                object or the component is unallocated (also NULL for an
 *                allocated but zero-sized component).
 * @param extents Receives the component's 1 extent(s); all zero when the
 *                component is unallocated.
 * @note The pointer is borrowed, not owned: it is invalidated by any fitting,
 *       assignment or destroy call on the object. Do not free it, and re-read
 *       it after every call that can reallocate the component.
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_getcomp_v(const fitpack_parametric_surface_c* self, double** ptr, int64_t extents[1]);

/**
 * @brief Borrow the storage of allocatable array component 'z'.
 *
 * @param self    Object handle
 * @param ptr     Receives the address of the first element, or NULL when the
 *                object or the component is unallocated (also NULL for an
 *                allocated but zero-sized component).
 * @param extents Receives the component's 3 extent(s); all zero when the
 *                component is unallocated.
 *
 * @note The elements are laid out COLUMN-MAJOR (Fortran order): extents[0] is
 *       the leading (fastest-varying) extent, so element (i, j) — zero-based —
 *       lives at ptr[i + j * extents[0]].
 * @note The pointer is borrowed, not owned: it is invalidated by any fitting,
 *       assignment or destroy call on the object. Do not free it, and re-read
 *       it after every call that can reallocate the component.
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_getcomp_z(const fitpack_parametric_surface_c* self, double** ptr, int64_t extents[3]);

/**
 * @brief Borrow the storage of allocatable array component 't'.
 *
 * @param self    Object handle
 * @param ptr     Receives the address of the first element, or NULL when the
 *                object or the component is unallocated (also NULL for an
 *                allocated but zero-sized component).
 * @param extents Receives the component's 2 extent(s); all zero when the
 *                component is unallocated.
 *
 * @note The elements are laid out COLUMN-MAJOR (Fortran order): extents[0] is
 *       the leading (fastest-varying) extent, so element (i, j) — zero-based —
 *       lives at ptr[i + j * extents[0]].
 * @note The pointer is borrowed, not owned: it is invalidated by any fitting,
 *       assignment or destroy call on the object. Do not free it, and re-read
 *       it after every call that can reallocate the component.
 */
FITPACK_CAPI_EXPORT void fitpack_parametric_surface_c_getcomp_t(const fitpack_parametric_surface_c* self, double** ptr, int64_t extents[2]);

/* ===========================================================================================
 * Scalar Property Accessors
 * =========================================================================================== */

FITPACK_CAPI_EXPORT int32_t* fitpack_parametric_surface_c_ref_idim(const fitpack_parametric_surface_c* self);

FITPACK_CAPI_EXPORT int32_t* fitpack_parametric_surface_c_ref_nmax(const fitpack_parametric_surface_c* self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_PARAMETRIC_SURFACE_C_H_INCLUDED */
