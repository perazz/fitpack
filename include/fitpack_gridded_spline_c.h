/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fitpack_gridded_spline_c.h                                                                          **
 *   ** @brief Standalone C interface to fitpack_gridded_spline (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_GRIDDED_SPLINE_C_H_INCLUDED
#define FITPACK_GRIDDED_SPLINE_C_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

/* Minimal fx_status definition for standalone use.
 * Layout contract: this struct is the C half of type(fx_status) in the
 * generated <project>_fx_status module (templates/fortran_fx_status.f90.jinja2),
 * which in turn mirrors fortran-arrays' arrays_c. Field order, widths and
 * FX_LEN_STATUS_MSG must change in all three places at once.
 * ARRAYS_C_H_INCLUDED: a gated build pulls the real arrays_c.h into the
 * translation unit ahead of this header (via fxArrays.hpp), and it defines the
 * same struct and macro behind no FX_STATUS_DEFINED guard — so where both are
 * present the real definition wins and this copy stands down. */
#if !defined(FX_STATUS_DEFINED) && !defined(ARRAYS_C_H_INCLUDED)
#define FX_STATUS_DEFINED
#define FX_LEN_STATUS_MSG 248
typedef struct fx_status { bool ok; int code; char message[FX_LEN_STATUS_MSG]; } fx_status;
#endif

#include "fitpack_capi_export.h"
#include "fitpack_gridded_splines_c_types.h"  /* For fitpack_gridded_spline_c, fitpack_gridded_spline_c_null */
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
 * @brief Allocate new fitpack_gridded_spline object
 * @param self Pointer to wrapper (will be initialized)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_allocate(fitpack_gridded_spline_c* self, fx_status* status);

/**
 * @brief Deallocate fitpack_gridded_spline object
 * @param self Pointer to wrapper (will be nullified)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_destroy(fitpack_gridded_spline_c* self, fx_status* status);

/**
 * @brief Copy fitpack_gridded_spline object.
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
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_copy(fitpack_gridded_spline_c* self, const fitpack_gridded_spline_c* other, bool deep_copy, fx_status* status);

/**
 * @brief Shallow copy (pointer semantics — Fortran "associate" construct)
 * @param self Destination wrapper (will point to same object as 'other')
 * @param other Source wrapper (read-only)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_associate(fitpack_gridded_spline_c* self, const fitpack_gridded_spline_c* other, fx_status* status);

/**
 * @brief Move allocation (transfer ownership)
 * @param to Destination wrapper (receives ownership)
 * @param from Source wrapper (becomes null)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_move_alloc(fitpack_gridded_spline_c* to, fitpack_gridded_spline_c* from, fx_status* status);

/* ===========================================================================================
 * Method Wrappers (standalone-compatible only)
 * =========================================================================================== */

/**
 * @brief fit
 * @param smoothing 
 * @param order 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_gridded_spline_c_fit(const fitpack_gridded_spline_c* self, double* smoothing, int32_t* order);

/**
 * @brief least_squares
 * @param smoothing 
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_gridded_spline_c_least_squares(const fitpack_gridded_spline_c* self, double* smoothing, bool* reset_knots);

/**
 * @brief interpolate
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_gridded_spline_c_interpolate(const fitpack_gridded_spline_c* self);

/**
 * @brief cross_section
 * @param ax 
 * @param u 
 * @param ierr 
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_cross_section(fitpack_gridded_spline_c* self, int32_t ax, double u, int32_t* ierr, fitpack_gridded_spline_c* result);

/**
 * @brief comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_gridded_spline_c_comm_size(const fitpack_gridded_spline_c* self);

/**
 * @brief mse
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_gridded_spline_c_mse(const fitpack_gridded_spline_c* self);

/**
 * @brief core_comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_gridded_spline_c_core_comm_size(const fitpack_gridded_spline_c* self);

/**
 * @brief destroy_base
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_destroy_base(fitpack_gridded_spline_c* self);

/**
 * @brief new_points
 * @param xg column-major (Fortran order) buffer of xg_n1 x xg_n2 elements; xg_n1 is the leading (first) extent
 * @param z flat payload of n_z elements (binds rank-1 to the assumed-rank dummy; length independent of sibling arrays)
 * @param m 
 * @param row_major 
 * @param order 
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_new_points(fitpack_gridded_spline_c* self, int32_t xg_n1, int32_t xg_n2, const double* xg, int32_t n_z, const double* z, int32_t n, int32_t* m, bool* row_major, int32_t* order);

/**
 * @brief new_fit
 * @param xg column-major (Fortran order) buffer of xg_n1 x xg_n2 elements; xg_n1 is the leading (first) extent
 * @param z flat payload of n_z elements (binds rank-1 to the assumed-rank dummy; length independent of sibling arrays)
 * @param m 
 * @param row_major 
 * @param order 
 * @param smoothing 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_gridded_spline_c_new_fit(const fitpack_gridded_spline_c* self, int32_t xg_n1, int32_t xg_n2, const double* xg, int32_t n_z, const double* z, int32_t n, int32_t* m, bool* row_major, int32_t* order, double* smoothing);

/**
 * @brief eval_ongrid
 * @param xg column-major (Fortran order) buffer of xg_n1 x xg_n2 elements; xg_n1 is the leading (first) extent
 * @param m 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of (PRODUCT(m)) — known to the caller from object state / inputs
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_eval_ongrid(fitpack_gridded_spline_c* self, int32_t xg_n1, int32_t xg_n2, const double* xg, int32_t n, const int32_t* m, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief grid_eval_one
 * @param x 
 * @param ierr 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_gridded_spline_c_grid_eval_one(const fitpack_gridded_spline_c* self, int32_t n, const double* x, int32_t* ierr);

/**
 * @brief grid_eval_many
 * @param xp column-major (Fortran order) buffer of xp_n1 x xp_n2 elements; xp_n1 is the leading (first) extent
 * @param ierr 
 * @note result: caller-allocated buffer of max_size elements; n_result receives the actual count (= product of (SIZE(xp, 2)))
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_grid_eval_many(fitpack_gridded_spline_c* self, int32_t xp_n1, int32_t xp_n2, const double* xp, int32_t* ierr, double* result, int32_t* n_result, int32_t max_size);

/**
 * @brief dfdx_ongrid
 * @param xg column-major (Fortran order) buffer of xg_n1 x xg_n2 elements; xg_n1 is the leading (first) extent
 * @param m 
 * @param nu 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of (PRODUCT(m)) — known to the caller from object state / inputs
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_dfdx_ongrid(fitpack_gridded_spline_c* self, int32_t xg_n1, int32_t xg_n2, const double* xg, int32_t n, const int32_t* m, const int32_t* nu, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief grid_derivatives_one
 * @param x 
 * @param nu 
 * @param ierr 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_gridded_spline_c_grid_derivatives_one(const fitpack_gridded_spline_c* self, int32_t n, const double* x, const int32_t* nu, int32_t* ierr);

/**
 * @brief grid_derivatives_many
 * @param xp column-major (Fortran order) buffer of xp_n1 x xp_n2 elements; xp_n1 is the leading (first) extent
 * @param nu 
 * @param ierr 
 * @note result: caller-allocated buffer of max_size elements; n_result receives the actual count (= product of (SIZE(xp, 2)))
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_grid_derivatives_many(fitpack_gridded_spline_c* self, int32_t xp_n1, int32_t xp_n2, const double* xp, int32_t n, const int32_t* nu, int32_t* ierr, double* result, int32_t* n_result, int32_t max_size);

/**
 * @brief integral
 * @param lower 
 * @param upper 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_gridded_spline_c_integral(const fitpack_gridded_spline_c* self, int32_t n, const double* lower, const double* upper);

/**
 * @brief derivative_spline
 * @param nu 
 * @param ierr 
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_derivative_spline(fitpack_gridded_spline_c* self, int32_t n, const int32_t* nu, int32_t* ierr, fitpack_gridded_spline_c* result);

/**
 * @brief comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_comm_pack(fitpack_gridded_spline_c* self, int32_t n, double* buffer);

/**
 * @brief comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_comm_expand(fitpack_gridded_spline_c* self, int32_t n, const double* buffer);

/**
 * @brief core_comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_core_comm_pack(fitpack_gridded_spline_c* self, int32_t n, double* buffer);

/**
 * @brief core_comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_core_comm_expand(fitpack_gridded_spline_c* self, int32_t n, const double* buffer);

/* ===========================================================================================
 * Component Array Accessors (raw pointer + extents)
 * =========================================================================================== */

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
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_getcomp_t(const fitpack_gridded_spline_c* self, double** ptr, int64_t extents[2]);

/**
 * @brief Borrow the storage of allocatable array component 'xg'.
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
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_getcomp_xg(const fitpack_gridded_spline_c* self, double** ptr, int64_t extents[2]);

/**
 * @brief Borrow the storage of allocatable array component 'z'.
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
FITPACK_CAPI_EXPORT void fitpack_gridded_spline_c_getcomp_z(const fitpack_gridded_spline_c* self, double** ptr, int64_t extents[1]);

/* ===========================================================================================
 * Scalar Property Accessors
 * =========================================================================================== */

FITPACK_CAPI_EXPORT int32_t* fitpack_gridded_spline_c_ref_dims(const fitpack_gridded_spline_c* self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_GRIDDED_SPLINE_C_H_INCLUDED */
