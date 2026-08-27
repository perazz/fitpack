/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fitpack_grid_sphere_c.h                                                                          **
 *   ** @brief Standalone C interface to fitpack_grid_sphere (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_GRID_SPHERE_C_H_INCLUDED
#define FITPACK_GRID_SPHERE_C_H_INCLUDED

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
#include "fitpack_gridded_sphere_c_types.h"  /* For fitpack_grid_sphere_c, fitpack_grid_sphere_c_null */
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
 * @brief Allocate new fitpack_grid_sphere object
 * @param self Pointer to wrapper (will be initialized)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_allocate(fitpack_grid_sphere_c* self, fx_status* status);

/**
 * @brief Deallocate fitpack_grid_sphere object
 * @param self Pointer to wrapper (will be nullified)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_destroy(fitpack_grid_sphere_c* self, fx_status* status);

/**
 * @brief Copy fitpack_grid_sphere object.
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
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_copy(fitpack_grid_sphere_c* self, const fitpack_grid_sphere_c* other, bool deep_copy, fx_status* status);

/**
 * @brief Shallow copy (pointer semantics — Fortran "associate" construct)
 * @param self Destination wrapper (will point to same object as 'other')
 * @param other Source wrapper (read-only)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_associate(fitpack_grid_sphere_c* self, const fitpack_grid_sphere_c* other, fx_status* status);

/**
 * @brief Move allocation (transfer ownership)
 * @param to Destination wrapper (receives ownership)
 * @param from Source wrapper (becomes null)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_move_alloc(fitpack_grid_sphere_c* to, fitpack_grid_sphere_c* from, fx_status* status);

/* ===========================================================================================
 * Method Wrappers (standalone-compatible only)
 * =========================================================================================== */

/**
 * @brief BC_north_pole
 * @param z0 
 * @param exact 
 * @param differentiable 
 * @param zero_grad 
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_BC_north_pole(fitpack_grid_sphere_c* self, double* z0, bool* exact, bool* differentiable, bool* zero_grad);

/**
 * @brief BC_south_pole
 * @param z0 
 * @param exact 
 * @param differentiable 
 * @param zero_grad 
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_BC_south_pole(fitpack_grid_sphere_c* self, double* z0, bool* exact, bool* differentiable, bool* zero_grad);

/**
 * @brief fit
 * @param smoothing 
 * @param keep_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_grid_sphere_c_fit(const fitpack_grid_sphere_c* self, double* smoothing, bool* keep_knots);

/**
 * @brief least_squares
 * @param smoothing 
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_grid_sphere_c_least_squares(const fitpack_grid_sphere_c* self, double* smoothing, bool* reset_knots);

/**
 * @brief interpolate
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_grid_sphere_c_interpolate(const fitpack_grid_sphere_c* self, bool* reset_knots);

/**
 * @brief gridded_eval_one
 * @param u 
 * @param v 
 * @param ierr 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_grid_sphere_c_gridded_eval_one(const fitpack_grid_sphere_c* self, double u, double v, int32_t* ierr);

/**
 * @brief comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_grid_sphere_c_comm_size(const fitpack_grid_sphere_c* self);

/**
 * @brief write
 * @param filename 
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_write(fitpack_grid_sphere_c* self, const char* filename);

/**
 * @brief mse
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_grid_sphere_c_mse(const fitpack_grid_sphere_c* self);

/**
 * @brief core_comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_grid_sphere_c_core_comm_size(const fitpack_grid_sphere_c* self);

/**
 * @brief destroy_base
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_destroy_base(fitpack_grid_sphere_c* self);

/**
 * @brief new_points
 * @param u 
 * @param v 
 * @param z column-major (Fortran order) buffer of z_n1 x z_n2 elements; z_n1 is the leading (first) extent
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_new_points(fitpack_grid_sphere_c* self, int32_t n, const double* u, const double* v, int32_t z_n1, int32_t z_n2, const double* z);

/**
 * @brief new_fit
 * @param u 
 * @param v 
 * @param z column-major (Fortran order) buffer of z_n1 x z_n2 elements; z_n1 is the leading (first) extent
 * @param smoothing 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_grid_sphere_c_new_fit(const fitpack_grid_sphere_c* self, int32_t n, const double* u, const double* v, int32_t z_n1, int32_t z_n2, const double* z, double* smoothing);

/**
 * @brief gridded_eval_many
 * @param u 
 * @param v 
 * @param ierr 
 * @note result: caller-allocated buffer of max_size elements; n_result receives the actual count (= product of (SIZE(v), SIZE(u)))
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_gridded_eval_many(fitpack_grid_sphere_c* self, int32_t n, const double* u, const double* v, int32_t* ierr, double* result, int32_t* n_result, int32_t max_size);

/**
 * @brief comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_comm_pack(fitpack_grid_sphere_c* self, int32_t n, double* buffer);

/**
 * @brief comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_comm_expand(fitpack_grid_sphere_c* self, int32_t n, const double* buffer);

/**
 * @brief core_comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_core_comm_pack(fitpack_grid_sphere_c* self, int32_t n, double* buffer);

/**
 * @brief core_comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_core_comm_expand(fitpack_grid_sphere_c* self, int32_t n, const double* buffer);

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
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_getcomp_u(const fitpack_grid_sphere_c* self, double** ptr, int64_t extents[1]);

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
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_getcomp_v(const fitpack_grid_sphere_c* self, double** ptr, int64_t extents[1]);

/**
 * @brief Borrow the storage of allocatable array component 'z'.
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
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_getcomp_z(const fitpack_grid_sphere_c* self, double** ptr, int64_t extents[2]);

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
FITPACK_CAPI_EXPORT void fitpack_grid_sphere_c_getcomp_t(const fitpack_grid_sphere_c* self, double** ptr, int64_t extents[2]);

/* ===========================================================================================
 * Scalar Property Accessors
 * =========================================================================================== */

FITPACK_CAPI_EXPORT int32_t* fitpack_grid_sphere_c_ref_nmax(const fitpack_grid_sphere_c* self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_GRID_SPHERE_C_H_INCLUDED */
