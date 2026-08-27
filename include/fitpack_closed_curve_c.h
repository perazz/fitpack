/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fitpack_closed_curve_c.h                                                                          **
 *   ** @brief Standalone C interface to fitpack_closed_curve (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_CLOSED_CURVE_C_H_INCLUDED
#define FITPACK_CLOSED_CURVE_C_H_INCLUDED

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
#include "fitpack_parametric_curves_c_types.h"  /* For fitpack_closed_curve_c, fitpack_closed_curve_c_null */
#include "fitpack_parametric_curve_c.h"

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
 * @brief Allocate new fitpack_closed_curve object
 * @param self Pointer to wrapper (will be initialized)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_allocate(fitpack_closed_curve_c* self, fx_status* status);

/**
 * @brief Deallocate fitpack_closed_curve object
 * @param self Pointer to wrapper (will be nullified)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_destroy(fitpack_closed_curve_c* self, fx_status* status);

/**
 * @brief Copy fitpack_closed_curve object.
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
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_copy(fitpack_closed_curve_c* self, const fitpack_closed_curve_c* other, bool deep_copy, fx_status* status);

/**
 * @brief Shallow copy (pointer semantics — Fortran "associate" construct)
 * @param self Destination wrapper (will point to same object as 'other')
 * @param other Source wrapper (read-only)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_associate(fitpack_closed_curve_c* self, const fitpack_closed_curve_c* other, fx_status* status);

/**
 * @brief Move allocation (transfer ownership)
 * @param to Destination wrapper (receives ownership)
 * @param from Source wrapper (becomes null)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_move_alloc(fitpack_closed_curve_c* to, fitpack_closed_curve_c* from, fx_status* status);

/* ===========================================================================================
 * Method Wrappers (standalone-compatible only)
 * =========================================================================================== */

/**
 * @brief set_default_parameters
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_set_default_parameters(fitpack_closed_curve_c* self);

/**
 * @brief fit
 * @param smoothing 
 * @param order 
 * @param keep_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_closed_curve_c_fit(const fitpack_closed_curve_c* self, double* smoothing, int32_t* order, bool* keep_knots);

/**
 * @brief interpolate
 * @param order 
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_closed_curve_c_interpolate(const fitpack_closed_curve_c* self, int32_t* order, bool* reset_knots);

/**
 * @brief least_squares
 * @param smoothing 
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_closed_curve_c_least_squares(const fitpack_closed_curve_c* self, double* smoothing, bool* reset_knots);

/**
 * @brief comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_closed_curve_c_comm_size(const fitpack_closed_curve_c* self);

/**
 * @brief mse
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_closed_curve_c_mse(const fitpack_closed_curve_c* self);

/**
 * @brief core_comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_closed_curve_c_core_comm_size(const fitpack_closed_curve_c* self);

/**
 * @brief destroy_base
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_destroy_base(fitpack_closed_curve_c* self);

/**
 * @brief new_points
 * @param x column-major (Fortran order) buffer of x_n1 x x_n2 elements; x_n1 is the leading (first) extent
 * @param u length = x_n2
 * @param w length = x_n2
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_new_points(fitpack_closed_curve_c* self, int32_t x_n1, int32_t x_n2, const double* x, double* u, double* w);

/**
 * @brief new_fit
 * @param x column-major (Fortran order) buffer of x_n1 x x_n2 elements; x_n1 is the leading (first) extent
 * @param u length = x_n2
 * @param w length = x_n2
 * @param smoothing 
 * @param order 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_closed_curve_c_new_fit(const fitpack_closed_curve_c* self, int32_t x_n1, int32_t x_n2, const double* x, double* u, double* w, double* smoothing, int32_t* order);

/**
 * @brief eval_one
 * @param u 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of (this % idim) — known to the caller from object state / inputs
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_eval_one(fitpack_closed_curve_c* self, double u, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief eval_many
 * @param u 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of (this % idim, SIZE(u)) — known to the caller from object state / inputs
 * @note result is the rank-2 Fortran result flattened column-major
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_eval_many(fitpack_closed_curve_c* self, int32_t n, const double* u, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief curve_derivative
 * @param u 
 * @param order 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of ((THIS % IDIM)) — known to the caller from object state / inputs
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_curve_derivative(fitpack_closed_curve_c* self, double u, int32_t order, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief curve_derivatives
 * @param u 
 * @param order 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of ((THIS % IDIM, SIZE(U))) — known to the caller from object state / inputs
 * @note result is the rank-2 Fortran result flattened column-major
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_curve_derivatives(fitpack_closed_curve_c* self, int32_t n, const double* u, int32_t order, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief curve_all_derivatives
 * @param u 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of ((THIS % IDIM, 0 : THIS % ORDER)) — known to the caller from object state / inputs
 * @note result is the rank-2 Fortran result flattened column-major
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_curve_all_derivatives(fitpack_closed_curve_c* self, double u, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_comm_pack(fitpack_closed_curve_c* self, int32_t n, double* buffer);

/**
 * @brief comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_comm_expand(fitpack_closed_curve_c* self, int32_t n, const double* buffer);

/**
 * @brief core_comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_core_comm_pack(fitpack_closed_curve_c* self, int32_t n, double* buffer);

/**
 * @brief core_comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_closed_curve_c_core_comm_expand(fitpack_closed_curve_c* self, int32_t n, const double* buffer);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_CLOSED_CURVE_C_H_INCLUDED */
