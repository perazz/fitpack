/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fitpack_periodic_curve_c.h                                                                          **
 *   ** @brief Standalone C interface to fitpack_periodic_curve (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_PERIODIC_CURVE_C_H_INCLUDED
#define FITPACK_PERIODIC_CURVE_C_H_INCLUDED

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
#include "fitpack_curves_c_types.h"  /* For fitpack_periodic_curve_c, fitpack_periodic_curve_c_null */
#include "fitpack_curve_c.h"

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
 * @brief Allocate new fitpack_periodic_curve object
 * @param self Pointer to wrapper (will be initialized)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_allocate(fitpack_periodic_curve_c* self, fx_status* status);

/**
 * @brief Deallocate fitpack_periodic_curve object
 * @param self Pointer to wrapper (will be nullified)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_destroy(fitpack_periodic_curve_c* self, fx_status* status);

/**
 * @brief Copy fitpack_periodic_curve object.
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
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_copy(fitpack_periodic_curve_c* self, const fitpack_periodic_curve_c* other, bool deep_copy, fx_status* status);

/**
 * @brief Shallow copy (pointer semantics — Fortran "associate" construct)
 * @param self Destination wrapper (will point to same object as 'other')
 * @param other Source wrapper (read-only)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_associate(fitpack_periodic_curve_c* self, const fitpack_periodic_curve_c* other, fx_status* status);

/**
 * @brief Move allocation (transfer ownership)
 * @param to Destination wrapper (receives ownership)
 * @param from Source wrapper (becomes null)
 * @param status Optional error status (NULL = error stop on failure)
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_move_alloc(fitpack_periodic_curve_c* to, fitpack_periodic_curve_c* from, fx_status* status);

/* ===========================================================================================
 * Method Wrappers (standalone-compatible only)
 * =========================================================================================== */

/**
 * @brief fit
 * @param smoothing 
 * @param order 
 * @param keep_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_periodic_curve_c_fit(const fitpack_periodic_curve_c* self, double* smoothing, int32_t* order, bool* keep_knots);

/**
 * @brief interpolate
 * @param order 
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_periodic_curve_c_interpolate(const fitpack_periodic_curve_c* self, int32_t* order, bool* reset_knots);

/**
 * @brief least_squares
 * @param smoothing 
 * @param reset_knots 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_periodic_curve_c_least_squares(const fitpack_periodic_curve_c* self, double* smoothing, bool* reset_knots);

/**
 * @brief curve_eval_one
 * @param x 
 * @param ierr 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_periodic_curve_c_curve_eval_one(const fitpack_periodic_curve_c* self, double x, int32_t* ierr);

/**
 * @brief curve_eval_one_noerr
 * @param x 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_periodic_curve_c_curve_eval_one_noerr(const fitpack_periodic_curve_c* self, double x);

/**
 * @brief integral
 * @param from 
 * @param to 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_periodic_curve_c_integral(const fitpack_periodic_curve_c* self, double from, double to);

/**
 * @brief curve_derivative
 * @param x 
 * @param order 
 * @param ierr 
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_periodic_curve_c_curve_derivative(const fitpack_periodic_curve_c* self, double x, int32_t order, int32_t* ierr);

/**
 * @brief curve_insert_knot_one
 * @param x 
 * @param ierr 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_curve_insert_knot_one(fitpack_periodic_curve_c* self, double x, int32_t* ierr);

/**
 * @brief comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_periodic_curve_c_comm_size(const fitpack_periodic_curve_c* self);

/**
 * @brief mse
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_periodic_curve_c_mse(const fitpack_periodic_curve_c* self);

/**
 * @brief core_comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_periodic_curve_c_core_comm_size(const fitpack_periodic_curve_c* self);

/**
 * @brief destroy_base
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_destroy_base(fitpack_periodic_curve_c* self);

/**
 * @brief new_points
 * @param x 
 * @param y 
 * @param w 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_new_points(fitpack_periodic_curve_c* self, int32_t n, const double* x, const double* y, double* w);

/**
 * @brief new_fit
 * @param x 
 * @param y 
 * @param w 
 * @param smoothing 
 * @param order 
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_periodic_curve_c_new_fit(const fitpack_periodic_curve_c* self, int32_t n, const double* x, const double* y, double* w, double* smoothing, int32_t* order);

/**
 * @brief curve_eval_many
 * @param x 
 * @param ierr 
 * @note result: caller-allocated buffer of max_size elements; n_result receives the actual count (= product of (SIZE(x)))
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_curve_eval_many(fitpack_periodic_curve_c* self, int32_t n, const double* x, int32_t* ierr, double* result, int32_t* n_result, int32_t max_size);

/**
 * @brief curve_eval_many_pure
 * @param x 
 * @note result: caller-allocated buffer of max_size elements; n_result receives the actual count (= product of (SIZE(x)))
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_curve_eval_many_pure(fitpack_periodic_curve_c* self, int32_t n, const double* x, double* result, int32_t* n_result, int32_t max_size);

/**
 * @brief fourier_coefficients
 * @param alpha 
 * @param a 
 * @param b 
 * @param ierr 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_fourier_coefficients(fitpack_periodic_curve_c* self, int32_t n, const double* alpha, double* a, double* b, int32_t* ierr);

/**
 * @brief zeros
 * @param ierr 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_zeros(fitpack_periodic_curve_c* self, int32_t* ierr, double* result, int32_t* n_result, int32_t max_size);

/**
 * @brief curve_derivatives
 * @param x 
 * @param order 
 * @param ierr 
 * @note result: caller-allocated buffer of max_size elements; n_result receives the actual count (= product of ((SIZE(X))))
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_curve_derivatives(fitpack_periodic_curve_c* self, int32_t n, const double* x, int32_t order, int32_t* ierr, double* result, int32_t* n_result, int32_t max_size);

/**
 * @brief curve_all_derivatives
 * @param x 
 * @param ierr 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of ((THIS % ORDER + 1)) — known to the caller from object state / inputs
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_curve_all_derivatives(fitpack_periodic_curve_c* self, double x, int32_t* ierr, double* result, int32_t n_result);

/**
 * @brief curve_all_derivatives_pure
 * @param x 
 * @note result: caller-allocated buffer of n_result elements; n_result = total element count of ((THIS % ORDER + 1)) — known to the caller from object state / inputs
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_curve_all_derivatives_pure(fitpack_periodic_curve_c* self, double x, double* result, int32_t n_result);

/**
 * @brief curve_insert_knot_many
 * @param x 
 * @param ierr 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_curve_insert_knot_many(fitpack_periodic_curve_c* self, int32_t n, const double* x, int32_t* ierr);

/**
 * @brief comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_comm_pack(fitpack_periodic_curve_c* self, int32_t n, double* buffer);

/**
 * @brief comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_comm_expand(fitpack_periodic_curve_c* self, int32_t n, const double* buffer);

/**
 * @brief core_comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_core_comm_pack(fitpack_periodic_curve_c* self, int32_t n, double* buffer);

/**
 * @brief core_comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_periodic_curve_c_core_comm_expand(fitpack_periodic_curve_c* self, int32_t n, const double* buffer);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_PERIODIC_CURVE_C_H_INCLUDED */
