/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_fitter_c.h (module fitpack_fitters)
!> @brief Standalone C interface to fitpack_fitter (no fortran-arrays dependency)
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

#ifndef FITPACK_FITTER_C_H_INCLUDED
#define FITPACK_FITTER_C_H_INCLUDED

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
#include "fitpack_fitters_c_types.h"  /* For fitpack_fitter_c, fitpack_fitter_c_null */

#ifdef __cplusplus
extern "C" {
#endif


/* ===========================================================================================
 * Method Wrappers (standalone-compatible only)
 * =========================================================================================== */

/**
 * @brief mse
 * @return Result value
 */
FITPACK_CAPI_EXPORT double fitpack_fitter_c_mse(const fitpack_fitter_c* self);

/**
 * @brief core_comm_size
 * @return Result value
 */
FITPACK_CAPI_EXPORT int32_t fitpack_fitter_c_core_comm_size(const fitpack_fitter_c* self);

/**
 * @brief destroy_base
 */
FITPACK_CAPI_EXPORT void fitpack_fitter_c_destroy_base(fitpack_fitter_c* self);

/**
 * @brief core_comm_pack
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_fitter_c_core_comm_pack(const fitpack_fitter_c* self, int32_t n, double* buffer);

/**
 * @brief core_comm_expand
 * @param buffer 
 */
FITPACK_CAPI_EXPORT void fitpack_fitter_c_core_comm_expand(fitpack_fitter_c* self, int32_t n, const double* buffer);

/* ===========================================================================================
 * Component Array Accessors (raw pointer + extents)
 * =========================================================================================== */

/**
 * @brief Borrow the storage of allocatable array component 'c'.
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
FITPACK_CAPI_EXPORT void fitpack_fitter_c_getcomp_c(const fitpack_fitter_c* self, double** ptr, int64_t extents[1]);

/**
 * @brief Borrow the storage of allocatable array component 'wrk'.
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
FITPACK_CAPI_EXPORT void fitpack_fitter_c_getcomp_wrk(const fitpack_fitter_c* self, double** ptr, int64_t extents[1]);

/**
 * @brief Borrow the storage of allocatable array component 'iwrk'.
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
FITPACK_CAPI_EXPORT void fitpack_fitter_c_getcomp_iwrk(const fitpack_fitter_c* self, int32_t** ptr, int64_t extents[1]);

/* ===========================================================================================
 * Scalar Property Accessors
 * =========================================================================================== */

FITPACK_CAPI_EXPORT int32_t* fitpack_fitter_c_ref_iopt(const fitpack_fitter_c* self);

FITPACK_CAPI_EXPORT double* fitpack_fitter_c_ref_smoothing(const fitpack_fitter_c* self);

FITPACK_CAPI_EXPORT double* fitpack_fitter_c_ref_fp(const fitpack_fitter_c* self);

FITPACK_CAPI_EXPORT int32_t* fitpack_fitter_c_ref_lwrk(const fitpack_fitter_c* self);

FITPACK_CAPI_EXPORT int32_t* fitpack_fitter_c_ref_liwrk(const fitpack_fitter_c* self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_FITTER_C_H_INCLUDED */
