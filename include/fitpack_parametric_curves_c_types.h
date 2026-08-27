/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_parametric_curves_c_types.h (module fitpack_parametric_curves)
!> @brief C type definitions for fitpack_parametric_curves
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

#ifndef FITPACK_PARAMETRIC_CURVES_C_TYPES_H_INCLUDED
#define FITPACK_PARAMETRIC_CURVES_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_parametric_curve. */
typedef struct fitpack_parametric_curve_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_parametric_curve_c;

extern const char fitpack_parametric_curve_c_typename[];

static const fitpack_parametric_curve_c fitpack_parametric_curve_c_null = { NULL, false, NULL };

const char* fitpack_parametric_curve_c_fortran_type_name(fitpack_parametric_curve_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_parametric_curve_c"). */
const char* fitpack_parametric_curve_c_c_type_name(fitpack_parametric_curve_c self);

/* Opaque C wrapper for fitpack_closed_curve. */
typedef struct fitpack_closed_curve_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_closed_curve_c;

extern const char fitpack_closed_curve_c_typename[];

static const fitpack_closed_curve_c fitpack_closed_curve_c_null = { NULL, false, NULL };

const char* fitpack_closed_curve_c_fortran_type_name(fitpack_closed_curve_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_closed_curve_c"). */
const char* fitpack_closed_curve_c_c_type_name(fitpack_closed_curve_c self);

/* Opaque C wrapper for fitpack_constrained_curve. */
typedef struct fitpack_constrained_curve_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_constrained_curve_c;

extern const char fitpack_constrained_curve_c_typename[];

static const fitpack_constrained_curve_c fitpack_constrained_curve_c_null = { NULL, false, NULL };

const char* fitpack_constrained_curve_c_fortran_type_name(fitpack_constrained_curve_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_constrained_curve_c"). */
const char* fitpack_constrained_curve_c_c_type_name(fitpack_constrained_curve_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_PARAMETRIC_CURVES_C_TYPES_H_INCLUDED */
