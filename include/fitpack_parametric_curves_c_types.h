/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_parametric_curves_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_parametric_curves
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

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
