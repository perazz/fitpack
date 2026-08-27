/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_curves_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_curves
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_CURVES_C_TYPES_H_INCLUDED
#define FITPACK_CURVES_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_curve. */
typedef struct fitpack_curve_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_curve_c;

extern const char fitpack_curve_c_typename[];

static const fitpack_curve_c fitpack_curve_c_null = { NULL, false, NULL };

const char* fitpack_curve_c_fortran_type_name(fitpack_curve_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_curve_c"). */
const char* fitpack_curve_c_c_type_name(fitpack_curve_c self);

/* Opaque C wrapper for fitpack_periodic_curve. */
typedef struct fitpack_periodic_curve_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_periodic_curve_c;

extern const char fitpack_periodic_curve_c_typename[];

static const fitpack_periodic_curve_c fitpack_periodic_curve_c_null = { NULL, false, NULL };

const char* fitpack_periodic_curve_c_fortran_type_name(fitpack_periodic_curve_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_periodic_curve_c"). */
const char* fitpack_periodic_curve_c_c_type_name(fitpack_periodic_curve_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_CURVES_C_TYPES_H_INCLUDED */
