/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_convex_curves_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_convex_curves
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_CONVEX_CURVES_C_TYPES_H_INCLUDED
#define FITPACK_CONVEX_CURVES_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_convex_curve. */
typedef struct fitpack_convex_curve_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_convex_curve_c;

extern const char fitpack_convex_curve_c_typename[];

static const fitpack_convex_curve_c fitpack_convex_curve_c_null = { NULL, false, NULL };

const char* fitpack_convex_curve_c_fortran_type_name(fitpack_convex_curve_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_convex_curve_c"). */
const char* fitpack_convex_curve_c_c_type_name(fitpack_convex_curve_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_CONVEX_CURVES_C_TYPES_H_INCLUDED */
