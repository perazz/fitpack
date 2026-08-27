/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_gridded_splines_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_gridded_splines
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_GRIDDED_SPLINES_C_TYPES_H_INCLUDED
#define FITPACK_GRIDDED_SPLINES_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_gridded_spline. */
typedef struct fitpack_gridded_spline_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_gridded_spline_c;

extern const char fitpack_gridded_spline_c_typename[];

static const fitpack_gridded_spline_c fitpack_gridded_spline_c_null = { NULL, false, NULL };

const char* fitpack_gridded_spline_c_fortran_type_name(fitpack_gridded_spline_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_gridded_spline_c"). */
const char* fitpack_gridded_spline_c_c_type_name(fitpack_gridded_spline_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_GRIDDED_SPLINES_C_TYPES_H_INCLUDED */
