/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_gridded_polar_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_gridded_polar
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_GRIDDED_POLAR_C_TYPES_H_INCLUDED
#define FITPACK_GRIDDED_POLAR_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_grid_polar. */
typedef struct fitpack_grid_polar_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_grid_polar_c;

extern const char fitpack_grid_polar_c_typename[];

static const fitpack_grid_polar_c fitpack_grid_polar_c_null = { NULL, false, NULL };

const char* fitpack_grid_polar_c_fortran_type_name(fitpack_grid_polar_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_grid_polar_c"). */
const char* fitpack_grid_polar_c_c_type_name(fitpack_grid_polar_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_GRIDDED_POLAR_C_TYPES_H_INCLUDED */
