/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_parametric_surfaces_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_parametric_surfaces
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_PARAMETRIC_SURFACES_C_TYPES_H_INCLUDED
#define FITPACK_PARAMETRIC_SURFACES_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_parametric_surface. */
typedef struct fitpack_parametric_surface_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_parametric_surface_c;

extern const char fitpack_parametric_surface_c_typename[];

static const fitpack_parametric_surface_c fitpack_parametric_surface_c_null = { NULL, false, NULL };

const char* fitpack_parametric_surface_c_fortran_type_name(fitpack_parametric_surface_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_parametric_surface_c"). */
const char* fitpack_parametric_surface_c_c_type_name(fitpack_parametric_surface_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_PARAMETRIC_SURFACES_C_TYPES_H_INCLUDED */
