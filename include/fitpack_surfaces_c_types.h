/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_surfaces_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_surfaces
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_SURFACES_C_TYPES_H_INCLUDED
#define FITPACK_SURFACES_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_surface. */
typedef struct fitpack_surface_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_surface_c;

extern const char fitpack_surface_c_typename[];

static const fitpack_surface_c fitpack_surface_c_null = { NULL, false, NULL };

const char* fitpack_surface_c_fortran_type_name(fitpack_surface_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_surface_c"). */
const char* fitpack_surface_c_c_type_name(fitpack_surface_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_SURFACES_C_TYPES_H_INCLUDED */
