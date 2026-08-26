/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_sphere_domains_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_sphere_domains
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_SPHERE_DOMAINS_C_TYPES_H_INCLUDED
#define FITPACK_SPHERE_DOMAINS_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_sphere. */
typedef struct fitpack_sphere_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_sphere_c;

extern const char fitpack_sphere_c_typename[];

static const fitpack_sphere_c fitpack_sphere_c_null = { NULL, false, NULL };

const char* fitpack_sphere_c_fortran_type_name(fitpack_sphere_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_sphere_c"). */
const char* fitpack_sphere_c_c_type_name(fitpack_sphere_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_SPHERE_DOMAINS_C_TYPES_H_INCLUDED */
