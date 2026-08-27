/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_polar_domains_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_polar_domains
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_POLAR_DOMAINS_C_TYPES_H_INCLUDED
#define FITPACK_POLAR_DOMAINS_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_polar. */
typedef struct fitpack_polar_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_polar_c;

extern const char fitpack_polar_c_typename[];

static const fitpack_polar_c fitpack_polar_c_null = { NULL, false, NULL };

const char* fitpack_polar_c_fortran_type_name(fitpack_polar_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_polar_c"). */
const char* fitpack_polar_c_c_type_name(fitpack_polar_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_POLAR_DOMAINS_C_TYPES_H_INCLUDED */
