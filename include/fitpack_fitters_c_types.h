/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fitpack_fitters_c_types.h                                                               **
 *   ** @brief C type definitions for fitpack_fitters
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FITPACK_FITTERS_C_TYPES_H_INCLUDED
#define FITPACK_FITTERS_C_TYPES_H_INCLUDED

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque C wrapper for fitpack_fitter. */
typedef struct fitpack_fitter_c {
    void* cptr;
    bool is_pointer;
    const char* name_cptr;
} fitpack_fitter_c;

extern const char fitpack_fitter_c_typename[];

static const fitpack_fitter_c fitpack_fitter_c_null = { NULL, false, NULL };

const char* fitpack_fitter_c_fortran_type_name(fitpack_fitter_c self);

/* Return the static null-terminated C wrapper type name (e.g. "fitpack_fitter_c"). */
const char* fitpack_fitter_c_c_type_name(fitpack_fitter_c self);

#ifdef __cplusplus
}
#endif

#endif /* FITPACK_FITTERS_C_TYPES_H_INCLUDED */
