/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_sphere_domains_c_types.h (module fitpack_sphere_domains)
!> @brief C type definitions for fitpack_sphere_domains
!
!   @author Federico Perini
!   @date   2026-08-27
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************/

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
