#ifndef FITPACK_CURVES_C_H_INCLUDED
#define FITPACK_CURVES_C_H_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_curves_c
!> @brief C interface to fitpack_curves
!
!   Author: (C) Federico Perini
!> @since     12/09/2023
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************/

#include <stddef.h>
#include <stdbool.h>

#include "fitpack_core_c.h"

#ifdef __cplusplus
extern "C" {
#endif

// Define C opaque pointer struct
typedef struct fitpack_curve_c {
	void* cptr;
} fitpack_curve_c;

// Null-initialized object
static const fitpack_curve_c fitpack_curve_c_null = { .cptr = NULL };

// Methods
void    fitpack_curve_c_allocate       (fitpack_curve_c *self);
void    fitpack_curve_c_destroy        (fitpack_curve_c *self);
void    fitpack_curve_c_copy           (fitpack_curve_c *self, const fitpack_curve_c* other);
void    fitpack_curve_c_move_alloc     (fitpack_curve_c *to, fitpack_curve_c* from);
void    fitpack_curve_c_new_points     (fitpack_curve_c *self, FP_SIZE npts, FP_REAL* x, FP_REAL* y, FP_REAL* w);
FP_FLAG fitpack_curve_c_new_fit        (fitpack_curve_c *self, FP_SIZE npts, FP_REAL* x, FP_REAL* y, FP_REAL* w, FP_REAL* smoothing);
FP_FLAG fitpack_curve_c_fit            (fitpack_curve_c *self, FP_REAL* smoothing, FP_SIZE* order);
FP_FLAG fitpack_curve_c_interpolating  (fitpack_curve_c *self, FP_SIZE* order);
FP_REAL fitpack_curve_c_eval_one       (fitpack_curve_c *self, FP_REAL x, FP_FLAG* ierr);
void    fitpack_curve_c_eval_many      (fitpack_curve_c *self, FP_SIZE npts, FP_REAL* x, FP_REAL* y, FP_FLAG* ierr);
FP_REAL fitpack_curve_c_integral       (fitpack_curve_c *self, FP_REAL from, FP_REAL to);
void    fitpack_curve_c_fourier        (fitpack_curve_c *self, FP_SIZE nparm, const FP_REAL* alpha, FP_REAL* A, FP_REAL* B, FP_FLAG* ierr);
FP_REAL fitpack_curve_c_derivative     (fitpack_curve_c *self, FP_REAL x, FP_SIZE order, FP_SIZE* ierr);
FP_FLAG  fitpack_curve_c_all_derivatives(fitpack_curve_c *self, FP_REAL x, FP_REAL* ddx);
FP_REAL fitpack_curve_c_smoothing      (fitpack_curve_c *self);
FP_REAL fitpack_curve_c_mse            (fitpack_curve_c *self);
FP_SIZE  fitpack_curve_c_degree         (fitpack_curve_c *self);
void    fitpack_curve_c_set_bc         (fitpack_curve_c *self, FP_FLAG bc);
FP_FLAG fitpack_curve_c_get_bc         (fitpack_curve_c *self);

#ifdef __cplusplus
}
#endif

#endif // FITPACK_CURVES_C_H_INCLUDED

