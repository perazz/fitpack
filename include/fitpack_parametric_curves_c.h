#ifndef FITPACK_PARAMETRIC_CURVES_C_H_INCLUDED
#define FITPACK_PARAMETRIC_CURVES_C_H_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fitpack_parametric_curves_c
!> @brief C interface to fitpack_parametric_curves
!
!   Author: (C) Federico Perini
!> @since     01/06/2024
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
typedef struct fitpack_parametric_curve_c {
	void* cptr;
} fitpack_parametric_curve_c;

// Null-initialized object
static const fitpack_parametric_curve_c fitpack_parametric_curve_c_null = { .cptr = NULL };

// Methods
void    fitpack_parametric_curve_c_allocate              (fitpack_parametric_curve_c *self);
void    fitpack_parametric_curve_c_destroy               (fitpack_parametric_curve_c *self);
void    fitpack_parametric_curve_c_copy                  (fitpack_parametric_curve_c *self, const fitpack_parametric_curve_c* other);
void    fitpack_parametric_curve_c_move_alloc            (fitpack_parametric_curve_c *to, fitpack_parametric_curve_c* from);
void    fitpack_parametric_curve_c_new_points            (fitpack_parametric_curve_c *self, FP_SIZE ndim, FP_SIZE npts, FP_REAL* x, FP_REAL* y, FP_REAL* w);
void    fitpack_parametric_curve_c_set_default_parameters(fitpack_parametric_curve_c *self);
FP_FLAG fitpack_parametric_curve_c_new_fit               (fitpack_parametric_curve_c *self, FP_SIZE ndim, FP_SIZE npts, FP_REAL* x, FP_REAL* u, FP_REAL* w, FP_REAL* smoothing, FP_SIZE* order);
FP_FLAG fitpack_parametric_curve_c_fit                   (fitpack_parametric_curve_c *self, FP_REAL* smoothing, FP_SIZE* order);
FP_FLAG fitpack_parametric_curve_c_interpolating         (fitpack_parametric_curve_c *self, FP_SIZE* order);
FP_FLAG fitpack_parametric_curve_c_eval_one              (fitpack_parametric_curve_c *self, FP_REAL u, FP_REAL* y);
//FP_REAL fitpack_parametric_curve_c_integral       (fitpack_parametric_curve_c *self, FP_REAL from, FP_REAL to);
//void    fitpack_parametric_curve_c_fourier        (fitpack_parametric_curve_c *self, FP_SIZE nparm, const FP_REAL* alpha, FP_REAL* A, FP_REAL* B, FP_FLAG* ierr);
FP_FLAG fitpack_parametric_curve_c_derivative      (fitpack_parametric_curve_c *self, FP_REAL u, FP_SIZE order, FP_REAL* ddx);
FP_REAL fitpack_parametric_curve_c_smoothing       (fitpack_parametric_curve_c *self);
FP_REAL fitpack_parametric_curve_c_mse             (fitpack_parametric_curve_c *self);
FP_SIZE  fitpack_parametric_curve_c_degree         (fitpack_parametric_curve_c *self);
FP_SIZE  fitpack_parametric_curve_c_idim           (fitpack_parametric_curve_c *self);
FP_REAL  fitpack_parametric_curve_c_ubegin         (const fitpack_parametric_curve_c *self);
FP_REAL  fitpack_parametric_curve_c_uend           (const fitpack_parametric_curve_c *self);
FP_REAL* fitpack_parametric_curve_c_ubegin_ref    (fitpack_parametric_curve_c *self);
FP_REAL* fitpack_parametric_curve_c_uend_ref      (fitpack_parametric_curve_c *self);

#ifdef __cplusplus
}
#endif

#endif // FITPACK_PARAMETRIC_CURVES_C_H_INCLUDED

