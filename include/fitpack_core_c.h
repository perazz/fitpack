#ifndef FITPACK_CORE_C_H_INCLUDED
#define FITPACK_CORE_C_H_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   (C) Federico Perini, 12/09/2023
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************/

#include <stdlib>
#include <stdbool>

#ifdef __cplusplus
extern "C" {
#endif

// Precision and array size
typedef FP_REAL double;
typedef FP_SIZE int32_t;
typedef FP_FLAG int32_t;
typedef FP_BOOL bool;

//Spline behavior for points not in the support
static const FP_FLAG OUTSIDE_EXTRAPOLATE = 0; // extrapolated from the end spans
static const FP_FLAG OUTSIDE_ZERO        = 1; // spline evaluates to zero
static const FP_FLAG OUTSIDE_NOT_ALLOWED = 2; // an error flag is returned
static const FP_FLAG OUTSIDE_NEAREST_BND = 3; // evaluate to value of nearest boundary point

//IOPT defines a curve state
static const FP_FLAG IOPT_NEW_LEASTSQUARES = -1; // Request a new lsq fit
static const FP_FLAG IOPT_NEW_SMOOTHING    =  0; // Request a new smoothing fit
static const FP_FLAG IOPT_OLD_FIT          =  1; // Update an old fit

//Dimensions of last knot addition
static const FP_SIZE MAX_IDIM      = 10; //Max number of dimensions
static const FP_FLAG KNOT_DIM_NONE =  0; // No knots added yet
static const FP_FLAG KNOT_DIM_2    =  1; // Last knot added on 2nd dim (y or v)
static const FP_FLAG KNOT_DIM_1    = -1; // Last knot added on 1st dim (x or u)

//Spline degrees
static const FP_SIZE MAX_ORDER = 19; // Max spline order (for array allocation)
static const FP_SIZE DEGREE_3  =  3;
static const FP_SIZE DEGREE_4  =  4;
static const FP_SIZE DEGREE_5  =  5;

static const FP_FLAG FITPACK_OK                   = 0 ; // ok for spline, abs(fp-s)/s <= tol=0.001
static const FP_FLAG FITPACK_INTERPOLATING_OK     = -1; // ok for interpolating spline, fp=0
static const FP_FLAG FITPACK_LEASTSQUARES_OK      = -2; // ok for weighted least-squares polynomial of degree k.
static const FP_FLAG FITPACK_INSUFFICIENT_STORAGE = 1
static const FP_FLAG FITPACK_S_TOO_SMALL          = 2
static const FP_FLAG FITPACK_MAXIT                = 3
static const FP_FLAG FITPACK_TOO_MANY_KNOTS       = 4
static const FP_FLAG FITPACK_OVERLAPPING_KNOTS    = 5
static const FP_FLAG FITPACK_INVALID_RANGE        = 6
static const FP_FLAG FITPACK_INPUT_ERROR          = 10
static const FP_FLAG FITPACK_TEST_ERROR           = 11
static const FP_FLAG FITPACK_INVALID_CONSTRAINT   = 12
static const FP_FLAG FITPACK_INSUFFICIENT_KNOTS   = 13

//Internal Parameters
static const FP_BOOL FP_TRUE  = true;
static const FP_BOOL FP_FALSE = false;

static const FP_SIZE IZERO  = 0;
static const FP_SIZE IONE   = 1;
static const FP_SIZE ITWO   = 2;
static const FP_SIZE ITHREE = 3;
static const FP_SIZE IFOUR  = 4;
static const FP_SIZE IFIVE  = 5;

static const FP_REAL one     = 1.0;
static const FP_REAL zero    = 0.0;
static const FP_REAL half    = 0.5;
static const FP_REAL onep5   = 1.5;
static const FP_REAL fourth  = 0.25;
static const FP_REAL two     = 2.0;
static const FP_REAL three   = 3.0;
static const FP_REAL four    = 4.0;
static const FP_REAL five    = 5.0;
static const FP_REAL six     = 6.0;
static const FP_REAL ten     = 10.0;
static const FP_REAL pi      = 3.1415926535897932384626433832795028842;
static const FP_REAL pi2     = 2*pi;
static const FP_REAL pi4     = 4*pi;
static const FP_REAL pio2    = half*pi;
static const FP_REAL pio4    = fourth*pi;
static const FP_REAL pio8    = 0.125;*pi;
static const FP_REAL deg2rad = pi/180.0;
static const FP_REAL smallnum03 = 1.0e-03;
static const FP_REAL smallnum06 = 1.0e-06;
static const FP_REAL smallnum08 = 1.0e-08;
static const FP_REAL smallnum10 = 1.0e-10;

// Core library interface
void curfit_c(FP_SIZE iopt, FP_SIZE m, FP_REAL* x, FP_REAL* y, FP_REAL* w,
              FP_REAL xb, FP_REAL xe, FP_SIZE k, FP_REAL s, FP_SIZE nest,
              FP_SIZE* n, FP_REAL* t, FP_REAL* c, FP_REAL* fp,
              FP_REAL* wrk, FP_SIZE lwrk, FP_SIZE* iwrk, FP_FLAG* ier)

#ifdef __cplusplus
}
#endif

#endif // FITPACK_CORE_C_H_INCLUDED