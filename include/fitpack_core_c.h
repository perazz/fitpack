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

#include <stdlib.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Precision and array size
typedef double  FP_REAL;
typedef int32_t FP_SIZE;
typedef int32_t FP_FLAG;
typedef bool    FP_BOOL;

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

// Function interface to define a polar boundary shape: radius as a function of angle
// rad = boundary(theta);
typedef FP_REAL (*FP_POLAR_BC)(FP_REAL); 

// Error codes
             FP_BOOL FITPACK_SUCCESS_c(FP_FLAG ierr);
static const FP_FLAG FITPACK_OK                   = 0 ; // ok for spline, abs(fp-s)/s <= tol=0.001
static const FP_FLAG FITPACK_INTERPOLATING_OK     = -1; // ok for interpolating spline, fp=0
static const FP_FLAG FITPACK_LEASTSQUARES_OK      = -2; // ok for weighted least-squares polynomial of degree k.
static const FP_FLAG FITPACK_INSUFFICIENT_STORAGE = 1;
static const FP_FLAG FITPACK_S_TOO_SMALL          = 2;
static const FP_FLAG FITPACK_MAXIT                = 3;
static const FP_FLAG FITPACK_TOO_MANY_KNOTS       = 4;
static const FP_FLAG FITPACK_OVERLAPPING_KNOTS    = 5;
static const FP_FLAG FITPACK_INVALID_RANGE        = 6;
static const FP_FLAG FITPACK_INPUT_ERROR          = 10;
static const FP_FLAG FITPACK_TEST_ERROR           = 11;
static const FP_FLAG FITPACK_INVALID_CONSTRAINT   = 12;
static const FP_FLAG FITPACK_INSUFFICIENT_KNOTS   = 13;

void fitpack_message_c(FP_FLAG ierr, char* msg);

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
static const FP_REAL pio8    = 0.125*pi;
static const FP_REAL deg2rad = pi/180.0;
static const FP_REAL third   = one/three;
static const FP_REAL smallnum03 = 1.0e-03;
static const FP_REAL smallnum06 = 1.0e-06;
static const FP_REAL smallnum08 = 1.0e-08;
static const FP_REAL smallnum10 = 1.0e-10;


// Core library interface
void curfit_c(FP_SIZE iopt, FP_SIZE m, const FP_REAL* x, const FP_REAL* y, const FP_REAL* w,
              FP_REAL xb, FP_REAL xe, FP_SIZE k, FP_REAL s, FP_SIZE nest,
              FP_SIZE* n, FP_REAL* t, FP_REAL* c, FP_REAL* fp,
              FP_REAL* wrk, FP_SIZE lwrk, FP_SIZE* iwrk, FP_FLAG* ier);


void percur_c(FP_SIZE iopt, FP_SIZE m, const FP_REAL* x, const FP_REAL* y, const FP_REAL* w,
              FP_SIZE k, FP_REAL s, FP_SIZE nest, FP_SIZE* n, FP_REAL* t, FP_REAL* c,
              FP_REAL* fp, FP_REAL* wrk, FP_SIZE lwrk, FP_SIZE* iwrk, FP_FLAG* ier);


void parcur_c(FP_SIZE iopt, FP_SIZE ipar, FP_SIZE idim, FP_SIZE m, FP_REAL* u, FP_SIZE mx,
              const FP_REAL* x,  FP_REAL* w, FP_REAL* ub, FP_REAL* ue, FP_SIZE k,
              FP_REAL* s, FP_SIZE nest, FP_SIZE* n, FP_REAL *t, FP_SIZE nc, FP_REAL* c,
              FP_REAL* fp, FP_REAL *wrk, FP_SIZE lwrk, FP_SIZE *iwrk, FP_FLAG* ier);


void clocur_c(FP_SIZE iopt, FP_SIZE ipar, FP_SIZE idim, FP_SIZE m, FP_REAL* u, FP_SIZE mx,
              const FP_REAL* x, const FP_REAL* w, FP_SIZE k, FP_REAL s, FP_SIZE nest, FP_SIZE *n, FP_REAL *t,
              FP_SIZE nc, FP_REAL* c, FP_REAL* fp, FP_REAL* wrk, FP_SIZE lwrk, FP_SIZE *iwrk,
              FP_FLAG* ier);


void cocosp_c(FP_SIZE m, const FP_REAL* x, const FP_REAL* y, const FP_REAL* w, FP_SIZE n,
              const FP_REAL* t,FP_REAL* e,FP_SIZE maxtr,FP_SIZE maxbin, FP_REAL* c, FP_REAL* sq,
              FP_REAL* sx, FP_BOOL* bind, FP_REAL* wrk, FP_SIZE lwrk, FP_SIZE* iwrk, FP_SIZE kwrk,
              FP_FLAG* ier);


void concon_c(FP_SIZE iopt,FP_SIZE m,const FP_REAL* x,const FP_REAL* y,const FP_REAL* w,FP_REAL* v,FP_REAL s,
              FP_SIZE nest,FP_SIZE maxtr,FP_SIZE maxbin, FP_SIZE* n, FP_REAL* t,FP_REAL* c,FP_REAL* sq,
              FP_REAL* sx,FP_BOOL* bind, FP_REAL* wrk,FP_SIZE lwrk,FP_SIZE* iwrk,FP_SIZE kwrk,FP_FLAG* ier);

void splev_c(const FP_REAL* t, FP_SIZE n, const FP_REAL* c, FP_SIZE k, const FP_REAL* x, 
             FP_REAL* y, FP_SIZE m, FP_FLAG e,FP_FLAG* ier);
              
void splder_c(const FP_REAL* t,FP_SIZE n,const FP_REAL* c,FP_SIZE k,FP_SIZE nu,const FP_REAL* x,
              FP_REAL* y, FP_SIZE m, FP_FLAG e, FP_REAL* wrk, FP_FLAG*ier);

void spalde_c(const FP_REAL* t,FP_SIZE n,const FP_REAL* c,FP_SIZE k1,
              FP_REAL x,FP_REAL* d, FP_FLAG* ier);

void curev_c(FP_SIZE idim,const FP_REAL* t,FP_SIZE n, const FP_REAL* c,FP_SIZE nc, FP_SIZE k,
             const FP_REAL* u,FP_SIZE m, FP_REAL* x,FP_SIZE mx, FP_FLAG* ier);

void cualde_c(FP_SIZE idim, const FP_REAL* t,FP_SIZE n,const FP_REAL* c,FP_SIZE nc,FP_SIZE k1, 
              FP_REAL u, FP_REAL* d, FP_SIZE nd, FP_FLAG* ier);
  
void insert_c(FP_SIZE iopt, FP_REAL* t, FP_SIZE* n, FP_REAL* c, FP_SIZE k, FP_REAL x, 
              FP_SIZE nest, FP_FLAG* ier);  
           
FP_REAL splint_c(const FP_REAL* t, FP_SIZE n, const FP_REAL* c, FP_SIZE k, FP_REAL a, 
              FP_REAL b, FP_REAL* wrk);            
    
void fourco_c(const FP_REAL* t, FP_SIZE n, const FP_REAL* c, const FP_REAL* alfa, FP_SIZE m, 
              FP_REAL* ress, FP_REAL* resc, FP_REAL* wrk1, FP_REAL* wrk2, FP_FLAG* ier);    
  
void sproot_c(const FP_REAL* t, FP_SIZE n, const FP_REAL* c, FP_REAL* zeros, FP_SIZE mest, 
              FP_SIZE* m, FP_FLAG* ier);
    
void surfit_c(FP_SIZE iopt, FP_SIZE m, FP_REAL* x, FP_REAL* y, const FP_REAL* z, const FP_REAL* w, 
              FP_REAL xb, FP_REAL xe, FP_REAL yb, FP_REAL ye, FP_SIZE kx, FP_SIZE ky, FP_REAL s, 
              FP_SIZE nxest, FP_SIZE nyest, FP_SIZE nmax, FP_REAL eps, FP_SIZE* nx, FP_REAL* tx, 
              FP_SIZE* ny, FP_REAL* ty, FP_REAL* c, FP_REAL* fp, FP_REAL* wrk1, FP_SIZE lwrk1, 
              FP_REAL* wrk2, FP_SIZE lwrk2, FP_SIZE* iwrk, FP_SIZE kwrk, FP_FLAG* ier); 

              
void regrid_c(FP_SIZE iopt, FP_SIZE mx, const FP_REAL* x, 
                            FP_SIZE my, const FP_REAL* y, const FP_REAL* z, 
                            FP_REAL xb, FP_REAL xe, FP_REAL yb, FP_REAL ye, 
                            FP_SIZE kx, FP_SIZE ky, FP_REAL s, 
                            FP_SIZE nxest, FP_SIZE nyest, FP_SIZE* nx, FP_REAL* tx, 
                            FP_SIZE* ny, FP_REAL* ty, FP_REAL* c, FP_REAL* fp, 
                            FP_REAL* wrk, FP_SIZE lwrk, FP_SIZE* iwrk, FP_SIZE kwrk, 
                            FP_FLAG* ier);
void polar_c(const FP_SIZE* iopt, FP_SIZE m, const FP_REAL* x, const FP_REAL* y, const FP_REAL* z,
             const FP_REAL* w, FP_POLAR_BC rad, FP_REAL s, FP_SIZE nuest, FP_SIZE nvest, FP_REAL eps,
             FP_SIZE* nu, FP_REAL* tu, FP_SIZE* nv, FP_REAL* tv, FP_REAL* u, FP_REAL* v, FP_REAL* c,
             FP_REAL* fp, FP_REAL* wrk1, FP_SIZE lwrk1, FP_REAL* wrk2, FP_SIZE lwrk2, FP_SIZE* iwrk,
             FP_SIZE kwrk, FP_FLAG* ier);
              
void pogrid_c(const FP_SIZE* iopt, const FP_SIZE* ider, FP_SIZE mu, const FP_REAL* u, FP_SIZE mv, 
              const FP_REAL* v, const FP_REAL* z, FP_REAL z0, FP_REAL r, FP_REAL s, FP_SIZE nuest, 
              FP_SIZE nvest, FP_SIZE* nu, FP_REAL* tu, FP_SIZE* nv, FP_REAL* tv, FP_REAL* c, 
              FP_REAL* fp, FP_REAL* wrk, FP_SIZE lwrk, FP_SIZE* iwrk, FP_SIZE kwrk, FP_FLAG* ier);

void sphere_c(FP_SIZE iopt, FP_SIZE m, const FP_REAL* teta, const FP_REAL* phi, const FP_REAL* r,
              const FP_REAL* w, FP_REAL s, FP_SIZE ntest, FP_SIZE npest, FP_REAL eps, FP_SIZE* nt,
              FP_REAL* tt, FP_SIZE* np, FP_REAL* tp, FP_REAL* c, FP_REAL* fp, FP_REAL* wrk1, 
              FP_SIZE lwrk1, FP_REAL* wrk2, FP_SIZE lwrk2, FP_SIZE* iwrk, FP_SIZE kwrk, FP_FLAG* ier);
#ifdef __cplusplus
}
#endif

#endif // FITPACK_CORE_C_H_INCLUDED
