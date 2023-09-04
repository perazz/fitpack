/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   Refactored by Federico Perini, 10/6/2022
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************/

#ifndef FITPACK_CURVES_C_H_INCLUDED
#define FITPACK_CURVES_C_H_INCLUDED

#include <stdlib>

#ifdef __cplusplus
extern "C" {
#endif

// Integer/real data types
typedef FP_INT  int;
typedef FP_REAL double;

// Interoperable opaque pointer type
struct fp_curve_c {
  void* cptr;
};

// C API
void    fp_curve_c_allocate       (fp_curve_c &self);
void    fp_curve_c_destroy        (fp_curve_c &self);
void    fp_curve_c_new_points     (fp_curve_c &self, FP_INT npts, FP_REAL* x, FP_REAL* y, FP_REAL* w);
FP_INT  fp_curve_c_new_fit        (fp_curve_c &self, FP_INT npts, FP_REAL* x, FP_REAL* y, FP_REAL* w, FP_REAL* smoothing);
FP_INT  fp_curve_c_fit            (fp_curve_c &self, FP_REAL* smoothing);
FP_INT  fp_curve_c_interpolating  (fp_curve_c &self);
FP_REAL fp_curve_c_eval_one       (fp_curve_c &self, FP_REAL x, FP_INT* ierr);
void    fp_curve_c_eval_many      (fp_curve_c &self, FP_INT npts, FP_REAL* x, FP_REAL* y, FP_INT* ierr);
FP_REAL fp_curve_c_integral       (fp_curve_c &self, FP_REAL from, FP_REAL to);
void    fp_curve_c_fourier        (fp_curve_c &self, FP_INT nparm, FP_REAL* alpha, FP_REAL* A, FP_REAL* B, FP_INT* ierr);
FP_REAL fp_curve_c_derivative     (fp_curve_c &self, FP_REAL x, FP_INT order, FP_INT* ierr);
FP_INT  fp_curve_c_all_derivatives(fp_curve_c &self, FP_REAL x, FP_REAL* ddx);
FP_REAL fp_curve_c_smoothing      (fp_curve_c &self);
FP_REAL fp_curve_c_mse            (fp_curve_c &self);
FP_INT  fp_curve_c_degree         (fp_curve_c &self);


#ifdef __cplusplus
}
#endif

// C++ class
#ifdef __cplusplus

class fpCurve
{
    public:
        fpCurve() { fp_curve_c_allocate(cptr); };
       ~fpCurve() { fp_curve_c_destroy (cptr); };

//       // Fitting
//        fpCurve() { fp_curve_c_new_points(); };
//
//       (fp_curve_c &self, FP_INT npts, FP_REAL* x, FP_REAL* y, FP_REAL* w);

       // Fit properties
       const FP_INT  degree   () { return fp_curve_c_degree(cptr); };
       const FP_REAL smoothing() { return fp_curve_c_smoothing(cptr); };
       const FP_REAL mse      () { return fp_curve_c_mse(cptr); };

       // Opaque C type
       fp_curve_c cptr;

    protected:

    private:

};
#endif

#endif // FITPACK_CURVES_C_H_INCLUDED
