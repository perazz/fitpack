#ifndef FPPERIODICCURVE_HPP_INCLUDED
#define FPPERIODICCURVE_HPP_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpPeriodicCurve
!> @brief C++ interface to fitpack_periodic_curve
!
!   Author: (C) Federico Perini
!> @since     12/30/2023
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************/

// Import Fortran-C interface
#include "fitpack_periodic_curves_c.h"
#include <vector>
using std::vector;

class fpPeriodicCurve
{
    public:

        // Constructors/destructors
         fpPeriodicCurve() { fitpack_periodic_curve_c_allocate(&cptr); };
        ~fpPeriodicCurve() { fitpack_periodic_curve_c_destroy(&cptr); };
        fpPeriodicCurve(const fpPeriodicCurve &rhs) { fitpack_periodic_curve_c_copy(&cptr, &rhs.cptr); };
        fpPeriodicCurve(fitpack_periodic_curve_c & rhs, const bool move_alloc=false)
            { move_alloc ? fitpack_periodic_curve_c_move_alloc(&cptr, &rhs)
                         : fitpack_periodic_curve_c_copy(&cptr, &rhs); };
        void operator= ( const fpPeriodicCurve &rhs) { fitpack_periodic_curve_c_copy(&cptr, &rhs.cptr); };
        void destroy() { fitpack_periodic_curve_c_destroy(&cptr); };

        // New curve from x,y only
        FP_FLAG new_fit(vector<FP_REAL> x, vector<FP_REAL> y, FP_REAL smoothing = 1000.0)
        {
            FP_SIZE npts = x.size();
            return fitpack_periodic_curve_c_new_fit(&cptr,npts,x.data(),y.data(),nullptr,&smoothing);
        }

        // New curve from x, y and weights w
        FP_FLAG new_fit(vector<FP_REAL> x, vector<FP_REAL> y, vector<FP_REAL> w, FP_REAL smoothing = 1000.0)
        {
            FP_SIZE npts = x.size();
            return fitpack_periodic_curve_c_new_fit(&cptr,npts,x.data(),y.data(),w.data(),&smoothing);
        }

        // Fit properties
        const FP_SIZE degree   () { return fitpack_periodic_curve_c_degree(&cptr); };
        const FP_REAL smoothing() { return fitpack_periodic_curve_c_smoothing(&cptr); };
        const FP_REAL mse      () { return fitpack_periodic_curve_c_mse(&cptr); };

        // Get value at x
        FP_REAL eval(FP_REAL x, FP_SIZE* ierr=nullptr)
        {
            return fitpack_periodic_curve_c_eval_one(&cptr,x,ierr);
        }

        // Get values at a vector of x coordinates
        vector<FP_REAL> eval(vector<FP_REAL> x, FP_SIZE* ierr=nullptr)
        {
           FP_SIZE npts = x.size();
           vector<FP_REAL> y;
           y.resize(npts);
           fitpack_periodic_curve_c_eval_many(&cptr,npts,x.data(),y.data(),ierr);
        }

        // Get single derivative at x
        FP_REAL ddx(FP_REAL x, FP_SIZE order, FP_SIZE* ierr=nullptr)
        {
            return fitpack_periodic_curve_c_derivative(&cptr,x,order,ierr);
        }

        // Get all derivatives at x
        vector<FP_REAL> ddx(FP_REAL x, FP_SIZE* ierr=nullptr)
        {
           vector<FP_REAL> deriv;
           deriv.resize(degree()+1);
           FP_SIZE ier = fitpack_periodic_curve_c_all_derivatives(&cptr,x,deriv.data());
           if (ierr) (*ierr)=ier;
           return deriv;
        }

        // Get integral in range
        FP_REAL integral(FP_REAL from, FP_REAL to)
        {
           return fitpack_periodic_curve_c_integral(&cptr, from, to);
        }

        // Get fourier coefficients
        FP_FLAG fourier(const vector<FP_REAL> &alpha,
                              vector<FP_REAL> &A, vector<FP_REAL> &B)
        {
           FP_FLAG ierr;
           FP_SIZE nparm = alpha.size();
           A.resize(nparm);
           B.resize(nparm);
           fitpack_periodic_curve_c_fourier(&cptr, nparm, alpha.data(), A.data(), B.data(), &ierr );
           return ierr;
        }


    protected:

    private:

        // Opaque C structure
        fitpack_periodic_curve_c cptr = fitpack_periodic_curve_c_null;

};

#endif // FPPERIODICCURVE_HPP_INCLUDED
