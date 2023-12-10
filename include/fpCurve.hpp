#ifndef FPCURVE_HPP_INCLUDED
#define FPCURVE_HPP_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpCurve
!> @brief C++ interface to fitpack_curve
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

// Import Fortran-C interface
#include "fitpack_curves_c.h"
#include <vector>

class fpCurve
{
    public:

        // Constructors/destructors
         fpCurve() { fitpack_curve_c_allocate(&cptr); };
        ~fpCurve() { fitpack_curve_c_destroy(&cptr); };
        fpCurve(const fpCurve &rhs) { fitpack_curve_c_copy(&cptr, &rhs.cptr); };
        fpCurve(fitpack_curve_c & rhs, const bool move_alloc=false)
            { move_alloc ? fitpack_curve_c_move_alloc(&cptr, &rhs)
                         : fitpack_curve_c_copy(&cptr, &rhs); };
        void operator= ( const fpCurve &rhs) { fitpack_curve_c_copy(&cptr, &rhs.cptr); };
        void destroy() { fitpack_curve_c_destroy(&cptr); };

        // New curve from x,y only
        FP_FLAG new_fit(std::vector<FP_REAL> x, std::vector<FP_REAL> y)
        {
            FP_SIZE npts = x.size();
            return fitpack_curve_c_new_fit(&cptr,npts,x.data(),y.data(),nullptr,nullptr);
        }

        // Fit properties
        const FP_SIZE degree   () { return fitpack_curve_c_degree(&cptr); };
        const FP_REAL smoothing() { return fitpack_curve_c_smoothing(&cptr); };
        const FP_REAL mse      () { return fitpack_curve_c_mse(&cptr); };

        // Get single derivative at x
        FP_REAL ddx(FP_REAL x, FP_SIZE order, FP_SIZE* ierr=nullptr)
        {
            return fitpack_curve_c_derivative(&cptr,x,order,ierr);
        }

        // Get all derivatives at x
        std::vector<FP_REAL> ddx(FP_REAL x, FP_SIZE* ierr=nullptr)
        {
           std::vector<FP_REAL> deriv;
           deriv.reserve(degree()+1);
           FP_SIZE ier = fitpack_curve_c_all_derivatives(&cptr,x,deriv.data());
           if (ierr) (*ierr)=ier;
           return deriv;
        }

        // Get integral in range
        FP_REAL integrate(FP_REAL from, FP_REAL to)
        {
           return fitpack_curve_c_integral(&cptr, from, to);
        }

        // Get fourier coefficients
        FP_FLAG fourier(const std::vector<FP_REAL> &alpha,
                              std::vector<FP_REAL> &A, std::vector<FP_REAL> &B)
        {
           FP_FLAG ierr;
           FP_SIZE nparm = alpha.size();
           A.reserve(nparm);
           B.reserve(nparm);
           fitpack_curve_c_fourier(&cptr, nparm, alpha.data(), A.data(), B.data(), &ierr );
           return ierr;
        }


    protected:

    private:

        // Opaque C structure
        fitpack_curve_c cptr = fitpack_curve_c_null;

};

#endif // FPCURVE_HPP_INCLUDED

