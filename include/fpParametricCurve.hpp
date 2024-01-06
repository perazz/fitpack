#ifndef FPPARAMETRICCURVE_HPP_INCLUDED
#define FPPARAMETRICCURVE_HPP_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpParametricCurve
!> @brief C++ interface to fitpack_parametric_curve
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

// Import Fortran-C interface
#include "fitpack_parametric_curves_c.h"
#include <vector>
using std::vector;

typedef vector<FP_REAL> fpPoint;

// Flatten a 2D vector
vector<FP_REAL> flatten_2d_vector(vector<fpPoint> x2d)
{

    // Get number of points
    FP_SIZE npts = x2d.size();

    // Get max vector size
    FP_SIZE ndim = 0;
    for (FP_SIZE i=0; i<x2d.size(); i++)
    {
        ndim = x2d[i].size()>ndim? x2d[i].size() : ndim;
    }

    // Allocate flattened vector
    FP_SIZE n = npts*ndim;

    vector<FP_REAL> x(n,0.0);

    for (FP_SIZE pt=0; pt<x2d.size(); pt++)
    {
        for (FP_SIZE idim=0; idim<x2d[pt].size(); idim++)
        {
            x[pt*ndim+idim] = x2d[pt][idim];
        }
    }

}

class fpParametricCurve
{
    public:

        // Constructors/destructors
         fpParametricCurve() { fitpack_parametric_curve_c_allocate(&cptr); };
        ~fpParametricCurve() { fitpack_parametric_curve_c_destroy(&cptr); };
        fpParametricCurve(const fpParametricCurve &rhs) { fitpack_parametric_curve_c_copy(&cptr, &rhs.cptr); };
        fpParametricCurve(fitpack_parametric_curve_c & rhs, const bool move_alloc=false)
            { move_alloc ? fitpack_parametric_curve_c_move_alloc(&cptr, &rhs)
                         : fitpack_parametric_curve_c_copy(&cptr, &rhs); };
        void operator= ( const fpParametricCurve &rhs) { fitpack_parametric_curve_c_copy(&cptr, &rhs.cptr); };
        void destroy() { fitpack_parametric_curve_c_destroy(&cptr); };

        // New curve from x,y only
        FP_FLAG new_fit(vector<fpPoint> x, vector<fpPoint> y, FP_REAL smoothing = 1000.0)
        {
            FP_SIZE npts = x.size();
            FP_SIZE ndim = npts>0? x[0].size() : 0;

            vector<FP_REAL> x1d = flatten_2d_vector(x);
            vector<FP_REAL> y1d = flatten_2d_vector(y);

            return fitpack_parametric_curve_c_new_fit(&cptr,ndim,npts,x1d.data(),y1d.data(),nullptr,&smoothing);
        }

        // New curve from x, y and weights w
        FP_FLAG new_fit(vector<fpPoint> x, vector<fpPoint> y, vector<fpPoint> w, FP_REAL smoothing = 1000.0)
        {
            FP_SIZE npts = x.size();
            FP_SIZE ndim = npts>0? x[0].size() : 0;

            vector<FP_REAL> x1d = flatten_2d_vector(x);
            vector<FP_REAL> y1d = flatten_2d_vector(y);
            vector<FP_REAL> w1d = flatten_2d_vector(w);

            return fitpack_parametric_curve_c_new_fit(&cptr,ndim,npts,x1d.data(),y1d.data(),w1d.data(),&smoothing);
        }
//
//        // Fit properties
//        const FP_SIZE degree   () { return fitpack_parametric_curve_c_degree(&cptr); };
//        const FP_REAL smoothing() { return fitpack_parametric_curve_c_smoothing(&cptr); };
//        const FP_REAL mse      () { return fitpack_parametric_curve_c_mse(&cptr); };
//
//        // Get value at x
//        FP_REAL eval(FP_REAL x, FP_SIZE* ierr=nullptr)
//        {
//            return fitpack_parametric_curve_c_eval_one(&cptr,x,ierr);
//        }
//
//        // Get values at a vector of x coordinates
//        vector<FP_REAL> eval(vector<FP_REAL> x, FP_SIZE* ierr=nullptr)
//        {
//           FP_SIZE npts = x.size();
//           vector<FP_REAL> y;
//           y.resize(npts);
//           fitpack_parametric_curve_c_eval_many(&cptr,npts,x.data(),y.data(),ierr);
//           return y;
//        }
//
//        // Get single derivative at x
//        FP_REAL ddx(FP_REAL x, FP_SIZE order, FP_SIZE* ierr=nullptr)
//        {
//            return fitpack_parametric_curve_c_derivative(&cptr,x,order,ierr);
//        }
//
//        // Get all derivatives at x
//        vector<FP_REAL> ddx(FP_REAL x, FP_SIZE* ierr=nullptr)
//        {
//           vector<FP_REAL> deriv;
//           deriv.resize(degree()+1);
//           FP_SIZE ier = fitpack_parametric_curve_c_all_derivatives(&cptr,x,deriv.data());
//           if (ierr) (*ierr)=ier;
//           return deriv;
//        }
//
//        // Get integral in range
//        FP_REAL integral(FP_REAL from, FP_REAL to)
//        {
//           return fitpack_parametric_curve_c_integral(&cptr, from, to);
//        }
//
//        // Get fourier coefficients
//        FP_FLAG fourier(const vector<FP_REAL> &alpha,
//                              vector<FP_REAL> &A, vector<FP_REAL> &B)
//        {
//           FP_FLAG ierr;
//           FP_SIZE nparm = alpha.size();
//           A.resize(nparm);
//           B.resize(nparm);
//           fitpack_parametric_curve_c_fourier(&cptr, nparm, alpha.data(), A.data(), B.data(), &ierr );
//           return ierr;
//        }


    protected:

    private:

        // Opaque C structure
        fitpack_parametric_curve_c cptr = fitpack_parametric_curve_c_null;

};

#endif // FPPARAMETRICCURVE_HPP_INCLUDED

