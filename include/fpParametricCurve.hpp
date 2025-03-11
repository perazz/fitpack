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

        // New curve from x (guess u with Euclidean distance)
        FP_FLAG new_fit(vector<fpPoint> x, FP_REAL smoothing = 1000.0, FP_SIZE order = 3)
        {
            FP_SIZE npts = x.size();
            FP_SIZE ndim = npts>0? x[0].size() : 0;
            vector<FP_REAL> x1d = flatten_2d_vector(x);

            return fitpack_parametric_curve_c_new_fit(&cptr,ndim,npts,x1d.data(),nullptr,nullptr,&smoothing,&order);
        }

        // New curve from x,u only
        FP_FLAG new_fit(vector<fpPoint> x, vector<FP_REAL> u, FP_REAL smoothing = 1000.0, FP_SIZE order = 3)
        {
            FP_SIZE npts = x.size();
            FP_SIZE ndim = npts>0? x[0].size() : 0;
            vector<FP_REAL> x1d = flatten_2d_vector(x);

            return fitpack_parametric_curve_c_new_fit(&cptr,ndim,npts,x1d.data(),u.data(),nullptr,&smoothing,&order);
        }

        // New curve from x, y and weights w
        FP_FLAG new_fit(vector<fpPoint> x, vector<FP_REAL> u, vector<FP_REAL> w, FP_REAL smoothing = 1000.0, FP_SIZE order = 3)
        {
            FP_SIZE npts = x.size();
            FP_SIZE ndim = npts>0? x[0].size() : 0;
            vector<FP_REAL> x1d = flatten_2d_vector(x);

            return fitpack_parametric_curve_c_new_fit(&cptr,ndim,npts,x1d.data(),u.data(),w.data(),&smoothing,&order);
        }

        // Update fit with new parameters
        FP_FLAG fit(FP_SIZE order)                    { return fitpack_parametric_curve_c_fit(&cptr,nullptr,&order); }
        FP_FLAG fit(FP_REAL smoothing)                { return fitpack_parametric_curve_c_fit(&cptr,&smoothing,nullptr); }
        FP_FLAG fit(FP_REAL smoothing, FP_SIZE order) { return fitpack_parametric_curve_c_fit(&cptr,&smoothing,&order); }

        // Get the interpolating fit
        FP_FLAG interpolate()                         { return fitpack_parametric_curve_c_interpolating(&cptr,nullptr); }
        FP_FLAG interpolate(FP_SIZE order)            { return fitpack_parametric_curve_c_interpolating(&cptr,nullptr); }

        // Fit properties
        const FP_SIZE degree   () { return fitpack_parametric_curve_c_degree(&cptr); };
        const FP_REAL smoothing() { return fitpack_parametric_curve_c_smoothing(&cptr); };
        const FP_REAL mse      () { return fitpack_parametric_curve_c_mse(&cptr); };
        const FP_SIZE ndim     () { return fitpack_parametric_curve_c_idim(&cptr); };
        const FP_REAL ubegin   () const { return fitpack_parametric_curve_c_ubegin(&cptr);} ;
        const FP_REAL uend     () const { return fitpack_parametric_curve_c_uend(&cptr); };
        
        FP_REAL&      ubegin() { return *fitpack_parametric_curve_c_ubegin_ref(&cptr); };
        FP_REAL&      uend()   { return *fitpack_parametric_curve_c_uend_ref(&cptr); };

        

        // Get value at u
        fpPoint eval(FP_REAL u, FP_FLAG* ierr=nullptr)
        {
            fpPoint y(fitpack_parametric_curve_c_idim(&cptr),0.0);
            FP_FLAG ierr0 = fitpack_parametric_curve_c_eval_one(&cptr, u, y.data());
            if (ierr) (*ierr) = ierr0;
            return y;
        }

        // Get value at many u
        vector<fpPoint> eval(vector<FP_REAL> u, FP_FLAG* ierr=nullptr)
        {
            FP_FLAG ierr0 = FITPACK_OK;
            fpPoint y1(fitpack_parametric_curve_c_idim(&cptr),0.0);
            vector<fpPoint> y(u.size(),y1);

            for (FP_SIZE i=0; i<u.size(); i++)
            {
                y[i] = eval(u[i],&ierr0);
                if (!FITPACK_SUCCESS_c(ierr0)) break;
            }

            // Return error flag
            if (ierr) (*ierr) = ierr0;
            return y;

        }

        // Get values at a range of parameter values
        vector<fpPoint> eval(FP_REAL umin, FP_REAL umax, FP_SIZE npts = 100, FP_SIZE* ierr=nullptr)
        {
           vector<FP_REAL> uu(npts);
           
           // Generate linearly spaced u values
           FP_REAL du = (umax - umin) / (npts - 1);
           for (FP_SIZE i = 0; i < npts; ++i) {
               uu[i] = umin + i * du;
           }
           
           return eval(uu,ierr);
           
        }             
        
        // Get single derivative at u
        fpPoint ddu(FP_REAL u, FP_SIZE order, FP_FLAG* ierr=nullptr)
        {
            fpPoint dx(fitpack_parametric_curve_c_idim(&cptr),0.0);
            FP_FLAG ierr0 = fitpack_parametric_curve_c_derivative(&cptr,u,order,dx.data());
            if (ierr) (*ierr) = ierr0;
            return dx;
        }

        // Get all derivatives at u
        vector<fpPoint> ddu_all(FP_REAL u, FP_FLAG* ierr=nullptr)
        {
           vector<fpPoint> deriv(degree()+1);
           FP_FLAG ierr0 = FITPACK_OK;
           for (FP_SIZE order = 0; order<degree()+1; order++) {
                deriv[order] = ddu(u,order,&ierr0);
                if (!FITPACK_SUCCESS_c(ierr0)) break;
           }
           if (ierr) (*ierr)=ierr0;
           return deriv;
        }

    protected:

    private:

        // Opaque C structure
        fitpack_parametric_curve_c cptr = fitpack_parametric_curve_c_null;

};

#endif // FPPARAMETRICCURVE_HPP_INCLUDED

