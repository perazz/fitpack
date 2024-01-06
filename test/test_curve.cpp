#include "../include/fitpack_core_c.h"
#include "../include/fpCurve.hpp"
#include "../include/fitpack_periodic_curves_c.h"
#include "../include/fpPeriodicCurve.hpp"
#include "../include/fitpack_parametric_curves_c.h"
#include "../include/fpParametricCurve.hpp"

#include <cmath>
#include <vector>
#include <iostream>
using std::vector;

extern "C" {

// Test curve
FP_BOOL test_cpp_sine_fit()
{
   static const int N = 201;
   FP_FLAG ierr = FITPACK_OK;

   // Create a sine function over 200 points in [0,2*pi]
   vector<FP_REAL> x; x.reserve(N);
   vector<FP_REAL> y; y.reserve(N);
   for (int i=0; i<N; i++)
   {
       x.push_back(pi2*i/(N-1));
       y.push_back(sin(x[i]));
   }

   // Create interpolating curve (smoothing=0)
   fpCurve sine;
   ierr = sine.new_fit(x,y,0.0);
   if (!FITPACK_SUCCESS_c(ierr)) return FP_FALSE;

   // Create 200 points halfway between the range
   vector<FP_REAL> xrand; xrand.reserve(N-1);
   for (int i=0; i<N-1; i++)
   {
       xrand.push_back(half*(x[i]+x[i+1]));
   }

   // Perform test far from the knots
   vector<FP_REAL> dfdx(4,0.0);
   for (int i=0; i<N-1; i++)
   {
       // Evaluate function
       FP_REAL yeval = sine.eval(xrand[i],&ierr);
       if (!FITPACK_SUCCESS_c(ierr)) return FP_FALSE;

       // Get exact function and derivatives here
       dfdx[0] = sin(xrand[i]);
       dfdx[1] = cos(xrand[i]);
       dfdx[2] = -dfdx[0];
       dfdx[3] = -dfdx[1];

       // Evaluate function and 3 derivatives
       for (FP_SIZE order = 0; order<4; order++)
       {
          FP_REAL yprime = sine.ddx(xrand[i],order,&ierr);
          if (!FITPACK_SUCCESS_c(ierr)) return FP_FALSE;

          // Check error
          if (abs(yprime-dfdx[order])>smallnum03*fmax(0.01,abs(dfdx[order]))) return FP_FALSE;

       }

   }

   // Evaluate integral, compare with exact
   FP_REAL fint = sine.integral(pi*third,pi2);
   FP_REAL eint = cos(pi*third)-cos(pi2);
   if (abs(eint-fint)>smallnum08) return FP_FALSE;

   // All checks passed: success!
   return FP_TRUE;

}

// ODE-style reciprocal error weight
FP_REAL rewt(FP_REAL RTOL, FP_REAL ATOL, FP_REAL x)
{
   return one/(RTOL*abs(x)+ATOL);
}

// Periodic curve test function
FP_REAL periodic_fun(FP_REAL x)
{
    return cos(x) + sin(2*x);
}

// Periodic curve analytical integral
FP_REAL intgl(FP_REAL x)
{
    return sin(x) - half*cos(2*x);
}

// Periodic test: fit a cosine function
FP_BOOL test_cpp_periodic_fit()
{
    static const int N = 200;
    static const FP_REAL RTOL = 1.0e-2;
    static const FP_REAL ATOL = 1.0e-4;
    FP_FLAG ierr = FITPACK_OK;

    // Generate a periodic function with 200 points
    vector<FP_REAL> x(N,0.0);
    vector<FP_REAL> y(N,0.0);
    for (int i=0; i<N; i++)
    {
        x[i] = (pi2*i/(N-1));
        y[i] = periodic_fun(x[i]);
    }

    // Create new INTERPOLATING curve
    fpPeriodicCurve curve;
    ierr = curve.new_fit(x,y,0.0);

    // Failed to create
    if (!FITPACK_SUCCESS_c(ierr)) return FP_FALSE;

    // Create 200 points halfway between the range
    vector<FP_REAL> xrand(N,0.0);
    for (int i=0; i<N-1; i++) xrand[i] = half*(x[i]+x[i+1]);

    // Also add the extremes
    xrand[0]   = x[0];
    xrand[N-1] = x[N-1];

    // Perform tests
    vector<FP_REAL> dfdx(4,0.0);
    for (int i=0; i<N; i++)
    {

        FP_REAL yeval = curve.eval(xrand[i],&ierr);
        if (!FITPACK_SUCCESS_c(ierr)) return FP_FALSE;


        // Get analytical function and derivatives
        dfdx[0] =  cos(xrand[i]) +   sin(2*xrand[i]); // function
        dfdx[1] = -sin(xrand[i]) + 2*cos(2*xrand[i]); // derivatives
        dfdx[2] = -cos(xrand[i]) - 4*sin(2*xrand[i]);
        dfdx[3] = +sin(xrand[i]) - 8*cos(2*xrand[i]);

        // Evaluate function error
        if (abs(yeval-dfdx[0])*rewt(RTOL,ATOL,dfdx[0])>one)
        {
            std::cout << "[periodic_fit] sine function error is too large: x="
                      << xrand[i] << " yspline=" << yeval << " analytical=" << dfdx[0] << std::endl;

            return FP_FALSE;
        }

        for (FP_SIZE order = 0; order<4; order++)
        {
          FP_REAL yprime = curve.ddx(xrand[i],order,&ierr);
          if (!FITPACK_SUCCESS_c(ierr)) return FP_FALSE;

          // Check error
          if (abs(yprime-dfdx[order])*rewt(RTOL,ATOL,dfdx[order])>one) return FP_FALSE;

        }

    }

    // Test integral calculation *across the periodic boundary*
    FP_REAL fint = curve.integral(half*pi,3*pi);
    FP_REAL eint = intgl(3*pi)-intgl(half*pi);
    FP_REAL error = abs(eint-fint)*rewt(RTOL,ATOL,eint);
    if (error>one)
    {
        std::cout << "[periodic_fit] integral error is too large: fit="
                  << fint << " analytical=" << eint << " error=" << error << std::endl;

        return FP_FALSE;
    }

    // All checks passed: success!
    return FP_TRUE;

}


// Test parametric curve fit
// Periodic test: fit a cosine function
FP_BOOL test_cpp_parametric_fit()
{

     char msg[256];

     // Normalized coordinates
     static const FP_REAL u[] = {120,128,133,136,138,141,144,146,149,151,154,161,170,180,190,
                                 200,210,220,230,240,250,262,269,273,278,282,287,291,295,299,305,315};

     // Data points
     static const FP_REAL x[][2] = {{-1.5141, 0.5150},{-2.0906, 1.3412},{-1.9253, 2.6094},{-0.8724, 3.2358},
                                    {-0.3074, 2.7401},{-0.5534, 2.7823},{ 0.0192, 3.5932},{ 1.2298, 3.8353},
                                    { 2.5479, 2.5863},{ 2.4710, 1.3105},{ 1.7063, 0.6841},{ 1.1183, 0.2575},
                                    { 0.5534, 0.2460},{ 0.4727, 0.3689},{ 0.3574, 0.2460},{ 0.1998, 0.2998},
                                    { 0.2882, 0.3651},{ 0.2613, 0.3343},{ 0.2652, 0.3881},{ 0.2805, 0.4573},
                                    { 0.4112, 0.5918},{ 0.9377, 0.7110},{ 1.3527, 0.4035},{ 1.5564, 0.0769},
                                    { 1.6141,-0.3920},{ 1.6333,-0.8570},{ 1.1567,-1.3412},{ 0.8109,-1.5641},
                                    { 0.2498,-1.7409},{-0.2306,-1.7178},{-0.7571,-1.2989},{-1.1222,-0.5572}};

     // Fit points into a vector
     vector<fpPoint> xv;
     for (FP_SIZE i=0; i<32; i++) xv.push_back(fpPoint(x[i],x[i]+2));

     vector<FP_REAL> uv(u, u+32);

     // Use flat weights
     vector<FP_REAL> wv(32,1.0);

     // Create curve object
     fpParametricCurve curve;

     // Run tests
     FP_FLAG ierr = FITPACK_OK;

     for (FP_SIZE loop=1; loop<=8; loop++)
     {
         switch (loop)
         {
         case 1: // start with a polynomial curve ( s very large)
            {
                 ierr = curve.new_fit(xv,uv,100.0);
                 break;
            }
         case 2: // Smaller values of s to get a tighter approximation
            {
                ierr = curve.fit(1.0);
                break;
            }
         case 3:
            {
                ierr = curve.fit(0.05);
                break;
            }
         case 4: // Larger values to get a smoother approximateion
            {
                ierr = curve.fit(0.25);
                break;
            }
         case 5: // New fit to get possibly fewer knots
            {
                ierr = curve.new_fit(xv,uv,wv,0.25);
                break;
            }
         case 6: // Let the program determine parameters u[i]
            {
                ierr = curve.new_fit(xv,0.25);
                break;
            }
         case 7: // 5-degree spline approximation
            {
                ierr = curve.new_fit(xv,0.25,5);
                break;
            }
         case 8:
            {
                // Interpolating curve
                ierr = curve.interpolate();
                break;
            }
         default:
            break;
         }

         if (!FITPACK_SUCCESS_c(ierr)) {
            fitpack_message_c(ierr,msg);
            std::cout << "[test_parametric_fit] test " << loop << " failed: " << msg << std::endl;
            break;
         }

          // Evaluate the spline at all nodes
          vector<fpPoint> y = curve.eval(uv,&ierr);

          if (!FITPACK_SUCCESS_c(ierr)) {
            fitpack_message_c(ierr,msg);
            std::cout << "[test_parametric_fit] point evaluation " << loop << " failed: " << msg << std::endl;
            break;
          }

          // Evaluate derivatives at a random point from the initial set
          if (loop<6) {
              vector<fpPoint> dy = curve.ddu_all(uv[12], &ierr);

              if (!FITPACK_SUCCESS_c(ierr)) {
                fitpack_message_c(ierr,msg);
                std::cout << "[test_parametric_fit] derivative evaluation " << loop << " failed: " << msg << std::endl;
                break;
              }
          }

     }

    // All checks passed: success!
    return FITPACK_SUCCESS_c(ierr);

}


}

