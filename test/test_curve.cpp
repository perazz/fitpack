#include "../include/fitpack_core_c.h"
#include "../include/fpCurve.hpp"
#include "../include/fitpack_periodic_curves_c.h"
#include "../include/fpPeriodicCurve.hpp"

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
   vector<FP_REAL> dfdx; dfdx.resize(4);
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



// Periodic test: fit a cosine function
FP_BOOL test_cpp_periodic_fit()
{

//
//
//   integer, parameter     :: N = 200
//   real(FP_REAL), parameter :: RTOL = 1.0e-2_FP_REAL
//   real(FP_REAL), parameter :: ATOL = 1.0e-4_FP_REAL
//   type(fitpack_periodic_curve) :: curve
//   real(FP_REAL) :: x(N),y(N),xrand(N),yeval,yprime,dfdx(0:3),fint,eint
//   integer :: ierr,i,order
//
//   success = .false.
//
//   ! Generate a periodic function with 200 points
//   x = linspace(zero,pi2,N)
//   y = fun(x)
//
//   ! Create INTERPOLATING
//   ierr = curve%new_fit(x,y,smoothing=zero)
//
//   ! Failed to create
//   if (.not.FITPACK_SUCCESS(ierr)) then
//      print *, '[periodic_fit] error generating sine function '
//      return
//   end if
//
//   ! Create 200 points in between the range. Include both extremes
//   xrand(1)     = x(1)
//   xrand(2:N-1) = half*(x(1:N-2)+x(2:N-1))
//   xrand(N)     = x(n)
//
//   do i=1,n
//
//      ! Evaluate curve
//      yeval = curve%eval(xrand(i),ierr);
//      if (.not.FITPACK_SUCCESS(ierr)) then
//         print *, '[periodic_fit] error evaluating sine function '
//         return
//      end if
//
//      ! Get analytical function and derivatives
//      dfdx(0) =  cos(xrand(i)) +   sin(2*xrand(i)) ! the function
//      dfdx(1) = -sin(xrand(i)) + 2*cos(2*xrand(i)) ! derivative
//      dfdx(2) = -cos(xrand(i)) - 4*sin(2*xrand(i))
//      dfdx(3) = +sin(xrand(i)) - 8*cos(2*xrand(i))
//
//      ! error
//      if (abs(yeval-dfdx(0))*rewt(RTOL,ATOL,dfdx(0))>one) then
//         print 1, xrand(i),yeval,dfdx(0)
//         return
//      end if
//
//      ! Evaluate first derivative
//      do order = 1,3
//         yprime = curve%dfdx(xrand(i),order=order,ierr=ierr)
//
//         ! Check evaluation
//         if (.not.FITPACK_SUCCESS(ierr)) then
//           print 2, order,xrand(i),FITPACK_MESSAGE(ierr)
//           return
//         end if
//
//         ! Check error
//         if (abs(yprime-dfdx(order))*rewt(RTOL,ATOL,dfdx(order))>one) then
//            print 3, order,xrand(i),yprime,dfdx(order)
//            return
//         end if
//      end do
//
//   end do
//
//   ! Test integral calculation *across the periodic boundary*
//   fint = curve%integral(from=half*pi,to=3*pi)
//   eint = intgl(3*pi)-intgl(half*pi)
//
//   ! Check error
//   if (abs(eint-fint)*rewt(RTOL,ATOL,eint)>one) then
//      print 4, pi/3,two*pi,fint,eint
//      return
//   end if
//
//   ! All checks passed: success!
     return FP_TRUE;
//
//   1 format('[periodic_fit] sine function error is too large: x=',f6.2,' yspline=',f6.2,' analytical=',f6.2)
//   2 format('[periodic_fit] cannot evaluate ',i0,'-th derivative at ',f6.2,': ',a)
//   3 format('[periodic_fit] ',i0,'-th derivative error is too large: x=',f6.2,' yp(spline)=',f6.2,' analytical=',f6.2)
//   4 format('[periodic_fit] integral error is too large: [a=',f6.2,',b=',f6.2,'] int(spline)=',f6.2,' analytical=',f6.2)
//
//   contains
//
//     elemental real(FP_REAL) function fun(x) result(y)
//        real(FP_REAL), intent(in) :: x
//        y = cos(x) + sin(2*x)
//     end function fun
//
//     elemental real(FP_REAL) function intgl(x) result(y)
//        real(FP_REAL), intent(in) :: x
//        y = sin(x) - half*cos(2*x)
//     end function intgl

}



}

