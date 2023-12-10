#include "../include/fitpack_core_c.h"
#include "../include/fpCurve.hpp"

#include <cmath>
#include <vector>
using std::vector;

extern "C" {

// Test curve
FP_BOOL test_cpp_sine_fit()
{
   static const int N = 201;

   vector<FP_REAL> x; x.reserve(N);
   vector<FP_REAL> y; y.reserve(N);

   for (int i=0; i<N; i++)
   {
       x.push_back(pi2*i/(N-1));
       y.push_back(sin(x[i]));
   }

   // Create interpolating curve
   fpCurve sine;

   sine.new_fit(x,y);

   return FP_TRUE;

}


}

