#ifndef FPCONSTRAINEDCURVETEMPLATE_HPP_INCLUDED
#define FPCONSTRAINEDCURVETEMPLATE_HPP_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpConstrainedCurveTemplate
!> @brief Dimension-templated C++ interface to fitpack_constrained_curve
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
#include <array>
#include <cstddef>
#include <vector>

#include "fitpack_constrained_curves_c.h"
#include "fpParametricCurve.hpp"

template <std::size_t N>
class fpConstrainedCurveTemplate {
 public:
  using Point = std::array<FP_REAL, N>;
  using Points = std::vector<Point>;

  // Constructors/destructors
  fpConstrainedCurveTemplate() { fitpack_constrained_curve_c_allocate(&cptr); };
  ~fpConstrainedCurveTemplate() { fitpack_constrained_curve_c_destroy(&cptr); };
  fpConstrainedCurveTemplate(const fpConstrainedCurveTemplate &rhs) { fitpack_constrained_curve_c_copy(&cptr, &rhs.cptr); };
  fpConstrainedCurveTemplate(fitpack_constrained_curve_c &rhs, const bool move_alloc = false) {
    move_alloc ? fitpack_constrained_curve_c_move_alloc(&cptr, &rhs) : fitpack_constrained_curve_c_copy(&cptr, &rhs);
  };
  void operator=(const fpConstrainedCurveTemplate &rhs) { fitpack_constrained_curve_c_copy(&cptr, &rhs.cptr); };
  void destroy() { fitpack_constrained_curve_c_destroy(&cptr); };

  // New curve from x (guess u with Euclidean distance)

  FP_FLAG new_fit(const Points &x, FP_REAL smoothing = 1000.0, FP_SIZE order = 3) {
    FP_SIZE npts = x.size();
    FP_SIZE ndim = N;
    auto p_x1d = flatten_pointer(x);
    return fitpack_constrained_curve_c_new_fit(&cptr, ndim, npts, p_x1d, nullptr, nullptr, &smoothing, &order);
  }

  // New curve from x,u only
  FP_FLAG new_fit(const Points &x, const std::vector<FP_REAL> &u, FP_REAL smoothing = 1000.0, FP_SIZE order = 3) {
    FP_SIZE npts = x.size();
    FP_SIZE ndim = N;
    auto p_x1d = flatten_pointer(x);
    return fitpack_constrained_curve_c_new_fit(&cptr, ndim, npts, p_x1d, mutable_data(u), nullptr, &smoothing, &order);
  }

  // New curve from x, y and weights w
  FP_FLAG new_fit(const Points &x,
                  const std::vector<FP_REAL> &u,
                  const std::vector<FP_REAL> &w,
                  FP_REAL smoothing = 1000.0,
                  FP_SIZE order = 3) {
    FP_SIZE npts = x.size();
    FP_SIZE ndim = N;
    auto p_x1d = flatten_pointer(x);

    return fitpack_constrained_curve_c_new_fit(&cptr, ndim, npts, p_x1d, mutable_data(u), mutable_data(w), &smoothing, &order);
  }

  // Update fit with new parameters
  FP_FLAG fit(FP_SIZE order) { return fitpack_constrained_curve_c_fit(&cptr, nullptr, &order); }
  FP_FLAG fit(FP_REAL smoothing) { return fitpack_constrained_curve_c_fit(&cptr, &smoothing, nullptr); }
  FP_FLAG fit(FP_REAL smoothing, FP_SIZE order) { return fitpack_constrained_curve_c_fit(&cptr, &smoothing, &order); }

  // Get the interpolating fit
  FP_FLAG interpolate() { return fitpack_constrained_curve_c_interpolating(&cptr); }

  // Fit properties
  const FP_SIZE degree() { return fitpack_constrained_curve_c_degree(&cptr); };
  const FP_REAL smoothing() { return fitpack_constrained_curve_c_smoothing(&cptr); };
  const FP_REAL mse() { return fitpack_constrained_curve_c_mse(&cptr); };
  const FP_SIZE ndim() { return fitpack_constrained_curve_c_idim(&cptr); };
  const FP_REAL ubegin() { return fitpack_constrained_curve_c_ubegin(&cptr); };
  const FP_REAL uend() { return fitpack_constrained_curve_c_uend(&cptr); };

  // Set constraints, begin point only

  FP_FLAG constrain_begin(const Points &ddx_begin) {
    FP_SIZE nbegin = ddx_begin.size();
    FP_SIZE nend = 0;
    auto pbegin_1d = flatten_pointer(ddx_begin);
    return fitpack_constrained_curve_c_set_constraints(&cptr, nbegin, nend, pbegin_1d, nullptr);
  }
  // Set constraints, both endpoints
  FP_FLAG constrain_both(const Points &ddx_begin, const Points &ddx_end) {
    FP_SIZE nbegin = ddx_begin.size();
    FP_SIZE nend = ddx_end.size();
    auto pbegin_1d = flatten_pointer(ddx_begin);
    auto pend_1d = flatten_pointer(ddx_end);
    return fitpack_constrained_curve_c_set_constraints(&cptr, nbegin, nend, pbegin_1d, pend_1d);
  }
  // Set constraints, end point only
  FP_FLAG constrain_end(const Points &ddx_end) {
    FP_SIZE nbegin = 0;
    FP_SIZE nend = ddx_end.size();
    auto pend_1d = flatten_pointer(ddx_end);
    return fitpack_constrained_curve_c_set_constraints(&cptr, nbegin, nend, nullptr, pend_1d);
  }

  // Clean constraints
  void clean_constraints() { fitpack_constrained_curve_c_clean_constraints(&cptr); }

  // Get value at u
  Point eval(FP_REAL u, FP_FLAG *ierr = nullptr) {
    Point y;
    FP_FLAG ierr0 = fitpack_constrained_curve_c_eval_one(&cptr, u, y.data());
    if (ierr) (*ierr) = ierr0;
    return y;
  }

  // Get value at many u
  Points eval(const std::vector<FP_REAL> &u, FP_FLAG *ierr = nullptr) {
    FP_FLAG ierr0 = FITPACK_OK;
    Point y1{};
    Points y(u.size(), y1);

    for (FP_SIZE i = 0; i < static_cast<FP_SIZE>(u.size()); i++) {
      y[i] = eval(u[i], &ierr0);
      if (!FITPACK_SUCCESS_c(ierr0)) break;
    }

    // Return error flag
    if (ierr) (*ierr) = ierr0;
    return y;
  }

  // Get single derivative at u
  Point ddu(FP_REAL u, FP_SIZE order, FP_FLAG *ierr = nullptr) {
    Point dx;
    dx.fill(0.0);
    FP_FLAG ierr0 = fitpack_constrained_curve_c_derivative(&cptr, u, order, dx.data());
    if (ierr) (*ierr) = ierr0;
    return dx;
  }

  // Get all derivatives at u
  Points ddu_all(FP_REAL u, FP_FLAG *ierr = nullptr) {
    Points deriv(degree() + 1);
    FP_FLAG ierr0 = FITPACK_OK;
    for (FP_SIZE order = 0; order < degree() + 1; order++) {
      deriv[order] = ddu(u, order, &ierr0);
      if (!FITPACK_SUCCESS_c(ierr0)) break;
    }
    if (ierr) (*ierr) = ierr0;
    return deriv;
  }

 protected:
 private:
  static FP_REAL *mutable_data(const std::vector<FP_REAL> &values) {
    return values.empty() ? nullptr : const_cast<FP_REAL *>(values.data());
  }

  static FP_REAL *flatten_pointer(const Points &pts) {
    return pts.empty() ? nullptr : const_cast<FP_REAL *>(reinterpret_cast<const FP_REAL *>(pts.front().data()));
  }

  // Opaque C structure
  fitpack_constrained_curve_c cptr = fitpack_constrained_curve_c_null;
};

#endif  // FPCONSTRAINEDCURVETEMPLATE_HPP_INCLUDED
