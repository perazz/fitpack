! **************************************************************************************************
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
! **************************************************************************************************
module fitpack
   use fitpack_core
   use fitpack_fitters
   use fitpack_curves
   use fitpack_surfaces
   use fitpack_grid_surfaces
   use fitpack_gridded_polar
   use fitpack_gridded_sphere
   use fitpack_polar_domains
   use fitpack_sphere_domains
   use fitpack_convex_curves
   use fitpack_parametric_curves
   use fitpack_parametric_surfaces

   implicit none
   private

   ! Public interface: kind parameters
   public :: FP_REAL
   public :: FP_SIZE
   public :: FP_FLAG
   public :: FP_BOOL

   ! Public interface: abstract base type
   public :: fitpack_fitter

   ! Public interface: types
   public :: fitpack_curve
   public :: fitpack_closed_curve
   public :: fitpack_periodic_curve
   public :: fitpack_parametric_curve
   public :: fitpack_constrained_curve
   public :: fitpack_convex_curve

   public :: fitpack_polar
   public :: fitpack_grid_polar

   public :: fitpack_sphere
   public :: fitpack_grid_sphere

   public :: fitpack_surface
   public :: fitpack_grid_surface
   public :: fitpack_parametric_surface

   ! Public interface: boundary behavior flags
   public :: OUTSIDE_EXTRAPOLATE, OUTSIDE_ZERO, OUTSIDE_NOT_ALLOWED, OUTSIDE_NEAREST_BND

   ! Public interface: fitting state flags
   public :: IOPT_NEW_LEASTSQUARES, IOPT_NEW_SMOOTHING, IOPT_OLD_FIT

   ! Public interface: error flags
   public :: FITPACK_OK, FITPACK_INTERPOLATING_OK, FITPACK_LEASTSQUARES_OK
   public :: FITPACK_INSUFFICIENT_STORAGE, FITPACK_S_TOO_SMALL, FITPACK_MAXIT
   public :: FITPACK_TOO_MANY_KNOTS, FITPACK_OVERLAPPING_KNOTS, FITPACK_INVALID_RANGE
   public :: FITPACK_INPUT_ERROR, FITPACK_TEST_ERROR
   public :: FITPACK_INVALID_CONSTRAINT, FITPACK_INSUFFICIENT_KNOTS
   public :: CONCON_MAXBIN, CONCON_MAXTR, CONCON_QP_FAIL

   ! Public interface: error utility functions
   public :: FITPACK_SUCCESS, FITPACK_MESSAGE, IOPT_MESSAGE

   ! Public interface: named constants
   public :: zero, one, half, pi


end module fitpack
