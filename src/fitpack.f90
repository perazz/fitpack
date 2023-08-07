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
   use fitpack_curves
   use fitpack_surfaces
   use fitpack_grid_surfaces
   use fitpack_polar_domains
   use fitpack_sphere_domains
   use fitpack_parametric_curves

   implicit none
   private

   ! Public interface
   public :: RKIND
   public :: fitpack_curve
   public :: fitpack_closed_curve
   public :: fitpack_periodic_curve
   public :: fitpack_parametric_curve
   public :: fitpack_constrained_curve

   public :: fitpack_polar
   public :: fitpack_sphere
   public :: fitpack_surface
   public :: fitpack_grid_surface




end module fitpack
