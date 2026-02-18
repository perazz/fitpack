# fitpack — Modern Fortran Spline Fitting {#mainpage}

**fitpack** is a Modern Fortran (2008+) library for curve and surface fitting with splines,
based on the FITPACK package by Paul Dierckx.

It provides an object-oriented API with derived types for all major fitting tasks:
smoothing, interpolation, least-squares approximation, and constrained fitting.
All routines support spline degrees 1 through 5 (cubic by default) and return
B-spline representations that can be evaluated, differentiated, integrated, and
root-found efficiently.

## Type Hierarchy

All fitter types extend the abstract base `fitpack_fitter`:

```
fitpack_fitter (abstract base)
├── fitpack_curve                   — 1D spline: y = s(x)
│   ├── fitpack_periodic_curve      — periodic boundary conditions
│   └── fitpack_convex_curve        — convexity-constrained fitting
├── fitpack_parametric_curve        — parametric curves: x_i = s_i(u)
│   ├── fitpack_closed_curve        — closed parametric curves
│   └── fitpack_constrained_curve   — endpoint derivative constraints
├── fitpack_surface                 — 2D spline from scattered data
├── fitpack_grid_surface            — 2D spline from gridded data
├── fitpack_parametric_surface      — parametric surfaces
├── fitpack_polar                   — polar-domain fitting (scattered)
├── fitpack_grid_polar              — polar-domain fitting (gridded)
├── fitpack_sphere                  — spherical-domain fitting (scattered)
└── fitpack_grid_sphere             — spherical-domain fitting (gridded)
```

## Quick Start

### Fit a curve

```fortran
use fitpack, only: fitpack_curve

type(fitpack_curve) :: curve
integer :: ierr

! Create a smoothing spline from data
curve = fitpack_curve(x, y, ierr=ierr)

! Evaluate at new points
y_new = curve%eval(x_new, ierr)

! Compute derivative, integral, zeros
dydx  = curve%dfdx(x_new, ierr)
area  = curve%integral(x(1), x(size(x)), ierr)
roots = curve%zeros(ierr)
```

### Fit a surface

```fortran
use fitpack, only: fitpack_grid_surface

type(fitpack_grid_surface) :: surf
integer :: ierr

! Create a smoothing spline from gridded data
surf = fitpack_grid_surface(x, y, z, ierr=ierr)

! Evaluate on a new grid
z_new = surf%eval(x_new, y_new, ierr)
```

## Building

fitpack uses [fpm](https://fpm.fortran-lang.org/) (Fortran Package Manager):

```bash
fpm build          # build the library
fpm test           # run all tests
```

## Reference

The algorithms are described in:

> P. Dierckx, *Curve and Surface Fitting with Splines*,
> Oxford University Press, 1993.

Each routine's documentation references the relevant book chapter, section,
and equation numbers.
