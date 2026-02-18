# Spherical Spline Fitting Tutorial {#tutorial_sphere}

This tutorial covers bivariate spline fitting on the unit sphere with FITPACK:
constructing a bicubic spline \f$ r = s(\theta, \phi) \f$ from data in spherical
coordinates, handling pole singularities, and evaluating the fitted surface.

## Overview

Many physical problems produce data on the surface of a sphere — geophysical
fields, satellite measurements, radiation patterns, or angular distributions.
FITPACK provides two sphere-specific types, selected by input data structure:

| Type | Input | Core routine | Use case |
|------|-------|-------------|----------|
| `fitpack_sphere` | Scattered \f$ (\theta_i, \phi_i, r_i) \f$ | `sphere` | Irregularly sampled observations |
| `fitpack_grid_sphere` | Grid \f$ z(v_j, u_i) \f$ | `spgrid` | Latitude-longitude grids |

**Always prefer `fitpack_grid_sphere` when data lies on a regular
latitude-longitude grid** — the `spgrid` algorithm exploits the grid structure
and is significantly more efficient.

## Mathematical Background

Points on the unit sphere are parameterized by colatitude
\f$ \theta \in [0, \pi] \f$ (measured from the north pole) and longitude
\f$ \phi \in [0, 2\pi) \f$. The Cartesian embedding is:

\f[
    x = \sin\theta \cos\phi, \quad
    y = \sin\theta \sin\phi, \quad
    z = \cos\theta
\f]

### The Pole Problem

At the north pole (\f$ \theta = 0 \f$) and south pole (\f$ \theta = \pi \f$),
all values of \f$ \phi \f$ correspond to the same geometric point. A naive
tensor-product spline in \f$ (\theta, \phi) \f$ would not enforce this
single-valuedness. FITPACK's sphere routines automatically impose pole
constraints so that the fitted spline is smooth and well-defined at both poles.

### Spherical Harmonics Connection

The real spherical harmonics \f$ Y_l^m(\theta, \phi) \f$ form a natural basis
for functions on the sphere. The first few are:

\f[
    Y_0^0 = 1, \quad
    Y_1^0 = \cos\theta, \quad
    Y_2^2 = \sin^2\!\theta \cos 2\phi
\f]

A smooth test function built from these harmonics provides a useful validation
target, since the spline should reproduce low-order harmonics to high accuracy
with only a moderate number of knots.

## Scattered Sphere Data — `fitpack_sphere`

Use `fitpack_sphere` when data points are irregularly distributed over the
sphere.

### Construction and Fitting

```fortran
use fitpack, only: fitpack_sphere, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half, pi
type(fitpack_sphere) :: sph
integer(FP_FLAG) :: ierr
real(FP_REAL) :: theta(m), phi(m), r(m)

! theta in [0, pi], phi in [0, 2*pi]
ierr = sph%new_fit(theta, phi, r, smoothing=real(m, FP_REAL))
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)
```

The `new_fit` method stores the data, allocates workspace, and performs a
smoothing fit in one call. The smoothing parameter \f$ s \f$ controls the
trade-off between fidelity and smoothness. A practical starting value is
\f$ s = m \f$ (the number of data points).

### Evaluation

The `eval` method returns values on a tensor-product grid:

```fortran
f = sph%eval(t_eval, p_eval, ierr)
! f has shape (size(p_eval), size(t_eval)) — f(j,i) = s(t_eval(i), p_eval(j))

! Single-point variant
val = sph%eval(theta_pt, phi_pt, ierr)
```

After fitting, the knot counts are available from `sph\%knots`:

```fortran
print *, 'Theta knots:', sph%knots(1), '  Phi knots:', sph%knots(2)
```

The residual (weighted sum of squared errors) is returned by `sph\%mse()`.

### Example: Spherical Harmonic Test Function

The following test function combines three harmonics:

\f[
    f(\theta, \phi) = Y_0^0 + Y_1^0 + \tfrac{1}{2}\, Y_2^2
                    = 1 + \cos\theta + \tfrac{1}{2}\sin^2\!\theta\,\cos 2\phi
\f]

```fortran
elemental real(FP_REAL) function test_sphere(theta, phi) result(f)
    real(FP_REAL), intent(in) :: theta, phi
    f = one + cos(theta) + half * sin(theta)**2 * cos(2.0_FP_REAL * phi)
end function

! Generate quasi-random scattered points on the sphere
do i = 1, m
    theta(i) = acos(one - 2.0_FP_REAL * mod(i * 0.618_FP_REAL, one))
    phi(i)   = 2.0_FP_REAL * pi * mod(i * 0.325_FP_REAL, one)
    r(i)     = test_sphere(theta(i), phi(i))
end do

! Fit and then refine with tighter smoothing
ierr = sph%new_fit(theta, phi, r, smoothing=real(m, FP_REAL))
ierr = sph%fit(smoothing=10.0_FP_REAL)
```

Tightening the smoothing parameter increases the number of knots and reduces
the residual, allowing the spline to capture finer angular structure.

## Gridded Sphere Data — `fitpack_grid_sphere`

Use `fitpack_grid_sphere` when data is sampled on a regular
colatitude-longitude grid — the typical case for reanalysis products, climate
model output, or any field stored on a lat-lon mesh.

### Construction and Fitting

```fortran
use fitpack, only: fitpack_grid_sphere, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, pi
type(fitpack_grid_sphere) :: gsph
integer(FP_FLAG) :: ierr
real(FP_REAL) :: u(nu), v(nv), zg(nv, nu)

! u (colatitude) strictly inside (0, pi)
! v (longitude) in [0, 2*pi)
! zg(j, i) = f(u(i), v(j))
ierr = gsph%new_fit(u, v, zg, smoothing=zero)
```

Note the array layout: `zg` has shape `(size(v), size(u))`, with the longitude
index first. The colatitude grid must satisfy \f$ 0 < u_1 < \cdots < u_{n_u} < \pi \f$
(strictly interior — the poles are not grid points). The longitude grid must
satisfy \f$ 0 \leq v_1 < v_2 < \cdots < v_{n_v} < 2\pi \f$.

### Pole Boundary Conditions

Because the grid excludes the poles, FITPACK allows you to supply pole values
and control the smoothness there:

```fortran
call gsph%BC_north_pole(z0=z_north, exact=.true., differentiable=.true.)
call gsph%BC_south_pole(z0=z_south, exact=.true., differentiable=.true.)
ierr = gsph%fit(smoothing=zero)  ! Re-fit to apply new BCs
```

| Argument | Type | Effect |
|----------|------|--------|
| `z0` | `real(FP_REAL)` | Prescribe the function value at the pole |
| `exact` | `logical` | `.true.`: spline passes exactly through `z0`; `.false.`: treated as data |
| `differentiable` | `logical` | `.true.`: enforce \f$ C^1 \f$ continuity at the pole |
| `zero_grad` | `logical` | `.true.`: enforce \f$ \nabla s = 0 \f$ at the pole |

When `z0` is omitted, no function-value constraint is imposed at that pole.

### Example: Function on a Lat-Lon Grid

```fortran
integer, parameter :: nu = 15, nv = 30
real(FP_REAL) :: u(nu), v(nv), zg(nv, nu), du, dv

! Interior colatitude grid, uniform longitude grid
du = pi / (nu + 1)
dv = 2.0_FP_REAL * pi / nv
do i = 1, nu
    u(i) = i * du
end do
do j = 1, nv
    v(j) = (j - 1) * dv
end do

! Fill data on the grid
do i = 1, nu
    do j = 1, nv
        zg(j, i) = test_sphere(u(i), v(j))
    end do
end do

! Fit with pole BCs
ierr = gsph%new_fit(u, v, zg, smoothing=zero)
call gsph%BC_north_pole(z0=test_sphere(zero, zero), exact=.true., differentiable=.true.)
call gsph%BC_south_pole(z0=test_sphere(pi, zero),   exact=.true., differentiable=.true.)
ierr = gsph%fit(smoothing=zero)

! Evaluate on a new grid
f = gsph%eval(u_eval, v_eval, ierr)
! f(j, i) = s(u_eval(i), v_eval(j))
```

## Complete Example

A full working program demonstrating both scattered fitting with progressive
refinement and gridded fitting with pole boundary conditions is provided in
`examples/example_sphere.f90`. Build and run it with:

```
fpm run --example example_sphere
```

@see @ref fitpack_sphere_domains, @ref fitpack_gridded_sphere,
     Dierckx Ch. 11, &sect;11.2 (pp. 263&ndash;269)
