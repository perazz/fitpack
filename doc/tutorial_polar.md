# Polar Domain Fitting Tutorial {#tutorial_polar}

This tutorial covers bivariate spline fitting on disc-shaped domains using
FITPACK's polar coordinate types: `fitpack_polar` for scattered data and
`fitpack_grid_polar` for gridded data.

## Overview

Many physical problems produce data on disc-shaped or circular regions.
FITPACK provides two types that handle the polar coordinate transformation
internally and enforce the continuity constraints required at the origin:

| Type | Input data | Boundary | Core routine |
|------|-----------|----------|-------------|
| `fitpack_polar` | Scattered \f$ (x_i, y_i, z_i) \f$ | General \f$ r(\theta) \f$ | `polar` |
| `fitpack_grid_polar` | Grid \f$ z(v_j, u_i) \f$ | Constant radius \f$ r \f$ | `pogrid` |

**Always prefer `fitpack_grid_polar` when data lies on a polar grid** -- the
`pogrid` algorithm exploits the grid structure and is significantly more
efficient.

## Mathematical Background

### Normalized polar coordinates

Given a disc bounded by \f$ r(\theta) \f$, FITPACK maps Cartesian
coordinates to normalized polar coordinates \f$ (u, v) \f$:

\f[
    x = u \, r(v) \cos v, \quad y = u \, r(v) \sin v, \quad
    0 \leq u \leq 1, \; -\pi \leq v \leq \pi
\f]

The normalized radial coordinate \f$ u = \sqrt{x^2 + y^2} / r(\theta) \f$
lies in \f$ [0, 1] \f$ regardless of boundary shape, and \f$ v = \theta \f$
is the polar angle. A bicubic spline \f$ s(u, v) \f$ is fitted in the
rectangular parameter domain \f$ [0, 1] \times [-\pi, \pi] \f$.

### Continuity at the origin

Because all angles \f$ v \f$ map to the same physical point when
\f$ u = 0 \f$, the spline must satisfy continuity constraints:

| `bc_continuity_origin` | Constraint |
|------------------------|------------|
| 0 | \f$ C^0 \f$: \f$ s(0, v) = \text{const} \f$ |
| 1 | \f$ C^1 \f$: additionally \f$ \partial s / \partial u \f$ consistent across angles |
| 2 | \f$ C^2 \f$: second-order consistency (default for `fitpack_polar`) |

@see Dierckx, Ch. 11, &sect;11.1 (pp. 255&ndash;263)

## Scattered Data: `fitpack_polar`

`fitpack_polar` accepts scattered Cartesian coordinates and a user-supplied
boundary function \f$ r(\theta) \f$:

```fortran
use fitpack, only: fitpack_polar, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, one
type(fitpack_polar) :: pol
integer(FP_FLAG) :: ierr

ierr = pol%new_fit(x, y, z, unit_disk, w, smoothing=real(m, FP_REAL))
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)
```

| Argument | Description |
|----------|-------------|
| `x(:)`, `y(:)` | Cartesian coordinates of the scattered data |
| `z(:)` | Function values at those points |
| `boundary` | Procedure matching `fitpack_polar_boundary` -- returns \f$ r(\theta) \f$ |
| `w(:)` | Optional per-point weights (default: unit weights) |
| `smoothing` | Optional smoothing parameter \f$ s \f$ |

After the initial fit, adjust smoothing without resupplying data:

```fortran
ierr = pol%fit(smoothing=10.0_FP_REAL)
```

Evaluate at Cartesian points -- the \f$ (u, v) \f$ transform is internal:

```fortran
z_val = pol%eval(x_val, y_val, ierr)        ! single point
z_arr = pol%eval(x_arr, y_arr, ierr)        ! array of points
```

### Example: smooth function on a unit disc

```fortran
integer, parameter :: m = 100
real(FP_REAL) :: x(m), y(m), z(m), w(m), r, theta
type(fitpack_polar) :: pol
integer(FP_FLAG) :: ierr
integer :: i

do i = 1, m
    r     = sqrt(mod(i * 0.618_FP_REAL, one))
    theta = 2*pi * mod(i * 0.325_FP_REAL, one) - pi
    x(i)  = r * cos(theta);  y(i) = r * sin(theta)
    z(i)  = (x(i)**2 + y(i)**2) / ((x(i) + y(i))**2 + half)
end do
w = one

ierr = pol%new_fit(x, y, z, unit_disk, w, smoothing=real(m, FP_REAL))
print *, pol%eval(zero, zero, ierr)
call pol%destroy()
```

## Gridded Data: `fitpack_grid_polar`

`fitpack_grid_polar` accepts data on a polar grid with a **constant**
boundary radius:

```fortran
use fitpack, only: fitpack_grid_polar, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, pi
type(fitpack_grid_polar) :: gpol
integer(FP_FLAG) :: ierr

ierr = gpol%new_fit(u, v, r, z, z0=z0)
```

| Argument | Description |
|----------|-------------|
| `u(:)` | Radial grid, \f$ 0 < u_1 < \cdots < u_{n_u} \leq 1 \f$ |
| `v(:)` | Angular grid, \f$ -\pi \leq v_1 < \cdots < v_{n_v} < \pi \f$ |
| `r` | Constant boundary radius |
| `z(nv, nu)` | Function values (angular index first) |
| `z0` | Optional function value at origin \f$ u = 0 \f$ |
| `smoothing` | Optional smoothing parameter \f$ s \f$ |

### Origin boundary conditions

Use `gpol\%set_origin_BC` to control spline behavior at the origin:

```fortran
call gpol%set_origin_BC(z0=z0, exact=.true., differentiable=.true.)
ierr = gpol%fit(smoothing=real(nu * nv, FP_REAL))
```

| Argument | Effect |
|----------|--------|
| `z0` | Prescribe the function value at the origin |
| `exact` | `.true.`: interpolate \f$ z_0 \f$ exactly; `.false.`: approximate |
| `differentiable` | `.true.`: enforce \f$ C^1 \f$ at origin; `.false.`: \f$ C^0 \f$ only |

When `z0` is omitted, no origin constraint is imposed.

Evaluate on \f$ (u, v) \f$ coordinates:

```fortran
z_val  = gpol%eval(u_val, v_val, ierr)       ! single point
z_grid = gpol%eval(u_new, v_new, ierr)       ! grid: shape [size(v), size(u)]
```

### Example: data on a polar grid

```fortran
integer, parameter :: nu = 10, nv = 20
real(FP_REAL) :: u(nu), v(nv), zg(nv, nu), z0
type(fitpack_grid_polar) :: gpol
integer(FP_FLAG) :: ierr
integer :: i, j

do i = 1, nu; u(i) = real(i, FP_REAL) / nu; end do
do j = 1, nv; v(j) = -pi + (j - 1) * 2*pi / nv; end do
do i = 1, nu
    do j = 1, nv
        zg(j, i) = (u(i)*cos(v(j)))**2 + (u(i)*sin(v(j)))**2
    end do
end do
z0 = zero

ierr = gpol%new_fit(u, v, one, zg, z0=z0)
call gpol%set_origin_BC(z0=z0, exact=.true., differentiable=.true.)
ierr = gpol%fit(smoothing=real(nu * nv, FP_REAL))
print *, 'Residual: ', gpol%mse()
call gpol%destroy()
```

## Boundary Function

The boundary function for `fitpack_polar` must match the abstract interface
`fitpack_polar_boundary` -- a `pure` function receiving
\f$ \theta \in [-\pi, \pi] \f$ and returning a positive radius:

```fortran
pure real(FP_REAL) function my_boundary(theta) result(r)
    real(FP_REAL), intent(in) :: theta
    r = ...
end function
```

**Unit disc** (constant radius):

```fortran
pure real(FP_REAL) function unit_disk(theta)
    real(FP_REAL), intent(in) :: theta
    unit_disk = one
end function
```

**Ellipse** with semi-axes \f$ a \f$ and \f$ b \f$:

```fortran
pure real(FP_REAL) function ellipse(theta)
    real(FP_REAL), intent(in) :: theta
    real(FP_REAL), parameter :: a = 2.0_FP_REAL, b = 1.0_FP_REAL
    ellipse = a * b / sqrt((b * cos(theta))**2 + (a * sin(theta))**2)
end function
```

## Complete Example

A fully worked program covering both scattered and gridded polar fitting is
provided in the source tree. Build and run with:

```
fpm run --example example_polar
```

@see `examples/example_polar.f90`

## Error Handling

All routines return an `integer(FP_FLAG)` error code:

```fortran
if (.not. FITPACK_SUCCESS(ierr)) then
    print *, 'Error: ', FITPACK_MESSAGE(ierr)
end if
```

| Flag | Meaning |
|------|---------|
| `FITPACK_OK` (0) | Success |
| `FITPACK_INTERPOLATING_OK` (-1) | Interpolating spline returned |
| `FITPACK_S_TOO_SMALL` (1) | Smoothing too small for the knot budget |
| `FITPACK_MAXIT` (2) | Maximum iterations reached |
| `FITPACK_INPUT_ERROR` (10) | Invalid input data |
