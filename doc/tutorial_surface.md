# Scattered Surface Fitting Tutorial {#tutorial_surface}

This tutorial covers bivariate spline fitting to scattered data with
`fitpack_surface`: constructing a tensor-product B-spline
\f$ z = s(x, y) \f$ from irregularly sampled \f$ (x_i, y_i, z_i) \f$
triples, controlling smoothness, and evaluating derived quantities.

Use `fitpack_surface` whenever your observations are not aligned on a
rectangular grid.  If the data **do** lie on a regular grid, prefer
`fitpack_grid_surface` instead (see @ref tutorial_grid_surface).

## Mathematical Background

Given \f$ m \f$ scattered observations \f$ (x_i, y_i, z_i) \f$ with positive
weights \f$ w_i \f$, FITPACK determines a tensor-product B-spline

\f[
    s(x, y) = \sum_{i=1}^{n_x - k_x - 1} \sum_{j=1}^{n_y - k_y - 1}
              c_{ij} \, B_i^{k_x}(x) \, B_j^{k_y}(y)
\f]

that minimizes the thin-plate smoothing energy subject to a residual
constraint:

\f[
    \sum_{i=1}^{m} \bigl( w_i \, (z_i - s(x_i, y_i)) \bigr)^2 \leq S
\f]

where \f$ S \f$ is the user-supplied smoothing parameter.  Knot positions
and the number of knots \f$ n_x, n_y \f$ are chosen automatically by the
`surfit` algorithm.

The default spline degree is bicubic (\f$ k_x = k_y = 3 \f$).  Degrees 1
through 5 are supported.

@see Dierckx, Ch. 5, &sect;5.3 (pp. 85&ndash;98)

## Basic Example

The code below fits a spline to 200 quasi-randomly distributed samples of
\f$ f(x,y) = \sin(\pi x) \cos(\pi y) \f$ on \f$ [0,1]^2 \f$.

```fortran
program surface_demo
    use fitpack, only: fitpack_surface, FP_REAL, FP_FLAG, &
                       FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half, pi
    implicit none

    integer, parameter :: m = 200
    real(FP_REAL) :: x(m), y(m), z(m)
    type(fitpack_surface) :: surf
    integer(FP_FLAG) :: ierr
    integer :: i

    ! Quasi-random sampling (golden-ratio sequences)
    do i = 1, m
        x(i) = mod(i * 0.6180339887_FP_REAL, one)
        y(i) = mod(i * 0.3247179572_FP_REAL, one)
        z(i) = sin(pi * x(i)) * cos(pi * y(i))
    end do

    ! Fit with smoothing = m (a practical starting value for unit weights)
    ierr = surf%new_fit(x, y, z, smoothing=real(m, FP_REAL))
    if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)

    print '(a,i0,a,i0)', 'Knots: ', surf%knots(1), ' x ', surf%knots(2)
    print '(a,es10.3)',   'Residual: ', surf%mse()

    call surf%destroy()
end program
```

`surf\%new_fit(x, y, z, w, smoothing, order)` stores the data, allocates
workspace, and performs a smoothing fit in a single call.  It returns an
`integer(FP_FLAG)` error code.

## Features

### Point evaluation

Evaluate the fitted surface at one or many scattered points:

```fortran
! Scalar
z_val = surf%eval(x_val, y_val, ierr)

! Array of scattered points
z_arr = surf%eval(x_arr, y_arr, ierr)   ! size(z_arr) == size(x_arr)
```

Both `surf\%eval` variants call the unscattered spline evaluator `bispeu`,
so the input points need not be sorted.

### Grid evaluation

Evaluate on every node of a rectangular output grid:

```fortran
! z_grid(j,i) = s(x_grid(i), y_grid(j))
z_grid = surf%eval_ongrid(x_grid, y_grid, ierr)
```

The result has shape `(size(y_grid), size(x_grid))`.

### Partial derivatives

Compute \f$ \partial^{n_x+n_y} s / \partial x^{n_x} \partial y^{n_y} \f$
at scattered points or on a grid:

```fortran
! ds/dx at a single point
dzdx = surf%dfdx(x0, y0, 1, 0, ierr)

! d2s/dxdy at many scattered points
d2 = surf%dfdx(x_arr, y_arr, 1, 1, ierr)

! ds/dy on a grid
dz_grid = surf%dfdx_ongrid(x_grid, y_grid, 0, 1, ierr)
```

The derivative order arguments `nux`, `nuy` must satisfy
\f$ 0 \leq n_x < k_x \f$ and \f$ 0 \leq n_y < k_y \f$.

### Integration

Integrate the surface over a rectangular sub-domain
\f$ [x_1, x_2] \times [y_1, y_2] \f$:

```fortran
vol = surf%integral([x1, y1], [x2, y2])
```

The two arguments are 2-element arrays `lower = [x_lo, y_lo]` and
`upper = [x_hi, y_hi]`.

### Cross-sections

Extract a one-dimensional profile by fixing one coordinate:

```fortran
type(fitpack_curve) :: profile

! Profile f(y) = s(x0, y)  (along_y = .true.)
profile = surf%cross_section(x0, .true., ierr)

! Profile g(x) = s(x, y0)  (along_y = .false.)
profile = surf%cross_section(y0, .false., ierr)
```

The returned `fitpack_curve` is a full 1D spline that can be evaluated,
differentiated, and integrated independently.

### Derivative spline

Build a new `fitpack_surface` whose evaluation gives a partial derivative
of the original:

```fortran
type(fitpack_surface) :: dsdx_surf

dsdx_surf = surf%derivative_spline(1, 0, ierr)   ! ds/dx
```

The resulting surface has reduced degree
\f$ (k_x - n_x, \; k_y - n_y) \f$.

## Smoothing Sweep

The smoothing parameter \f$ S \f$ controls the trade-off between fidelity
and smoothness.  Smaller \f$ S \f$ forces more knots and a closer fit;
larger \f$ S \f$ yields a smoother approximation.

```fortran
print '(a)', '     S       knots_x  knots_y  residual'
do i = 1, 4
    ierr = surf%new_fit(x, y, z, smoothing=10.0_FP_REAL**(3 - i))
    if (FITPACK_SUCCESS(ierr)) then
        print '(es10.1, 2i9, es12.3)', &
              surf%smoothing, surf%knots(1), surf%knots(2), surf%mse()
    end if
end do
```

Special cases:

- **\f$ S = 0 \f$** &mdash; interpolation (the spline passes through every
  data point). Use `surf\%interpolate()`.
- **Large \f$ S \f$** &mdash; very smooth, few knots, possibly large
  residual.
- **\f$ S = m \f$** (number of data points, unit weights) &mdash; a
  practical starting value.

After any fit, `surf\%mse()` returns the weighted sum of squared residuals
\f$ f_p \f$, and `surf\%knots` returns a 2-element integer array
\f$ [n_x, n_y] \f$.

## Fitting Modes

The same three strategies available for curves work for surfaces:

| Method | Call | Description |
|--------|------|-------------|
| Automatic knots | `surf\%fit(smoothing=s)` | Knots chosen to satisfy the smoothing constraint |
| Interpolation | `surf\%interpolate()` | Passes through all data points (\f$ S = 0 \f$) |
| Least-squares | `surf\%least_squares()` | User-supplied knots, no smoothing constraint |

## Error Handling

All routines return an `integer(FP_FLAG)` error code:

```fortran
use fitpack, only: FITPACK_OK, FITPACK_SUCCESS, FITPACK_MESSAGE

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

## Complete Example

A self-contained program demonstrating all of the features above is
provided in `examples/example_surface.f90`.  Build and run it with:

```
fpm run --example example_surface
```
