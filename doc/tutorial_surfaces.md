# Surface Fitting Tutorial {#tutorial_surfaces}

This tutorial covers bivariate spline fitting with FITPACK: constructing a
tensor-product B-spline \f$ z = s(x, y) \f$ from data and evaluating derived
quantities (partial derivatives, integrals, cross-sections).

## Scattered vs Gridded Data

FITPACK provides two surface types, selected by input data structure:

| Type | Input | Core routine | Efficiency |
|------|-------|-------------|------------|
| `fitpack_surface` | Scattered \f$ (x_i, y_i, z_i) \f$ | `surfit` | General |
| `fitpack_grid_surface` | Grid \f$ z(y_j, x_i) \f$ | `regrid` | Faster for gridded data |

**Always prefer `fitpack_grid_surface` when data lies on a rectangular grid** —
the `regrid` algorithm exploits the grid structure and is significantly more
efficient than `surfit`.

## Scattered Data

```fortran
use fitpack, only: fitpack_surface

real(FP_REAL) :: x(100), y(100), z(100)
type(fitpack_surface) :: surf
integer :: ierr

surf = fitpack_surface(x, y, z, ierr=ierr)
```

@see Dierckx, Ch. 5, &sect;5.3 (pp. 85&ndash;98)

## Gridded Data

```fortran
use fitpack, only: fitpack_grid_surface

real(FP_REAL) :: xg(20), yg(30), zg(30, 20)
type(fitpack_grid_surface) :: gsurf
integer :: ierr

gsurf = fitpack_grid_surface(xg, yg, zg, ierr=ierr)
```

Note that `zg` has shape `(size(yg), size(xg))` — the first index corresponds to
the \f$ y \f$-direction.

@see Dierckx, Ch. 5, &sect;5.4 (pp. 98&ndash;103)

## Smoothing, Interpolation, Least-Squares

The same three fitting modes are available as for curves:

```fortran
! Automatic knot placement with smoothing
call surf%fit(smoothing=s, ierr=ierr)

! Interpolation (s = 0)
call surf%interpolate(ierr=ierr)

! Least-squares with user-supplied knots
call surf%least_squares(tx, ty, ierr=ierr)
```

## Evaluation

Evaluate the fitted spline at scattered points or on a grid:

```fortran
! Single point
z_val = surf%eval(x_val, y_val, ierr)

! Array of scattered points
z_arr = surf%eval(x_arr, y_arr, ierr)

! Rectangular output grid
z_grid = surf%eval(x_new, y_new, ierr)  ! for fitpack_grid_surface
```

## Partial Derivatives

Compute partial derivatives \f$ \partial^{n_x+n_y} s / \partial x^{n_x} \partial y^{n_y} \f$:

```fortran
! Derivatives on a grid (for fitpack_grid_surface)
dz = gsurf%dfdx(x_eval, y_eval, order=[1, 0], ierr=ierr)   ! ds/dx
dz = gsurf%dfdx(x_eval, y_eval, order=[0, 1], ierr=ierr)   ! ds/dy
dz = gsurf%dfdx(x_eval, y_eval, order=[1, 1], ierr=ierr)   ! d2s/dxdy
```

## Integration

Compute the integral over a rectangular sub-domain:

\f[
    \iint_{[x_1, x_2] \times [y_1, y_2]} s(x, y) \, dx \, dy
\f]

```fortran
area = gsurf%integral(x1, x2, y1, y2, ierr)
```

## Cross-Sections

Extract a 1D profile (cross-section) at a fixed \f$ x \f$ or \f$ y \f$ value:

```fortran
type(fitpack_curve) :: profile

! Profile at fixed y = y0
profile = gsurf%cross_section(y0, axis=2, ierr=ierr)

! Evaluate the 1D profile
z_line = profile%eval(x_eval, ierr)
```

The returned `fitpack_curve` is a full 1D spline that can be evaluated,
differentiated, and integrated independently.

## Derivative Splines

Transform the B-spline coefficients to obtain a new surface whose evaluation
gives the partial derivative:

```fortran
type(fitpack_grid_surface) :: dsdx_surf

dsdx_surf = gsurf%derivative_spline(order=[1, 0], ierr=ierr)
```

This creates a new surface object whose `eval` method returns
\f$ \partial s / \partial x \f$.

## Parametric Surfaces

For surfaces in \f$ \mathbb{R}^d \f$ parameterized by \f$ (u, v) \f$, use
`fitpack_parametric_surface`:

```fortran
use fitpack, only: fitpack_parametric_surface

type(fitpack_parametric_surface) :: psurf

! z has shape [size(v), size(u), idim]
psurf = fitpack_parametric_surface(u, v, z, ierr=ierr)
```

Each dimension supports optional periodicity in \f$ u \f$ and/or \f$ v \f$.

@see Dierckx, Ch. 10, &sect;10.2 (pp. 241&ndash;254); parsur, surev
