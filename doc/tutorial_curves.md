# Curve Fitting Tutorial {#tutorial_curves}

This tutorial covers one-dimensional spline fitting with FITPACK: constructing
a B-spline \f$ y = s(x) \f$ from data, controlling smoothness, and evaluating
derived quantities (derivatives, integrals, roots).

## Basic Workflow

```fortran
use fitpack, only: fitpack_curve, FP_REAL, FITPACK_OK

real(FP_REAL) :: x(20), y(20)
type(fitpack_curve) :: curve
integer :: ierr

! 1. Prepare data (must be sorted by x)
x = [...]
y = [...]

! 2. Construct a smoothing spline (automatic knots)
curve = fitpack_curve(x, y, ierr=ierr)

! 3. Evaluate at new points
y_new = curve%eval(x_new, ierr)

! 4. Clean up
call curve%destroy()
```

The constructor calls `curfit` internally, which selects knots automatically
to balance closeness of fit against smoothness.

## Smoothing Parameter

The smoothing parameter \f$ s \f$ controls the trade-off between fidelity and
smoothness. Specifically, the fitted spline minimizes a roughness measure
subject to:

\f[
    \sum_{i=1}^{m} \left( w_i \, (y_i - s(x_i)) \right)^2 \leq s
\f]

- **\f$ s = 0 \f$**: Interpolation — the spline passes through every data point.
  Use `curve%interpolate()` or set `smoothing=0`.
- **Small \f$ s \f$**: Close to the data but may oscillate.
- **Large \f$ s \f$**: Very smooth but may miss data features.
- **Default**: The constructor uses \f$ s = 1000 \f$; call `curve%fit(smoothing=s)`
  to adjust.

A practical starting point is \f$ s = m \f$ (the number of data points) when weights
are unity, since the expected sum of squared residuals for a good fit is
\f$ \mathcal{O}(m) \f$.

@see Dierckx, Ch. 5, &sect;5.2.4 (pp. 78&ndash;84)

## Fitting Modes

FITPACK offers three fitting strategies:

| Method | Call | Description |
|--------|------|-------------|
| Automatic knots | `curve%fit(smoothing=s)` | Knots chosen to satisfy the smoothing constraint |
| Interpolation | `curve%interpolate()` | Passes through all data points (\f$ s = 0 \f$) |
| Least-squares | `curve%least_squares(knots)` | User-supplied knots, no smoothing constraint |

After any fit, `curve%mse()` returns the weighted sum of squared residuals \f$ f_p \f$.

## Spline Degree

The default spline degree is \f$ k = 3 \f$ (cubic). Degrees 1 through 5 are supported.
Set the degree at construction:

```fortran
curve = fitpack_curve(x, y, k=5, ierr=ierr)   ! quintic
```

## Weights

Optional per-point weights emphasize certain data points:

```fortran
curve = fitpack_curve(x, y, w=weights, ierr=ierr)
```

Points with larger weights are fitted more closely.

## Evaluation, Derivatives, and Integration

Once fitted, the spline supports several operations:

```fortran
! Evaluate at one or many points
y_val  = curve%eval(x_val, ierr)

! First derivative dy/dx
dydx   = curve%dfdx(x_val, order=1, ierr=ierr)

! n-th derivative (up to order k)
d2ydx2 = curve%dfdx(x_val, order=2, ierr=ierr)

! Definite integral over [a, b]
area   = curve%integral(a, b, ierr)

! All zeros (roots) in the data range
roots  = curve%zeros(ierr)
```

## Periodic Curves

For data with periodic boundary conditions — e.g. \f$ y(x_1) = y(x_m) \f$ — use
`fitpack_periodic_curve`:

```fortran
use fitpack, only: fitpack_periodic_curve

type(fitpack_periodic_curve) :: pcurve
pcurve = fitpack_periodic_curve(x, y, ierr=ierr)
```

The periodic variant enforces \f$ s^{(j)}(x_1) = s^{(j)}(x_m) \f$ for
\f$ j = 0, 1, \ldots, k-1 \f$.

## Parametric Curves

To fit a curve through points in \f$ \mathbb{R}^d \f$ (e.g. a 2D trajectory),
use `fitpack_parametric_curve`:

```fortran
use fitpack, only: fitpack_parametric_curve

real(FP_REAL) :: pts(2, 50)   ! 2D points
type(fitpack_parametric_curve) :: pcurve

pcurve = fitpack_parametric_curve(pts, ierr=ierr)
```

Each component \f$ s_j(u) \f$ shares a common knot vector. Parameter values
are computed automatically from cumulative chord lengths unless explicitly supplied.

Closed parametric curves use `fitpack_closed_curve`; curves with prescribed
endpoint derivatives use `fitpack_constrained_curve`.

@see @ref fitpack_parametric_curves, Dierckx Ch. 9 (pp. 199&ndash;228)

## Convexity-Constrained Curves

When the spline must be locally convex or concave, use `fitpack_convex_curve`:

```fortran
use fitpack, only: fitpack_convex_curve

type(fitpack_convex_curve) :: ccurve
ccurve = fitpack_convex_curve(x, y, ierr=ierr)

! Set per-point convexity flags: 1 = concave, -1 = convex, 0 = free
call ccurve%set_convexity(v)
call ccurve%fit(smoothing=s, ierr=ierr)
```

@see @ref fitpack_convex_curves, Dierckx Ch. 8, &sect;8.3&ndash;8.4 (pp. 173&ndash;196)

## Error Handling

All routines return an `integer(FP_FLAG)` error code:

```fortran
use fitpack, only: FITPACK_OK, FITPACK_SUCCESS, FITPACK_MESSAGE

if (.not. FITPACK_SUCCESS(ierr)) then
    print *, 'Error: ', FITPACK_MESSAGE(ierr)
end if
```

Common return values:

| Flag | Meaning |
|------|---------|
| `FITPACK_OK` (0) | Success |
| `FITPACK_INTERPOLATING_OK` (-1) | Interpolating spline returned |
| `FITPACK_S_TOO_SMALL` (1) | Smoothing parameter too small for the requested knot budget |
| `FITPACK_MAXIT` (2) | Maximum iterations reached |
| `FITPACK_INPUT_ERROR` (10) | Invalid input data |
