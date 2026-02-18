# Convexity-Constrained Curve Fitting {#tutorial_convex_curve}

This tutorial covers shape-preserving spline fitting with `fitpack_convex_curve`,
which enforces local convexity or concavity constraints on a cubic spline.
Use it when the data represents a physical quantity that must be convex or
concave — e.g. moisture content vs depth, dose-response curves, or
thermodynamic isotherms — and an unconstrained spline would introduce
spurious wiggles.

## Mathematical Background

A cubic B-spline \f$ s(x) = \sum_j c_j B_j(x) \f$ is locally convex
(\f$ s''(x) \geq 0 \f$) when the coefficient second differences satisfy

\f[
    c_{j-1} - 2\,c_j + c_{j+1} \geq 0
\f]

and locally concave (\f$ s''(x) \leq 0 \f$) when

\f[
    c_{j-1} - 2\,c_j + c_{j+1} \leq 0
\f]

FITPACK expresses these requirements through per-data-point flags
\f$ v_i \in \{-1, 0, +1\} \f$ with the constraint \f$ s''(x_i)\,v_i \leq 0 \f$:

| Flag | Constraint | Shape |
|------|------------|-------|
| \f$ v_i = +1 \f$ | \f$ s''(x_i) \leq 0 \f$ | Concave (curves downward) |
| \f$ v_i = -1 \f$ | \f$ s''(x_i) \geq 0 \f$ | Convex (curves upward) |
| \f$ v_i = 0 \f$ | none | Unconstrained |

The fitting always uses cubic splines (\f$ k = 3 \f$). The core solver
(`concon`) places knots automatically to satisfy the smoothing constraint
while respecting the requested convexity at every data point.

@see Dierckx, Ch. 8, &sect;8.3&ndash;8.4 (pp. 173&ndash;196)

## Basic Example

Consider moisture-content measurements as a function of depth. The data
rises quickly then levels off — a concave profile.

```fortran
use fitpack, only: fitpack_convex_curve, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, one
implicit none

integer, parameter :: m = 16
type(fitpack_convex_curve) :: curve
integer(FP_FLAG) :: ierr
real(FP_REAL)    :: x(m), y(m), v(m)

! Moisture content vs depth (16 points, concave trend)
x = [0.1_FP_REAL,  0.3_FP_REAL,  0.5_FP_REAL,  0.7_FP_REAL,  0.9_FP_REAL,  &
     1.25_FP_REAL, 1.75_FP_REAL, 2.25_FP_REAL, 2.75_FP_REAL, 3.5_FP_REAL,  &
     4.5_FP_REAL,  5.5_FP_REAL,  6.5_FP_REAL,  7.5_FP_REAL,  8.5_FP_REAL,  9.5_FP_REAL]
y = [0.124_FP_REAL, 0.234_FP_REAL, 0.256_FP_REAL, 0.277_FP_REAL, 0.278_FP_REAL, &
     0.291_FP_REAL, 0.308_FP_REAL, 0.311_FP_REAL, 0.315_FP_REAL, 0.322_FP_REAL, &
     0.317_FP_REAL, 0.326_FP_REAL, 0.323_FP_REAL, 0.321_FP_REAL, 0.322_FP_REAL, &
     0.328_FP_REAL]
```

### Step 1 — Load Data and Set Convexity Flags

```fortran
! Load data points
call curve%new_points(x, y)

! Concave everywhere: v = +1 means s''(x_i) <= 0
v = one
ierr = curve%set_convexity(v)
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)
```

### Step 2 — Fit with Constraints

Call `curve\%fit` with a smoothing parameter. The solver selects knots
automatically while enforcing \f$ s''(x_i) \leq 0 \f$ at every data point:

```fortran
ierr = curve%fit(smoothing=0.04_FP_REAL)
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)
```

### Step 3 — Evaluate and Inspect

```fortran
real(FP_REAL) :: y_eval, d2y
integer :: i

do i = 1, m
    y_eval = curve%eval(x(i), ierr)
    d2y    = curve%dfdx(x(i), order=2, ierr=ierr)
    print '(f8.3, f10.4, es12.3)', x(i), y_eval, d2y   ! d2y <= 0
end do

call curve%destroy()
```

## Comparing Constrained vs Unconstrained

Without constraints, a standard `fitpack_curve` may oscillate through noisy
data, creating regions where the second derivative changes sign. The
constrained fit eliminates these wiggles:

```fortran
use fitpack, only: fitpack_curve, fitpack_convex_curve, FP_REAL, FP_FLAG, one

type(fitpack_curve)        :: free_curve
type(fitpack_convex_curve) :: constrained
integer(FP_FLAG) :: ierr
real(FP_REAL)    :: v(m)

free_curve = fitpack_curve(x, y, ierr=ierr)          ! unconstrained

call constrained%new_points(x, y)                     ! constrained concave
v = one
ierr = constrained%set_convexity(v)
ierr = constrained%fit(smoothing=0.04_FP_REAL)

! Second derivative at an interior point
print *, free_curve%dfdx(3.0_FP_REAL, order=2, ierr=ierr)    ! may be > 0
print *, constrained%dfdx(3.0_FP_REAL, order=2, ierr=ierr)   ! guaranteed <= 0

call free_curve%destroy()
call constrained%destroy()
```

The constrained spline guarantees \f$ s''(x) \leq 0 \f$ throughout.

## Mixed Constraints

The constraint vector need not be uniform. Assign different flags per point
for data that changes curvature sign:

```fortran
v(1:6)  = -one   ! convex region: s''(x_i) >= 0
v(7:m)  =  one   ! concave region: s''(x_i) <= 0
ierr = curve%set_convexity(v)
```

Setting \f$ v_i = 0 \f$ leaves individual points unconstrained, which is
useful at transition regions where the curvature sign is unknown.

## Constructor Shortcut

A one-step constructor combines `new_points`, `set_convexity`, and `fit`:

```fortran
v = one
curve = fitpack_convex_curve(x, y, v, smoothing=0.04_FP_REAL, ierr=ierr)
```

## API Summary

`fitpack_convex_curve` extends `fitpack_curve` with convexity constraints.
The spline degree is always cubic (\f$ k = 3 \f$).

| Method | Description |
|--------|-------------|
| `curve\%new_points(x, y, w)` | Load data points and optional weights |
| `curve\%set_convexity(v)` | Set per-point convexity flags (\f$ +1 \f$, \f$ -1 \f$, or \f$ 0 \f$) |
| `curve\%fit(smoothing)` | Fit with automatic knots (calls `concon`) |
| `curve\%least_squares()` | Least-squares fit with current knots (calls `cocosp`) |
| `curve\%eval(x)` | Evaluate the spline |
| `curve\%dfdx(x, order)` | Evaluate derivatives |
| `curve\%integral(a, b)` | Definite integral over \f$ [a, b] \f$ |
| `curve\%destroy()` | Release allocated memory |

## Complete Example

See `examples/example_convex_curve.f90`
(`examples/example_convex_curve.f90`) for a full working program that loads
moisture-content data with endpoint-emphasized weights, fits at multiple
smoothing levels, and prints a comparison table of knot count, residual,
and maximum second derivative.

@see @ref fitpack_convex_curves, Dierckx Ch. 8 (pp. 173&ndash;196)
