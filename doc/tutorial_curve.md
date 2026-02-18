# Univariate Curve Fitting with fitpack_curve {#tutorial_curve}

This tutorial explains how to use the `fitpack_curve` type for univariate
spline smoothing and interpolation of the form \f$ y = s(x) \f$. It covers
the mathematical background, basic and advanced usage, and practical guidance
on choosing a smoothing parameter.

Use `fitpack_curve` when you have a single-valued relationship between an
independent variable \f$ x \f$ and a dependent variable \f$ y \f$ &mdash; that
is, each \f$ x_i \f$ maps to exactly one \f$ y_i \f$. The data must be sorted
in strictly increasing order of \f$ x \f$.

## Mathematical Background

Given \f$ m \f$ data points \f$ (x_i, y_i) \f$ with positive weights
\f$ w_i \f$, FITPACK constructs a B-spline \f$ s(x) \f$ of degree \f$ k \f$
(default \f$ k = 3 \f$, cubic) by solving the constrained optimization problem

\f[
    \min \int_{x_1}^{x_m} \left[ s''(x) \right]^2 dx
    \quad \text{subject to} \quad
    \sum_{i=1}^{m} w_i^2 \left( y_i - s(x_i) \right)^2 \leq S
\f]

The **smoothing parameter** \f$ S \geq 0 \f$ controls the trade-off between
roughness and fidelity to the data:

- **\f$ S = 0 \f$** &mdash; Interpolation. The spline passes through every
  data point exactly.
- **Small \f$ S \f$** &mdash; Close to the data but may oscillate between
  points.
- **Large \f$ S \f$** &mdash; Very smooth, but may miss genuine features.

### Choosing S

When the data satisfy \f$ y_i = g(x_i) + \varepsilon_i \f$ with independent
errors of variance \f$ \sigma^2 \f$ and unit weights, the expected residual
sum of squares for the true function is \f$ m \sigma^2 \f$. A natural starting
point is therefore

\f[
    S \approx m \, \sigma^2
\f]

If no noise estimate is available, \f$ S = m \f$ (assuming unit weights) is a
reasonable first guess. The algorithm selects knots automatically: fewer knots
for larger \f$ S \f$, more knots for smaller \f$ S \f$. After fitting, inspect
`curve\%mse()` and `curve\%knots` to judge the result.

## Basic Example

The following program fits a smoothing spline to noisy sine data, evaluates
it, and prints the result.

```fortran
program basic_curve
    use fitpack, only: fitpack_curve, FP_REAL, FP_FLAG, &
                       FITPACK_SUCCESS, FITPACK_MESSAGE, &
                       zero, one, half, pi
    implicit none

    integer, parameter :: m = 50
    real(FP_REAL), parameter :: two = 2.0_FP_REAL

    real(FP_REAL) :: x(m), y(m), dx
    type(fitpack_curve) :: curve
    integer(FP_FLAG) :: ierr
    real(FP_REAL) :: y_at_pi
    integer :: i

    ! Generate data on [0, 2*pi]
    dx = two * pi / (m - 1)
    do i = 1, m
        x(i) = (i - 1) * dx
        y(i) = sin(x(i)) + 0.05_FP_REAL * sin(17.0_FP_REAL * x(i))
    end do

    ! Fit a smoothing spline (S = m)
    ierr = curve%new_fit(x, y, smoothing=real(m, FP_REAL))
    if (.not. FITPACK_SUCCESS(ierr)) then
        print *, 'Fit failed: ', FITPACK_MESSAGE(ierr)
        error stop
    end if

    ! Evaluate at x = pi
    y_at_pi = curve%eval(pi, ierr)
    print '(a,f10.6)', 's(pi)  = ', y_at_pi
    print '(a,i0)',    'knots  = ', curve%knots
    print '(a,es10.3)', 'mse    = ', curve%mse()

    call curve%destroy()
end program basic_curve
```

### Step-by-step

1. **Import** &sect; The `use fitpack` statement exposes the curve type, kind
   parameters (`FP_REAL`, `FP_FLAG`), error-handling helpers
   (`FITPACK_SUCCESS`, `FITPACK_MESSAGE`), and named constants (`zero`, `one`,
   `half`, `pi`). Note that `two` is not a public export &mdash; define it
   locally when needed.

2. **Prepare data** &sect; The abscissae `x` must be strictly increasing. The
   corresponding values `y` may contain noise.

3. **Fit** &sect; `curve\%new_fit(x, y, smoothing=S)` constructs the spline.
   It returns an `integer(FP_FLAG)` error code. Check it with
   `FITPACK_SUCCESS(ierr)`.

4. **Evaluate** &sect; `curve\%eval(x, ierr)` accepts a scalar or an array and
   returns the spline value(s).

5. **Clean up** &sect; `curve\%destroy()` releases all internal allocations.

## Advanced Features

### Weights

Supply per-point weights via the `w` argument to emphasize or de-emphasize
individual observations. Points with larger weights are fitted more closely.

```fortran
real(FP_REAL) :: w(m)

! Inverse-variance weighting
w = one / sigma
ierr = curve%new_fit(x, y, w=w, smoothing=real(m, FP_REAL))
```

When weights are proportional to \f$ 1 / \sigma_i \f$, the smoothing
constraint becomes \f$ \sum (y_i - s(x_i))^2 / \sigma_i^2 \leq S \f$, and
\f$ S \approx m \f$ remains a natural choice.

### Derivatives

Evaluate derivatives of the fitted spline up to order \f$ k \f$ (the spline
degree, default 3):

```fortran
real(FP_REAL) :: dydx, d2ydx2

dydx   = curve%dfdx(pi / 4, order=1, ierr=ierr)   ! first derivative
d2ydx2 = curve%dfdx(pi / 4, order=2, ierr=ierr)   ! second derivative
```

The `order` parameter selects which derivative: 1 for \f$ s'(x) \f$,
2 for \f$ s''(x) \f$, 3 for \f$ s'''(x) \f$.

### Integration

Compute the definite integral of the spline over an arbitrary sub-interval
\f$ [a, b] \f$:

```fortran
real(FP_REAL) :: area

area = curve%integral(zero, pi)
```

This evaluates \f$ \int_a^b s(x) \, dx \f$ analytically from the B-spline
coefficients, so the result is exact for the fitted spline (no numerical
quadrature error).

### Zero Finding

Find all zeros (roots) of a cubic spline within the data range:

```fortran
real(FP_REAL), allocatable :: roots(:)

roots = curve%zeros(ierr)
print '(a,i0)', 'Number of zeros: ', size(roots)
```

The `curve\%zeros(ierr)` function returns an allocatable array containing
every root of \f$ s(x) = 0 \f$. This operation requires a cubic spline
(\f$ k = 3 \f$).

### Interpolation

For exact interpolation &mdash; a spline that passes through every data point
&mdash; set `smoothing=zero`:

```fortran
ierr = curve%new_fit(x, y, smoothing=zero)
```

Alternatively, call the dedicated interpolation method:

```fortran
ierr = curve%interpolate()
```

Both produce the smoothest spline (minimum roughness integral) that satisfies
\f$ s(x_i) = y_i \f$ for all \f$ i \f$.

### Smoothing Sweep

When the appropriate smoothing level is unknown, try a sequence of \f$ S \f$
values on a logarithmic scale and inspect the resulting number of knots and
residual:

```fortran
integer :: j

print '(a)', '     S       knots   residual'
do j = 1, 5
    ierr = curve%new_fit(x, y, smoothing=10.0_FP_REAL**(3 - j))
    if (FITPACK_SUCCESS(ierr)) then
        print '(es10.1,i8,es12.3)', &
            curve%smoothing, curve%knots, curve%mse()
    end if
end do
```

A good fit typically shows a "knee" where adding more knots (decreasing
\f$ S \f$) yields diminishing improvement in the residual. Choose the
\f$ S \f$ value at or just past this transition point.

## API Summary

| Method | Description |
|--------|-------------|
| `curve\%new_fit(x, y, w, smoothing, order)` | Fit a smoothing spline; returns `integer(FP_FLAG)` |
| `curve\%eval(x, ierr)` | Evaluate \f$ s(x) \f$ at a scalar or array of points |
| `curve\%dfdx(x, order, ierr)` | Evaluate derivative of order 1, 2, or 3 |
| `curve\%integral(from, to)` | Definite integral \f$ \int_a^b s(x) \, dx \f$ |
| `curve\%zeros(ierr)` | All zeros of \f$ s(x) \f$ (cubic splines only) |
| `curve\%interpolate()` | Exact interpolation (\f$ S = 0 \f$) |
| `curve\%knots` | Number of knots in the current fit |
| `curve\%mse()` | Weighted residual sum of squares \f$ f_p \f$ |
| `curve\%smoothing` | Current smoothing parameter \f$ S \f$ |
| `curve\%destroy()` | Release all allocated memory |

## Complete Example

A fully compilable program demonstrating smoothing, interpolation, derivatives,
integration, root finding, and a smoothing sweep is provided in
`examples/example_curve.f90`. Build and run it with:

```
fpm run --example example_curve
```

## Error Handling

All fitting and evaluation routines return an `integer(FP_FLAG)` error code.
Use `FITPACK_SUCCESS(ierr)` to test for success and `FITPACK_MESSAGE(ierr)` to
obtain a human-readable description:

```fortran
if (.not. FITPACK_SUCCESS(ierr)) then
    print *, 'Error: ', FITPACK_MESSAGE(ierr)
end if
```

Common return values:

| Flag | Meaning |
|------|---------|
| `FITPACK_OK` (0) | Success |
| `FITPACK_INTERPOLATING_OK` (\f$-1\f$) | Interpolating spline returned |
| `FITPACK_S_TOO_SMALL` (1) | Smoothing parameter too small for the knot budget |
| `FITPACK_MAXIT` (2) | Maximum iterations reached |
| `FITPACK_INPUT_ERROR` (10) | Invalid input data |

@see @ref theory_curve_fitting
@see @ref theory_bsplines
@see @ref tutorial_periodic_curve
