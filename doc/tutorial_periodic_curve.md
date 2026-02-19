# Periodic Curve Fitting {#tutorial_periodic_curve}

This tutorial covers periodic spline fitting with `fitpack_periodic_curve`:
constructing a B-spline \f$ y = s(x) \f$ that satisfies periodic boundary
conditions, evaluating it, and computing derivatives and integrals across
the period boundary.

## Overview

A periodic spline enforces continuity of the function and all its derivatives
up to order \f$ k-1 \f$ at the period boundary:

\f[
    s^{(j)}(x_1) = s^{(j)}(x_m), \quad j = 0, 1, \ldots, k-1
\f]

This is the natural choice whenever the underlying signal is inherently
periodic, for example:

- **Cam profiles** and rotary-engine port timing (functions of crank angle)
- **Fourier-like data** sampled over one full period
- **Angular measurements** such as wind direction, phase angles, or azimuth

The type `fitpack_periodic_curve` extends `fitpack_curve` and inherits all of
its methods (`eval`, `dfdx`, `integral`, `zeros`, `destroy`, etc.). The only
difference is that fitting calls the `percur` algorithm instead of `curfit`,
which constructs a periodic knot vector and ties the B-spline coefficients at
the boundary.

## Mathematical Background

Given \f$ m \f$ data points \f$ (x_i, y_i) \f$ with \f$ x_1 < x_2 < \cdots < x_m \f$,
the period is \f$ P = x_m - x_1 \f$. The `percur` routine places a **wrapped
knot vector**: interior knots \f$ t_{k+2}, \ldots, t_{n-k-1} \f$ lie inside
\f$ [x_1, x_m) \f$, while the boundary knots satisfy

\f[
    t_{j} = t_{n-2k+j} - P, \quad j = 1, \ldots, k, \qquad
    t_{n-k+j} = t_{k+1+j} + P, \quad j = 1, \ldots, k
\f]

The B-spline coefficients are tied so that \f$ c_{n-k-1+j} = c_j \f$ for
\f$ j = 1, \ldots, k \f$, which guarantees that the spline and its first
\f$ k-1 \f$ derivatives are continuous across the period boundary.

As with `fitpack_curve`, the smoothing parameter \f$ s \f$ controls the
trade-off between fidelity and smoothness:

\f[
    \sum_{i=1}^{m-1} \left( w_i \bigl(y_i - s(x_i)\bigr) \right)^2 \leq s
\f]

Note that the last data point \f$ (x_m, y_m) \f$ defines the period length
but is **not** included in the residual sum (only \f$ m-1 \f$ points are
fitted).

## Basic Example

Fit the periodic function \f$ f(x) = \cos(x) + \sin(2x) \f$ on
\f$ [0, 2\pi] \f$. The data **must** include the right endpoint
\f$ x_m = 2\pi \f$ so that `percur` can determine the period.

```fortran
program periodic_demo
    use fitpack, only: fitpack_periodic_curve, FP_REAL, FP_FLAG, &
                       FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half, pi
    implicit none

    integer, parameter :: m = 200
    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: pi2 = two * pi

    type(fitpack_periodic_curve) :: curve
    integer(FP_FLAG) :: ierr
    real(FP_REAL) :: x(m), y(m)
    integer :: i

    ! Generate data covering one full period [0, 2*pi]
    do i = 1, m
        x(i) = (i - 1) * pi2 / (m - 1)
        y(i) = cos(x(i)) + sin(two * x(i))
    end do

    ! Fit with smoothing = 1
    ierr = curve%new_fit(x, y, smoothing=one)
    if (.not. FITPACK_SUCCESS(ierr)) then
        print *, 'Fit failed: ', FITPACK_MESSAGE(ierr)
        error stop
    end if

    print '(a,i0)',     'Knots:    ', curve%knots
    print '(a,es10.3)', 'Residual: ', curve%mse()

    call curve%destroy()
end program periodic_demo
```

The call to `curve\%new_fit` loads the data and performs a smoothing fit in one
step. Because `fitpack_periodic_curve` extends `fitpack_curve`, all evaluation
and query methods work identically.

## Periodicity Verification

A key property of the periodic spline is that function values **and**
derivatives match at the endpoints. After fitting, verify:

```fortran
real(FP_REAL) :: y_left, y_right, dy_left, dy_right

y_left   = curve%eval(x(1), ierr)
y_right  = curve%eval(x(m), ierr)
dy_left  = curve%dfdx(x(1), order=1, ierr=ierr)
dy_right = curve%dfdx(x(m), order=1, ierr=ierr)

print '(a,f12.8)', 's(x_1)  = ', y_left
print '(a,f12.8)', 's(x_m)  = ', y_right
print '(a,f12.8)', "s'(x_1) = ", dy_left
print '(a,f12.8)', "s'(x_m) = ", dy_right
```

For a cubic spline (\f$ k = 3 \f$), the values, first derivatives, and second
derivatives all match to machine precision:

\f[
    s(x_1) = s(x_m), \quad
    s'(x_1) = s'(x_m), \quad
    s''(x_1) = s''(x_m)
\f]

This is **not** the case for a non-periodic `fitpack_curve`, where the
boundary derivatives are generally unrelated.

## Integration Across the Boundary

The `curve\%integral` method works with periodic splines exactly as with
non-periodic ones. The integration limits need not coincide with the data
range. For example, integrate from \f$ \pi/2 \f$ to \f$ 3\pi \f$
(which crosses the period boundary at \f$ 2\pi \f$):

```fortran
real(FP_REAL) :: area

area = curve%integral(half * pi, 3.0_FP_REAL * pi)
print '(a,f12.8)', 'integral [pi/2, 3*pi] = ', area
```

For the test function \f$ f(x) = \cos(x) + \sin(2x) \f$ the exact integral is

\f[
    \int_{\pi/2}^{3\pi} \bigl[\cos(x) + \sin(2x)\bigr] \, dx
    = \Bigl[\sin(x) - \frac{1}{2}\cos(2x)\Bigr]_{\pi/2}^{3\pi}
    = -\frac{3}{2}
\f]

## Smoothing and Interpolation

The periodic variant supports the same fitting modes as `fitpack_curve`:

| Method | Call | Description |
|--------|------|-------------|
| Smoothing | `curve\%new_fit(x, y, smoothing=s)` | Automatic knot placement with smoothing constraint |
| Re-fit | `curve\%fit(smoothing=s)` | Change smoothing on existing data |
| Interpolation | `curve\%interpolate()` | Pass through all data points (\f$ s = 0 \f$) |

Lowering the smoothing parameter increases the number of knots and brings the
spline closer to the data:

```fortran
integer :: i

print '(a)', '     S       knots   residual'
do i = 1, 4
    ierr = curve%new_fit(x, y, smoothing=10.0_FP_REAL**(2 - i))
    if (FITPACK_SUCCESS(ierr)) then
        print '(es10.1,i8,es12.3)', curve%smoothing, curve%knots, curve%mse()
    end if
end do
```

## Complete Example

A self-contained program demonstrating periodic fitting, evaluation,
periodicity verification, integration, interpolation, and a smoothing sweep
is provided in the repository:

@see `examples/example_periodic_curve.f90`

@see @ref tutorial_curve for the non-periodic curve tutorial
