# Parametric, Closed, and Constrained Curves {#tutorial_parametric_curves}

This tutorial covers parametric curve fitting in \f$ \mathbb{R}^d \f$ with
FITPACK: open parametric curves, closed (periodic) curves, and curves with
prescribed endpoint derivatives.

## Overview

A standard 1D spline \f$ y = s(x) \f$ cannot represent multi-valued or looping
paths. A **parametric curve** expresses each coordinate as a separate spline
over a common parameter \f$ u \f$:

\f[
    \mathbf{s}(u) = \bigl( s_1(u),\; s_2(u),\; \ldots,\; s_d(u) \bigr),
    \quad u \in [u_1, u_m]
\f]

All \f$ d \f$ component splines share a single knot vector and degree \f$ k \f$.
Use parametric curves for multi-dimensional paths, closed loops, or constrained
endpoints with prescribed positions and derivatives.

| Type | Description | Core routine |
|------|-------------|-------------|
| `fitpack_parametric_curve` | Open parametric spline | `parcur` |
| `fitpack_closed_curve` | Periodic parametric spline | `clocur` |
| `fitpack_constrained_curve` | Endpoint derivative constraints | `concur` |

@see Dierckx, Ch. 9 (pp. 199&ndash;228)

## Parametric Curves

### Chord-Length Parameterization

When no parameter values are supplied, FITPACK computes them from
**cumulative chord lengths**, normalized to \f$ [0, 1] \f$:

\f[
    u_1 = 0, \quad
    u_i = u_{i-1} + \| \mathbf{x}_i - \mathbf{x}_{i-1} \|, \quad
    u_m = 1
\f]

You may also supply your own strictly increasing parameter values (e.g. time
stamps) via the optional `u` argument.

### API

Construct and fit with `new_fit`:

```fortran
use fitpack, only: fitpack_parametric_curve, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero

type(fitpack_parametric_curve) :: curve
integer(FP_FLAG) :: ierr

! x has shape (idim, m): idim dimensions, m points
ierr = curve%new_fit(x, u, w, smoothing, order)
```

Arguments: **`x(idim, m)`** data points; optional **`u(m)`** parameter values;
optional **`w(m)`** weights; optional **`smoothing`** (\f$ s \f$, use `zero`
for interpolation); optional **`order`** (spline degree \f$ k \f$, default 3).

Evaluation returns arrays whose first dimension is `idim`:

- `curve\%eval_one(u, ierr)` &mdash; returns `array(idim)` at a single \f$ u \f$.
- `curve\%eval_many(u, ierr)` &mdash; returns `array(idim, m)` at multiple \f$ u \f$.
- The generic `curve\%eval` dispatches to either form.

Derivatives with respect to \f$ u \f$:

- `curve\%dfdx(u, order, ierr)` &mdash; returns `array(idim)`, the `order`-th derivative.
- `curve\%dfdx_all(u, ierr)` &mdash; returns `array(idim, 0:k)`, all derivatives at once.

### Example: Lissajous Figure

A Lissajous curve \f$ x = \sin(3t),\; y = \sin(2t) \f$ is multi-valued in
both coordinates, making it a natural test case.

```fortran
use fitpack, only: fitpack_parametric_curve, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, pi
implicit none

real(FP_REAL), parameter :: two = 2.0_FP_REAL, three = 3.0_FP_REAL
real(FP_REAL), parameter :: pi2 = two * pi
integer, parameter :: m = 60
real(FP_REAL) :: pts(2, m), u(m), du, y(2), dy(2)
type(fitpack_parametric_curve) :: curve
integer(FP_FLAG) :: ierr
integer :: i

du = pi2 / (m - 1)
do i = 1, m
    u(i) = (i - 1) * du
    pts(1, i) = sin(three * u(i))
    pts(2, i) = sin(two * u(i))
end do

! Interpolation (s = 0) with automatic chord-length parameters
ierr = curve%new_fit(pts, smoothing=zero)
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)

y  = curve%eval_one(zero, ierr)            ! position at u=0
dy = curve%dfdx(zero, order=1, ierr=ierr)  ! tangent at u=0

call curve%destroy()
```

## Closed Curves

`fitpack_closed_curve` extends `fitpack_parametric_curve` with periodic
boundary conditions, guaranteeing:

\f[
    s_j^{(\ell)}(u_1) = s_j^{(\ell)}(u_m), \quad
    \ell = 0, 1, \ldots, k-1, \quad j = 1, \ldots, d
\f]

**Do not** duplicate the first point as the last; periodicity is handled
internally by `clocur`.

### Example: Ellipse

```fortran
use fitpack, only: fitpack_closed_curve, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, pi
implicit none

real(FP_REAL), parameter :: two = 2.0_FP_REAL, pi2 = two * pi
integer, parameter :: m = 40
real(FP_REAL) :: pts(2, m), du, y_begin(2), y_end(2)
type(fitpack_closed_curve) :: curve
integer(FP_FLAG) :: ierr
integer :: i

! Sample ellipse x = 2*cos(t), y = sin(t); do NOT repeat first point
du = pi2 / m
do i = 1, m
    pts(1, i) = two * cos((i - 1) * du)
    pts(2, i) = sin((i - 1) * du)
end do

ierr = curve%new_fit(pts, smoothing=zero)
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)

! Verify closure at parameter boundaries
y_begin = curve%eval_one(curve%ubegin, ierr)
y_end   = curve%eval_one(curve%uend, ierr)
print '(a, 2f10.5)', 's(u_begin) = ', y_begin
print '(a, 2f10.5)', 's(u_end)   = ', y_end
! These should be nearly identical.

call curve%destroy()
```

The fields `curve\%ubegin` and `curve\%uend` give the parameter range. For a
closed curve, evaluating at either boundary returns the same point.

## Constrained Curves

`fitpack_constrained_curve` extends `fitpack_parametric_curve` to prescribe
function values and derivatives at the endpoints. This is useful when the curve
must join smoothly to another segment, or when physical boundary conditions
(position, velocity, acceleration) are known.

### Setting Constraints

```fortran
call curve%set_constraints(ddx_begin, ddx_end, ierr)
```

- **`ddx_begin(idim, 0:nb-1)`** &mdash; Left endpoint: column 0 is position
  \f$ \mathbf{s}(u_1) \f$, column \f$ \ell \f$ is \f$ \mathbf{s}^{(\ell)}(u_1) \f$.
- **`ddx_end(idim, 0:ne-1)`** &mdash; Right endpoint (same layout).
- Both arguments are optional; omitting one removes constraints at that endpoint.

Shape `(idim, 0:1)` prescribes position and first derivative; `(idim, 0:2)`
adds the second derivative, and so on.

### Example: Spiral with Prescribed Tangents

Fit \f$ x = u\cos u,\; y = u\sin u \f$ over \f$ [0, \pi] \f$ with prescribed
endpoint positions and tangent vectors.

```fortran
use fitpack, only: fitpack_constrained_curve, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, pi
implicit none

integer, parameter :: m = 25, idim = 2
real(FP_REAL) :: pts(idim, m), u(m), du
real(FP_REAL) :: ddx_begin(idim, 0:1), ddx_end(idim, 0:1)
real(FP_REAL), allocatable :: yall(:,:)
type(fitpack_constrained_curve) :: curve
integer(FP_FLAG) :: ierr
integer :: i

du = pi / (m - 1)
do i = 1, m
    u(i) = (i - 1) * du
    pts(1, i) = u(i) * cos(u(i))
    pts(2, i) = u(i) * sin(u(i))
end do

! Position and tangent at u = 0
ddx_begin(:, 0) = [zero, zero]    ! s(0) = (0, 0)
ddx_begin(:, 1) = [one, zero]     ! s'(0) = (1, 0)

! Position and tangent at u = pi
ddx_end(:, 0) = [-pi, zero]       ! s(pi) = (-pi, 0)
ddx_end(:, 1) = [-one, -pi]       ! s'(pi) = (-1, -pi)

! Step 1: initial fit with data and parameter values
ierr = curve%new_fit(pts, u=u, smoothing=real(m, FP_REAL))
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)

! Step 2: apply constraints and re-fit
call curve%set_constraints(ddx_begin, ddx_end, ierr)
ierr = curve%fit(smoothing=real(m, FP_REAL))
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)

! Verify constraints via dfdx_all (returns all derivatives 0..k)
yall = curve%dfdx_all(curve%ubegin, ierr)
print '(a, 2f10.5)', 's(0)   = ', yall(:, 0)
print '(a, 2f10.5)', 's''(0)  = ', yall(:, 1)

yall = curve%dfdx_all(curve%uend, ierr)
print '(a, 2f10.5)', 's(pi)  = ', yall(:, 0)
print '(a, 2f10.5)', 's''(pi) = ', yall(:, 1)

call curve%destroy()
```

Note the two-step workflow: `new_fit` loads data and produces an initial fit;
then `set_constraints` + `curve\%fit` re-fits with boundary conditions enforced.

## Complete Example

A full working program combining all three curve types is provided in
`examples/example_parametric_curves.f90`. Build and run with:

```
fpm run --example example_parametric_curves
```

## Error Handling

All routines return an `integer(FP_FLAG)` error code:

```fortran
use fitpack, only: FITPACK_OK, FITPACK_SUCCESS, FITPACK_MESSAGE

if (.not. FITPACK_SUCCESS(ierr)) then
    print *, 'Error: ', FITPACK_MESSAGE(ierr)
end if
```

When the optional `ierr` argument is omitted, the library halts on failure.

## Summary

| Type | Use case | Key methods |
|------|----------|-------------|
| `fitpack_parametric_curve` | Open paths in \f$ \mathbb{R}^d \f$ | `new_fit`, `eval`, `dfdx` |
| `fitpack_closed_curve` | Periodic loops | `new_fit`, `eval`, `ubegin`/`uend` |
| `fitpack_constrained_curve` | Prescribed endpoint derivatives | `set_constraints`, `fit`, `dfdx_all` |

@see @ref fitpack_parametric_curves, Dierckx Ch. 9 (pp. 199&ndash;228)
