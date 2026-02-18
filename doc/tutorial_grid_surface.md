# Grid and Parametric Surface Tutorial {#tutorial_grid_surface}

This tutorial covers two specialized surface types for data on rectangular
parameter grids: `fitpack_grid_surface` for scalar-valued surfaces
\f$ z = s(x, y) \f$, and `fitpack_parametric_surface` for multi-component
surfaces \f$ \mathbf{s}(u, v) \in \mathbb{R}^d \f$.

Both exploit the Kronecker product structure of tensor-product B-splines on
rectangular grids. When the data points \f$ (x_i, y_j) \f$ form a full grid,
the normal equations decouple into smaller systems along each axis, reducing the
computational cost from \f$ \mathcal{O}(n_x \, n_y) \f$ to
\f$ \mathcal{O}(n_x + n_y) \f$ per iteration. This makes gridded fitting
significantly faster than the general scattered-data algorithm (`surfit`).

## fitpack_grid_surface

### Creating a Fit

Construct a grid surface from coordinate vectors `x(nx)`, `y(ny)` and function
values `z(ny, nx)`. Note that **the first index of `z` corresponds to `y`** and
the second to `x`:

```fortran
use fitpack, only: fitpack_grid_surface, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half

type(fitpack_grid_surface) :: gsurf
integer(FP_FLAG) :: ierr

ierr = gsurf%new_fit(x, y, z, smoothing=s)
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)
```

The `smoothing` parameter controls the trade-off between fidelity and
smoothness, just as for curves. Setting `smoothing=zero` produces an
interpolating spline. An optional `order` argument sets the B-spline degree
(default 3, i.e. bicubic).

### Evaluation

Evaluate the fitted surface at scattered points or on a rectangular output grid:

```fortran
z_val = gsurf%eval(x0, y0, ierr)              ! single point
z_pts = gsurf%eval(x_pts, y_pts, ierr)        ! scattered points
z_out = gsurf%eval_ongrid(x_out, y_out, ierr)  ! grid -> z_out(my, mx)
```

The `gsurf\%eval` generic dispatches to scattered-point evaluation (`bispeu`),
while `gsurf\%eval_ongrid` uses the faster grid evaluator (`bispev`) and returns
an array with shape `(size(y_out), size(x_out))`.

### Partial Derivatives

Compute partial derivatives
\f$ \partial^{n_x + n_y} s / \partial x^{n_x} \partial y^{n_y} \f$:

```fortran
dzdx = gsurf%dfdx(x0, y0, 1, 0, ierr)              ! ds/dx at a point
d2z  = gsurf%dfdx(x_pts, y_pts, 1, 1, ierr)         ! d2s/dxdy scattered
dzdy = gsurf%dfdx_ongrid(x_out, y_out, 0, 1, ierr)  ! ds/dy on a grid
```

The integer arguments are the derivative orders \f$ n_x \f$ and \f$ n_y \f$
(both must satisfy \f$ 0 \leq n_x < k_x \f$ and \f$ 0 \leq n_y < k_y \f$).

### Integration

Compute the double integral over a rectangular sub-domain:

\f[
    I = \iint_{[x_1, x_2] \times [y_1, y_2]} s(x, y) \, dx \, dy
\f]

```fortran
vol = gsurf%integral([x1, y1], [x2, y2])
```

The two arguments are the lower and upper corners of the integration rectangle.

### Cross-Sections and Derivative Splines

Extract a 1D profile at a fixed coordinate value:

```fortran
type(fitpack_curve) :: profile
profile = gsurf%cross_section(x0, .true., ierr)  ! f(y) = s(x0, y)
```

The returned `fitpack_curve` is a full 1D spline that can be evaluated,
differentiated, and integrated independently.

To build a new surface object whose evaluation gives a partial derivative:

```fortran
type(fitpack_grid_surface) :: dsurf
dsurf = gsurf%derivative_spline(1, 0, ierr)   ! ds/dx
z_val = dsurf%eval(x0, y0, ierr)
```

### Example: Peaks Function on a Grid

```fortran
use fitpack, only: fitpack_grid_surface, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, half
implicit none

integer, parameter :: nx = 20, ny = 25
real(FP_REAL) :: xg(nx), yg(ny), zg(ny, nx)
type(fitpack_grid_surface) :: gsurf
integer(FP_FLAG) :: ierr
real(FP_REAL) :: z_val
integer :: i, j

! Build coordinate vectors: x, y in [-2, 2]
do i = 1, nx
    xg(i) = -2.0_FP_REAL + 4.0_FP_REAL * (i - 1) / (nx - 1)
end do
do j = 1, ny
    yg(j) = -2.0_FP_REAL + 4.0_FP_REAL * (j - 1) / (ny - 1)
end do

! Sample f(x,y) = exp(-(x^2+y^2)/2) * cos(x) * sin(y)
do i = 1, nx
    do j = 1, ny
        zg(j, i) = exp(-half * (xg(i)**2 + yg(j)**2)) &
                  * cos(xg(i)) * sin(yg(j))
    end do
end do

! Interpolating fit (smoothing = 0)
ierr = gsurf%new_fit(xg, yg, zg, smoothing=zero)
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)

! Evaluate at a test point
z_val = gsurf%eval(zero, zero, ierr)
print '(a,f12.6)', 's(0, 0) = ', z_val

call gsurf%destroy()
```

## fitpack_parametric_surface

### Multi-Component Surfaces

`fitpack_parametric_surface` fits a bicubic tensor-product spline surface
\f$ \mathbf{s}(u, v) = (s_1(u,v), \ldots, s_d(u,v)) \f$ through data given on
a rectangular parameter grid \f$ (u_i, v_j) \f$. This is the natural choice for
geometry embedded in \f$ \mathbb{R}^3 \f$ (or higher dimensions) where the
surface is described by a mapping from a 2D parameter domain.

### Creating a Fit

The data array `z` has shape `(nv, nu, idim)`, where `idim` is the number of
coordinate components:

```fortran
use fitpack, only: fitpack_parametric_surface, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, pi

type(fitpack_parametric_surface) :: psurf
integer(FP_FLAG) :: ierr

ierr = psurf%new_fit(u, v, z, smoothing=zero, &
                     periodic_BC=[.true., .true.])
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)
```

The `periodic_BC` argument is a 2-element logical array `[periodic_u, periodic_v]`.
When an element is `.true.`, the corresponding parameter direction is treated as
periodic, enforcing continuity of the spline and its derivatives across the
boundary.

### Evaluation

Evaluation returns an array with the surface coordinates:

```fortran
pt = psurf%eval(u0, v0, ierr)        ! single point -> pt(idim)
f  = psurf%eval(u_out, v_out, ierr)  ! grid -> f(mv, mu, idim)
```

### Example: Fitting a Torus

A torus with major radius \f$ R \f$ and tube radius \f$ r \f$ is periodic in
both \f$ u \f$ (around the hole) and \f$ v \f$ (around the tube):

\f[
    \mathbf{s}(u, v) = \begin{pmatrix}
        (R + r \cos v) \cos u \\
        (R + r \cos v) \sin u \\
        r \sin v
    \end{pmatrix},
    \quad u, v \in [0, 2\pi).
\f]

```fortran
use fitpack, only: fitpack_parametric_surface, FP_REAL, FP_FLAG, &
                   FITPACK_SUCCESS, FITPACK_MESSAGE, zero, pi
implicit none

integer, parameter :: nu = 20, nv = 30, idim = 3
real(FP_REAL), parameter :: two = 2.0_FP_REAL
real(FP_REAL), parameter :: R = two, r_tube = 0.5_FP_REAL
real(FP_REAL) :: u(nu), v(nv), z(nv, nu, idim)
type(fitpack_parametric_surface) :: psurf
integer(FP_FLAG) :: ierr
real(FP_REAL) :: pt(idim)
integer :: i, j

! Parameter grid in [0, 2*pi)
do i = 1, nu
    u(i) = two * pi * (i - 1) / nu
end do
do j = 1, nv
    v(j) = two * pi * (j - 1) / nv
end do

! Sample the torus
do i = 1, nu
    do j = 1, nv
        z(j, i, 1) = (R + r_tube * cos(v(j))) * cos(u(i))
        z(j, i, 2) = (R + r_tube * cos(v(j))) * sin(u(i))
        z(j, i, 3) = r_tube * sin(v(j))
    end do
end do

! Interpolating fit, periodic in both directions
ierr = psurf%new_fit(u, v, z, smoothing=zero, &
                     periodic_BC=[.true., .true.])
if (.not. FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)

! Evaluate: expect (R + r, 0, 0) at (0, 0)
pt = psurf%eval(zero, zero, ierr)
print '(a,3f10.4)', 's(0, 0) = ', pt

call psurf%destroy()
```

## Grid vs Scattered Fitting

Choose the surface type based on the structure of your input data:

| Criterion | `fitpack_grid_surface` | `fitpack_surface` |
|-----------|----------------------|-------------------|
| Input layout | Rectangular grid \f$ z(y_j, x_i) \f$ | Arbitrary \f$ (x_i, y_i, z_i) \f$ |
| Core routine | `regrid` | `surfit` |
| Speed | Fast (Kronecker structure) | General |
| Weights | Not supported | Per-point weights |
| Typical use | Tabulated functions, image data | Sensor readings, scattered measurements |

**Rule of thumb**: if the data can be arranged on a rectangular grid with no
missing values, `fitpack_grid_surface` is always the better choice. For
irregular point clouds, use `fitpack_surface`.

For surfaces embedded in \f$ \mathbb{R}^d \f$ (e.g. geometry defined by a
parameter mapping), use `fitpack_parametric_surface` instead.

## Complete Example

A full working program demonstrating both `fitpack_grid_surface` (peaks
function with evaluation, integration, cross-sections, and derivative splines)
and `fitpack_parametric_surface` (torus with periodic boundary conditions) is
available in:

```
examples/example_grid_surface.f90
```

Build and run it with:

```
fpm run --example example_grid_surface
```

@see @ref fitpack_grid_surfaces, @ref fitpack_parametric_surfaces,
     Dierckx Ch. 5 &sect;5.4 (pp. 98&ndash;103), Ch. 10 &sect;10.2 (pp. 241&ndash;254)
