# Polar and Spherical Domains {#tutorial_special}

This tutorial covers spline fitting on non-rectangular domains: polar coordinates
(disc-shaped regions) and spherical coordinates (the unit sphere).

## Polar Domain — Scattered Data

`fitpack_polar` fits a bicubic spline to data scattered over a general polar domain
\f$ x^2 + y^2 \leq r(\theta)^2 \f$, where \f$ r(\theta) \f$ is a user-supplied
boundary function.

The Cartesian coordinates are internally transformed to normalized polar
coordinates:

\f[
    x = u \, r(v) \cos v, \quad y = u \, r(v) \sin v, \quad
    0 \leq u \leq 1, \; -\pi \leq v \leq \pi
\f]

```fortran
use fitpack, only: fitpack_polar

type(fitpack_polar) :: pol
integer :: ierr

! rad_func is a function: real(FP_REAL) function rad_func(theta) result(r)
pol = fitpack_polar(x, y, z, rad_func, ierr=ierr)

! Evaluate
z_val = pol%eval(x_val, y_val, ierr)
```

### Continuity at the Origin

The continuity order at \f$ u = 0 \f$ (the origin) is configurable:

| `bc_continuity_origin` | Meaning |
|------------------------|---------|
| 0 | \f$ C^0 \f$ continuity |
| 1 | \f$ C^1 \f$ continuity |
| 2 | \f$ C^2 \f$ continuity (default) |

@see Dierckx, Ch. 11, &sect;11.1 (pp. 255&ndash;263)

## Polar Domain — Gridded Data

`fitpack_grid_polar` fits a bicubic spline to data on a polar grid with a
**constant** boundary radius \f$ r \f$:

```fortran
use fitpack, only: fitpack_grid_polar

type(fitpack_grid_polar) :: gpol
integer :: ierr

gpol = fitpack_grid_polar(u, v, z, r, ierr=ierr)
```

Additional options:
- **Origin value**: Set `z0` to prescribe the function value at the origin,
  with optional exact interpolation (`z0_exact = .true.`).
- **Zero gradient**: Set `z0_zero_gradient = .true.` to enforce
  \f$ \nabla s(0,0) = 0 \f$.

@see pogrid

## Spherical Domain — Scattered Data

`fitpack_sphere` fits a bicubic spline to data scattered on the unit sphere,
parameterized by colatitude \f$ \theta \in [0, \pi] \f$ and longitude
\f$ \phi \in [0, 2\pi] \f$:

\f[
    r = s(\theta, \phi)
\f]

```fortran
use fitpack, only: fitpack_sphere

type(fitpack_sphere) :: sph
integer :: ierr

sph = fitpack_sphere(theta, phi, r, ierr=ierr)

! Evaluate at new locations
r_val = sph%eval(theta_new, phi_new, ierr)
```

The fitted spline automatically satisfies pole constraints to ensure smoothness
at \f$ \theta = 0 \f$ (north pole) and \f$ \theta = \pi \f$ (south pole).

@see Dierckx, Ch. 11, &sect;11.2 (pp. 263&ndash;269)

## Spherical Domain — Gridded Data

`fitpack_grid_sphere` fits a bicubic spline to data on a latitude-longitude grid:

```fortran
use fitpack, only: fitpack_grid_sphere

type(fitpack_grid_sphere) :: gsph
integer :: ierr

gsph = fitpack_grid_sphere(u, v, z, ierr=ierr)
```

### Pole Boundary Conditions

Each pole (north and south) has independently configurable boundary conditions:

```fortran
! North pole: prescribed value z0 = 1.0, exact, with C1 continuity
call gsph%BC_north_pole(z0=1.0_FP_REAL, exact=.true., continuity=1)

! South pole: prescribed value z0 = -1.0, approximate
call gsph%BC_south_pole(z0=-1.0_FP_REAL, exact=.false.)
```

| Setting | Description |
|---------|-------------|
| `z0` | Function value at the pole |
| `exact` | If `.true.`, the spline passes exactly through `z0` |
| `continuity` | Continuity order at the pole (0 = \f$ C^0 \f$, 1 = \f$ C^1 \f$) |
| `zero_gradient` | If `.true.`, enforce vanishing gradient at the pole |

@see spgrid

## Summary

| Domain | Scattered | Gridded |
|--------|-----------|---------|
| Polar disc | `fitpack_polar` (polar) | `fitpack_grid_polar` (pogrid) |
| Sphere | `fitpack_sphere` (sphere) | `fitpack_grid_sphere` (spgrid) |
