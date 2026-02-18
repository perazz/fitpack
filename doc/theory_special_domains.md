# Polar and Spherical Domains {#theory_special_domains}

Fitting on non-rectangular domains&mdash;disc-shaped regions in the plane and the
surface of a sphere&mdash;requires coordinate transformations that map the physical
domain to a rectangular parameter space where tensor-product B-splines can be
applied. These transformations introduce singularities (the origin of a disc,
the poles of a sphere) where special continuity conditions must be imposed on
the B-spline coefficients to ensure a smooth, single-valued result.

The presentation follows Dierckx (1993), Chapter 11.

## Polar Coordinates

### The polar domain

A general polar domain is the set of points \f$ (x, y) \f$ satisfying
\f$ x^2 + y^2 \leq R(\theta)^2 \f$, where \f$ R(\theta) \f$ is a positive
boundary function and \f$ \theta = \mathrm{atan2}(y, x) \f$. When
\f$ R(\theta) \f$ is constant the domain is a circular disc.

### Normalized polar transform

FITPACK maps the physical coordinates \f$ (x, y) \f$ to normalized polar
coordinates \f$ (u, v) \f$ via

\f[
    x = u \, R(v) \cos v, \qquad y = u \, R(v) \sin v
\f]

with \f$ u \in [0, 1] \f$ and \f$ v \in [-\pi, \pi] \f$. The radial variable
\f$ u \f$ runs from 0 (the origin) to 1 (the boundary), and \f$ v \f$ is the
polar angle. In this parameter space the fitting problem becomes rectangular:
find a tensor-product B-spline

\f[
    s(u, v) = \sum_{i} \sum_{j} c_{ij} \, N_{i,k}(u) \, M_{j,k}(v)
\f]

that approximates the data in a weighted least-squares sense, subject to a
smoothing condition or prescribed knot placement.

### FITPACK routines

| Routine | Input | Boundary |
|---------|-------|----------|
| `polar` (scattered) | Arbitrary \f$ (x_i, y_i, z_i) \f$ | General \f$ R(\theta) \f$ via function pointer |
| `pogrid` (gridded) | Data on a polar grid \f$ (u_i, v_j) \f$ | Constant radius \f$ R \f$ |

## Continuity at the Origin

### The problem

At \f$ u = 0 \f$ all values of the angle \f$ v \f$ map to the same physical
point&mdash;the origin. A naive tensor-product spline in \f$ (u, v) \f$ would
assign potentially different values to \f$ s(0, v) \f$ for different \f$ v \f$,
producing a discontinuous or non-smooth surface at the centre of the disc.

### Continuity orders

To obtain a physically meaningful result the B-spline coefficients near
\f$ u = 0 \f$ must satisfy algebraic constraints that enforce continuity to the
desired order.

**\f$ C^0 \f$ continuity.** The function value at the origin must be
single-valued:

\f[
    s(0, v) = c_0 \quad \text{(constant for all } v\text{)}
\f]

This is achieved by constraining all coefficients associated with the first
radial B-spline to be equal.

**\f$ C^1 \f$ continuity.** In addition to \f$ C^0 \f$, the first radial
derivative \f$ \partial s / \partial u \f$ at \f$ u = 0 \f$ must be consistent
with a well-defined gradient vector in the \f$ (x, y) \f$ plane. Without this
constraint, the directional derivative at the origin could depend on the
approach angle in an unphysical way. The constraint takes the form of a linear
relation among the coefficients of the first two radial B-splines.

**\f$ C^2 \f$ continuity.** The second radial derivatives must likewise be
consistent with a well-defined Hessian at the origin. This imposes additional
linear relations on the first three rows of coefficients.

### Implementation in FITPACK

The `bc_continuity_origin` parameter selects the continuity order:

| `bc_continuity_origin` | Continuity | Available for |
|------------------------|------------|---------------|
| 0 | \f$ C^0 \f$ | Scattered and gridded polar |
| 1 | \f$ C^1 \f$ | Scattered and gridded polar |
| 2 | \f$ C^2 \f$ (default) | Scattered polar only |

The constraints are incorporated into the least-squares system during fitting,
reducing the effective number of free coefficients near \f$ u = 0 \f$.

## Boundary Conditions on the Disc Edge

At the outer boundary \f$ u = 1 \f$ the following options are available:

- **Extrapolation (default)**: No constraint is imposed at the boundary. The
  spline is free to take any value at \f$ u = 1 \f$.
- **Zero boundary**: The spline is forced to vanish on the boundary,
  \f$ s(1, v) = 0 \f$ for all \f$ v \f$. This is useful for problems where
  the function is known to be zero at the edge of the disc.

For gridded polar fitting (`pogrid`), additional options control the behaviour
at the origin:

- **Prescribed origin value**: The value \f$ z_0 = s(0, v) \f$ can be supplied
  explicitly. It can be enforced exactly (`z0_exact = .true.`) or treated as an
  additional data point to be approximated.
- **Zero gradient at the origin**: Setting `z0_zero_gradient = .true.` enforces
  \f$ \nabla s(0, 0) = 0 \f$, creating a flat spot at the centre of the disc.

## Spherical Coordinates

### The spherical domain

The unit sphere is parameterized by colatitude \f$ \theta \in [0, \pi] \f$ and
longitude \f$ \phi \in [0, 2\pi] \f$. A scalar function on the sphere is
represented as a tensor-product B-spline:

\f[
    s(\theta, \phi) = \sum_{i} \sum_{j} c_{ij} \, N_{i,k}(\theta) \, M_{j,k}(\phi)
\f]

### Periodicity in longitude

Because the sphere is periodic in \f$ \phi \f$, the spline must satisfy

\f[
    s(\theta, 0) = s(\theta, 2\pi)
\f]

together with matching of all derivatives at \f$ \phi = 0 \f$ and
\f$ \phi = 2\pi \f$. FITPACK enforces this by using a periodic knot vector and
periodic B-spline basis in the \f$ \phi \f$ direction.

### FITPACK routines

| Routine | Input |
|---------|-------|
| `sphere` (scattered) | Arbitrary \f$ (\theta_i, \phi_i, r_i) \f$ on the sphere |
| `spgrid` (gridded) | Data on a \f$ (\theta_i, \phi_j) \f$ grid |

## The Pole Problem

### Analogy with the origin problem

The pole problem on the sphere is the direct analogue of the origin problem in
polar coordinates. At the **north pole** (\f$ \theta = 0 \f$) and the
**south pole** (\f$ \theta = \pi \f$), all longitudes \f$ \phi \f$ correspond
to the same physical point. A naive tensor-product spline would allow different
values of \f$ s(0, \phi) \f$ at different longitudes, producing a discontinuity
at the pole.

### Continuity conditions

**\f$ C^0 \f$ continuity.** The function value at each pole must be
single-valued:

\f[
    s(0, \phi) = c_{\mathrm{NP}} \quad \text{(constant for all } \phi\text{)},
    \qquad
    s(\pi, \phi) = c_{\mathrm{SP}} \quad \text{(constant for all } \phi\text{)}
\f]

**\f$ C^1 \f$ continuity.** Additionally, the angular derivatives at each pole
must be consistent with a well-defined gradient on the sphere. Without this
constraint the directional derivative would depend on the approach azimuth in a
manner incompatible with the smoothness of the underlying surface.

FITPACK enforces these conditions by imposing linear constraints on the B-spline
coefficients near \f$ \theta = 0 \f$ and \f$ \theta = \pi \f$, reducing the
number of free parameters at each pole.

### Pole boundary conditions for gridded spheres

For gridded spherical fitting (`spgrid`), the north and south poles have
**independently configurable** boundary conditions:

| Setting | Description |
|---------|-------------|
| `z0` | Prescribed function value at the pole |
| `exact` | If `.true.`, the spline passes exactly through `z0` |
| `continuity` | Continuity order at the pole (0 = \f$ C^0 \f$, 1 = \f$ C^1 \f$) |
| `zero_gradient` | If `.true.`, enforce vanishing gradient at the pole |

Each pole is configured independently via `BC_north_pole` and `BC_south_pole`,
allowing asymmetric physical conditions (e.g., a known temperature at the north
pole with a free south pole).

## Practical Considerations

- **Arbitrary boundary shape**: The scattered polar fitter (`polar`) accepts a
  general boundary function \f$ R(\theta) \f$ through a function pointer,
  allowing fitting on non-circular disc-shaped domains. The gridded polar fitter
  (`pogrid`) requires a constant radius \f$ R \f$.

- **Evaluation grid layout**: The gridded sphere evaluation method `eval_many`
  returns an output array of shape \f$ (n_\phi, n_\theta) \f$, not paired
  points. This matches the tensor-product structure and allows efficient grid
  evaluation.

- **Independent pole settings**: The gridded sphere supports independent north
  and south pole boundary conditions, making it possible to impose different
  physical constraints at each pole without affecting the other.

- **Smoothing vs. interpolation**: All four routines (`polar`, `pogrid`,
  `sphere`, `spgrid`) support both smoothing and interpolation modes. In
  smoothing mode the continuity and boundary constraints are incorporated into
  the penalized least-squares problem; in interpolation mode they appear as
  hard constraints.

@see @ref theory_surface_fitting
@see @ref tutorial_polar
@see @ref tutorial_sphere
@see @ref book_reference

> P. Dierckx, *Curve and Surface Fitting with Splines*,
> Oxford University Press, 1993, Chapter 11.
