# Surface Fitting Theory {#theory_surface_fitting}

This page describes the mathematical foundations for bivariate spline fitting
in FITPACK: extending the one-dimensional B-spline framework to approximate
surfaces \f$ z = s(x, y) \f$ defined over a rectangular domain.

## Tensor Product Splines

A bivariate spline on the rectangular domain \f$ [a, b] \times [c, d] \f$ is
defined as a tensor product of univariate B-splines:

\f[
    s(x, y) = \sum_{i=1}^{n_x - k_x - 1} \sum_{j=1}^{n_y - k_y - 1}
    c_{ij} \, N_{i, k_x}(x) \, M_{j, k_y}(y)
\f]

where:

- \f$ N_{i, k_x}(x) \f$ are B-splines of degree \f$ k_x \f$ on the knot
  vector \f$ \mathbf{t}_x = (t_{x,1}, \ldots, t_{x, n_x}) \f$,
- \f$ M_{j, k_y}(y) \f$ are B-splines of degree \f$ k_y \f$ on the knot
  vector \f$ \mathbf{t}_y = (t_{y,1}, \ldots, t_{y, n_y}) \f$,
- \f$ c_{ij} \f$ are the \f$ (n_x - k_x - 1)(n_y - k_y - 1) \f$ unknown
  coefficients.

The two knot vectors are independent: the number of knots and their placement
can differ between the \f$ x \f$- and \f$ y \f$-directions. In FITPACK, both
degrees default to \f$ k_x = k_y = 3 \f$ (bicubic).

@see @ref theory_bsplines for a detailed treatment of univariate B-spline
basis functions and the de Boor--Cox recurrence.

## Scattered Data Fitting (surfit)

Given \f$ m \f$ scattered observations \f$ (x_i, y_i, z_i) \f$ with weights
\f$ w_i > 0 \f$, the goal is to find a tensor product spline \f$ s(x, y) \f$
that approximates the data while remaining smooth.

### The Observation Matrix

Each data point \f$ (x_i, y_i) \f$ lies in a unique cell of the knot grid.
For a given point, let \f$ \mu \f$ and \f$ \nu \f$ index the non-zero basis
functions in \f$ x \f$ and \f$ y \f$, respectively. The observation matrix
\f$ \mathbf{A} \f$ has entries:

\f[
    A_{i, \ell} = N_{\mu, k_x}(x_i) \, M_{\nu, k_y}(y_i),
    \qquad \ell = (\mu - 1)(n_y - k_y - 1) + \nu
\f]

where the two-dimensional index pair \f$ (\mu, \nu) \f$ is mapped to a single
column index \f$ \ell \f$. Each row of \f$ \mathbf{A} \f$ has at most
\f$ (k_x + 1)(k_y + 1) \f$ non-zero entries, so the matrix is sparse.

### Least-Squares Formulation

The weighted least-squares problem is:

\f[
    \min_{\mathbf{c}} \sum_{i=1}^{m} w_i^2
    \left( z_i - s(x_i, y_i) \right)^2
\f]

which in matrix form becomes \f$ \min \| \mathbf{W}(\mathbf{z} -
\mathbf{A} \mathbf{c}) \|_2^2 \f$ with \f$ \mathbf{W} =
\mathrm{diag}(w_1, \ldots, w_m) \f$.

### Smoothing Formulation

The smoothing spline minimizes a roughness functional subject to a data
fidelity constraint:

\f[
    \min \iint \left[
        \left( \frac{\partial^2 s}{\partial x^2} \right)^2
        + 2 \left( \frac{\partial^2 s}{\partial x \, \partial y} \right)^2
        + \left( \frac{\partial^2 s}{\partial y^2} \right)^2
    \right] dx \, dy
    \quad \text{subject to} \quad
    \sum_{i=1}^{m} w_i^2 (z_i - s(x_i, y_i))^2 \leq S
\f]

where \f$ S \f$ is the smoothing factor. This is the bivariate analogue of
the univariate smoothing spline problem.

### Rank Deficiency

When few data points fall within certain cells of the knot grid, the
observation matrix \f$ \mathbf{A} \f$ can become rank deficient. This is more
common in the bivariate case than in one dimension, because adding knots in
both directions creates cells that may contain no data at all.

FITPACK handles rank deficiency using Householder QR factorization with
column pivoting. When the numerical rank \f$ r \f$ of the system is less than
the number of unknowns, the algorithm computes a minimum-norm solution among
all least-squares solutions.

### Adaptive Knot Placement

Starting from a minimal knot set, `surfit` iteratively refines the
approximation:

1. Fit the current knot configuration by solving the least-squares problem.
2. If the residual sum of squares exceeds \f$ S \f$, identify regions with
   large residuals.
3. Insert new knots alternately in the \f$ x \f$- and \f$ y \f$-directions,
   placing each knot where the local residual is largest.
4. Repeat until \f$ \sum w_i^2 (z_i - s(x_i, y_i))^2 \leq S \f$ is
   satisfied or a maximum knot count is reached.

The alternating strategy ensures balanced refinement across both coordinate
directions.

@see Dierckx, Ch. 5, &sect;5.3 (pp. 85&ndash;98)

## Gridded Data Fitting (regrid)

When the data lie on a rectangular grid \f$ z_{ji} = z(y_j, x_i) \f$ with
\f$ i = 1, \ldots, m_x \f$ and \f$ j = 1, \ldots, m_y \f$, the problem has
a Kronecker product structure that can be exploited for dramatic computational
savings.

### Kronecker Product Structure

On a grid, the observation matrix factors as a Kronecker product:

\f[
    \mathbf{A}_{\text{grid}} = \mathbf{B}_x \otimes \mathbf{B}_y
\f]

where \f$ \mathbf{B}_x \f$ is the \f$ m_x \times (n_x - k_x - 1) \f$ matrix
of univariate B-spline values at the \f$ x \f$-grid points, and
\f$ \mathbf{B}_y \f$ is the \f$ m_y \times (n_y - k_y - 1) \f$ matrix at
the \f$ y \f$-grid points. The coefficient matrix conceptually satisfies:

\f[
    \mathbf{C} = \mathbf{B}_x^{-1} \, \mathbf{Z} \, (\mathbf{B}_y^T)^{-1}
\f]

where \f$ \mathbf{Z} \f$ is the \f$ m_x \times m_y \f$ matrix of observed
values and \f$ \mathbf{C} \f$ holds the spline coefficients.

### Computational Advantage

The factored form means the normal equations can be solved dimension by
dimension rather than simultaneously:

| Approach | Complexity |
|----------|-----------|
| Scattered (`surfit`) | \f$ \mathcal{O}(m \cdot n_x \cdot n_y) \f$ |
| Gridded (`regrid`) | \f$ \mathcal{O}(m_x \cdot m_y \cdot (k_x + k_y)) \f$ |

For a 100 x 100 grid with 20 knots in each direction, the gridded approach is
roughly two orders of magnitude faster. **Always prefer `fitpack_grid_surface`
when data lie on a rectangular grid.**

The smoothing framework (roughness minimization subject to a fidelity
constraint) is identical to the scattered case; only the numerical solution
method differs.

@see Dierckx, Ch. 5, &sect;5.4 (pp. 98&ndash;103)

## The Smoothing Norm

For bivariate splines, the roughness measure used in the smoothing formulation
is the thin-plate energy functional:

\f[
    R[s] = \iint \left[
        \left( \frac{\partial^2 s}{\partial x^2} \right)^2
        + 2 \left( \frac{\partial^2 s}{\partial x \, \partial y} \right)^2
        + \left( \frac{\partial^2 s}{\partial y^2} \right)^2
    \right] dx \, dy
\f]

This is the natural extension of \f$ \int (s'')^2 \, dx \f$ to two
dimensions. It penalizes curvature equally in all directions and is invariant
under rotation of the coordinate axes. The cross-derivative term
\f$ (\partial^2 s / \partial x \, \partial y)^2 \f$ ensures that diagonal
oscillations are penalized as well.

For a tensor product spline, this integral can be expressed as a quadratic
form in the coefficients \f$ c_{ij} \f$, involving integrals of products of
B-spline second derivatives. FITPACK computes these integrals analytically
using the knot vectors.

## Choosing the Smoothing Factor

The smoothing factor \f$ S \f$ controls the balance between fidelity and
smoothness:

- **\f$ S = 0 \f$**: Interpolation. The spline passes through every data
  point (requires enough knots to satisfy the interpolation conditions).
- **Small \f$ S \f$**: Close fit to the data, but the surface may exhibit
  spurious oscillations.
- **Large \f$ S \f$**: Very smooth surface that may underfit the data.

### Practical Guidelines

For **scattered data** with unit weights (\f$ w_i = 1 \f$), a reasonable
starting point is:

\f[
    S \approx m
\f]

where \f$ m \f$ is the number of data points. This follows from the
expectation that, for a good fit with normally distributed residuals, the sum
of squared residuals should be \f$ \mathcal{O}(m) \f$.

For **gridded data** with unit weights, the analogous rule is:

\f[
    S \approx m_x \cdot m_y
\f]

After fitting, inspect the mean squared error via the `mse()` method to assess
fit quality. If the residual is too large, decrease \f$ S \f$; if the surface
oscillates, increase \f$ S \f$.

When data have non-unit weights, interpret \f$ S \f$ relative to the weighted
residual. Setting \f$ w_i = 1 / \sigma_i \f$ (where \f$ \sigma_i \f$ is the
measurement standard deviation) makes the residual sum an approximate
\f$ \chi^2 \f$ statistic, and \f$ S \approx m \f$ corresponds to a
reduced \f$ \chi^2 \f$ of unity.

## Parametric Surfaces

For surfaces embedded in \f$ \mathbb{R}^d \f$ and parameterized by
\f$ (u, v) \f$, each coordinate component is represented as a separate
bivariate spline sharing common knot vectors:

\f[
    s_l(u, v) = \sum_{i=1}^{n_u - k_u - 1} \sum_{j=1}^{n_v - k_v - 1}
    c_{ij}^{(l)} \, N_{i, k_u}(u) \, M_{j, k_v}(v),
    \qquad l = 1, \ldots, d
\f]

All \f$ d \f$ components share the same knot vectors
\f$ \mathbf{t}_u \f$ and \f$ \mathbf{t}_v \f$, and the same spline degrees
\f$ k_u \f$ and \f$ k_v \f$. Only the coefficient arrays
\f$ c_{ij}^{(l)} \f$ differ between components.

### Periodicity

Either parameter direction (or both) may be periodic:

- **Periodic in \f$ u \f$**: enforces
  \f$ s_l^{(p)}(u_{\min}, v) = s_l^{(p)}(u_{\max}, v) \f$ for
  \f$ p = 0, 1, \ldots, k_u - 1 \f$.
- **Periodic in \f$ v \f$**: enforces
  \f$ s_l^{(p)}(u, v_{\min}) = s_l^{(p)}(u, v_{\max}) \f$ for
  \f$ p = 0, 1, \ldots, k_v - 1 \f$.

Periodicity is useful for closed surfaces such as cylinders and tori. The
FITPACK routine `parsur` handles the general case.

@see Dierckx, Ch. 10, &sect;10.2 (pp. 241&ndash;254)

## Summary

| Problem | FITPACK Type | Core Routine | Key Advantage |
|---------|-------------|-------------|---------------|
| Scattered \f$ (x_i, y_i, z_i) \f$ | `fitpack_surface` | `surfit` | Handles arbitrary point distributions |
| Gridded \f$ z(y_j, x_i) \f$ | `fitpack_grid_surface` | `regrid` | Kronecker structure, much faster |
| Parametric \f$ \mathbf{r}(u, v) \f$ | `fitpack_parametric_surface` | `parsur` | Surfaces in \f$ \mathbb{R}^d \f$ |

@see @ref tutorial_surface for code examples and practical usage
@see @ref tutorial_polar for polar and spherical domain fitting
@see @ref book_reference for a complete mapping of routines to book chapters
