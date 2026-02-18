# Curve Fitting Theory {#theory_curve_fitting}

This page describes the mathematical foundations of the curve fitting algorithms
in FITPACK: given data points \f$ (x_i, y_i) \f$ for \f$ i = 1, \ldots, m \f$,
find a spline \f$ s(x) \f$ in the B-spline space \f$ S(k, \mathbf{t}) \f$ that
represents the data well. The library supports several formulations of this
problem, from simple least-squares with user-supplied knots to fully automatic
smoothing with adaptive knot placement.

## The Approximation Problem

Given \f$ m \f$ data points \f$ (x_i, y_i) \f$ with strictly increasing abscissae

\f[
    x_1 < x_2 < \cdots < x_m
\f]

and optional positive weights \f$ w_i \f$, the goal is to find a spline
\f$ s \in S(k, \mathbf{t}) \f$ of degree \f$ k \f$ on a knot vector
\f$ \mathbf{t} \f$ that approximates the data. Three common formulations arise:

- **Interpolation.** Require \f$ s(x_i) = y_i \f$ for every data point. This
  determines the spline exactly when \f$ n = m \f$ (number of knots minus
  degree minus one equals number of data points). No smoothing is applied.

- **Least-squares with fixed knots.** Given a prescribed knot vector, find
  the B-spline coefficients \f$ c_j \f$ that minimize the weighted residual
  sum of squares. The number and placement of knots are chosen by the user.

- **Smoothing.** Let the algorithm choose the knot vector automatically,
  balancing closeness-of-fit against smoothness of the spline. This is
  FITPACK's primary mode of operation.

@see Dierckx, Ch. 4&ndash;5 (pp. 49&ndash;103)

## Least-Squares with Fixed Knots

Suppose a knot vector \f$ \mathbf{t} = (t_1, \ldots, t_{n+k+1}) \f$ is given.
The spline has the B-spline representation

\f[
    s(x) = \sum_{j=1}^{n} c_j \, B_j(x; k, \mathbf{t})
\f]

where \f$ n \f$ is the number of B-spline basis functions and \f$ c_j \f$ are
the unknown coefficients. The least-squares problem is to minimize

\f[
    \sum_{i=1}^{m} w_i^2 \left( y_i - \sum_{j=1}^{n} c_j \, B_j(x_i) \right)^2
\f]

over all coefficient vectors \f$ \mathbf{c} \in \mathbb{R}^n \f$.

### Normal Equations and Banded Structure

Setting the gradient to zero yields the normal equations
\f$ A^T W^2 A \, \mathbf{c} = A^T W^2 \mathbf{y} \f$, where \f$ A_{ij} = B_j(x_i) \f$
is the \f$ m \times n \f$ collocation matrix and \f$ W = \mathrm{diag}(w_i) \f$.
Because each B-spline \f$ B_j \f$ is nonzero on at most \f$ k+1 \f$ knot spans,
the matrix \f$ A \f$ has at most \f$ k+1 \f$ nonzero entries per row. The resulting
normal equations matrix \f$ A^T W^2 A \f$ is symmetric, positive semi-definite,
and banded with bandwidth \f$ k+1 \f$.

### QR Factorization via Givens Rotations

Rather than forming and solving the normal equations directly (which squares the
condition number), FITPACK solves the least-squares problem by computing a QR
factorization of the weighted matrix \f$ W A \f$ using Givens rotations. Each
data point \f$ (x_i, y_i) \f$ contributes one row to the system; Givens rotations
annihilate entries below the diagonal one at a time, maintaining the banded
structure throughout.

This approach has several advantages:

- **Numerical stability.** The condition number of the triangular factor
  \f$ R \f$ equals \f$ \kappa(WA) \f$, not \f$ \kappa(WA)^2 \f$.
- **Efficiency.** Each rotation touches only \f$ \mathcal{O}(k) \f$ elements,
  giving \f$ \mathcal{O}(m k^2) \f$ total work.
- **Incremental updates.** Adding or removing knots requires only local
  updates to the factorization.

The back-substitution step solves the upper triangular system \f$ R \mathbf{c} = Q^T W \mathbf{y} \f$
for the coefficients.

@see Dierckx, Ch. 4, &sect;4.2 (pp. 53&ndash;62)

## The Smoothing Problem

Instead of fixing the knot vector a priori, the smoothing formulation treats
the number and positions of knots as unknowns. The problem is to find a spline
\f$ s \in S(k, \mathbf{t}) \f$ that minimizes a roughness measure subject to
a constraint on the residual:

\f[
    \min_{s \in S(k, \mathbf{t})} \int_{x_1}^{x_m} \left[ s''(x) \right]^2 \, dx
    \quad \text{subject to} \quad
    \sum_{i=1}^{m} w_i^2 \left( y_i - s(x_i) \right)^2 \leq S
\f]

Here \f$ S \geq 0 \f$ is the **smoothing parameter** that controls the trade-off:

- **\f$ S = 0 \f$**: The constraint forces interpolation. The spline passes
  through every data point, and the roughness integral is minimized among all
  interpolating splines.
- **Small \f$ S \f$**: The spline stays close to the data but may oscillate
  between data points, requiring many knots.
- **Large \f$ S \f$**: The spline is very smooth (few knots) but may miss
  important features in the data.

For a given \f$ S \f$, the theoretical solution is a natural spline of degree
\f$ 2k + 1 \f$ with knots at the data sites. FITPACK approximates this solution
using a spline of degree \f$ k \f$ on a much smaller knot set, determined
adaptively.

@see Dierckx, Ch. 5, &sect;5.2 (pp. 67&ndash;84)

## Adaptive Knot Placement

FITPACK's `curfit` routine implements an iterative strategy to determine both
the number and positions of knots:

1. **Initialize.** Start with the minimal knot vector containing \f$ 2(k+1) \f$
   knots: \f$ k+1 \f$ knots stacked at each end of the data range.

2. **Solve.** Compute the least-squares spline for the current knot vector
   using Givens rotations. Let \f$ f_p = \sum w_i^2 (y_i - s(x_i))^2 \f$ be the
   weighted residual sum of squares.

3. **Test convergence.** If \f$ f_p \leq S \f$, the smoothing constraint is
   satisfied and the algorithm terminates.

4. **Add knots.** Identify the knot interval \f$ [t_j, t_{j+1}) \f$ with the
   largest contribution to the residual. Insert a new knot at the data abscissa
   in that interval where the local residual is greatest.

5. **Iterate.** Return to step 2 with the augmented knot vector.

6. **Smoothing parameter update.** Once enough knots are available, the
   algorithm uses rational interpolation (Dierckx's `fprati`) to refine the
   Lagrange multiplier \f$ p \f$ that connects the unconstrained minimization

   \f[
       \eta(p) = p \sum_{i=1}^{m} w_i^2 (y_i - s_p(x_i))^2
       + (1 - p) \int [s_p''(x)]^2 \, dx
   \f]

   to the constrained formulation. This avoids solving many full least-squares
   problems and typically converges in a few iterations.

The procedure is designed to produce a spline with as few knots as possible
while satisfying the smoothing constraint \f$ f_p \leq S \f$.

@see Dierckx, Ch. 5, &sect;5.2.4 (pp. 78&ndash;84)

## Choosing the Smoothing Factor

Selecting an appropriate value of \f$ S \f$ is the most important practical
decision in spline smoothing. The following guidelines apply:

- **Known noise variance.** If the data satisfy \f$ y_i = g(x_i) + \varepsilon_i \f$
  where \f$ \varepsilon_i \f$ are independent errors with variance \f$ \sigma^2 \f$
  and the weights are \f$ w_i = 1 \f$, the expected value of the residual sum
  of squares for the true function is \f$ m \sigma^2 \f$. A natural choice is

  \f[
      S \approx m \, \sigma^2
  \f]

  More generally, for non-uniform weights, \f$ S \approx \sum w_i^2 \sigma_i^2 \f$.

- **Unit weights, well-scaled data.** When \f$ w_i = 1 \f$ and no noise estimate
  is available, \f$ S = m \f$ is a reasonable starting point.

- **Over-fitting (\f$ S \f$ too small).** The spline chases noise in the data,
  introducing spurious oscillations and using too many knots. The number of
  knots grows until the maximum is reached, and the fit may become numerically
  ill-conditioned.

- **Under-fitting (\f$ S \f$ too large).** The spline is too smooth and fails
  to capture genuine features. The knot count stays near the minimum \f$ 2(k+1) \f$,
  and the residuals show systematic patterns.

- **Diagnostic.** After fitting, inspect the mean squared error via
  `curve\%mse()`. If the residual is much larger than the expected noise level,
  decrease \f$ S \f$; if the spline oscillates visibly, increase \f$ S \f$.

A practical workflow is to try a sequence of \f$ S \f$ values (e.g., on a
logarithmic scale) and select the one that gives the best balance between
smoothness and fidelity, or use cross-validation techniques.

@see Dierckx, Ch. 5, &sect;5.2.4 (pp. 78&ndash;84)

## Periodic Splines

When the data represent a function that is periodic on \f$ [x_1, x_m] \f$ (for
example, angular measurements as a function of time), the fitted spline should
satisfy periodicity constraints:

\f[
    s^{(j)}(x_1) = s^{(j)}(x_m), \quad j = 0, 1, \ldots, k-1
\f]

These \f$ k \f$ conditions ensure that the spline and its first \f$ k-1 \f$
derivatives match at the endpoints, producing a globally smooth periodic
function.

In terms of the B-spline representation, periodicity is enforced by wrapping
the knot vector and tying the first and last \f$ k \f$ coefficients:

\f[
    c_j = c_{n - k + j}, \quad j = 1, \ldots, k
\f]

This reduces the number of free coefficients from \f$ n \f$ to \f$ n - k \f$.
The modified least-squares problem is solved with the same Givens rotation
approach.

FITPACK's `percur` routine implements automatic-knot periodic curve fitting.
The corresponding type is `fitpack_periodic_curve`.

@see Dierckx, Ch. 6, &sect;6.1 (pp. 105&ndash;114)

## Parametric Curves

A parametric curve in \f$ \mathbb{R}^d \f$ is defined by \f$ d \f$ coordinate
functions of a common parameter \f$ u \f$:

\f[
    \mathbf{x}(u) = \bigl( s_1(u), \, s_2(u), \, \ldots, \, s_d(u) \bigr)
\f]

Given data points \f$ \mathbf{x}_i \in \mathbb{R}^d \f$ for \f$ i = 1, \ldots, m \f$,
the first step is to assign parameter values \f$ u_i \f$.

### Chord-Length Parameterization

The default choice is cumulative chord length:

\f[
    u_1 = 0, \quad u_{i+1} = u_i + \| \mathbf{x}_{i+1} - \mathbf{x}_i \|,
    \quad i = 1, \ldots, m-1
\f]

This ensures that parameter increments are proportional to the spacing between
successive data points, preventing bunching in regions of high curvature.

### Shared Knot Vector

All \f$ d \f$ component splines \f$ s_j(u) \f$ share the same knot vector
\f$ \mathbf{t} \f$ and degree \f$ k \f$. The smoothing problem becomes: minimize

\f[
    \sum_{j=1}^{d} \int \left[ s_j''(u) \right]^2 \, du
    \quad \text{subject to} \quad
    \sum_{i=1}^{m} \sum_{j=1}^{d} w_i^2 \left( x_{ij} - s_j(u_i) \right)^2 \leq S
\f]

FITPACK's `parcur` routine handles this by solving \f$ d \f$ coupled
least-squares problems with a common knot vector.

### Closed Curves

For closed parametric curves, the spline must satisfy periodic boundary
conditions in every component:

\f[
    s_j^{(l)}(u_1) = s_j^{(l)}(u_m), \quad l = 0, 1, \ldots, k-1, \quad j = 1, \ldots, d
\f]

This is implemented by `clocur` (type `fitpack_closed_curve`).

### Endpoint Constraints

In some applications, the tangent vector or higher derivatives at the curve
endpoints are prescribed. FITPACK's `concur` routine supports derivative
constraints of the form:

\f[
    s_j^{(l)}(u_1) = \alpha_{jl}, \quad s_j^{(l)}(u_m) = \beta_{jl}
\f]

for selected orders \f$ l \f$ and components \f$ j \f$, implemented by
`fitpack_constrained_curve`.

@see Dierckx, Ch. 9 (pp. 199&ndash;228)

## Convexity-Constrained Fitting

In many applications the fitted curve must be shape-preserving: convex data
should produce a convex spline, concave data a concave spline. FITPACK supports
local convexity constraints via the second derivative.

### B-Spline Coefficient Conditions

A spline \f$ s(x) \f$ of degree \f$ k \geq 2 \f$ is convex on an interval if
and only if \f$ s''(x) \geq 0 \f$ there. For B-splines, this translates into
conditions on the second-order divided differences of the coefficients:

\f[
    \Delta^2 c_j = c_{j+2} - 2 c_{j+1} + c_j \geq 0
\f]

for all \f$ j \f$ whose corresponding B-splines overlap the interval.
Concavity requires \f$ \Delta^2 c_j \leq 0 \f$.

### Per-Interval Convexity Flags

FITPACK allows the user to specify convexity constraints independently for each
data interval via a flag vector \f$ v_i \f$:

| Flag value | Meaning |
|------------|---------|
| \f$ +1 \f$ | The spline must be concave on \f$ [x_i, x_{i+1}] \f$ |
| \f$ -1 \f$ | The spline must be convex on \f$ [x_i, x_{i+1}] \f$ |
| \f$ 0 \f$  | No constraint (the spline is free) |

### Quadratic Programming Formulation

The constrained fitting problem becomes a quadratic program (QP): minimize the
roughness integral (or weighted residual) subject to the linear inequality
constraints on \f$ \Delta^2 c_j \f$ and the smoothing constraint
\f$ f_p \leq S \f$. FITPACK solves this QP iteratively, combining the adaptive
knot strategy with active-set management for the convexity constraints.

Two routines are provided:

- **`concon`** (type `fitpack_convex_curve`): automatic knot placement with
  convexity constraints. The algorithm adds knots as in `curfit` but enforces
  the shape constraints at each step.
- **`cocosp`**: fixed knot vector with convexity constraints. The user supplies
  the knots and the routine solves the constrained least-squares problem.

@see Dierckx, Ch. 8, &sect;8.3&ndash;8.4 (pp. 173&ndash;196)

## Summary

| Formulation | Routine | Type | Key parameter |
|-------------|---------|------|---------------|
| Automatic smoothing | `curfit` | `fitpack_curve` | Smoothing factor \f$ S \f$ |
| Interpolation | `curfit` | `fitpack_curve` | \f$ S = 0 \f$ |
| Fixed-knot least-squares | `curfit` | `fitpack_curve` | User knot vector |
| Periodic smoothing | `percur` | `fitpack_periodic_curve` | Smoothing factor \f$ S \f$ |
| Open parametric | `parcur` | `fitpack_parametric_curve` | Smoothing factor \f$ S \f$ |
| Closed parametric | `clocur` | `fitpack_closed_curve` | Smoothing factor \f$ S \f$ |
| Endpoint-constrained | `concur` | `fitpack_constrained_curve` | Derivative values |
| Convex (auto knots) | `concon` | `fitpack_convex_curve` | Convexity flags \f$ v_i \f$ |
| Convex (fixed knots) | `cocosp` | &mdash; | Convexity flags \f$ v_i \f$ |

@see @ref theory_bsplines
@see @ref tutorial_curve
@see @ref book_reference
