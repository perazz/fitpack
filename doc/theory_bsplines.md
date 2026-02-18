# B-Spline Foundations {#theory_bsplines}

This page introduces the mathematical foundations of B-splines as used throughout
FITPACK. All curve and surface fitting routines in the library represent their
results as linear combinations of B-spline basis functions, so understanding these
foundations is essential for interpreting knot vectors, coefficients, and the
various evaluation, differentiation, and integration operations.

The presentation follows Dierckx (1993), Chapters 1 and 6.

## Piecewise Polynomials and Knot Vectors

A **spline** of degree \f$ k \f$ is a piecewise polynomial function whose pieces
join with \f$ k - 1 \f$ continuous derivatives at designated breakpoints called
**knots**. The knots are collected into a non-decreasing **knot vector**

\f[
    t_0 \leq t_1 \leq \cdots \leq t_{n+k+1}
\f]

which partitions the approximation interval \f$ [a, b] \f$ into sub-intervals.
The number of distinct interior knots determines the flexibility of the spline:
more knots allow closer approximation to complicated shapes, while fewer knots
enforce smoothness.

The set of all splines of degree \f$ k \f$ on a given knot vector \f$ \mathbf{t} \f$
forms a linear space \f$ \mathcal{S}(k, \mathbf{t}) \f$ of dimension \f$ n + 1 \f$,
where \f$ n + 1 \f$ is the number of B-spline basis functions (and therefore the
number of coefficients needed to represent any element of the space).

In FITPACK the first and last \f$ k + 1 \f$ knots are clamped to the boundary:

\f[
    t_0 = t_1 = \cdots = t_k = a, \qquad
    t_{n+1} = t_{n+2} = \cdots = t_{n+k+1} = b
\f]

so that the spline interpolates (rather than merely approximates) its coefficient
values at the endpoints.

## B-Spline Basis Functions

Every element of \f$ \mathcal{S}(k, \mathbf{t}) \f$ can be written uniquely as a
linear combination of **B-spline basis functions** \f$ N_{i,k}(x) \f$, which are
defined recursively by the **Cox&ndash;de Boor recurrence**.

### Degree zero

For \f$ k = 0 \f$ the basis functions are simple indicator functions:

\f[
    N_{i,0}(x) =
    \begin{cases}
        1 & \text{if } t_i \leq x < t_{i+1}, \\
        0 & \text{otherwise}.
    \end{cases}
\f]

### Recurrence for degree \f$ k \geq 1 \f$

Higher-degree basis functions are built from lower-degree ones:

\f[
    N_{i,k}(x) =
        \frac{x - t_i}{t_{i+k} - t_i} \, N_{i,k-1}(x)
      + \frac{t_{i+k+1} - x}{t_{i+k+1} - t_{i+1}} \, N_{i+1,k-1}(x)
\f]

where any fraction with a zero denominator is taken to be zero (this arises when
consecutive knots coincide). The routine `fpbspl` implements this recurrence.

### Key properties

The B-spline basis functions enjoy several properties that make them ideal for
numerical computation:

1. **Non-negativity**: \f$ N_{i,k}(x) \geq 0 \f$ for all \f$ x \f$.
2. **Local support**: \f$ N_{i,k}(x) = 0 \f$ outside the interval
   \f$ [t_i,\, t_{i+k+1}] \f$. Each basis function is nonzero on at most
   \f$ k + 1 \f$ knot spans, so spline evaluation is a local operation.
3. **Partition of unity**: At every point in \f$ [a, b] \f$ the basis functions
   sum to one: \f$ \sum_{i} N_{i,k}(x) = 1 \f$. This ensures that any convex
   combination of coefficients produces a spline value within the convex hull of
   those coefficients.
4. **Smoothness**: At a simple knot (multiplicity 1), \f$ N_{i,k} \f$ has
   \f$ k - 1 \f$ continuous derivatives. Each additional coincident knot reduces
   smoothness by one order.

@see Dierckx, Ch. 1, &sect;1.1&ndash;1.2

## Spline as a Linear Combination

Any spline \f$ s \in \mathcal{S}(k, \mathbf{t}) \f$ is represented as

\f[
    s(x) = \sum_{i=0}^{n} c_i \, N_{i,k}(x)
\f]

where the \f$ c_i \f$ are the **B-spline coefficients** (also called control
points or de Boor points). The entire fitting problem therefore reduces to
determining the \f$ n + 1 \f$ coefficients for a given knot vector.

In FITPACK the coefficient array is stored alongside the knot vector in each
derived type (e.g., `fitpack_curve`, `fitpack_surface`). Once computed by the
fitting routines, the same coefficients are used by all evaluation, derivative,
and integration routines.

## Evaluation: the de Boor Algorithm

Given a point \f$ x \in [t_j,\, t_{j+1}) \f$, the spline value \f$ s(x) \f$ can
be computed without explicitly forming all basis functions. The **de Boor
algorithm** is a triangular recurrence analogous to Horner's method for ordinary
polynomials.

Define the initial values

\f[
    d_i^{[0]} = c_i, \qquad i = j-k, \ldots, j
\f]

and compute for \f$ r = 1, \ldots, k \f$:

\f[
    d_i^{[r]} =
        \frac{x - t_i}{t_{i+k+1-r} - t_i} \, d_i^{[r-1]}
      + \frac{t_{i+k+1-r} - x}{t_{i+k+1-r} - t_i} \, d_{i-1}^{[r-1]},
    \qquad i = j-k+r, \ldots, j
\f]

After \f$ k \f$ stages the result is \f$ s(x) = d_j^{[k]} \f$. The cost is
\f$ \mathcal{O}(k^2) \f$ per evaluation point, independent of the total number
of knots. This is the core loop in `splev` and underlies every `eval` method
in the library.

@see Dierckx, Ch. 1, Eq. 1.11

## Derivatives

Differentiation of a B-spline representation yields another B-spline of one
degree lower. The derivative of \f$ s(x) = \sum_i c_i \, N_{i,k}(x) \f$ is

\f[
    s'(x) = \sum_{i=1}^{n} c_i^{(1)} \, N_{i,k-1}(x)
\f]

where the new coefficients are obtained by **differencing**:

\f[
    c_i^{(1)} = \frac{k}{t_{i+k} - t_i} \, (c_i - c_{i-1})
\f]

Higher-order derivatives follow by repeated application: the \f$ \nu \f$-th
derivative is a spline of degree \f$ k - \nu \f$ with coefficients obtained from
\f$ \nu \f$ successive differencing steps. This structure is exploited by `splder`,
which computes derivatives of any order up to \f$ k \f$.

Because each differentiation reduces the degree by one, at most \f$ k \f$
derivatives are nonzero. A cubic spline (\f$ k = 3 \f$) therefore has a
piecewise-constant third derivative and a zero fourth derivative, regardless of
the coefficients.

@see Dierckx, Ch. 1, &sect;1.3

## Integration

The definite integral of a spline over an interval can be expressed directly in
terms of its coefficients, without numerical quadrature. For a spline
\f$ s(x) = \sum_i c_i \, N_{i,k}(x) \f$ on the interval \f$ [a, b] \f$:

\f[
    \int_a^b s(x) \, dx = \sum_{i=0}^{n} c_i \, \frac{t_{i+k+1} - t_i}{k + 1}
\f]

This formula follows from the property that the integral of each basis function
satisfies

\f[
    \int_a^b N_{i,k}(x) \, dx = \frac{t_{i+k+1} - t_i}{k + 1}
\f]

when the knot vector is clamped at both ends. The routine `splint` implements
this formula, giving exact results (up to floating-point arithmetic) regardless
of the complexity of the spline.

@see Dierckx, Ch. 1, &sect;1.4

## Knot Insertion (Oslo Algorithm)

A fundamental operation on B-spline representations is **knot insertion**:
adding a new knot \f$ \bar{t} \f$ to the knot vector without changing the spline
function itself. The resulting representation has one more knot, one more
coefficient, and the same function values everywhere.

Suppose a knot \f$ \bar{t} \in [t_j,\, t_{j+1}) \f$ is inserted. The new
coefficients \f$ \bar{c}_i \f$ are obtained from the old ones by a **convex
combination**:

\f[
    \bar{c}_i =
    \begin{cases}
        c_i & \text{if } i \leq j - k, \\[4pt]
        \alpha_i \, c_i + (1 - \alpha_i) \, c_{i-1}
            & \text{if } j - k < i \leq j, \\[4pt]
        c_{i-1} & \text{if } i > j,
    \end{cases}
\f]

with weights

\f[
    \alpha_i = \frac{\bar{t} - t_i}{t_{i+k} - t_i}.
\f]

Because \f$ 0 \leq \alpha_i \leq 1 \f$, the new coefficients lie in the convex
hull of their neighbours&mdash;a key property that preserves shape.

Knot insertion is the basis of the **Oslo algorithm** (Cohen, Lyche, and
Schumaker, 1986), which generalizes the operation to insert several knots at
once. In FITPACK the `insert` routine performs single knot insertion and is used
both for refinement and for the zero-finding algorithm in `sproot`.

@see Dierckx, Ch. 6, &sect;6.2 (pp. 114&ndash;118)

## Tensor Products for Bivariate Splines

The B-spline framework extends naturally to two dimensions via **tensor
products**. Given knot vectors \f$ \mathbf{t}_x \f$ and \f$ \mathbf{t}_y \f$
with associated basis functions \f$ N_{i,k_x}(x) \f$ and \f$ M_{j,k_y}(y) \f$,
a bivariate spline is

\f[
    s(x, y) = \sum_{i=0}^{n_x} \sum_{j=0}^{n_y}
        c_{ij} \, N_{i,k_x}(x) \, M_{j,k_y}(y)
\f]

The coefficient array \f$ c_{ij} \f$ has \f$ (n_x + 1)(n_y + 1) \f$ entries and
is stored in column-major order in FITPACK.

Evaluation of a tensor-product spline decomposes into two sequences of
univariate evaluations: first evaluate in \f$ x \f$ for each row of coefficients,
then evaluate the resulting intermediate values in \f$ y \f$. This separability
is the reason tensor-product splines scale well to two dimensions.

The routines `bispev` (grid evaluation) and `bispeu` (scattered evaluation)
implement this strategy. Partial derivatives (`parder`, `pardeu`) and double
integrals (`dblint`) follow the same tensor-product decomposition.

@see Dierckx, Ch. 1, &sect;1.5

## Further Reading

The foundations described on this page underpin all fitting algorithms in the
library. The fitting problem itself&mdash;choosing knots and coefficients to
approximate given data&mdash;is covered separately.

@see @ref theory_curve_fitting
@see @ref book_reference

> P. Dierckx, *Curve and Surface Fitting with Splines*,
> Oxford University Press, 1993.
