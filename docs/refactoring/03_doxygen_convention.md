# FITPACK Doxygen + MathJax Documentation Convention

> Adapted from the [fortran-lapack](https://github.com/perazz/fortran-lapack) Doxygen style.

## Overview

All routines in `fitpack_core.F90` should have structured documentation headers using
Doxygen-compatible Fortran comment markers (`!>` / `!!`), with mathematical content rendered
via MathJax. Each routine's documentation should link back to the relevant equations and sections
in the book (*Curve and Surface Fitting with Splines*, P. Dierckx, Oxford Univ. Press, 1993).

**Key principle**: The original FITPACK routines contain extensive documentation (up to 200+ lines
for public API routines like `curfit` and `surfit`). This documentation must be **preserved in
full** — not shortened or discarded. The refactoring converts it from plain-text F77 comment style
to structured Doxygen format, rewriting mathematical expressions with MathJax and adding book
equation references, but retaining all original content: parameter descriptions, error codes,
usage guidance, restrictions, and references.

---

## Comment Syntax

- `!>` — First line of a documentation block (triggers Doxygen parsing)
- `!!` — Continuation lines within the same block
- Place the block **immediately before** the `subroutine` / `function` statement

---

## Document Structure

The overall structure for a routine's documentation block:

```fortran
!> @brief One-line summary.
!!
!! ### Description
!! Multi-paragraph description with math, preserving all original content.
!!
!! ### Matrix structure
!! Band structure, storage layout (when applicable).
!!
!! ### Parameters
!! @param[in]    x  Description of input parameter
!! @param[out]   y  Description of output parameter
!! @param[in,out] a Description of modified parameter
!!
!! ### Return value
!! @return Description of result (for functions only)
!!
!! ### Error flags
!! Detailed description of all error codes (for public API routines).
!!
!! ### Usage guidance
!! Original "further comments" section with advice on parameter choices.
!!
!! ### Restrictions
!! Input validity conditions.
!!
!! @note  Important usage notes
!! @warning Warnings about destructive operations
!!
!! ### References
!! @see Dierckx, Ch. N, §N.M (pp. XX-YY), Eq. N.M
!! @see Original paper references
pure subroutine my_routine(x, y, a)
```

Not all sections are required for every routine — use what is appropriate:

| Section               | Public API | Internal core | Low-level utility |
|-----------------------|:----------:|:-------------:|:-----------------:|
| `@brief`              | Required   | Required      | Required          |
| Description + math    | Required   | Required      | Required          |
| Matrix structure      | If applicable | If applicable | If applicable  |
| `@param`              | Required   | Required      | Required          |
| Error flags           | Required   | If applicable | —                 |
| Usage guidance        | Required   | —             | —                 |
| Restrictions          | Required   | If applicable | —                 |
| `@see` book refs      | Required   | Required      | Required          |

---

## Math Syntax

**Inline math**: `\f$ expression \f$`
```fortran
!! Computes the rotation parameters \f$ c \f$ and \f$ s \f$ such that
!! \f$ c^2 + s^2 = 1 \f$.
```

**Display equations**: `\f[ expression \f]`
```fortran
!! The Givens rotation eliminates element \f$ e_i \f$ by computing:
!! \f[
!!     r'_i = \sqrt{r_i^2 + e_i^2}, \quad
!!     c = \frac{r_i}{r'_i}, \quad
!!     s = \frac{e_i}{r'_i}
!! \f]
```

**Equation tags** (matching book numbering):
```fortran
!! \f[
!!     F_g(p) = \sum_{i=1}^{m} w_i \left( y_i - s_g(x_i; p) \right)^2
!!     \tag{5.1}
!! \f]
```

---

## Complete Examples

### Example 1: Public API routine (`curfit`) — full-length documentation

This shows how the original ~200-line F77 comment block gets refactored into Doxygen format.
All content is preserved; math is rewritten with MathJax; book equation references are added.

```fortran
!> @brief Smoothing spline curve fitting.
!!
!! ### Description
!!
!! Given data points \f$ (x_i, y_i) \f$ and positive weights \f$ w_i \f$,
!! \f$ i = 1, \ldots, m \f$, determines a smooth spline approximation
!! \f$ s(x) \f$ of degree \f$ k \f$ on the interval \f$ [x_b, x_e] \f$.
!!
!! If `iopt = -1`, computes the weighted least-squares spline for a given
!! set of knots. If `iopt >= 0`, the number and position of knots
!! \f$ \lambda_j \f$, \f$ j = 1, \ldots, n \f$, are chosen automatically.
!! Smoothness is achieved by minimizing the discontinuity jumps of the
!! \f$ k \f$-th derivative of \f$ s(x) \f$ at the interior knots, subject to
!! the constraint:
!!
!! \f[
!!     F(p) = \sum_{i=1}^{m} \left( w_i (y_i - s(x_i)) \right)^2 \leq S
!!     \tag{5.1}
!! \f]
!!
!! where \f$ S \geq 0 \f$ is a given smoothing factor.
!!
!! The fit is returned in B-spline representation (coefficients
!! \f$ c_j \f$, \f$ j = 1, \ldots, n-k-1 \f$) and can be evaluated using
!! subroutine `splev`.
!!
!! ### Parameters
!!
!! @param[in] iopt  Fit mode flag:
!!                  - `iopt = -1`: weighted least-squares spline with given knots.
!!                  - `iopt = 0`: smoothing spline; starts with initial knot set
!!                    \f$ \lambda_i = x_b \f$, \f$ \lambda_{i+k+1} = x_e \f$,
!!                    \f$ i = 1, \ldots, k+1 \f$.
!!                  - `iopt = 1`: continue with knots from last call.
!!                    **Must** be immediately preceded by a call with `iopt = 0` or `1`.
!! @param[in] m     Number of data points. Must satisfy \f$ m > k \f$.
!! @param[in] x     Abscissae \f$ x_i \f$, strictly ascending:
!!                  \f$ x_b \leq x_1 < x_2 < \cdots < x_m \leq x_e \f$.
!! @param[in] y     Ordinates \f$ y_i \f$, \f$ i = 1, \ldots, m \f$.
!! @param[in] w     Weights \f$ w_i > 0 \f$, \f$ i = 1, \ldots, m \f$.
!!                  See *Usage guidance* below for recommended values.
!! @param[in] xb,xe Boundaries of the approximation interval.
!!                  \f$ x_b \leq x_1 \f$, \f$ x_e \geq x_m \f$.
!! @param[in] k     Spline degree, \f$ 1 \leq k \leq 5 \f$. Cubic splines
!!                  (\f$ k = 3 \f$) are recommended. Even \f$ k \f$ with small
!!                  \f$ S \f$ is discouraged.
!! @param[in] s     Smoothing factor \f$ S \geq 0 \f$ (only used when `iopt >= 0`).
!!                  See *Usage guidance* below.
!! @param[in] nest  Over-estimate of the total number of knots. Must satisfy
!!                  `nest >= 2*k+2`. In practice, `nest = m/2` is usually sufficient;
!!                  `nest = m+k+1` is always sufficient (interpolation case \f$ S = 0 \f$).
!! @param[in,out] n Total number of knots. On exit (unless `ier = 10`), contains
!!                  the number of knots of the returned spline. If `iopt = 1`, the
!!                  value must be left unchanged between calls. If `iopt = -1`,
!!                  must be specified on entry.
!! @param[in,out] t Knot vector of dimension `nest`. On successful exit, contains
!!                  the knots: interior knots \f$ \lambda_{k+2}, \ldots, \lambda_{n-k-1} \f$
!!                  and boundary knots \f$ \lambda_1 = \cdots = \lambda_{k+1} = x_b \f$,
!!                  \f$ \lambda_{n-k} = \cdots = \lambda_n = x_e \f$.
!!                  If `iopt = -1`, interior knots must be supplied on entry.
!!                  If `iopt = 1`, must be left unchanged between calls.
!! @param[out]   c  B-spline coefficients \f$ c_1, \ldots, c_{n-k-1} \f$.
!! @param[out]   fp Weighted sum of squared residuals \f$ F(p) \f$.
!! @param[in,out] wrk  Workspace array of dimension `lwrk`.
!!                     If `iopt = 1`, `wrk(1:n)` must be left unchanged between calls.
!! @param[in] lwrk Dimension of `wrk`. Must satisfy
!!                 `lwrk >= (k+1)*m + nest*(7+3*k)`.
!! @param[in,out] iwrk Integer workspace of dimension `nest`.
!!                     If `iopt = 1`, `iwrk(1:n)` must be left unchanged between calls.
!! @param[out]   ier Error flag. See *Error flags* below.
!!
!! ### Error flags
!!
!! - `ier = 0`: Normal return. The spline satisfies
!!   \f$ |F(p) - S| / S \leq \tau \f$ with \f$ \tau = 0.001 \f$.
!! - `ier = -1`: Normal return. Interpolating spline (\f$ F(p) = 0 \f$).
!! - `ier = -2`: Normal return. The spline is the weighted least-squares
!!   polynomial of degree \f$ k \f$. In this case, `fp` gives the upper
!!   bound \f$ F_0 \f$ for the smoothing factor.
!! - `ier = 1`: Error. Required storage exceeds `nest`. Probable cause:
!!   `nest` too small, or \f$ S \f$ too small. The returned spline uses the
!!   current knot set; `fp` gives the residual sum (\f$ F(p) > S \f$).
!! - `ier = 2`: Error. Theoretically impossible result during iteration
!!   for \f$ F(p) = S \f$. Probable cause: \f$ S \f$ too small.
!! - `ier = 3`: Error. Maximum iterations (`maxit = 20`) reached.
!!   Probable cause: \f$ S \f$ too small.
!! - `ier = 10`: Error. Invalid input. See *Restrictions* below.
!!
!! ### Usage guidance
!!
!! The parameter \f$ S \f$ controls the trade-off between closeness of fit
!! and smoothness. If \f$ S \f$ is too large, signal will be lost; if too
!! small, the spline will pick up noise. In the extreme cases:
!! - \f$ S = 0 \f$: interpolating spline
!! - \f$ S \to \infty \f$: weighted least-squares polynomial of degree \f$ k \f$
!!
!! If the weights are taken as \f$ w_i = 1/d_i \f$ where \f$ d_i \f$ is an
!! estimate of the standard deviation of \f$ y_i \f$, a good \f$ S \f$ should
!! lie in the range \f$ (m - \sqrt{2m},\; m + \sqrt{2m}) \f$.
!!
!! If nothing is known about the statistical error, set all \f$ w_i = 1 \f$
!! and determine \f$ S \f$ by trial and error: start with a very large
!! \f$ S \f$ (to obtain the least-squares polynomial and its upper bound
!! \f$ F_0 \f$), then progressively decrease (e.g., \f$ S = F_0/10, F_0/100,
!! \ldots \f$).
!!
!! To economize the search, use `iopt = 0` for the first call and
!! `iopt = 1` for subsequent calls with different \f$ S \f$ values. Once a
!! satisfactory fit is found, a final call with `iopt = 0` may return
!! an equally good fit with fewer knots.
!!
!! ### Restrictions
!!
!! - `-1 <= iopt <= 1`, `1 <= k <= 5`, `m > k`, `nest >= 2*k+2`
!! - \f$ w_i > 0 \f$, \f$ x_b \leq x_1 < x_2 < \cdots < x_m \leq x_e \f$
!! - `lwrk >= (k+1)*m + nest*(7+3*k)`
!! - If `iopt = -1`: `2*k+2 <= n <= min(nest, m+k+1)`,
!!   \f$ x_b < \lambda_{k+2} < \cdots < \lambda_{n-k-1} < x_e \f$,
!!   and the Schoenberg-Whitney conditions must hold (Eq. 4.5).
!! - If `iopt >= 0`: \f$ S \geq 0 \f$; if \f$ S = 0 \f$ then `nest >= m+k+1`.
!!
!! ### References
!!
!! @see Dierckx, Ch. 4, §4.1-4.2 (pp. 53-73): least-squares spline fitting
!! @see Dierckx, Ch. 5, §5.1-5.3 (pp. 75-94): smoothing spline theory
!! @see Dierckx, Ch. 13, §13.1 (pp. 249-255): `curfit` routine description
!! @see P. Dierckx, "An algorithm for smoothing, differentiation and integration
!!      of experimental data using spline functions",
!!      J. Comp. Appl. Math. 1 (1975), 165-184.
!! @see P. Dierckx, "An improved algorithm for curve fitting with spline
!!      functions", Report TW54, K.U. Leuven, 1981.
pure subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
```

### Example 2: Low-level utility (`fpgivs`)

```fortran
!> @brief Compute parameters of a Givens plane rotation.
!!
!! Given a pivot element \f$ e_i \f$ and a diagonal element \f$ r_i \f$, computes
!! the cosine \f$ c \f$ and sine \f$ s \f$ of a Givens rotation that eliminates
!! \f$ e_i \f$:
!!
!! \f[
!!     r'_i = \sqrt{r_i^2 + e_i^2}, \quad
!!     c = \frac{r_i}{r'_i}, \quad
!!     s = \frac{e_i}{r'_i}
!!     \tag{4.15}
!! \f]
!!
!! The rotation preserves orthogonality (\f$ c^2 + s^2 = 1 \f$) and is used
!! to triangularize the observation matrix \f$ E \f$ row by row during the
!! QR decomposition of the least-squares system.
!!
!! @param[in]     piv  Pivot element \f$ e_i \f$ to be eliminated
!! @param[in,out] ww   On entry, diagonal element \f$ r_i \f$.
!!                     On exit, updated value \f$ r'_i = \sqrt{r_i^2 + e_i^2} \f$.
!! @param[out]    cos  Cosine of rotation angle
!! @param[out]    sin  Sine of rotation angle
!!
!! @see Dierckx, Ch. 4, §4.1.2 (pp. 55-58), Eq. 4.15
elemental subroutine fpgivs(piv,ww,cos,sin)
```

### Example 3: Back-substitution (`fpback`)

```fortran
!> @brief Solve an upper triangular banded system by back-substitution.
!!
!! Computes the solution \f$ c \f$ of the triangular system:
!!
!! \f[
!!     R_1 \, c = z_1
!!     \tag{4.14}
!! \f]
!!
!! where \f$ R_1 \f$ is an \f$ n \times n \f$ upper triangular matrix with
!! bandwidth \f$ k \f$, obtained from the QR factorization of the observation
!! matrix \f$ E \f$.
!!
!! ### Matrix structure
!!
!! \f$ R_1 \f$ is stored as `a(nest, k)` where `a(i,1)` holds the diagonal
!! element \f$ r_{i,i} \f$ and `a(i,j)` for \f$ j > 1 \f$ holds \f$ r_{i,i+j-1} \f$.
!! See Dierckx, Ch. 4, Fig. 4.2.
!!
!! @param[in] a     Upper triangular band matrix \f$ R_1 \f$, stored as `a(nest, k)`
!! @param[in] z     Right-hand side vector \f$ z_1 \f$ of length \f$ n \f$
!! @param[in] n     Number of equations (= number of B-spline coefficients)
!! @param[in] k     Bandwidth of \f$ R_1 \f$ (= spline order \f$ k+1 \f$)
!! @param[in] nest  Leading dimension of array `a`
!!
!! @return Solution vector \f$ c \f$ of length \f$ n \f$
!!
!! @see Dierckx, Ch. 4, §4.1.2 (pp. 55-58), Eq. 4.14
!! @see Dierckx, Ch. 5, §5.2.2 (pp. 76-79), Eq. 5.15
pure function fpback(a,z,n,k,nest) result(c)
```

### Example 4: Refactored helper (`fp_rotate_row`)

```fortran
!> @brief Rotate a row into an upper triangular band matrix using Givens rotations.
!!
!! Eliminates all non-zero elements of row \f$ h(1:\text{band}) \f$ by applying
!! successive Givens rotations against the diagonal of \f$ R_1 \f$. Each rotation
!! also updates the right-hand side vector \f$ z \f$ with the corresponding
!! scalar contribution \f$ y_i \f$.
!!
!! This is the inner loop of the QR factorization used throughout FITPACK
!! (Eq. 4.15). It encapsulates the pattern:
!!
!! ```
!! do i = 1, band
!!     call fpgivs(h(i), a(j,1), c, s)
!!     call fprota(c, s, yi, z(j))
!!     call fprota(c, s, h(i+1:), a(j,2:))
!! end do
!! ```
!!
!! @param[in,out] h        Row to rotate, length `band`. Modified in place.
!! @param[in]     band     Number of elements in the row (typically \f$ k+1 \f$)
!! @param[in,out] a        Upper triangular band matrix \f$ R_1 \f$
!! @param[in,out] yi       Scalar RHS contribution. Modified by rotations.
!! @param[in,out] z        RHS vector. Elements `z(j_start+1 : j_start+band)` modified.
!! @param[in]     j_start  Starting row index minus 1 (j increments before use)
!!
!! @see Dierckx, Ch. 4, §4.1.2 (pp. 55-58), Eq. 4.15
pure subroutine fp_rotate_row(h, band, a, yi, z, j_start)
```

---

## Doxygen Build Configuration

Use the same setup as [fortran-lapack](https://github.com/perazz/fortran-lapack):

- **Theme**: [doxygen-awesome-css](https://github.com/jothepro/doxygen-awesome-css) with dark mode toggle
- **MathJax**: Version 2.7.7 from CDN
- **Key Doxyfile settings**:

```
OPTIMIZE_FOR_FORTRAN   = YES
USE_MATHJAX            = YES
MATHJAX_RELPATH        = https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/
FORMULA_FONTSIZE       = 16
HTML_EXTRA_STYLESHEET  = doxygen-awesome.css
HTML_EXTRA_FILES       = doxygen-awesome-darkmode-toggle.js
GENERATE_TREEVIEW      = YES
```

---

## Checklist

When documenting a routine, verify:

- [ ] `@brief` is concise (one line)
- [ ] All original documentation content is preserved (not shortened)
- [ ] Mathematical expressions rewritten with `\f$ \f$` / `\f[ \f]`
- [ ] Book equation numbers added as `\tag{N.M}` to display equations
- [ ] Matrix structure described if routine operates on band matrices
- [ ] All parameters documented with `@param[direction]`
- [ ] Error flags fully described (public API routines)
- [ ] Usage guidance preserved (public API routines)
- [ ] `@see` references specific book chapter, section, page range, and equation numbers
- [ ] Math renders correctly (no unescaped underscores in non-math context)
