# Book Reference Index {#book_reference}

Quick-reference table mapping FITPACK routines to chapters, sections, and equations
in the primary reference:

> P. Dierckx, *Curve and Surface Fitting with Splines*,
> Oxford University Press, 1993.

## Core Algorithms

| Routine | Book Reference | Description |
|---------|---------------|-------------|
| `fpbspl` | Ch. 1, Eq. 1.11 | B-spline basis evaluation (de Boor&ndash;Cox recurrence) |
| `fprati` | Ch. 5, &sect;5.2.4 | Rational interpolation for smoothing parameter update |
| `fpdisc` | Ch. 4, &sect;4.3 | Discontinuity jumps of B-spline derivatives at knots |
| `fpgivs` | Ch. 4, &sect;4.2 | Givens rotation computation |
| `fprota` | Ch. 4, &sect;4.2 | Apply Givens rotation to matrix rows |
| `fpback` | Ch. 4, &sect;4.2 | Back-substitution for upper triangular systems |
| `fprank` | Ch. 4, &sect;4.3 | Rank-deficient least-squares via Householder |
| `fporde` | Ch. 5, &sect;5.3 | Assign scattered data to knot intervals |
| `fpchec` | Ch. 4, &sect;4.1 | Check knot sequence validity |
| `fpchep` | Ch. 4, &sect;4.1 | Check periodic knot sequence |
| `fppara` | Ch. 9, &sect;9.1 | Compute chord-length parameterization |

## Curve Fitting (1D)

| Routine | Book Reference | Description |
|---------|---------------|-------------|
| `curfit` | Ch. 5, &sect;5.2 (pp. 67&ndash;84) | Automatic-knot curve fitting |
| `percur` | Ch. 6, &sect;6.1 (pp. 105&ndash;114) | Periodic curve fitting |
| `parcur` | Ch. 9, &sect;9.1 (pp. 199&ndash;212) | Open parametric curve fitting |
| `clocur` | Ch. 9, &sect;9.2 (pp. 213&ndash;216) | Closed parametric curve fitting |
| `concur` | Ch. 9, &sect;9.3 (pp. 217&ndash;228) | Constrained parametric curve fitting |
| `concon` | Ch. 8, &sect;8.3 (pp. 173&ndash;188) | Convexity-constrained fitting (auto knots) |
| `cocosp` | Ch. 8, &sect;8.4 (pp. 189&ndash;196) | Convexity-constrained fitting (given knots) |
| `insert` | Ch. 6, &sect;6.2 (pp. 114&ndash;118) | Knot insertion (Oslo algorithm) |

## Curve Evaluation

| Routine | Book Reference | Description |
|---------|---------------|-------------|
| `splev` | Ch. 1, Eq. 1.11 | Evaluate spline at given points |
| `splder` | Ch. 1, &sect;1.3 | Evaluate spline derivatives |
| `splint` | Ch. 1, &sect;1.4 | Definite integral of spline |
| `sproot` | Ch. 6, &sect;6.3 (pp. 118&ndash;123) | Zeros of a cubic spline |
| `fourco` | Ch. 7, &sect;7.3 (pp. 153&ndash;165) | Fourier coefficients of spline |
| `curev` | Ch. 9 | Evaluate parametric spline curve |
| `cualde` | Ch. 9 | All derivatives of parametric curve |

## Surface Fitting (2D)

| Routine | Book Reference | Description |
|---------|---------------|-------------|
| `surfit` | Ch. 5, &sect;5.3 (pp. 85&ndash;98) | Scattered bivariate fitting |
| `regrid` | Ch. 5, &sect;5.4 (pp. 98&ndash;103) | Gridded bivariate fitting |
| `parsur` | Ch. 10, &sect;10.2 (pp. 241&ndash;254) | Parametric surface fitting |

## Surface Evaluation

| Routine | Book Reference | Description |
|---------|---------------|-------------|
| `bispev` | Ch. 1, &sect;1.5 | Evaluate bivariate spline on grid |
| `bispeu` | Ch. 1, &sect;1.5 | Evaluate bivariate spline at scattered points |
| `parder` | Ch. 1, &sect;1.5 | Partial derivatives on grid |
| `pardeu` | Ch. 1, &sect;1.5 | Partial derivatives at scattered points |
| `pardtc` | &mdash; | Transform coefficients for derivative spline |
| `dblint` | Ch. 1, &sect;1.5 | Double integral over rectangle |
| `profil` | &mdash; | Cross-section at fixed x or y |
| `surev`  | Ch. 10 | Evaluate parametric surface on grid |
| `evapol` | &mdash; | Evaluate polar-domain spline |

## Polar and Spherical Domains

| Routine | Book Reference | Description |
|---------|---------------|-------------|
| `polar` | Ch. 11, &sect;11.1 (pp. 255&ndash;263) | Scattered polar fitting |
| `pogrid` | Ch. 11, &sect;11.1 (pp. 255&ndash;263) | Gridded polar fitting |
| `sphere` | Ch. 11, &sect;11.2 (pp. 263&ndash;269) | Scattered spherical fitting |
| `spgrid` | Ch. 11, &sect;11.2 (pp. 263&ndash;269) | Gridded spherical fitting |
