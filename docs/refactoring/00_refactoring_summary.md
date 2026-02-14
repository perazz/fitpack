# Fitpack Core Refactoring Summary

> **Style Guide**: Follow coding conventions in [CLAUDE.md](../../CLAUDE.md)

## Overview

The `src/fitpack_core.F90` file (~19,000 lines) was converted from legacy F77 code and contains
repeated patterns that can be extracted into reusable routines. The full mathematical context is
now available from the transcribed book (*Curve and Surface Fitting with Splines*, P. Dierckx, 1993)
in `docs/book/` and papers in `docs/papers/`.

**Important constraints**:
- Public API/interfaces cannot change
- Each refactoring step = one PR, fully tested against previous version
- Keep the original routine structure; no merging/consolidation of routines

---

## Completed Work

### PR 1: Knot Interval Search (DONE)

**Occurrences**: ~20 locations replaced
**New routine**: `fp_knot_interval` (hybrid linear/binary search)

### PR 2: Named Constants (DONE)

**Occurrences**: ~10 replacements
**Changes**: Magic numbers → `MAX_ORDER`, `DEGREE_3`, `DEGREE_5`, etc.

### PR 3: QR Row Rotation — Step 1 (DONE)

**Changes**: Added `fp_rotate_row` subroutine; replaced 1 occurrence in `fpcurf`

### PR 4: Complete Variant A Givens Rotations (DONE)

**Changes**: Replaced all ~15 scalar-RHS rotation patterns with `fp_rotate_row`
in `fpcurf`, `fpcons`, `fppara`, `fprank`

### PR 5: Variant B Standard Givens Rotations (DONE)

**Changes**: Added `fp_rotate_row_vec` subroutine for vector-RHS rotations
with explicit-shape arrays `h(band)`, `a(nest,*)`, `xi(idim)`, `z(n,idim)`.
Replaced 3 standard Variant B call sites in `fpclos`, `fpcons`, `fppara`.
Removed unused variables `i1`, `i3` from `fpcons` and `fppara`.

### PR 6: Shifting-Pivot and Two-Matrix Givens Rotations (DONE)

**Changes**: Added `fp_rotate_shifted_vec` (single-matrix shifted-pivot rotation)
and `fp_rotate_row_2mat_vec` (two-matrix rotation for periodic splines).
Replaced 4 remaining shifting-pivot call sites: 2 in `fpclos` (two-matrix periodic),
1 in `fpcons`, 1 in `fppara` (single-matrix smoothing iterations).
Removed unused variables: `piv`, `ij` from `fpclos`; `cos`, `sin`, `piv`, `i2`
from `fpcons` and `fppara`.

### PR 7: Scalar Shifted-Pivot and Two-Matrix Givens Rotations (DONE)

**Changes**: Added `fp_rotate_shifted` (scalar-RHS shifted-pivot rotation) and
`fp_rotate_row_2mat` (scalar-RHS two-matrix rotation for periodic splines).
Replaced 1 standard-walk site in `fpsurf` with `fp_rotate_row` and 7 shifted-pivot
sites: 1 in `fpcurf`, 2 in `fpsurf`, 2 in `fppola`, 2 in `fpsphe`.
Replaced 2 two-matrix sites in `fpperi` (observation and smoothing matrices).
Removed unused variables: `cos`, `sin`, `piv`, `i2` from `fpcurf`;
`cos`, `sin`, `piv`, `irot`, `i2` from `fpsurf`; `i2` from `fppola`;
`piv`, `ij` from `fpperi`.

---

## Remaining Plan

Each item below is a separate PR. PRs are ordered by dependency; documentation (PR 9) can
proceed in parallel with any code PR.

### PR 8: Variant D Givens Rotations

**Scope**: Grid/surface fitting rotation patterns
**New routine**: `fp_rotate_row_grid` (2D RHS)
**Occurrences**: ~15 locations in `fpgrdi`, `fpgrre`, `fpgrsp`, `fptrnp`, `fptrpe`,
`fpsurf`, `fpsphe`, `fppola`
**Difficulty**: Medium — 2D indexing and varying loop structures
**See**: [01_qr_row_rotation_refactor.md](01_qr_row_rotation_refactor.md), Variant D

---

### PR 9: Back-Substitution Interface

**Scope**: Unify `fpback` and `fpbacp` under a generic `fp_backsolve` interface
**Occurrences**: ~25 call sites
**Difficulty**: Low — clean abstraction, no logic changes
**Book references**: Ch. 4 Eq. 4.14 (R₁·c = z₁), Ch. 5 Eq. 5.15 (R₁*·c = z₁*)

---

### PR 10+: Doxygen + MathJax Documentation

**Scope**: Add structured documentation headers to all core routines in `fitpack_core.F90`
**Format**: See [03_doxygen_convention.md](03_doxygen_convention.md)
**Content sources**: Book chapters 1-13, papers in `docs/papers/`
**Approach**: Incremental — document routines as they are touched by code PRs,
plus dedicated documentation-only PRs for untouched routines
**Key routines** (with book references):

| Routine    | Book Section               | Key Equations     | Description                              |
|------------|----------------------------|-------------------|------------------------------------------|
| `fpgivs`   | §4.1.2, pp. 55-58          | 4.15, 4.16        | Givens rotation parameters               |
| `fprota`   | §4.1.2, pp. 55-58          | 4.15              | Apply Givens rotation                    |
| `fpback`   | §4.1.2, pp. 55-58          | 4.14              | Back-substitution, bandwidth k+1         |
| `fpbacp`   | §6.1, pp. 95-100           | 6.6               | Back-substitution, periodic (a1+a2)      |
| `fpbspl`   | §1.3, pp. 8-11             | 1.22, 1.30        | B-spline evaluation (de Boor-Cox)        |
| `fpcurf`   | §4.1-5.3, pp. 53-94        | 4.12, 5.10, 5.37  | Core curve fitting algorithm             |
| `fpdisc`   | §5.2.2, pp. 76-79          | 5.5, 5.6          | Discontinuity jumps of k-th derivative   |
| `fpknot`   | §5.3, pp. 87-94            | 5.37-5.43         | Adaptive knot placement                  |
| `fprati`   | §5.2.4, pp. 83-86          | 5.30-5.32         | Rational interpolation for F(p) = S      |
| `fprank`   | §9.1.2, pp. 150-152        | 9.8-9.10          | Rank-deficient system solution           |
| `fpchec`   | §4.1.1, pp. 53-55          | 4.5, 4.17         | Schoenberg-Whitney condition check       |
| `fpgrdi`   | §10.2, pp. 170-172         | 10.4-10.8         | Grid fitting (Kronecker decomposition)   |
| `fpgrre`   | §10.2, pp. 170-172         | 10.4-10.8         | Grid fitting (rectangular)               |
| `fptrnp`   | §10.2, pp. 170-172         | 10.8              | Triangularize observation matrix         |
| `fpsurf`   | §9.1-9.2, pp. 147-167      | 9.2-9.6           | Surface fitting to scattered data        |
| `fpsphe`   | §11.2, pp. 205-213         | 11.12-11.16       | Spherical coordinate fitting             |
| `fppola`   | §11.1, pp. 197-205         | 11.1-11.9         | Polar coordinate fitting                 |

---

## Testing Strategy

All PRs must:
1. Pass the full test suite (`fpm test`) — 49 tests, 0 failures
2. Maintain `pure` procedure attributes where present
3. Preserve numerical results to machine precision
4. Not change any public API signatures

---

## Workflow

1. Each PR branches from `main` (or the previous PR's merge)
2. Run `fpm test` before and after each change
3. Small, focused commits within each PR
4. Squash-merge to `main` for clean history
