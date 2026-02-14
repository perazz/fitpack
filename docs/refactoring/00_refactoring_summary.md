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

---

## Remaining Plan

Each item below is a separate PR. PRs are ordered by dependency; documentation (PR 8) can
proceed in parallel with any code PR.

### PR 4: Complete Variant A Givens Rotations

**Scope**: Replace all scalar-RHS rotation patterns with `fp_rotate_row`
**Occurrences**: ~15 locations in `fpcurf`, `fpcons`, `fppara`, `fprank`
**Difficulty**: Low — mechanical replacement, same signature as existing helper
**Variant pattern**:
```fortran
do i=1,k1
    j = j+1
    piv = h(i); if (equal(piv,zero)) cycle
    call fpgivs(piv,a(j,1),cos,sin)
    call fprota(cos,sin,yi,z(j))
    if (i<k1) call fprota(cos,sin,h(i+1:k1),a(j,2:k1-i+1))
end do
```
**See**: [01_qr_row_rotation_refactor.md](01_qr_row_rotation_refactor.md), Variant A

---

### PR 5: Reshape Array Indexing in `fpcurf`

**Scope**: Clean up 1D array stride patterns in the main curve fitting routine
**Goal**: Make band matrix access and multi-dimensional RHS indexing explicit,
preparing the ground for cleaner Variant B/C/D rotation replacements
**Difficulty**: Medium — requires careful index verification
**Key patterns to address**:
- Strided RHS access: `z(j:j+(idim-1)*n:n)` → clearer views
- Any remaining manual offset arithmetic for band matrix columns

**Book references**: Ch. 4 Fig. 4.2 (band structure of E and R₁),
Ch. 5 Fig. 5.1 (band structure of P and R₁*)

---

### PR 6: Variant B/C Givens Rotations

**Scope**: Vector RHS and two-matrix periodic variants
**New routines**: `fp_rotate_row_vec` (vector RHS with stride),
`fp_rotate_row_2mat` (two-matrix rotation for periodic splines)
**Occurrences**: ~15 locations in `fpclos`, `fpcons`, `fppara`, `fpperi`
**Difficulty**: Medium — strided access and secondary matrix add complexity
**See**: [01_qr_row_rotation_refactor.md](01_qr_row_rotation_refactor.md), Variants B & C

---

### PR 7: Variant D Givens Rotations

**Scope**: Grid/surface fitting rotation patterns
**New routine**: `fp_rotate_row_grid` (2D RHS)
**Occurrences**: ~15 locations in `fpgrdi`, `fpgrre`, `fpgrsp`, `fptrnp`, `fptrpe`,
`fpsurf`, `fpsphe`, `fppola`
**Difficulty**: Medium — 2D indexing and varying loop structures
**See**: [01_qr_row_rotation_refactor.md](01_qr_row_rotation_refactor.md), Variant D

---

### PR 8: Back-Substitution Interface

**Scope**: Unify `fpback` and `fpbacp` under a generic `fp_backsolve` interface
**Occurrences**: ~25 call sites
**Difficulty**: Low — clean abstraction, no logic changes
**Book references**: Ch. 4 Eq. 4.14 (R₁·c = z₁), Ch. 5 Eq. 5.15 (R₁*·c = z₁*)

---

### PR 9+: Doxygen + MathJax Documentation

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
