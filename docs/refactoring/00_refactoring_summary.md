# Fitpack Core Refactoring Summary

> **Style Guide**: Follow coding conventions in [CLAUDE.md](../../CLAUDE.md)

## Overview

The `src/fitpack_core.f90` file (~19,000 lines) was converted from legacy F77 code and contains several repeated patterns that could be extracted into reusable routines. This document summarizes the identified refactoring opportunities.

**Important constraint**: Public API/interfaces cannot be changed.

---

## Identified Refactoring Opportunities

### 1. QR Row Rotation Pattern (High Priority)

**Occurrences**: ~45 locations
**Difficulty**: Low
**Impact**: High

Repeated pattern for rotating a row into an upper triangular matrix using Givens transformations:

```fortran
do i=1,k1
    piv = h(i); if (equal(piv,zero)) cycle
    call fpgivs(piv,a(j,1),cos,sin)
    call fprota(cos,sin,yi,z(j))
    if (i<k1) call fprota(cos,sin,h(i+1:k1),a(j,2:k1-i+1))
end do
```

**Affected routines**: `fpcurf`, `fpclos`, `fpcons`, `fppara`, `fpperi`, `fpgrre`, `fpgrdi`, `fpgrsp`, `fptrnp`, `fptrpe`, `fprank`, `fpsurf`, `fpsphe`, `fppola`, and others.

**See**: `01_qr_row_rotation_refactor.md` for detailed plan.

---

### 2. Knot Interval Search (Medium-High Priority)

**Occurrences**: ~20 locations
**Difficulty**: Low
**Impact**: Medium

Repeated pattern for finding the knot interval containing a point:

```fortran
! search for knot interval t(l) <= x < t(l+1).
do while (xi>=t(l+1) .and. l/=nk1)
    l = l+1
end do
```

**Refactor proposal**: Create a utility function:

```fortran
pure function fp_find_knot_interval(t, n, x, l_start, l_max) result(l)
    real(FP_REAL), intent(in) :: t(n), x
    integer(FP_SIZE), intent(in) :: n, l_start, l_max
    integer(FP_SIZE) :: l
end function
```

**Locations**: Lines 3193, 4059, 4415, 4726, 4938, 8312, 9169, 14261, 17819, 17938, etc.

---

### 3. Consolidate `fptrnp` and `fptrpe` (Medium Priority)

**Occurrences**: 2 routines
**Difficulty**: Medium
**Impact**: Medium

These routines at lines 13846 and 13959 are nearly identical:
- `fptrnp`: Reduces (m+n-7) x (n-4) matrix to upper triangular form
- `fptrpe`: Same for cyclic bandmatrix, with extra `aa` array for periodic handling

**Refactor proposal**: Merge into single routine with `periodic` flag:

```fortran
pure subroutine fp_triangularize(m,mm,idim,n,nr,sp,p,b,z,a,q,right,periodic,aa)
    logical(FP_BOOL), intent(in) :: periodic
    real(FP_REAL), intent(out), optional :: aa(n,4)  ! only for periodic
```

---

### 4. Back-Substitution Wrapper (Medium Priority)

**Occurrences**: ~25 calls
**Difficulty**: Low
**Impact**: Medium

Two back-substitution functions are used based on matrix structure:
- `fpback`: Standard upper triangular solve
- `fpbacp`: Periodic/cyclic matrix solve (uses both `a1` and `a2` blocks)

**Refactor proposal**: Generic interface or wrapper:

```fortran
interface fp_backsolve
    module procedure fpback
    module procedure fpbacp
end interface
```

Or a unified routine that handles both cases internally.

---

### 5. B-spline Evaluation + Weight Pattern (Medium Priority)

**Occurrences**: ~10 locations
**Difficulty**: Low
**Impact**: Low-Medium

Repeated pattern:

```fortran
h = fpbspl(t,n,k,xi,l)
q(it,1:k1) = h(1:k1)
h(:k1) = wi*h(:k1)
```

**Refactor proposal**: Either a combined function `fpbspl_weighted` or a helper that stores and weights in one call.

---

### 6. Consolidate Similar Fitting Routines (Lower Priority)

**Occurrences**: 6+ major routines
**Difficulty**: High
**Impact**: High (long-term maintainability)

Several routines share 80%+ identical logic:
- `fpcurf` (line 4794) and `fpperi` (line 8964) - 1D curve fitting
- `fpcons` (line 3899) and `fppara` (line 8176) - constrained/parametric curves
- `fpgrdi`, `fpgrre`, `fpgrsp` - grid fitting variants

**Refactor proposal**: Extract common core with configuration parameters. This is a larger undertaking requiring careful analysis of the differences.

---

### 7. Named Constants for Loop Limits (Quick Win)

**Occurrences**: Many
**Difficulty**: Trivial
**Impact**: Low (readability)

Some magic numbers could be named constants:
- Bandwidth values `4`, `5` → `BAND_CUBIC`, `BAND_QUINTIC`
- `MAX_ORDER+1` → `BSPLINE_BUFFER_SIZE`

---

### 8. Simplify Strided Array Access (Medium Priority)

**Occurrences**: ~30 locations
**Difficulty**: Medium
**Impact**: Medium (readability)

Patterns like:

```fortran
z(j:j+(idim-1)*n:n)  ! strided access across dimensions
```

Could use helper subroutines for clearer index computation.

---

## Recommended Execution Order

1. **Knot Interval Search** - Quick win, isolated change, enables binary search optimization
2. **QR Row Rotation** - Highest value, lowest risk
3. **Named Constants** - Trivial, improves readability
4. **Back-Substitution Wrapper** - Clean abstraction
5. **B-spline + Weight Pattern** - Minor improvement
6. **Merge `fptrnp`/`fptrpe`** - Moderate effort
7. **Similar Fitting Routines** - Major refactor, do last

**See**: `02_knot_interval_search_refactor.md` for detailed plan on item 1.

---

## Testing Strategy

All refactors must:
1. Pass the existing test suite (`fpm test`)
2. Maintain `pure` procedure attributes where present
3. Preserve numerical results to machine precision
4. Not change any public API signatures
