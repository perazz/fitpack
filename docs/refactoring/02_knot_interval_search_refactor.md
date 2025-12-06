# Refactor Plan: Knot Interval Search

> **Style Guide**: Follow coding conventions in [CLAUDE.md](../../CLAUDE.md)

## Summary

Extract the repeated knot interval search pattern into a reusable `pure` function. This pattern finds the knot interval `t(l) <= x < t(l+1)` and appears ~20 times in `fitpack_core.f90`. Currently uses linear search; can be optimized with binary search for large knot vectors.

**Priority**: TOP - Quick win with potential performance benefit.

---

## The Pattern

### Current Implementation (Linear Search)

```fortran
! search for knot interval t(l) <= x < t(l+1)
do while (xi >= t(l+1) .and. l /= nk1)
    l = l + 1
end do
```

Or equivalently:

```fortran
do while (.not.(arg < t(l1) .or. l == nk1))
    l = l1
    l1 = l + 1
end do
```

### Locations in `fitpack_core.f90`

| Line | Routine | Context |
|------|---------|---------|
| 1495-1499 | `cualde` | Derivative evaluation |
| 1602-1606 | `curev` | Curve evaluation |
| 3193 | `fpclos` | Closed curve fitting |
| 4059 | `fpcons` | Constrained curve |
| 4415 | `fpcurf` | Standard curve fitting |
| 4726 | `fpcurf` | Knot insertion |
| 4938 | `fpcurf` | Main fitting loop |
| 8312 | `fppara` | Parametric curve |
| 9169 | `fpperi` | Periodic curve |
| 14261 | `splev` | Spline evaluation |
| 17819 | `splev` (pure) | Pure evaluation variant |
| 17938 | `splder` | Spline derivative |

Additional locations in surface fitting routines (`fpsurf`, `fpsphe`, `fppola`, etc.).

---

## Proposed Refactored Routine

### Option A: Simple Linear Search (drop-in replacement)

```fortran
!> Find knot interval index l such that t(l) <= x < t(l+1)
!> Starts search from l_start and stops at l_max.
pure function fp_find_knot_interval(t, x, l_start, l_max) result(l)
    real(FP_REAL), intent(in) :: t(:)       ! Knot vector
    real(FP_REAL), intent(in) :: x          ! Point to locate
    integer(FP_SIZE), intent(in) :: l_start ! Starting index (typically k+1)
    integer(FP_SIZE), intent(in) :: l_max   ! Maximum index (typically n-k-1)
    integer(FP_SIZE) :: l

    l = l_start
    do while (x >= t(l+1) .and. l /= l_max)
        l = l + 1
    end do
end function fp_find_knot_interval
```

### Option B: Binary Search (O(log n) for large knot vectors)

Based on `sorted_find` from FRESCO's sorting module:

```fortran
!> Find knot interval using binary search. O(log n) complexity.
!> Returns l such that t(l) <= x < t(l+1), clamped to [l_min, l_max].
pure function fp_find_knot_interval_binary(t, x, l_min, l_max) result(l)
    real(FP_REAL), intent(in) :: t(:)
    real(FP_REAL), intent(in) :: x
    integer(FP_SIZE), intent(in) :: l_min, l_max
    integer(FP_SIZE) :: l

    integer(FP_SIZE) :: lo, hi, mid

    ! Handle boundary cases
    if (x < t(l_min+1)) then
        l = l_min
        return
    end if
    if (x >= t(l_max+1)) then
        l = l_max
        return
    end if

    ! Binary search
    lo = l_min
    hi = l_max
    do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (x < t(mid+1)) then
            hi = mid
        else
            lo = mid
        end if
    end do

    ! Final check
    if (x >= t(hi)) then
        l = hi
    else
        l = lo
    end if
end function fp_find_knot_interval_binary
```

### Option C: Hybrid (linear for small, binary for large)

```fortran
!> Hybrid knot interval search: linear for small n, binary for large n.
pure function fp_find_knot_interval(t, x, l_start, l_max) result(l)
    real(FP_REAL), intent(in) :: t(:)
    real(FP_REAL), intent(in) :: x
    integer(FP_SIZE), intent(in) :: l_start, l_max
    integer(FP_SIZE) :: l

    integer(FP_SIZE), parameter :: LINEAR_THRESHOLD = 16

    if (l_max - l_start <= LINEAR_THRESHOLD) then
        l = fp_find_knot_interval_linear(t, x, l_start, l_max)
    else
        l = fp_find_knot_interval_binary(t, x, l_start, l_max)
    end if
end function fp_find_knot_interval
```

---

## Reference: FRESCO `sorted_find`

From `/Users/federico/fresco/fresco/src/core/sorting/sorting.fypp`, the `quickfind` routine uses:

```fortran
pure recursive function quickfind(list, x, bounds) result(it)
    ! Binary search with fallback to linear for small arrays (n <= 20)
    ! Returns positive index if found, negative insertion point if not
```

Key features to adapt:
- Recursive binary search for large arrays
- Linear search fallback for `n <= SMALL_SIZE` (20)
- Returns signed index (positive = exact match, negative = insertion point)

For FITPACK, we only need the interval (not exact match), so the logic simplifies.

---

## Implementation Steps

### Step 1: Add utility function

Add `fp_find_knot_interval` near other utility functions in `fitpack_core.f90` (around line 500, near `fpbspl`).

```fortran
!> Find knot interval index l such that t(l) <= x < t(l+1).
!> Uses linear search from l_start, stopping at l_max.
pure function fp_find_knot_interval(t, x, l_start, l_max) result(l)
    real(FP_REAL), intent(in) :: t(:)
    real(FP_REAL), intent(in) :: x
    integer(FP_SIZE), intent(in) :: l_start
    integer(FP_SIZE), intent(in) :: l_max
    integer(FP_SIZE) :: l

    l = l_start
    do while (x >= t(l+1) .and. l < l_max)
        l = l + 1
    end do
end function fp_find_knot_interval
```

### Step 2: Replace occurrences

Replace each occurrence with a call:

**Before:**
```fortran
l = k1
do while (.not.(xi < t(l+1) .or. l == nk1))
    l = l + 1
end do
```

**After:**
```fortran
l = fp_find_knot_interval(t, xi, k1, nk1)
```

### Step 3: Add binary search variant (optional optimization)

After basic replacement works, add `fp_find_knot_interval_binary` and benchmark.

### Step 4: Test thoroughly

Run `fpm test` after each batch of replacements.

---

## Testing Plan

1. **Unit tests**
   - Test with known knot vectors and query points
   - Test boundary cases: `x = t(k+1)`, `x = t(n-k)`, `x` outside range
   - Compare linear vs binary search results

2. **Regression tests**
   - All existing tests must pass
   - Numerical results must match exactly

3. **Performance benchmark** (optional)
   - Compare timing for spline evaluation with 100, 1000, 10000 knots
   - Binary search should show improvement for large n

---

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Off-by-one errors | Careful comparison with original code |
| Boundary conditions differ | Test edge cases explicitly |
| `l /= nk1` vs `l < l_max` | Verify equivalence in all contexts |

---

## Success Criteria

- [ ] All ~20 occurrences replaced with `fp_find_knot_interval`
- [ ] All existing tests pass
- [ ] No numerical differences
- [ ] Code reduction of ~40-60 lines
- [ ] (Optional) Binary search variant shows speedup for n > 100
