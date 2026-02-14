# Refactor Plan: QR Row Rotation Pattern

> **Style Guide**: Follow coding conventions in [CLAUDE.md](../../CLAUDE.md)

## Summary

Extract the repeated Givens rotation pattern used to rotate a row into an upper triangular matrix. This pattern appears ~45 times throughout `fitpack_core.f90`.

---

## The Pattern

### Variant A: Scalar RHS (most common)

```fortran
j = l-k1
do i=1,k1
    j = j+1
    piv = h(i); if (equal(piv,zero)) cycle
    call fpgivs(piv,a(j,1),cos,sin)
    call fprota(cos,sin,yi,z(j))
    if (i<k1) call fprota(cos,sin,h(i+1:k1),a(j,2:k1-i+1))
end do
```

**Found in**: `fpcurf` (line 4950), `fpperi` (line 9281), `fppara` (line 8329), `fpcons` (line 4086), etc.

### Variant B: Vector RHS (multi-dimensional)

```fortran
do i=1,k1
    j = j+1
    piv = h(i); if (equal(piv,zero)) cycle
    call fpgivs(piv,a(j,1),cos,sin)
    call fprota(cos,sin,xi(1:idim),z(j:j+(idim-1)*n:n))  ! strided
    if (i<k1) call fprota(cos,sin,h(i+1:k1),a(j,2:k1-i+1))
end do
```

**Found in**: `fpclos` (line 3306), `fpcons` (line 4086), `fppara` (line 8330), etc.

### Variant C: With secondary matrix rotation

```fortran
do i=1,kk1
    j = j+1
    piv = h(i); if (equal(piv,zero)) cycle
    call fpgivs(piv,a1(j,1),cos,sin)
    call fprota(cos,sin,yi,z(j))
    call fprota(cos,sin,h2(1:kk),a2(j,1:kk))  ! extra matrix
    if (i<kk1) call fprota(cos,sin,h1(2:i2),a1(j,2:i2))
end do
```

**Found in**: `fpclos` (line 3257), `fpperi` (line 9234), periodic spline routines.

### Variant D: Grid/surface fitting (2D indices)

```fortran
do i=1,iband
    irot = irot+1
    piv = h(i); if (equal(piv,zero)) cycle
    call fpgivs(piv,a(irot,1),cos,sin)
    call fprota(cos,sin,right(1:my),q(iq+1:iq+my))
    if (i<iband) then
        do j=i+1,iband
            call fprota(cos,sin,h(j),a(irot,j-i+1))
        end do
    endif
end do
```

**Found in**: `fpgrre` (line 6665), `fpgrdi` (line 5874), `fptrnp` (line 13924), `fptrpe`, `fpsurf`, etc.

---

## Proposed Refactored Routines

### Core routine: `fp_rotate_row`

```fortran
!> Rotate a row h(1:band) into upper triangular matrix A using Givens rotations.
!> Also applies the same rotations to scalar RHS z.
pure subroutine fp_rotate_row(h, band, a, na, z, j_start, n)
    real(FP_REAL), intent(inout) :: h(:)      ! Row to rotate (modified in place)
    integer(FP_SIZE), intent(in) :: band      ! Bandwidth (k1 or iband)
    real(FP_REAL), intent(inout) :: a(na,:)   ! Upper triangular matrix
    integer(FP_SIZE), intent(in) :: na        ! Leading dimension of a
    real(FP_REAL), intent(inout) :: z(:)      ! RHS vector (scalar case: size 1)
    integer(FP_SIZE), intent(in) :: j_start   ! Starting row index in A
    integer(FP_SIZE), intent(in) :: n         ! Stride for z (1 for scalar, n for vector)

    real(FP_REAL) :: cos, sin, piv
    integer(FP_SIZE) :: i, j

    j = j_start
    do i = 1, band
        j = j + 1
        piv = h(i)
        if (equal(piv, zero)) cycle

        call fpgivs(piv, a(j,1), cos, sin)
        call fprota(cos, sin, z(1), z(j))  ! or strided access

        if (i < band) then
            call fprota(cos, sin, h(i+1:band), a(j, 2:band-i+1))
        end if
    end do
end subroutine fp_rotate_row
```

### Extended routine: `fp_rotate_row_2mat`

For periodic splines that require rotating into two matrices:

```fortran
!> Rotate row into two matrices (A1, A2) for periodic spline handling
pure subroutine fp_rotate_row_2mat(h1, h2, band1, band2, a1, a2, na, z, j_start, n)
    real(FP_REAL), intent(inout) :: h1(:), h2(:)
    integer(FP_SIZE), intent(in) :: band1, band2
    real(FP_REAL), intent(inout) :: a1(na,:), a2(na,:)
    integer(FP_SIZE), intent(in) :: na
    real(FP_REAL), intent(inout) :: z(:)
    integer(FP_SIZE), intent(in) :: j_start, n
    ! ... implementation
end subroutine
```

### Grid variant: `fp_rotate_row_grid`

For surface fitting with 2D right-hand side:

```fortran
!> Rotate row for grid fitting (RHS is a matrix section)
pure subroutine fp_rotate_row_grid(h, band, a, na, rhs, q, irot_start, nrhs)
    real(FP_REAL), intent(inout) :: h(:)
    integer(FP_SIZE), intent(in) :: band
    real(FP_REAL), intent(inout) :: a(na,:)
    integer(FP_SIZE), intent(in) :: na
    real(FP_REAL), intent(inout) :: rhs(:), q(:)
    integer(FP_SIZE), intent(in) :: irot_start, nrhs
    ! ... implementation
end subroutine
```

---

## Implementation Steps

### Step 1: Identify all occurrences

Search patterns:
```
grep -n "fpgivs.*fprota" src/fitpack_core.f90
grep -n "rotate.*row\|rotate.*triangle" src/fitpack_core.f90
```

Key locations (line numbers are approximate):
- `fpcurf`: 4950-4965
- `fpclos`: 3075-3079, 3257-3269, 3306-3311, 3552-3564, 3576-3584
- `fpcons`: 4086-4098, 4285-4293
- `fppara`: 8329-8341, 8516-8524
- `fpperi`: 9064-9065, 9234-9245, 9259-9267, 9281-9292, 9517-9529, 9541-9550
- `fpgrdi`: 5874-5885, 5998-6013, 6028-6038, 6053-6068
- `fpgrre`: 6665-6676, 6733-6746
- `fpgrsp`: 7117-7131, 7207-7216, 7272-7281, 7299-7305
- `fptrnp`: 13924-13943
- `fptrpe`: 14082-14100, 14114-14128, 14143-14156
- `fprank`: 11008-11013, 11089-11095
- `fpsurf`: 13316-13322, 13593-13601, 13644-13652
- `fpsphe`: 12578-12584, 12871-12879, 12935-12943
- `fppola`: 10184-10191, 10334-10338, 10451-10457, 10748-10755, 10814-10821

### Step 2: Categorize by variant

Group occurrences by which variant they match (A, B, C, or D).

### Step 3: Implement base routine

1. Add `fp_rotate_row` to the module (private, near `fpgivs`/`fprota`)
2. Ensure it is `pure` and `elemental` where possible
3. Write unit tests for the isolated routine

### Step 4: Replace Variant A occurrences

Start with the simplest cases (scalar RHS, single matrix):
- `fpcurf` line ~4950
- Test after each replacement

### Step 5: Replace Variant B occurrences

Multi-dimensional RHS with strided access:
- May need an additional parameter for stride
- `fpclos`, `fpcons`, `fppara`

### Step 6: Implement and replace Variant C

Two-matrix rotation for periodic splines:
- Add `fp_rotate_row_2mat`
- Replace in `fpclos`, `fpperi`

### Step 7: Implement and replace Variant D

Grid fitting variant:
- Add `fp_rotate_row_grid`
- Replace in `fpgrdi`, `fpgrre`, `fpgrsp`, `fptrnp`, `fptrpe`, `fpsurf`, `fpsphe`, `fppola`

### Step 8: Final cleanup

- Remove any dead code
- Verify all tests pass
- Check for performance regression (should be none)

---

## Testing Plan

1. **Unit tests for new routines**
   - Test `fp_rotate_row` with known input/output
   - Verify Givens rotation properties (orthogonality)

2. **Regression tests**
   - Run full `fpm test` suite after each step
   - Compare numerical outputs to baseline (should match exactly)

3. **Performance check**
   - Time a representative fitting operation before/after
   - Subroutine call overhead should be negligible

---

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Subtle differences in loop bounds | Careful side-by-side comparison before each replacement |
| Strided array access complexity | Keep original pattern as fallback; test thoroughly |
| Breaking `pure` attribute | Ensure new routines are `pure` |
| Performance regression | Inline hints if needed; benchmark |

---

## Estimated Effort

- Step 1-2: 1-2 hours (analysis)
- Step 3: 1 hour (implement base routine + tests)
- Step 4: 2-3 hours (replace ~15 simple cases)
- Step 5-6: 2-3 hours (replace ~15 vector/periodic cases)
- Step 7: 3-4 hours (replace ~15 grid cases)
- Step 8: 1 hour (cleanup)

**Total**: ~10-14 hours of focused work

---

## Success Criteria

- [ ] All ~45 occurrences replaced with calls to 3-4 helper routines
- [ ] All existing tests pass
- [ ] No numerical differences in outputs
- [ ] Code reduction of ~300-400 lines
- [ ] Improved readability of fitting routines
