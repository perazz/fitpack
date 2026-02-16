# OOP Interface Analysis

## Current Type Hierarchy

```
fitpack_curve
  fitpack_periodic_curve          (extends, no new members)

fitpack_parametric_curve
  fitpack_closed_curve            (extends, no new members)
  fitpack_constrained_curve       (extends, adds derivative constraints)

fitpack_surface                   (scattered bivariate)
fitpack_grid_surface              (gridded bivariate)

fitpack_parametric_surface        (gridded parametric)

fitpack_polar                     (scattered polar)
fitpack_grid_polar                (gridded polar)

fitpack_sphere                    (scattered spherical)
fitpack_grid_sphere               (gridded spherical)
```

11 types total, no shared abstract base. Each type independently declares the same data
components and method patterns.

---

## 1. Duplicated Data Components

Every type carries these components (names vary slightly):

| Component | Curve | Surface | Grid Surf | Param Curve | Param Surf | Polar | Grid Polar | Sphere | Grid Sphere |
|-----------|:-----:|:-------:|:---------:|:-----------:|:----------:|:-----:|:----------:|:------:|:-----------:|
| knots / knots(2) | x | x | x | x | x | x | x | x | x |
| t(:) / t(:,:) | x | x | x | x | x | x | x | x | x |
| c(:) | x | x | x | x | x | x | x | x | x |
| order / order(2) | x | x | x | x | x | x | x | x | x |
| smoothing | x | x | x | x | x | x | x | x | x |
| fp | x | x | x | x | x | x | x | x | x |
| iopt | x | x | x | x | x | x | x | x | x |
| nest / nest(2) | x | x | x | x | x | x | x | x | x |
| wrk(:) | x | - | x | x | x | - | x | - | x |
| wrk1(:),wrk2(:) | - | x | - | - | - | x | - | x | - |
| iwrk(:) | x | x | x | x | x | x | x | x | x |
| bc | x | x | x | - | - | x | x | - | - |
| w(:) | x | x | - | x | - | x | - | x | - |

9 types x ~12 identical fields = significant duplication.

---

## 2. Duplicated Methods

| Method | Purpose | Present in |
|--------|---------|-----------|
| `destroy` | Deallocate, reset state | All 9 types |
| `new_points` | Store data, allocate workspace | All 9 types |
| `new_fit` | `new_points` + `fit` combined | All 9 types |
| `fit` | Automatic knot smoothing fit | All 9 types |
| `eval` | Evaluate spline at point(s) | All 9 types |
| `mse` | Return `fp` | Curve, Param Curve, Param Surf |
| `interpolate` | Fit with s=0 | Curve, Grid Surf, Param Curve, Param Surf, Polar, Grid Polar, Sphere, Grid Sphere |
| `least_squares` | Fit with fixed knots | Grid Surf, Param Surf, Polar, Grid Polar, Sphere, Grid Sphere |
| `dfdx` | Derivative(s) | Curve, Grid Surf, Surface, Param Curve |
| `integral` | Integration | Curve only |
| `zeros` | Root finding | Curve only |
| `fourier_coefficients` | Fourier analysis | Curve only |

The `destroy`, `new_points`, `new_fit`, `fit`, `eval` cycle is universal. An abstract
base could define the interface pattern, even if implementations differ.

---

## 3. Settings Management — Current State

Settings are bare public components with no accessors:

```fortran
type(fitpack_curve) :: c
c%smoothing = 0.5     ! direct access, no validation
c%order = 3           ! can set to invalid values
c%bc = OUTSIDE_ZERO   ! works but not discoverable
```

**Problems:**
- No validation on set (order must be 1-5, smoothing must be >= 0)
- No way to query valid ranges
- Boundary conditions use integer flags — not self-documenting
- Some settings only make sense at specific times (e.g., can't change `order` after fit
  without invalidating workspace sizes)
- No `get_knots()` / `set_knots()` — knots accessed directly via `t(:)` and `knots`

**Comparison across types:**

| Setting | Curve | Surface | Polar | Sphere |
|---------|-------|---------|-------|--------|
| Smoothing | `smoothing` | `smoothing` | `smoothing` | `smoothing` |
| Degree | `order` | `order(2)` | implicit (3,3) | implicit (3,3) |
| Boundary | `bc` flag | `bc` flag | `bc_boundary` + `bc_continuity_origin` | `pole_continuity(2)` + `pole_zero_grad(2)` |
| Domain | `xleft,xright` | `left(2),right(2)` | `rad` function ptr | fixed [0,pi]x[0,2pi] |
| Weights | `w(:)` | `w(:)` | `w(:)` | `w(:)` |

---

## 4. Orphaned fitpack_core Routines (Not in OOP Layer)

### High-value — should be wrapped

| Routine | Book Section | What it does | Suggested home |
|---------|-------------|--------------|----------------|
| `cocosp` | §7.1 | Convexity-constrained least-squares with given knots | New `fitpack_constrained_spline` or extend `fitpack_curve` |
| `concon` | §7.2 | Convexity-constrained smoothing (automatic knots) | Same |
| `dblint` | §8.3 | Double integration of bivariate spline | `fitpack_surface%integral`, `fitpack_grid_surface%integral` |
| `profil` | §8.2 | Cross-section of bivariate spline → univariate spline | `fitpack_surface%cross_section` → returns `fitpack_curve` |
| `insert` | §4.2 | Knot insertion (out-of-place) | All curve/surface types: `insert_knot` method |
| `pardtc` | §8.1 | Create derivative spline of bivariate spline | `fitpack_surface%derivative_spline` → returns new surface |

### Low-value — internal utilities, no wrapping needed

| Routine | Why not needed |
|---------|---------------|
| `bispeu` | Duplicate of `bispev` for unstructured points; already handled internally |
| `insert_inplace` | Used internally by `concur` path |
| `equal`, `fitpack_swap`, `fitpack_argsort` | Internal helpers |
| `FP_COMM_*` | Specialized MPI communication — already used by `fitpack_curve` |
| `evapol` | Called internally by `fitpack_polar%eval` |

---

## 5. Proposed Abstract Base Type

```fortran
type, abstract :: fitpack_fitter
    ! --- Spline representation (universal) ---
    integer(FP_SIZE)              :: iopt = IOPT_NEW_SMOOTHING
    real(FP_REAL)                 :: smoothing = 1000.0_FP_REAL
    real(FP_REAL)                 :: fp = zero          ! actual MSE after fit
    real(FP_REAL),    allocatable :: c(:)                ! B-spline coefficients
    integer(FP_SIZE), allocatable :: iwrk(:)             ! integer workspace
contains
    ! --- Deferred (each subtype implements) ---
    procedure(destroy_if),    deferred :: destroy
    procedure(fit_if),        deferred :: fit

    ! --- Non-overridable shared logic ---
    procedure, non_overridable :: mse              ! returns fp
    procedure, non_overridable :: set_smoothing    ! with validation
    procedure, non_overridable :: smoothing_status  ! ok/interpolating/lsq/error
end type
```

**What goes in the base vs. subtypes:**

| In abstract base | In subtypes |
|-----------------|-------------|
| `smoothing`, `fp`, `iopt` | knot arrays (1D vs 2D) |
| `c(:)`, `iwrk(:)` | data arrays (x,y vs x,y,z vs theta,phi,r) |
| `mse()`, `set_smoothing()` | `eval()` (return type varies) |
| `smoothing_status()` | `new_points()` (signature varies) |
| `destroy` interface | `fit()` implementation |
| `fit` interface | domain-specific settings |

**Key design question**: The eval/dfdx signatures differ too much across types (scalar vs
vector return, 1D vs 2D input) to put in the abstract base. The base would contain only the
fitting lifecycle and shared state management.

---

## 6. Settings Improvement Ideas

### Validated setters

```fortran
! Instead of: curve%order = 7  (silently invalid)
call curve%set_order(7, ierr)  ! returns FITPACK_INPUT_ERROR

! Instead of: curve%bc = 42  (meaningless)
call curve%set_boundary(OUTSIDE_EXTRAPOLATE)
```

### Query methods

```fortran
n = curve%num_knots()          ! instead of curve%knots
t = curve%get_knots()          ! returns copy of t(1:knots)
k = curve%get_order()          ! returns order
s = curve%get_smoothing()      ! returns current smoothing
```

### State-aware validation

```fortran
call curve%set_order(5)        ! invalidates workspace if different from current
call curve%fit()               ! reallocates workspace if needed
```

---

## 7. Missing OOP Functionality Summary

| Feature | Book ref | Difficulty | Impact |
|---------|---------|------------|--------|
| Convexity constraints (`concon`/`cocosp`) | §7.1-7.2 | Medium | Unique capability, no other library does this easily |
| Surface integration (`dblint`) | §8.3 | Low | Natural extension of `splint` to 2D |
| Surface cross-section (`profil`) | §8.2 | Low | Returns 1D spline from 2D — very useful |
| Derivative spline (`pardtc`) | §8.1 | Low | Creates new spline = d/dx of original |
| Knot insertion (`insert`) | §4.2 | Low | Useful for refinement and h-adaptivity |
| `bispeu` for scattered eval | — | Low | Currently surfaces eval on grids only; scattered eval missing from OOP |

---

## 8. Communication Interface Gap

Only `fitpack_curve` has `comm_size`/`comm_pack`/`comm_expand` for MPI serialization.
None of the surface types support this. If the communication pattern is worth keeping, it
should be extended to all types (or moved to the abstract base with a default implementation
that serializes all allocatable components).

---
---

# PR Plans

## PR A: Public API Cleanup & Umbrella Module (DONE)

**Goal**: Make the `fitpack` umbrella module self-sufficient — users should never need
`use fitpack_core` directly.

**Problem**: `fitpack.f90` exports the 11 types but not the flag constants they depend on.
A user writing `curve%bc = OUTSIDE_EXTRAPOLATE` must also `use fitpack_core`. Similarly,
`fitpack_surface` has no `eval` method (only `dfdx`), and `mse` is only on 3 of 9 types.

### Steps

1. **Re-export constants from `fitpack.f90`**:
   - Boundary flags: `OUTSIDE_EXTRAPOLATE`, `OUTSIDE_ZERO`, `OUTSIDE_NOT_ALLOWED`,
     `OUTSIDE_NEAREST_BND`
   - Fitting flags: `IOPT_NEW_LEASTSQUARES`, `IOPT_NEW_SMOOTHING`, `IOPT_OLD_FIT`
   - Error flags: `FITPACK_OK`, `FITPACK_SUCCESS`, `FITPACK_MESSAGE`,
     `FITPACK_INPUT_ERROR`, etc.
   - Kind parameters: already exported (`FP_REAL`, `FP_SIZE`, `FP_FLAG`, `FP_BOOL`)
   - Named constants: `zero`, `one`, `half`, `pi`, etc.

2. **Add `eval` to `fitpack_surface`**: It currently only has `dfdx`/`dfdx_ongrid`.
   Add `eval` generic (one-point and many-point) wrapping `bispev`/`bispeu`.

3. **Add `mse` to all types**: Currently only `fitpack_curve`, `fitpack_parametric_curve`,
   and `fitpack_parametric_surface` have it. It's a trivial one-liner (`mse = this%fp`).
   Add to: `fitpack_surface`, `fitpack_grid_surface`, `fitpack_polar`, `fitpack_grid_polar`,
   `fitpack_sphere`, `fitpack_grid_sphere`.

4. **Add `interpolate` and `least_squares` consistently**: `interpolate` is missing from
   `fitpack_surface`; `least_squares` is missing from `fitpack_curve`,
   `fitpack_parametric_curve`, `fitpack_surface`. Add where they make sense.

### Files touched
- `src/fitpack.f90` (re-exports)
- `src/fitpack_surfaces.f90` (add `eval`, `mse`, `interpolate`, `least_squares`)
- `src/fitpack_grid_surfaces.f90` (add `mse`)
- `src/fitpack_polar.f90` (add `mse`)
- `src/fitpack_gridded_polar.f90` (add `mse`)
- `src/fitpack_spheres.f90` (add `mse`)
- `src/fitpack_gridded_sphere.f90` (add `mse`)

### Testing
- Existing tests must pass unchanged
- Add a test that `use fitpack` alone is sufficient (no `use fitpack_core` needed)

---

## PR B: Abstract Base Type (DONE — PR #41)

**Goal**: Extract shared state and logic into an abstract `fitpack_fitter` base type that
all 9 concrete types extend.

### Design

The base contains **only** what is truly universal across all 9 types — the fitting state
machine and spline coefficient storage:

```fortran
!> Abstract base for all fitpack spline fitters.
type, abstract :: fitpack_spline

    !> Fitting state: IOPT_NEW_SMOOTHING, IOPT_NEW_LEASTSQUARES, IOPT_OLD_FIT
    integer(FP_SIZE) :: iopt = IOPT_NEW_SMOOTHING

    !> Smoothing parameter (controls fit-vs-data tradeoff)
    real(FP_REAL) :: smoothing = 1000.0_FP_REAL

    !> Weighted sum of squared residuals after last fit
    real(FP_REAL) :: fp = zero

    !> B-spline coefficients
    real(FP_REAL), allocatable :: c(:)

    !> Integer workspace
    integer(FP_SIZE), allocatable :: iwrk(:)

contains

    !> Weighted MSE of last fit
    procedure, non_overridable :: mse

    !> Fitting status after last call
    procedure, non_overridable :: fitting_status

    !> Parallel communication interface (deferred — each subtype knows its own layout)
    procedure(comm_size_if),   deferred :: comm_size    ! buffer size needed
    procedure(comm_pack_if),   deferred :: comm_pack    ! serialize state to buffer
    procedure(comm_expand_if), deferred :: comm_expand  ! deserialize state from buffer

end type fitpack_spline
```

The communication interface is **deferred** because each subtype has different scalar
fields and allocatable arrays to serialize. The abstract interfaces are:

```fortran
abstract interface
    elemental integer(FP_SIZE) function comm_size_if(this)
        import fitpack_spline
        class(fitpack_spline), intent(in) :: this
    end function

    pure subroutine comm_pack_if(this, buffer)
        import fitpack_spline, FP_COMM
        class(fitpack_spline), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)
    end subroutine

    pure subroutine comm_expand_if(this, buffer)
        import fitpack_spline, FP_COMM
        class(fitpack_spline), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)
    end subroutine
end interface
```

Currently only `fitpack_curve` implements these. All other types will gain communication
support as they migrate to the base type.

**What stays OUT of the base** (and why):

| Component | Why not in base |
|-----------|----------------|
| `knots` / `knots(2)` | Scalar vs 2-element array |
| `t(:)` / `t(:,:)` | 1D vs 2D |
| `order` / `order(2)` | Scalar vs 2-element array |
| `nest` / `nest(2)` | Same |
| `wrk(:)` / `wrk1(:),wrk2(:)` | 1 vs 2 workspace arrays |
| `bc` / `bc_boundary` / `pole_*` | Domain-specific meaning |
| Data arrays | Completely different per type |
| `destroy`, `fit`, `eval` | Signatures differ too much |

`fit()` and `eval()` are NOT deferred — there is no single abstract interface that fits
all 9 signatures (different return types, different argument lists). The base provides
shared data, shared non-overridable methods, and the deferred communication interface.

### Steps

1. **Define `fitpack_spline`** in a new `src/fitpack_types.f90` module (or in
   `fitpack_core.F90` near the type definitions section). Include the abstract
   communication interfaces (`comm_size_if`, `comm_pack_if`, `comm_expand_if`).

2. **Migrate `fitpack_curve`**: Change to `type, extends(fitpack_spline) :: fitpack_curve`.
   Remove `iopt`, `smoothing`, `fp`, `c(:)`, `iwrk(:)` from the type body (inherited).
   Remove the local `mse` implementation. The existing `comm_size`/`comm_pack`/`comm_expand`
   implementations already satisfy the deferred interfaces — just verify they compile.

3. **Migrate remaining types** one at a time, same mechanical process:
   - `fitpack_parametric_curve` (and its extensions)
   - `fitpack_surface`, `fitpack_grid_surface`
   - `fitpack_parametric_surface`
   - `fitpack_polar`, `fitpack_grid_polar`
   - `fitpack_sphere`, `fitpack_grid_sphere`
   Each migration requires implementing `comm_size`/`comm_pack`/`comm_expand` for that
   type. Follow the same pattern as `fitpack_curve`: pack scalar fields into FP_REAL slots,
   then pack each allocatable array using `FP_COMM_PACK`.

4. **Fix `destroy` implementations**: Each subtype's `destroy` must also reset the
   inherited base fields. Add a `reset_base` helper or just reset them explicitly.

5. **Update `fitpack.f90`**: Re-export `fitpack_spline` if users need to declare
   `class(fitpack_spline)` polymorphic containers, or keep it private if not needed yet.

### Risks
- `elemental` on `destroy` may conflict with `class` dummy arguments if the base is
  abstract (Fortran requires `impure elemental` for polymorphic dummy args in some
  compilers). Test with gfortran early.
- `non_overridable` on `mse` in the base means subtypes cannot have their own — verify
  no subtype has a different `mse` semantic (they don't, all return `fp`).

### Files touched
- New: `src/fitpack_types.f90` (or section in `fitpack_core.F90`)
- Modified: all 9 OOP wrapper modules
- Modified: `src/fitpack.f90`

### Testing
- All 49 existing tests must pass
- No public API change — this is purely internal refactoring

---

## PR C: Missing Curve Functionality

**Goal**: Wrap `concon`/`cocosp` (convexity constraints) and `insert` (knot insertion).

### Convexity-constrained fitting

**Book reference**: §7.1-7.2 — unique FITPACK capability for fitting cubic splines under
local convexity/concavity constraints.

**Design options**:
1. New type `fitpack_convex_curve` extending `fitpack_curve` — adds `v(:)` (constraint
   signs), `maxtr`, `maxbin`, `bind(:)` (active constraint flags)
2. Add a `fit_convex` method to `fitpack_curve` — simpler but mixes concerns

Option 1 is cleaner because convexity fitting has different workspace requirements,
different error semantics, and different output (`bind` array, `sq` instead of `fp`).

```fortran
type, extends(fitpack_curve) :: fitpack_convex_curve
    real(FP_REAL), allocatable    :: v(:)     ! Constraint signs at data points
    logical(FP_BOOL), allocatable :: bind(:)  ! Active constraints after fit
    integer(FP_SIZE) :: maxtr  = 100          ! Tree storage estimate
    integer(FP_SIZE) :: maxbin = 10           ! Max zero-curvature knots
contains
    procedure :: fit => convex_fit_automatic_knots   ! calls concon
    procedure :: fit_fixed_knots                     ! calls cocosp
end type
```

### Knot insertion

Two interfaces — in-place (primary) and out-of-place (returns new curve):

```fortran
!> Insert a knot into this curve in-place. The curve's knot vector, coefficients,
!> and workspace are updated. Calls insert_inplace internally.
subroutine insert_knot(this, x, ierr)
    class(fitpack_curve), intent(inout) :: this
    real(FP_REAL), intent(in) :: x             ! knot to insert
    integer(FP_FLAG), intent(out) :: ierr
end subroutine

!> Return a new curve with an additional knot inserted. The original is unchanged.
!> Calls insert (out-of-place) internally.
function curve_with_knot(this, x, ierr) result(refined)
    class(fitpack_curve), intent(in) :: this
    real(FP_REAL), intent(in) :: x             ! knot to insert
    integer(FP_FLAG), intent(out) :: ierr
    type(fitpack_curve) :: refined             ! new curve with extra knot
end function
```

Both should be overloaded under a single `generic :: insert => insert_knot, curve_with_knot`
if the signatures are distinguishable (they are: `intent(inout)` subroutine vs `intent(in)`
function). Otherwise keep them as separate named methods.

### Steps

1. Define `fitpack_convex_curve` type in `fitpack_curves.f90` (or new module).
2. Implement `convex_fit_automatic_knots` wrapping `concon`.
3. Implement `fit_fixed_knots` wrapping `cocosp`.
4. Add `insert_knot` (in-place) and `curve_with_knot` (out-of-place) to `fitpack_curve`.
5. Export from `fitpack.f90`.
6. Add tests based on the existing `concon`/`cocosp` test cases.

### Files touched
- `src/fitpack_curves.f90` (or new `src/fitpack_convex_curves.f90`)
- `src/fitpack.f90`
- `test/` — new test functions

---

## PR D: Missing Surface Functionality

**Goal**: Wrap `dblint` (integration), `profil` (cross-section), `pardtc` (derivative
spline), and add scattered-point evaluation to surface types.

### Surface integration (`dblint`)

Add `integral` method to `fitpack_surface` and `fitpack_grid_surface`:
```fortran
real(FP_REAL) function integral(this, xb, xe, yb, ye)
    ! Integrates the bivariate spline over [xb,xe] x [yb,ye]
    ! Calls dblint internally
end function
```

### Cross-section (`profil`)

Add `cross_section` method to surface types:
```fortran
type(fitpack_curve) function cross_section(this, iopt, u, ierr)
    ! iopt=0: x=u cross-section (returns curve in y)
    ! iopt=1: y=u cross-section (returns curve in x)
    ! Calls profil, returns a fully initialized fitpack_curve
end function
```

This is particularly useful — it creates a bridge between the 2D and 1D types.

### Derivative spline (`pardtc`)

Add `derivative_spline` method:
```fortran
subroutine derivative_spline(this, nux, nuy, deriv, ierr)
    ! Creates a new surface whose values are d^(nux+nuy)/dx^nux dy^nuy of this
    ! Calls pardtc to compute new coefficients
    class(fitpack_surface), intent(in) :: this
    type(fitpack_surface), intent(out) :: deriv
end function
```

### Scattered-point evaluation

`fitpack_surface` currently has no `eval` — only `dfdx`/`dfdx_ongrid`. Add:
```fortran
procedure, private :: surf_eval_one    ! single (x,y) -> z
procedure, private :: surf_eval_many   ! arrays x(:),y(:) -> z(:)
generic :: eval => surf_eval_one, surf_eval_many
```

These would call `bispev` (for gridded evaluation internally, building a 1x1 or mx1 grid)
or `bispeu` (for unstructured point lists).

### Steps

1. Add `eval` to `fitpack_surface` (wrapping `bispev`/`bispeu`).
2. Add `integral` to `fitpack_surface` and `fitpack_grid_surface` (wrapping `dblint`).
3. Add `cross_section` to `fitpack_surface` and `fitpack_grid_surface` (wrapping `profil`).
4. Add `derivative_spline` to `fitpack_surface` and `fitpack_grid_surface`
   (wrapping `pardtc`).
5. Add tests for each new method.

### Files touched
- `src/fitpack_surfaces.f90`
- `src/fitpack_grid_surfaces.f90`
- `test/` — new test functions

---

## PR Dependency Order

```
PR A  (API cleanup)          — no dependencies, goes first
  |
PR B  (abstract base type)   — benefits from PR A's consistent interface
  |
PR C  (convexity + insert)   — independent of B, depends on A
PR D  (surface features)     — independent of B, depends on A
```

PRs C and D can proceed in parallel with B. PR A should go first to establish the
consistent baseline (especially adding `mse` and `eval` everywhere, which PR B then
consolidates into the abstract base).
