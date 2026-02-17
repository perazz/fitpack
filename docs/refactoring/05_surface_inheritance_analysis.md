# Surface Type Inheritance: Feasibility Analysis

Could `fitpack_grid_surface` and/or `fitpack_parametric_surface` extend `fitpack_surface`,
analogous to how `fitpack_periodic_curve` extends `fitpack_curve`?

**Conclusion: Not recommended.** The surface types fail the key prerequisite for Fortran
type extension — a shared data model.

---

## Why Curve Inheritance Works

The curve hierarchy succeeds because subtypes share **the same data model**:

```
fitpack_curve                  — x(:), y(:), w(:), order, knots, t(:), c(:), ...
  └── fitpack_periodic_curve   — empty extension (marker type, same fields)
  └── fitpack_convex_curve     — adds constraint fields, overrides fit
```

- `fitpack_periodic_curve` is literally an empty `type, extends(fitpack_curve)` — same
  fields, same eval, same destroy. The **only** difference is a `select type` in
  `curve_fit_automatic_knots` that calls `percur` instead of `curfit`.
- `fitpack_convex_curve` adds fields (`v`, `sx`, `bind`) and overrides `fit`/`destroy`/
  `new_points`, but the base data model (scattered x/y points, 1D knots, 1D coefficients)
  is identical.

---

## Why Surface Inheritance Does NOT Work

### fitpack_grid_surface extending fitpack_surface

**Data model conflict:**

| Field | `fitpack_surface` | `fitpack_grid_surface` |
|-------|-------------------|------------------------|
| Input data | `x(m), y(m), z(m)` scattered | `x(nx), y(ny), z(ny,nx)` grid |
| Weights | `w(m)` | *(none)* |
| Boundary mode | `bc` | *(none)* |
| 2nd workspace | `wrk2(:), lwrk2` | *(none)* |
| Fit routine | `surfit` | `regrid` |
| Eval (scattered) | `bispeu` | *(N/A)* |
| Eval (grid) | `bispev` | `bispev` |

The `z` arrays have **different ranks** (`z(:)` vs `z(:,:)`). Fortran cannot redefine a
parent field's shape. If grid extended surface, it would carry `m`, `w(:)`, `bc`, `wrk2(:)`,
`lwrk2` — all meaningless ballast. Every method (`fit`, `eval`, `destroy`, `new_points`,
all 3 `comm_*`) would need to be overridden, so there is zero code reuse from the parent.

### fitpack_parametric_surface extending fitpack_surface

Even less feasible. The parametric surface is a fundamentally different beast:
- Uses `u(:), v(:)` parameters, not `x(:), y(:)` coordinates
- Has `idim` (multi-dimensional output) and `periodic_dim(2)` — no analogs in `fitpack_surface`
- No `order` field (always bicubic, hardcoded in `parsur`)
- No `left`/`right` boundaries
- `z(:,:,:)` is rank-3
- `eval` returns `f(mv,mu,idim)` — multi-dimensional, not scalar
- Uses completely different core routines: `parsur`/`surev`

---

## What IS Shared

Between `fitpack_surface` and `fitpack_grid_surface`, the **2D spline representation**
(the output of fitting) is identical:

```fortran
integer :: order(2), nest(2), nmax, knots(2)
real(FP_REAL) :: left(2), right(2)
real(FP_REAL), allocatable :: t(:,:), c(:)
```

These fields power the evaluation-side methods that are currently duplicated: `integral`,
`cross_section`, `derivative_spline`, `dfdx`/`dfdx_ongrid`.

`fitpack_parametric_surface` shares only `nest(2)`, `nmax`, `knots(2)`, `t(:,:)` — it lacks
`order` and `left/right`.

---

## Summary

Unlike `fitpack_periodic_curve` (same struct, different fitting dispatch), the surface types
differ in input data layout, array ranks, supported features, and underlying core routines.
Forcing inheritance would create types carrying unused fields, overriding every method,
gaining no polymorphism (since `eval` return types differ). The ~100 lines of duplicated
code in `integral`/`cross_section`/`derivative_spline` is the cost of keeping clean, honest
types — and it is manageable.
