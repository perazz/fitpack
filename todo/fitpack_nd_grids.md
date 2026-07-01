# Extending FITPACK to N-Dimensional Gridded Splines — Feasibility & Design

**Status:** IMPLEMENTED. The dimension-generic core (`regrid → fpregr → fpgrre` fit, `ndspev → fpndsp` eval) replaced the legacy 2-D trio and now backs both the 2-D `fitpack_grid_surface` and the new runtime-`dims` `fitpack_gridded_spline` class (3D+). Note: the endgame differed from §5.3 below — the team chose **one runtime-`dims` public type** (not fypp-generated concrete faces), with an F2018 `select rank` + `row_major` face for rank-natural input, since a `select rank` construct sidesteps the "rank problem" of §5.2 without preprocessing. The **peripherals are now generalized too** — scattered evaluation (`ndspeu`), partial derivatives (`pardtc`/`parder`/`pardeu`), the box integral (`dblint`), and the cross-section (`profil`), wired as `fitpack_gridded_spline` methods (`eval`/`dfdx`/`dfdx_ongrid`/`integral`/`cross_section`/`derivative_spline`). The legacy 2-D peripherals were **collapsed into these dimension-generic routines** (the classic `(tx,nx,ty,ny,…)` bodies were removed; the `_nd` suffix was dropped, and the scattered evaluator was renamed `bispeu_nd → ndspeu` to pair with `ndspev`). The 2-D surface classes and the `fp_bispeu_c`/`fp_parder_c` C ABI wrappers marshal their 2-D arguments into these routines; the classic-signature regression tests were migrated likewise (Gates I–M retain the dims=3 analytic/separable oracles). **C bindings (§ slice 10) are the sole remaining future work.**
**Scope:** tensor-product *gridded* smoothing splines over a `d`-dimensional domain (`d = 3 … ~6`), generalizing the current 1-D (`fitpack_curve`) and 2-D (`fitpack_grid_surface`) gridded fitters.
**Out of scope:** N-D *scattered* fitting (see §3).

---

## Remaining work (handoff for the next agent)

The dimension-generic **fit** (`regrid`/`fpregr`/`fpgrre`) and **eval** (`ndspev`/`fpndsp`)
cores are done, the legacy 2-D fit trio is retired, and the public runtime-`dims` class
`fitpack_gridded_spline` (fit + eval + comm serialization) is implemented and tested — Gate H
in `test/fitpack_grid_nd_tests.f90` asserts the class equals a direct `regrid`+`ndspev` call
bit-for-bit at dims=3 and dims=5, plus the `row_major` flag and comm round-trip. The peripherals
below are **now DONE** (Gates I–M); only the C bindings remain:

| Item | Dimension-generic routine (classic 2-D removed) | Status |
|------|------------------------------|--------|
| **Scattered evaluation** | `ndspeu` (was `bispeu`) | **DONE** — mixed-radix odometer, self-managed basis. Class `%eval`. Gate I. |
| **Partial derivatives** | `parder`/`pardeu`/`pardtc` (classic 2-D bodies removed) | **DONE** — `pardtc` (per-axis derivative recurrence + repack), `parder` (= `pardtc` + `ndspev`), `pardeu` (= `pardtc` + `ndspeu`). Class `%dfdx` / `%dfdx_ongrid` / `%derivative_spline`. Gate J/M. |
| **Domain integral** | `dblint` (classic 2-D body removed) | **DONE** — per-axis `fpintb` vectors contracted against the coefficient tensor (fpndsp odometer). Class `%integral`. Gate K. |
| **Cross-section / slice** | `profil` (classic 2-D body removed) | **DONE** — fix one axis → `dims-1` spline. Class `%cross_section` returns a `fitpack_gridded_spline` of `dims-1`. Gate L/M. |
| **C bindings (slice 10)** | `regrid_c` / `bispev_c` template | **REMAINING.** First surface-family opaque-pointer C binding. Follow `fitpack_curve_c` (`src/fitpack_curves_c.f90`): opaque `type, bind(C)` + `_allocate`/`_destroy`/`_get_pointer` + flat `fp_*_c` method functions (scalars by value, arrays by pointer). New file `src/fitpack_gridded_splines_c.f90`. |

**Where things live**
- Core N-D routines (all in `src/fitpack_core.F90`): `regrid`, `fpregr`, `fpgrre` (fit);
  `ndspev`, `fpndsp` (eval); `new_knot_dimension_nd` (knot-direction arbiter); peripherals
  `ndspeu`, `pardtc`, `parder`, `pardeu`, `dblint`, `profil` (dimension-generic; the classic 2-D
  bodies were removed and their public names/signatures now belong to these routines).
- Public class: `src/fitpack_gridded_splines.f90` (`fitpack_gridded_spline`, extends
  `fitpack_fitter`; runtime `dims`, fixed-`MAX_IDIM` metadata, `select rank` + `row_major` face;
  methods `eval`/`eval_ongrid`/`dfdx`/`dfdx_ongrid`/`integral`/`cross_section`/`derivative_spline`).
- Tests: `test/fitpack_grid_nd_tests.f90` — Gates A/D/D'/E/F/G (cores) + H (class) + I/J/K/L
  (peripheral cores, bit-for-bit at dims=2 & independent oracle at dims=3) + M (class peripherals).

**Method to follow** — same as slices 1–11: generalize each peripheral behind an
independent-oracle gate (the pattern of Gates D/D'/E/F/G/H, e.g. `fp==SSR`, separable-product,
closed-form), verify `fpm test` green on normal + `-fcheck=bounds` + `-finit-real=inf` poison,
then wire a type-bound method on `fitpack_gridded_spline`.

---

## 0. Executive summary

| Path | Core routines | N-D feasibility | Call |
|---|---|---|---|
| **Gridded** (`fitpack_grid_surface`) | `regrid → fpregr → fpgrre` | **High** — already a loop-over-2-directions of independent 1-D solves | **Do it.** One dimension-generic core + thin generated typed faces |
| **Scattered** (`fitpack_surface`) | `surfit → fpsurf` | **Low** — genuinely 2-D-coupled (fused banded system, panel indexing, rank handling) | **Leave out.** Same boundary bspline-fortran draws (delegates to SPLPAK) |

**Recommended architecture:** one *dimension-generic procedural core* (runtime `d`, flat arrays) + **fypp-generated concrete faces** `fitpack_gridded_3d … fitpack_gridded_Nd` (rank-natural `z(:,:,:…)`, compile-time dimension safety). **Avoid** parameterized derived types (PDT) and **avoid** a single runtime-`ndim` type for the public API (rank problem, §5.2).

**Practical ceiling:** memory/eval cost is `O(∏ nᵢ)` / `O(kᵈ)` — the curse of dimensionality, not the code, is the limit. Realistic sweet spot **d = 3–4 (maybe 5)**. `MAX_IDIM = 10` is aspirational.

---

## 1. Background: two orthogonal "dimensions"

FITPACK already has a notion of "N-D", but it is **not** the one we want to raise. Keep these distinct in code and docs:

- **Codomain / parametric dimension `idim`** — a curve `s(u) = (s₁(u), …, s_idim(u))` whose *output* lives in ℝ^idim. Already general up to `MAX_IDIM = 10` (`fitpack_core.F90:105`); used by parametric/closed curves (`clocur`, etc.; `idim` appears ~243× in the core).
- **Domain / tensor dimension `d`** — number of *independent variables*, i.e. the rank of the tensor-product B-spline basis. Currently **1** (curve) or **2** (surface). **This is what we raise.**

They compose eventually (a vector field over a `d`-D grid: domain `d`, codomain `idim`), but the first target is **scalar codomain** (`idim = 1`). Suggested vocabulary: `idim` = codomain, `ndim`/`d` = domain.

---

## 2. Why the gridded path generalizes (the decisive fact)

FITPACK's gridded solver **never forms a `d`-dimensional system**. It solves the tensor-product least-squares problem as a *sequence of independent 1-D banded solves*, one per axis (de Boor's alternating-direction method). This is dimension-agnostic by construction.

In `fpgrre` (`fitpack_core.F90:6940`):

1. Build one observation matrix **per axis** — `spx(mx,kx1)`, `spy(my,ky1)` — each from the same 1-D kernel `fpbspl`.
2. QR-factorize each axis independently into a band matrix (`ax`, `ay`) via Givens rotations.
3. Solve `(ry)·c·(rx)' = h` as **two passes of 1-D back-substitution** (`fpback`), one per tensor axis (`fitpack_core.F90:7120-7141`).

For `d` dimensions: build `d` observation matrices, `d` independent QR factorizations, a `d`-pass alternating-direction back-substitution. **No new linear algebra.** It is the existing 1-D machinery inside a `do idim = 1, d` loop.

Crucially, the per-axis primitives are **already refactored into dimension-agnostic 1-D form** (these were modernized away from the original Dierckx code):

| Primitive | `fitpack_core.F90` | Generic? | Role |
|---|---|---|---|
| `fp_rotate_row_block(h,band,a,na,right,q,nrhs,irot)` | 5859 | ✅ | single-axis QR row insert, `nrhs` configurable |
| `fp_rotate_row_stride(...)` | 5901 | ✅ | strided variant for the 2nd axis |
| `fp_rotate_2mat_stride(...)` | 5946 | ✅ | periodic two-matrix variant |
| `fpback` | 1763 | ✅ | banded back-substitution |
| `fpbspl` | 2219 | ✅ | Cox–de Boor basis values |
| `fpknot` | 8035 | ✅ | single-axis knot insertion |
| `fpdisc` | 5284 | ✅ | discontinuity/smoothing matrix |

`fpregr`'s knot-placement loop (`fitpack_core.F90:11907`) is **already per-axis**: `reducx/reducy`, `fpintx/fpinty`, `nplusx/nplusy`, `nrdatx/nrdaty` are paired, and the "which direction gets the next knot" decision is already abstracted into a pure arbiter:

```fortran
elemental integer function new_knot_dimension(n1,n1add,n1max,n2,n2add,n2max,last)  ! :12305
```

Binary today → becomes an **argmax over `d`** candidate directions. The rest of `fpregr` (the `p`-iteration / rational-interpolation root find for `f(p)=s`, the interpolation/least-squares bootstrap, even/odd-order knot init) is scalar or trivially looped.

**What is hard-coded to 2 (all mechanical, all the actual work):**
- `wrk`/`iwrk` partition pointer arithmetic in `regrid` (`fitpack_core.F90:17133-17148`: `lsx,lsy,lri,lq,lax,lay,lbx,lby`) and `fpregr` → becomes an offset loop over `d`.
- The intermediate **"ping-pong" RHS tensor** `q` (`mynx = nxest*my` in 2-D) — in N-D this buffer alternates between "axes already reduced to coefficients" × "axes still in data space". **Single most error-prone piece** (sizing + stride bookkeeping). Study bspline-fortran's `dbtpcf` for the pattern — it does exactly this.
- The residual-accounting contraction in `fpgrre` (`fitpack_core.F90:7157-7194`) and the eval contraction in `fpbisp` (`2168-2182`): nested `do i=1,mx; do j=1,my` → `d`-nested / recursive contraction, cost `O(kᵈ)` per point. Care needed with the half-splitting boundary attribution (`fac = term*half`) for the per-axis residual sums driving knot placement.

---

## 3. Why the scattered path does *not* generalize

`fpsurf` (`fitpack_core.F90:14021`) assembles a **single fused observation matrix** `a(nc, ib1)` with `nc = nk1x·nk1y` (full tensor coefficient count) and coupled bandwidths `ib1`/`ib3`, plus:
- spatial **panel/cell indexing** (`index(nrest)`, `nummer(m)`, `nreg = nxx·nyy`) mapping data points to knot cells,
- **rank-deficiency handling** (`rank`, `fprank`, Schoenberg–Whitney checks),
- a **direction-interchange** optimization to minimize bandwidth.

None of this decomposes into independent 1-D solves; bandwidth and panel structure are intrinsically `d`-dimensional and grow combinatorially. Extending it to N-D is research-grade, not templating. **Explicitly scope out** (document the boundary, as bspline-fortran does by pointing scattered/LSQ users to SPLPAK/FINTERP).

---

## 4. bspline-fortran comparison (reference: github.com/jacobwilliams/bspline-fortran)

| Aspect | bspline-fortran | FITPACK (this repo) |
|---|---|---|
| Dimensions | **1D–6D** (`db1ink…db6ink`, `db1val…db6val`, `bspline_1d…6d`) | 1D + 2D today |
| Problem | **pure tensor-product interpolation** (passes through every grid node) | **smoothing** with `s` tradeoff + **automatic knot insertion** |
| Knots | data-defined or user-supplied; cannot place *fewer* knots than data | automatic, residual-driven (the valuable part) |
| Scattered | **none** (grid-only) | `surfit` (2-D) |
| Core math | de Boor alternating-direction (`dbtpcf` setup, nested `dbvalu` eval) — **same factoring as `fpgrre`** | same de Boor tensor product |
| **Templating** | **none** — hand-written ladder, 6 copies of each routine; no fypp/codegen. *That* is why modules are 140–250 KB and why it stops at 6D | — |

**Takeaways:**
1. The *numerical* scaling is a solved problem — the de Boor factoring scales cleanly to 6D. What blew up there was **source duplication** (inherited from F77, never re-architected). Design that out.
2. bspline-fortran has **none** of FITPACK's value-add in N-D (no smoothing, no auto-knots, no scattered). So an N-D FITPACK *gridded smoother* would be something that **does not currently exist** in the Fortran ecosystem — a strong reason to build it, and to build it on the gridded path.

---

## 5. Architecture decision

Separate **logic** from **face**:

| Layer | Form | Count |
|---|---|---|
| **Core logic** | dimension-generic *procedures* — `regrid_nd(d, m, x_flat, k, …)`, `fpregr_nd`, `fpgrre_nd`, `fpbisp_nd` — runtime `d`, flat arrays | **one** implementation |
| **Public face** | concrete `fitpack_gridded_3d … _Nd`, rank-natural `z(:,:,:…)`, type-safe `eval_ongrid` | N **types**, **one** fypp template |

### 5.1 Avoid PDT (`fitpack_grid(d)`)
- **Compiler reality.** PDT `len`/`kind` parameters + allocatable components + type-bound procedures are exactly gfortran's fragile corner (long-standing ICEs/miscompiles). This repo already works around gfortran-13 codegen bugs (`surf_new_points`: *"Do not use sum() or it will segfault gfortran 13"*, `fitpack_gridded_surfaces.f90:296`). Don't bet the architecture on the toolchain's weakest feature.
- **Ergonomic friction.** A `len` parameter propagates into every declaration and complicates the C-interop layer.

### 5.2 Avoid a single runtime-`ndim` type for the public API — the rank problem
A spline library's UX is "give me rank-`d` data, get rank-`d` values." Today: `z(:,:)` on `fitpack_grid_surface`. Natural extension: `z(:,:,:)` (3-D), `z(:,:,:,:)` (4-D) … (Fortran allows rank ≤ 15, so all of `MAX_IDIM=10` is legal).

A **single** `fitpack_ndgrid` cannot hold an allocatable whose rank varies at runtime → it must store `z` **flattened** (`z(:)` + shape vector `m(d)`). Internally fine (the core already flattens: `regrid` takes `z(mx*my)`), but at the **API boundary** it regresses:
- user must flatten/reshape and carry shape metadata by hand,
- `eval_ongrid` returns flat `f(:)` not `f(:,:,:)`,
- inconsistent with the shipped 2-D `z(:,:)` type,
- loses compile-time safety (4-D data into a 3-D fit).

`select rank` (F2018) does **not** rescue it — it is *itself* a per-rank ladder (`rank(3); rank(4); …`) and gfortran support is shaky.

### 5.3 Recommended: generic core + fypp-generated concrete faces
The concrete faces are **pure marshalling**: pack rank-`d` arrays → flat → call the generic core → unpack. Identical except for rank → exactly what fypp emits from one loop:

```fypp
#:set DIMS = range(3, MAX_IDIM + 1)
#:for d in DIMS
  ! type fitpack_gridded_${d}$d with z(${ ':,' * (d-1) }$:), order(${d}$), knots(${d}$), t(:,${d}$), left(${d}$), right(${d}$)
  ! new_points / fit / eval_ongrid / dfdx … all delegate to regrid_nd(d=${d}$, …)
#:endfor
```

So "many classes" costs **one template**, not N hand-written copies (the bspline-fortran trap). You keep:
- rank-natural ergonomic I/O,
- compile-time dimension safety,
- consistency with the existing pattern (`order(2)`, `knots(2)`, `t(:,2)`, `left(2)/right(2)`, base `fitpack_fitter` carrying `wrk/lwrk/iwrk/c/fp/iopt/smoothing`, `fitpack_fitters.f90:37`).

Keep `fitpack_grid_surface` as the 2-D name (back-compat); add `fitpack_gridded_3d …`.

**When to flip to the single type instead:** only if (a) `d` is genuinely runtime-variable (rare for a fitter — users know their dimensionality at author time), or (b) the "no preprocessor" rule is non-negotiable, in which case a single flat-`z` runtime type is the only way to ship N-D without N hand-written types (trading ergonomics for zero codegen). A middle path — internal `fitpack_ndgrid_core` (flat, runtime `d`) *composed into* generated rank-typed faces — is possible but mild over-engineering; faces calling the generic procedures directly is simpler.

> **Note on CLAUDE.md:** the project currently states *"No preprocessor directives."* fypp is a build-time *source generator*, not `#ifdef` conditional compilation, but introducing it is a deliberate departure — decide consciously and update the doc. fpm supports an fypp preprocessing step.

---

## 6. What must change (2-D-coupling inventory)

The paired-scalar pattern is replicated widely; mechanical but it is the bulk of the line count (~4× the 1-D code). fypp earns its keep here.

**OOP types** (all share `order(2)`, `knots(2)`, `t(:,2)`, `nest(2)`, `left(2)/right(2)`):
- `fitpack_surface` (`fitpack_surfaces.f90:45`) — scattered, **not** generalized here
- `fitpack_grid_surface` (`fitpack_gridded_surfaces.f90:50`) — **the template seed**, `z(iy,ix)`
- `fitpack_parametric_surface`, `fitpack_polar`, `fitpack_grid_polar`, `fitpack_sphere`, `fitpack_grid_sphere` — specialized 2-D domains, separate effort

**Core evaluation/analysis (gridded), all hard-wired to 2 directions:**
| Routine | `fitpack_core.F90` | Generalize to |
|---|---|---|
| `bispev` / `fpbisp` | 400 / 2119 | `bispev_nd` / `fpbisp_nd` (per-axis weight tables → `d`-contraction) |
| `bispeu` | 342 | scattered-point eval over `d` axes |
| `parder` / `pardeu` / `pardtc` | 15275 / 15413 / 15549 | per-axis derivative passes → loop |
| `profil` | 16823 | cross-section: fix one of `d` axes |
| `dblint` | 1294 | `nd`-box integral |

**C bindings** (`fitpack_core_c.f90`): ~10 surface-related entry points (`regrid_c`, `bispev_c`, `bispeu_c`, `parder_c`, `surfit_c`, `parsur_c`, `polar_c`, `pogrid_c`, `sphere_c`, `spgrid_c`) — opaque-pointer types assume 2-D; generate the N-D ones.

**Fixed-size work arrays** sized by `MAX_ORDER = 19` (`fitpack_core.F90:111`): `h(MAX_ORDER+1)` and *paired* `hx/hy`, `hu/hv`, `hp/ht` (`:2130, :10753, :13178, :14050`) → per-axis `h(:, d)` or computed in a loop.

---

## 7. Caveats

1. **Curse of dimensionality is the real ceiling.** Coefficient count `∏(nᵢ − kᵢ − 1)`; eval `O(kᵈ)`/point. Memory-bound well before `d = 10` (why bspline-fortran stops at 6D, and even 6D is rarely usable). Realistic: **d = 3–4, maybe 5.** State this in docs.
2. **`idim` ≠ `ndim`** (§1). Don't conflate the parametric codomain dimension with the new domain dimension.
3. **Breadth.** ~7 OOP types + ~10 C entry points + the eval routines carry the 2-D pattern; the marshalling layer is most of the work even though the algorithm is the interesting part.

---

## 8. Refactor strategy — parallel implementation gated by a d=2 oracle

The method: build N-D versions *alongside* the 2-D ones, parameterized by `dims`, and gate every step on **bit-for-bit equivalence with the existing 2-D code at dims=2**. This separates *mechanical correctness* (the gate) from *new capability* (dims>2).

### The four moves
1. **Duplicate** the 2-D gridded vertical slice into `*_nd` routines (`regrid → regrid_nd`, `fpregr → fpregr_nd`, `fpgrre → fpgrre_nd`, `fpbisp → fpbisp_nd`). Keep the originals untouched.
2. **Dimensionalize at dims=2.** Turn paired scalars/arrays into dimensioned ones (`nx,ny → n(dims)`, `kx,ky → k(dims)`, `tx,ty → t(:,dims)`, offsets → a cumulative loop), add `integer(FP_DIM), intent(in) :: dims`, size arrays by `dims` (`nc = product(nest-k-1)`, `mz = product(m)`). **Run it at dims=2.** This is *not* "stay 2-D" — it is "parameterize, exercise at 2."
   - **Introduce the general index helpers here, at dims=2** (linear↔multi-index, e.g. `fp_grid_index`). If step 2 keeps 2-D-special indexing (`z((ix-1)*my+iy)`, hardcoded `lsx/lsy`), the gate validates code that step 4 must rip out and the N-D path is never tested. Generalize *forward* — only the three kernels in §2 (alternating-direction solve, the two contractions) stay literally-2 until step 4.
3. **Equivalence gate.** Standalone harness `test_nd_equivalence`: run `regrid` vs `regrid_nd(dims=2)` over the existing 2-D surface battery, **all `iopt` modes (−1/0/1) and s=0**, assert equality of `tx,ty,c,fp,nx,ny`.
   - **Aim bit-for-bit**, not tolerance: it catches floating-point *reassociation* bugs a tolerance test sleeps through. The cost is constraining collapsed loops to iterate axes / accumulate in the *same order* as the 2-D code at dims=2. Use tolerance only where reassociation is genuinely unavoidable. This harness is the backbone of the whole refactor.
4. **Extend to dims>2.** Three categories of code, three treatments:
   - *Per-axis-independent machinery* (observation matrices, knot init/placement, validation, workspace offsets): a **runtime `do idim=1,dims` collapse** of the two existing copies. Trivial.
   - *The three genuinely d-dimensional kernels* (alternating-direction solve + `d−1` transpose reshuffles, residual contraction, eval contraction): **runtime multi-index addressing**, still no fypp. The 2-D `q(mynx)` *is* the single transpose buffer; at dims>2 it generalizes to the ping-pong buffer of §9.
   - *Array-rank-dependent code* (the public faces `z(:,:,:…)`, rank-d↔flat marshalling): **fypp** — rank is compile-time and cannot be a runtime loop.
   Then `new_knot_dimension` (binary) → **argmax over `dims`**; confirm the `s`-tradeoff converges in 3-D.

### fypp boundary
**Core = runtime `dims` loop (no fypp)** — this is the bspline-fortran-beating property; everything including the hard kernels is runtime. **fypp only at the boundary**: concrete faces `fitpack_gridded_3d…_Nd` and rank-d marshalling. Defer fypp'd unrolled contraction kernels unless profiling later demands them.

### Sequencing — vertical slice, not breadth-first
Slice 1: `regrid_nd → fpregr_nd → fpgrre_nd → fpbisp_nd` (eval is needed *by* the gate). Green gate at dims=2, then dims=3 interpolation+smoothing test. **Only then** duplicate the peripherals: `parder_nd`, `dblint_nd` (N-D box integral), `profil`/cross-section, comm pack/expand, C bindings, and the fypp faces.

### Endgame (optional consolidation)
Once the dims=2 gate is permanently green, re-point `fitpack_grid_surface` at `regrid_nd(dims=2)` and **retire `regrid`/`fpregr`/`fpgrre`** — net code reduction. The gate is what licenses deleting the originals.

### Notes
- `FP_DIM`: alias to `c_int32_t`/`FP_SIZE` for C-interop; validate `1 ≤ dims ≤ MAX_IDIM`.
- Leave scattered N-D (`surfit`/`fpsurf`) explicitly unsupported (SPLPAK-style scope boundary).

---

## 9. Hardest part to get right

The **alternating-direction ping-pong buffer** in `fpgrre_nd`: the intermediate RHS tensor that, after solving axis `i`, holds (coefficients along axes `1..i`) × (data along axes `i+1..d`). Sizing is `max over i of [∏_{j≤i}(nⱼ−kⱼ−1) · ∏_{j>i} mⱼ]` and the stride pattern flips each pass. The 2-D code hides this as the single `q(mynx)` with `mynx = nxest*my`. Reference implementation to mirror: bspline-fortran `dbtpcf`.

---

## References
- Dierckx, *Curve and Surface Fitting with Splines*, Oxford UP, 1993 — Ch. 5 §5.4 (regrid).
- Dierckx, "A fast algorithm for smoothing data on a rectangular grid…", SIAM J. Numer. Anal. 19 (1982) 1286–1304.
- bspline-fortran (Jacob Williams) — de Boor tensor-product interpolation, 1D–6D, hand-written ladder; NIST CMLIB DBSPLIN/DTENSBS provenance.
- This repo: `regrid` `:16958`, `fpregr` `:11907`, `fpgrre` `:6940`, `fpbisp` `:2119`, `fpsurf` `:14021`, refactored rotations `:5859`/`:5901`/`:5946`.
