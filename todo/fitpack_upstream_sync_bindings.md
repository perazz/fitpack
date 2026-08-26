# Upstream Sync, SciPy-Ported Fixes, and Generated C/C++ Bindings — Implementation Plan

**Status:** PLANNED (2026-08-26).
**Scope:** (a) port the concrete FITPACK algorithm fixes SciPy accumulated during its C rewrite; (b) replace the hand-written C/C++ curve bindings with output from the fortran-arrays bindings generator, extended with a standalone/gated mode so this repository stays free of fortran-arrays dependencies; (c) close "slice 10" of [fitpack_nd_grids.md](fitpack_nd_grids.md) — C bindings for `fitpack_gridded_spline`; (d) supersede PR #61 with a dimension-templated `fpPoint<dims>` point type across the whole parametric family.
**Out of scope (deferred):** SciPy's architectural ideas — knot selection as a decoupled generator (would retire `iopt` cold/warm starts), 5-way granular `fpchec` error codes, `nest` as a first-class stopping criterion. Worth a later arc of their own.

---

## 1. Background audit (2026-08-26)

- `main` @ `f4503d8` (2026-07-01, N-D gridded splines). All remote branches are fully merged (verified 0-ahead) except `tvdp-tree-node-type` (PR #39, closed unmerged — derived-type tree refactor for `fpadno`/`fpdeno`/`fpfrno`/`fpseno`; keep for a future decision) and the trivial `docs/doxygen-mathjax` fix.
- Open items: **PR #61** (`fpConstrainedCurveTemplate<N>`, `std::array`-based fixed-dim points, constrained-curve only), **issue #30** (const-correctness of the C++ wrappers), **issue #11** (no golden-reference comparison against the F77 library — every fix ported here must carry its own regression test).
- The C/C++ interface covers only the 5 curve types (91 `fitpack_<type>_c_*` method symbols + 25 `fp_*_c` core procedural symbols + 2 unprefixed helpers = 118). Surfaces, polar, sphere, gridded, and the new N-D `fitpack_gridded_spline` have no C ABI.
- Known C++-layer defects that retire with the hand-written layer: `fpParametricCurve::interpolate(order)` ignores its argument (`fpParametricCurve.hpp:83`); `typedef std::vector<FP_REAL> fpPoint` + `using std::vector` leak C++ into `fitpack_core_c.h` (a nominally C header) and cost a heap allocation per point.
- SciPy completed its Fortran-FITPACK removal in 1.18.0 (PR scipy#24022). Its issue tracker yields the concrete portable fixes in §2. Notably, SciPy has **not** fixed the zero-residual `fpknot` degeneracy this repo fixed in PR #57 — worth reporting to them.
- The fortran-arrays bindings generator already has a `--standalone` mode and a 13-type fitpack prototype (`fitpack_test/`, 89/89 tests at the time), but the prototype is stale and standalone mode has three gaps (§3).

## 2. Phase 1 — SciPy-ported algorithm fixes (PR train; one test per fix)

All in `src/fitpack_core.F90` unless noted:

| # | Fix | Source |
|---|---|---|
| 1 | `pardtc`: x-derivative order compared against `ky` instead of `kx` — re-check in the dimension-generic rewrite | scipy#24084 |
| 2 | `regrid`/`fpregr`: NaN output when knots saturate in one direction; make `maxit` (hardcoded 20) configurable — re-verify on the new N-D `regrid`/`fpgrre` path | scipy#22433 |
| 3 | `fpsurf`: audit for uninitialized `acc` (same class as the PR #57 `fpknot` fix) | scipy#24270 |
| 4 | `splder`/`spalde`/`parder`: guard derivative order `nu > k` (return zeros) | scipy#24076 |
| 5 | Precondition `count_nonzero(w) > k`, not `m > k`, in `curfit`/`percur`/`parcur` (`w=0` is the documented point-exclusion mechanism) | scipy#23755 |
| 6 | `fpback`/`fpbacp`: zero-pivot guard on rank-deficient systems (extend the PR #48 `percur` pattern) | scipy#22704 (open everywhere, incl. the F77 original) |
| 7 | Clamp the initial-knot estimate against integer overflow (`nest`/`nxest`/`nyest` sizing); coherent `k < 1` rejection | scipy#25726, scipy#25381 |

## 3. Phase 2 — Generator work (fortran-arrays repo)

Design decision: **one standalone-pure C ABI; fxArray support gated C++-side only.** The generated `.f90`/`.h` never reference `array_c`. A generated `include/fitpack_config.h` sets `HAVE_FXARRAY` via `__has_include("fxArrays.hpp")` with a `FITPACK_NO_FXARRAY` opt-out and a `constexpr fitpackHasFxArray()` query (the wxPlotLib pattern). Zero-copy `fxArray<T>` views are built client-side from raw `ptr + extents` C getters — no fortran-arrays symbol in this repo's ABI.

Gaps to close in the generator (all changes gated so non-standalone output stays byte-identical):
1. `standalone:` YAML key (`false | true | "gated"`) — today CLI-only.
2. Self-contained standalone `.f90`: swap the unconditional `use arrays_c` import for a generated `fitpack_fx_status` module (distinct name — downstream projects compile these sources next to the real `arrays_c`).
3. Component array accessors (today dropped in standalone): raw `getcomp_<name>_raw(this, ptr, extents)` symbols; C++ emits `std::vector<T> <name>_vector()` always and zero-copy `fxArray<T> <name>()` under `#if HAVE_FXARRAY`.
4. Rank-2 raw argument variants, assumed-rank→rank-1 mapping, and object-state-sized result buffers (`y(idim)`, `zeval(product(m))`) — **the parity blocker**: without these, parametric `new_fit`/`eval` are silently skipped.
5. Plain-name raw symbols in standalone mode (no `_raw` suffix when the descriptor primary is filtered) — restores most of the 91 legacy method names.
6. Export-header branding parametrized from `project_name` (`FITPACK_CAPI_EXPORT`).
7. Refresh `fitpack_test/` against current sources; its test suite is the regression gate, plus a compile-only TU with `HAVE_FXARRAY` defined.

## 4. Phase 3 — Bindings adoption in this repository

- Commit `fitpack_bindings.yaml` + `extra_methods/` fragments + `generate_bindings.sh`; generate into `src/` + `include/`. **Delete** the 5 hand-written type shims (`fitpack_curves_c.f90`, `fitpack_periodic_curves_c.f90`, `fitpack_parametric_c.f90`, `fitpack_closed_c.f90`, `fitpack_constrained_c.f90`) and their headers. **Keep `fitpack_core_c.f90`/`fitpack_core_c.h` hand-written** — the 25 `fp_*_c` + 2 helper symbols are preserved exactly; regenerating them is all risk, no benefit. Remove the `fpPoint` typedef + `flatten_2d_vector` from `fitpack_core_c.h`.
- **`fitpack_gridded_spline` bindings (slice 10)**: extends the bound abstract-fitter base; assumed-rank `z(..)` uses the flat entry path with explicit `m`; `cross_section`/`derivative_spline` derived returns are supported; fixed-size `MAX_IDIM` metadata component accessors deferred (known generator gap).
- **`fpPoint<dims>`** (`include/fpPoint.hpp`, hand-maintained): `std::array<FP_REAL,dims>`-backed standard-layout POD, `static_assert`ed byte-compatible with `real(FP_REAL) :: x(dims)`; `std::vector<fpPoint<dims>>` is bit-identical to Fortran `x(dims,m)` and feeds rank-2 raw calls with zero copies. Dimension-templated method splices (runtime `idim() == dims` check) on `fpParametricCurve`/`fpClosedCurve`/`fpConstrainedCurve` (`new_fit`/`eval`/`ddu`/`ddu_all`/`constrain_*`) and an N-D analog on `fpGriddedSpline`. This supersedes PR #61, generalizing its idea (credit to @illionj) from one parallel class to the whole family with no type duplication.
- **Ergonomic splices** restore the established call shapes (`new_fit(x, y[, w], smoothing = 1000)`, the `fit` trio, `interpolate([order])`, range-`eval`, `ddx`, `integral`) so `test/test_curve.cpp` compiles with minimal edits — that compile is the source-compat test.
- **ABI parity report**: a committed script diffs `bind(c, name=...)` symbols old→new (+ `nm` on the built library). Expected: 27 core symbols identical; lifecycle/`new_fit`/`new_points`/`fit` names identical with signature deltas (`fx_status*` trailing args; `_copy` gains `deep_copy`); ~1/3 generic-specific renames (e.g. `_eval_one` → `_curve_eval_one`); ~60 new symbols (8 newly bound types + accessors). The report is the release-notes break list.
- Issue #30: verify const-correctness of the generated headers; close if satisfied.
- Docs: README bindings table (add surface/N-D rows, `fpPoint` migration note), version bump to 0.3.0 carrying the documented breaks (opaque handle 8→24 bytes; status/deep-copy arguments; `fpPoint` type change; per-type header filenames). Verify the build sets C++17 for the generated headers.

## 5. Phase 4 — Downstream adoption (FRESCO)

Consumers sync the repo (subtree), drop their vendored copy, and compile these sources directly; the generated headers auto-enable their fxArray layer via `fitpack_config.h`. The one breaking OO change since their last sync (`fitpack_grid_surface%eval` → `eval_ongrid`) has no consumers there. Details tracked downstream.

## 6. Verification

- `fpm build` + `fpm test` on gcc (normal, `-fcheck=bounds`, `-finit-real=inf` poison — the nd_grids gate discipline), `test_curve.cpp`, examples, CI green on every phase PR.
- Generator: its pytest suite (non-standalone output byte-identical) + the refreshed `fitpack_test` gate.
- ABI parity script output committed with the adoption PR.

## 7. Risks

1. Rank-2 raw / assumed-rank / state-sized results are the largest new generator surface; hard-skips are silent — diff the verbose skip report against an expected-methods checklist (parametric `new_fit`/`eval` are parity blockers).
2. The ABI break is real (24-byte handle, status args, `deep_copy`, renamed generic specifics). No known binary-only consumers; loud header comments + 0.3.0 release notes.
3. Gated `fxArray` component views alias allocatables that refitting reallocates — document "views are invalid after any fitting call".
4. Old C callers of retained plain-name symbols would get garbage `present()` checks rather than compile errors if they pass the old argument counts — the signature break must be loud in the headers.
5. Branch hygiene: prune the ~30 fully-merged remote branches; keep `tvdp-tree-node-type` and `gh-pages`.
