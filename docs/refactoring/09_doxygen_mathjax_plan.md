# PR 9: Doxygen + MathJax Documentation

## Context

All 8 code refactoring PRs are complete. The next step per `docs/refactoring/00_refactoring_summary.md`
is to add structured Doxygen documentation to the entire library. The documentation convention is
defined in `docs/refactoring/03_doxygen_convention.md`. The goal is a professional GitHub Pages site
(like [fortran-lapack](https://perazz.github.io/fortran-lapack/)) with dark mode, MathJax equations,
and cross-referenced API docs.

---

## Phase A: Infrastructure Setup — COMPLETE

All infrastructure is in place:
- `project/doxygen/` — Doxyfile, theme CSS, dark mode JS, header
- `doc/mainpage.md` — Landing page
- `.github/workflows/deploy-docs.yml` — GitHub Pages deployment
- `.gitignore` — excludes `project/doxygen/html/`

---

## Phase B: Core Routine Documentation (`fitpack_core.F90`) — COMPLETE

Each routine gets a new Doxygen comment block placed **immediately before** the
`subroutine`/`function` statement. The existing F77-style comments inside the routine body
are left as-is (they serve as inline implementation notes). The new Doxygen block provides:

1. `@brief` — one-line summary
2. **Description** — what the routine computes, written clearly with MathJax equations from the book
3. `@param` — all parameters with `[in]`/`[out]`/`[in,out]` direction
4. **Error flags** — for public API routines
5. **Usage guidance** — for public API routines (from original F77 comment blocks)
6. `@see` — book chapter, section, page range, equation numbers

Math syntax: `\f$ inline \f$`, `\f[ display \f]`, `\tag{N.M}` for book equation numbers.

### B1. Low-level utilities (10 routines)

Small, self-contained routines — good warm-up, establish the pattern.

| Routine | Book ref | Key equations |
|---------|----------|---------------|
| `fpgivs` | §4.1.2, pp.55-58 | 4.15, 4.16 |
| `fprota` | §4.1.2, pp.55-58 | 4.15 |
| `fpback` | §4.1.2, pp.55-58 | 4.14 |
| `fpbacp` | §6.1, pp.95-100 | 6.6 |
| `fpbspl` | §1.3, pp.8-11 | 1.22, 1.30 |
| `fpdisc` | §5.2.2, pp.76-79 | 5.5, 5.6 |
| `fprati` | §5.2.4, pp.83-86 | 5.30-5.32 |
| `fpchec` | §4.1.1, pp.53-55 | 4.5, 4.17 |
| `fpbisp` | §2.1.2, pp.28-30 | 2.14-2.17 |
| `fporde` | §9.1, pp.147-150 | — |

### B2. Refactored Givens rotation helpers (9 routines)

These already have Doxygen comments — review and enhance with book refs:

| Routine | Book ref |
|---------|----------|
| `fp_rotate_row` | §4.1.2, Eq.4.15 |
| `fp_rotate_row_vec` | §4.1.2, Eq.4.15 |
| `fp_rotate_shifted` | §5.2.2, Eq.5.15 |
| `fp_rotate_shifted_vec` | §5.2.2, Eq.5.15 |
| `fp_rotate_row_2mat_vec` | §6.1, Eq.6.4 |
| `fp_rotate_row_2mat` | §6.1, Eq.6.4 |
| `fp_rotate_row_block` | §10.2, Eq.10.4-10.8 |
| `fp_rotate_row_stride` | §10.2, Eq.10.8 |
| `fp_rotate_2mat_stride` | §10.2, Eq.10.8 |

### B3. Knot/tree/polynomial utilities (9 routines)

| Routine | Book ref | Key equations |
|---------|----------|---------------|
| `fpknot` | §5.3, pp.87-94 | 5.37-5.43 |
| `fprank` | §9.1.2, pp.150-152 | 9.8-9.10 |
| `fpinst` | §4.2, pp.63-65 | — |
| `fpintb` | §3.2, pp.44-46 | 3.8-3.10 |
| `fpadno` | §7.2, pp.125-130 | — |
| `fpdeno` | §7.2, pp.125-130 | — |
| `fpfrno` | §7.2, pp.125-130 | — |
| `fpseno` | §7.2, pp.125-130 | — |
| `fpcuro` | — | cubic polynomial roots |

### B4. Cyclic/polynomial helpers (5 routines)

| Routine | Book ref |
|---------|----------|
| `fpcyt1` | §6.1, pp.95-100, Eq.6.5-6.6 |
| `fpcyt2` | §6.1, pp.95-100 |
| `fppocu` | §7.1, pp.115-120 |
| `fpcoco` | §7.2, pp.125-130 |
| `fpcosp` | §7.1, pp.115-120 |

### B5. Core fitting algorithms (9 routines)

The big ones — these need the most detailed documentation with full equation chains.

| Routine | Book ref | Key equations |
|---------|----------|---------------|
| `fpcurf` | §4.1-5.3, pp.53-94 | 4.12, 5.10, 5.37 |
| `fpcons` | §8.2, pp.141-146 | — |
| `fpclos` | §6.1-6.2, pp.95-112 | 6.1-6.8 |
| `fppara` | §6.3, pp.112-114 | 6.9 |
| `fpperi` | §6.1-6.2, pp.95-112 | 6.1-6.8 |
| `fpsurf` | §9.1-9.2, pp.147-167 | 9.2-9.6 |
| `fppola` | §11.1, pp.197-205 | 11.1-11.9 |
| `fpsphe` | §11.2, pp.205-213 | 11.12-11.16 |
| `fppasu` | §10.3, pp.173-178 | 10.9-10.12 |

### B6. Grid fitting algorithms (6 routines)

| Routine | Book ref | Key equations |
|---------|----------|---------------|
| `fpgrre` | §10.2, pp.170-172 | 10.4-10.8 |
| `fpgrdi` | §10.2, pp.170-172 | 10.4-10.8 |
| `fptrnp` | §10.2, pp.170-172 | 10.8 |
| `fptrpe` | §10.2, pp.170-172 | 10.8 (periodic) |
| `fpgrsp` | §11.2, pp.205-213 | — |
| `fpgrpa` | §10.3, pp.173-178 | — |

### B7. Specialized helpers (12 routines)

| Routine | Book ref |
|---------|----------|
| `fpcsin` | §3.3, Fourier integrals |
| `fpbfou` | §3.3, Fourier coefficients |
| `fpadpo` | §7.1, polynomial spline |
| `fprppo` | §11.1, polar representation |
| `fprpsp` | §11.2, spherical representation |
| `fpopdi` | §11.1, polar derivatives |
| `fpopsp` | §11.2, spherical operations |
| `fpsysy` | — symmetric system solve |
| `sort_xy` | — reorder scattered data |
| `fppogr` | §11.1, polar grid fitting |
| `fpspgr` | §11.2, spherical grid setup |
| `fpregr` | §10.2, regression residuals |

### B8. Public API wrappers (30 routines)

These have the longest existing F77 comment blocks — full conversion to Doxygen format.
Preserve ALL original content (parameter descriptions, error codes, usage guidance).

**Curve routines:**
`curfit`, `splev`, `splder`, `splint`, `sproot`, `spalde`, `cualde`, `curev`,
`percur`, `clocur`, `parcur`, `concur`, `fourco`, `insert`, `insert_inplace`,
`cocosp`, `concon`

**Surface routines:**
`surfit`, `bispev`, `bispeu`, `parder`, `pardeu`, `pardtc`, `dblint`, `profil`,
`regrid`, `parsur`, `surev`

**Polar/sphere routines:**
`polar`, `pogrid`, `sphere`, `spgrid`, `evapol`

### B9. Utility/error/comm routines (~15 routines)

| Routine | Notes |
|---------|-------|
| `fitpack_error_handling` | Error handler |
| `FITPACK_MESSAGE` | Error code → string |
| `FITPACK_SUCCESS` | Check success |
| `get_smoothing` | Smoothing parameter logic |
| `IOPT_MESSAGE` | iopt → string |
| `fp_knot_interval` | Already has Doxygen — review |
| `fitpack_argsort` | Index sort |
| `equal`, `not_equal` | Float comparison |
| `swap_data`, `swap_size` | Element swap |
| `FP_*_COMM_*` | Already have Doxygen — review |

---

## Phase C: OOP API Documentation — PARTIALLY COMPLETE

All OOP modules already have module-level `!> @brief` documentation. Remaining work:
type-level and method-level Doxygen blocks for each public type-bound procedure.

### C1–C9: Type/method documentation

| Module | Status |
|--------|--------|
| `fitpack_fitters.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_curves.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_convex_curves.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_parametric.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_surfaces.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_grid_surfaces.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_parametric_surfaces.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_polar.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_gridded_polar.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_spheres.f90` | Module `@brief` done; type/method docs needed |
| `fitpack_gridded_sphere.f90` | Module `@brief` done; type/method docs needed |
| `fitpack.f90` | Module `@brief` done |

---

## Phase D: Tutorial / Guide Pages — COMPLETE

All tutorial and reference pages exist:
- `doc/mainpage.md` — Landing page
- `doc/tutorial_curve.md` — Curve fitting
- `doc/tutorial_periodic_curve.md` — Periodic curves
- `doc/tutorial_parametric_curves.md` — Parametric curves
- `doc/tutorial_convex_curve.md` — Convex curves
- `doc/tutorial_surface.md` — Surface fitting
- `doc/tutorial_grid_surface.md` — Grid surface fitting
- `doc/tutorial_polar.md` — Polar domains
- `doc/tutorial_sphere.md` — Spherical domains
- `doc/book_reference.md` — Book equation index

---

## Phase E: Verification

1. `fpm build && fpm test` — all tests pass (no code changes, only comments)
2. `cd project/doxygen && doxygen` — builds without errors/warnings
3. Open `project/doxygen/html/index.html` — verify rendering
4. Push to branch, verify GitHub Actions deploys to gh-pages
