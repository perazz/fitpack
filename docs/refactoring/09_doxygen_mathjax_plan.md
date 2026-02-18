# PR 9: Doxygen + MathJax Documentation

## Context

All 8 code refactoring PRs are complete. The next step per `docs/refactoring/00_refactoring_summary.md`
is to add structured Doxygen documentation to the entire library. The documentation convention is
defined in `docs/refactoring/03_doxygen_convention.md`. The goal is a professional GitHub Pages site
(like [fortran-lapack](https://perazz.github.io/fortran-lapack/)) with dark mode, MathJax equations,
and cross-referenced API docs.

This is a large PR — it touches every source file. The plan breaks it into sequential,
independently-testable steps.

---

## Phase A: Infrastructure Setup

### A1. Create `project/doxygen/` directory with theme files

Copy from `~/code/fortran-lapack/project/doxygen/`:
- `Doxyfile` — adapt for fitpack (project name, input paths, warn log path)
- `header.html` — unchanged (dark mode toggle init)
- `doxygen-awesome.css` — unchanged (2,681 lines)
- `doxygen-awesome-darkmode-toggle.js` — unchanged (258 lines)

Key Doxyfile changes from fortran-lapack:
```
PROJECT_NAME           = fitpack
INPUT                  = "../../src/"
INPUT                 += "../../doc/mainpage.md"
USE_MDFILE_AS_MAINPAGE = "../../doc/mainpage.md"
EXTENSION_MAPPING      = .F90=FortranFree
EXAMPLE_PATH           = "../../example/"
WARN_LOGFILE           = ""
EXTRACT_PRIVATE        = NO
CALL_GRAPH             = YES
CALLER_GRAPH           = NO
```

### A2. Create `doc/mainpage.md` — Doxygen landing page

High-level overview page (NOT the README — Doxygen main page is separate):
- Project description and purpose
- Links to type hierarchy overview
- Quick-start code examples
- Links to book reference

### A3. Create `.github/workflows/deploy-docs.yml`

Copy from fortran-lapack, same structure:
- Trigger on push to `main`
- Install doxygen 1.13.2, texlive-full, ghostscript, graphviz
- Run `cd project/doxygen && doxygen`
- Create `.nojekyll`
- Deploy to `gh-pages` via `peaceiris/actions-gh-pages@v3`

### A4. Add `project/doxygen/html/` to `.gitignore`

### A5. Verify: run `doxygen` locally, open HTML, confirm theme renders

---

## Phase B: Core Routine Documentation (`fitpack_core.F90`)

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

| Routine | Line | Book ref | Key equations |
|---------|------|----------|---------------|
| `fpgivs` | 5500 | §4.1.2, pp.55-58 | 4.15, 4.16 |
| `fprota` | 11402 | §4.1.2, pp.55-58 | 4.15 |
| `fpback` | 2368 | §4.1.2, pp.55-58 | 4.14 |
| `fpbacp` | 2406 | §6.1, pp.95-100 | 6.6 |
| `fpbspl` | 2721 | §1.3, pp.8-11 | 1.22, 1.30 |
| `fpdisc` | 5356 | §5.2.2, pp.76-79 | 5.5, 5.6 |
| `fprati` | 10973 | §5.2.4, pp.83-86 | 5.30-5.32 |
| `fpchec` | (in fpcons area) | §4.1.1, pp.53-55 | 4.5, 4.17 |
| `fpbisp` | 2648 | §2.1.2, pp.28-30 | 2.14-2.17 |
| `fporde` | 8121 | §9.1, pp.147-150 | — |

### B2. Refactored Givens rotation helpers (9 routines)

These already have Doxygen comments — review and enhance with book refs:

| Routine | Line | Book ref |
|---------|------|----------|
| `fp_rotate_row` | 5537 | §4.1.2, Eq.4.15 |
| `fp_rotate_row_vec` | 5569 | §4.1.2, Eq.4.15 |
| `fp_rotate_shifted` | 5608 | §5.2.2, Eq.5.15 |
| `fp_rotate_shifted_vec` | 5642 | §5.2.2, Eq.5.15 |
| `fp_rotate_row_2mat_vec` | 5675 | §6.1, Eq.6.4 |
| `fp_rotate_row_2mat` | 5716 | §6.1, Eq.6.4 |
| `fp_rotate_row_block` | 5758 | §10.2, Eq.10.4-10.8 |
| `fp_rotate_row_stride` | 5783 | §10.2, Eq.10.8 |
| `fp_rotate_2mat_stride` | 5808 | §10.2, Eq.10.8 |

### B3. Knot/tree/polynomial utilities (9 routines)

| Routine | Line | Book ref | Key equations |
|---------|------|----------|---------------|
| `fpknot` | 7632 | §5.3, pp.87-94 | 5.37-5.43 |
| `fprank` | 10738 | §9.1.2, pp.150-152 | 9.8-9.10 |
| `fpinst` | 7426 | §4.2, pp.63-65 | — |
| `fpintb` | 7494 | §3.2, pp.44-46 | 3.8-3.10 |
| `fpadno` | 2229 | §7.2, pp.125-130 | — |
| `fpdeno` | 5266 | §7.2, pp.125-130 | — |
| `fpfrno` | 5405 | §7.2, pp.125-130 | — |
| `fpseno` | 11560 | §7.2, pp.125-130 | — |
| `fpcuro` | 5071 | — | cubic polynomial roots |

### B4. Cyclic/polynomial helpers (5 routines)

| Routine | Line | Book ref |
|---------|------|----------|
| `fpcyt1` | 5183 | §6.1, pp.95-100, Eq.6.5-6.6 |
| `fpcyt2` | 5232 | §6.1, pp.95-100 |
| `fppocu` | 9478 | §7.1, pp.115-120 |
| `fpcoco` | 3670 | §7.2, pp.125-130 |
| `fpcosp` | 4275 | §7.1, pp.115-120 |

### B5. Core fitting algorithms (9 routines)

The big ones — these need the most detailed documentation with full equation chains.

| Routine | Line | Book ref | Key equations |
|---------|------|----------|---------------|
| `fpcurf` | 4717 | §4.1-5.3, pp.53-94 | 4.12, 5.10, 5.37 |
| `fpcons` | 3864 | §8.2, pp.141-146 | — |
| `fpclos` | 3031 | §6.1-6.2, pp.95-112 | 6.1-6.8 |
| `fppara` | 8156 | §6.3, pp.112-114 | 6.9 |
| `fpperi` | 8906 | §6.1-6.2, pp.95-112 | 6.1-6.8 |
| `fpsurf` | 12808 | §9.1-9.2, pp.147-167 | 9.2-9.6 |
| `fppola` | 9954 | §11.1, pp.197-205 | 11.1-11.9 |
| `fpsphe` | 12044 | §11.2, pp.205-213 | 11.12-11.16 |
| `fppasu` | 8514 | §10.3, pp.173-178 | 10.9-10.12 |

### B6. Grid fitting algorithms (6 routines)

| Routine | Line | Book ref | Key equations |
|---------|------|----------|---------------|
| `fpgrre` | 6655 | §10.2, pp.170-172 | 10.4-10.8 |
| `fpgrdi` | 5844 | §10.2, pp.170-172 | 10.4-10.8 |
| `fptrnp` | 13488 | §10.2, pp.170-172 | 10.8 |
| `fptrpe` | 13574 | §10.2, pp.170-172 | 10.8 (periodic) |
| `fpgrsp` | 6913 | §11.2, pp.205-213 | — |
| `fpgrpa` | 6316 | §10.3, pp.173-178 | — |

### B7. Specialized helpers (7 routines)

| Routine | Line | Book ref |
|---------|------|----------|
| `fpcsin` | 4664 | §3.3, Fourier integrals |
| `fpbfou` | 2469 | §3.3, Fourier coefficients |
| `fpadpo` | 2301 | §7.1, polynomial spline |
| `fprppo` | 11423 | §11.1, polar representation |
| `fprpsp` | 11498 | §11.2, spherical representation |
| `fpopdi` | 7740 | §11.1, polar derivatives |
| `fpopsp` | 7890 | §11.2, spherical operations |
| `fpsysy` | 13429 | — symmetric system solve |
| `sort_xy` | 13386 | — reorder scattered data |
| `fppogr` | 9548 | §11.1, polar grid fitting |
| `fpspgr` | 11596 | §11.2, spherical grid setup |
| `fpregr` | 11005 | §10.2, regression residuals |

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

| Routine | Line | Notes |
|---------|------|-------|
| `fitpack_error_handling` | 214 | Error handler |
| `FITPACK_MESSAGE` | 228 | Error code → string |
| `FITPACK_SUCCESS` | 292 | Check success |
| `get_smoothing` | 258 | Smoothing parameter logic |
| `IOPT_MESSAGE` | 280 | iopt → string |
| `fp_knot_interval` | 2756 | Already has Doxygen — review |
| `fitpack_argsort` | 18151 | Index sort |
| `equal`, `not_equal` | 18294 | Float comparison |
| `swap_data`, `swap_size` | 18255 | Element swap |
| `FP_*_COMM_*` | 18308+ | Already have Doxygen — review |

---

## Phase C: OOP API Documentation

Each OOP type gets a Doxygen comment block on the `type` declaration and on each
public method. These reference the book but focus on **usage**: what the type represents,
how to use each method, what the parameters mean.

### C1. Abstract base: `fitpack_fitter` (`src/fitpack_fitters.f90`)

Document: type overview, `mse()`, `fitting_status()`, `destroy_base()`, comm interfaces.

### C2. Curve types (`src/fitpack_curves.f90`)

- `fitpack_curve` — full documentation: constructor, `new_fit`, `fit`, `interpolate`,
  `least_squares`, `eval`, `dfdx`, `dfdx_all`, `integral`, `zeros`,
  `fourier_coefficients`, `insert_knot`, `destroy`, comm methods
- `fitpack_periodic_curve` — brief (marker type, inherits everything)

### C3. Convex curve (`src/fitpack_convex_curves.f90`)

- `fitpack_convex_curve` — `set_convexity`, `fit` (concon), `least_squares` (cocosp)

### C4. Parametric curves (`src/fitpack_parametric_curves.f90`)

- `fitpack_parametric_curve` — constructor, `eval` (returns multi-dim), `dfdx`, `dfdx_all`
- `fitpack_closed_curve` — brief
- `fitpack_constrained_curve` — `set_constraints`, derivative boundary conditions

### C5. Surface types (`src/fitpack_surfaces.f90`, `src/fitpack_grid_surfaces.f90`)

- `fitpack_surface` — full: `eval` (scattered, bispeu), `eval_ongrid` (grid, bispev),
  `dfdx`, `dfdx_ongrid`, `integral`, `cross_section`, `derivative_spline`
- `fitpack_grid_surface` — same methods, note eval/eval_ongrid distinction

### C6. Parametric surface (`src/fitpack_parametric_surfaces.f90`)

- `fitpack_parametric_surface` — `eval`, periodicity, `idim`

### C7. Polar types (`src/fitpack_polar.f90`, `src/fitpack_gridded_polar.f90`)

- `fitpack_polar` — boundary function pointer, origin continuity, eval in Cartesian
- `fitpack_grid_polar` — origin BC, grid structure

### C8. Sphere types (`src/fitpack_spheres.f90`, `src/fitpack_gridded_sphere.f90`)

- `fitpack_sphere` — theta/phi coordinates, pole handling
- `fitpack_grid_sphere` — pole BC, grid structure

### C9. Umbrella module (`src/fitpack.f90`)

Document the module itself: what it exports, how to use it.

---

## Phase D: Tutorial / Guide Pages

### D1. `doc/mainpage.md` — Landing page (created in A2)

- What is FITPACK?
- Type hierarchy diagram (text-based)
- Quick-start: fit a curve in 5 lines
- Link to tutorials below

### D2. `doc/tutorial_curves.md` — Curve fitting tutorial

- Basic curve fitting workflow
- Smoothing parameter guidance (from book §5.2.4)
- Periodic curves
- Parametric curves
- Convexity constraints
- Evaluation, derivatives, integration, roots

### D3. `doc/tutorial_surfaces.md` — Surface fitting tutorial

- Scattered vs gridded data
- Smoothing, interpolation, least-squares
- Evaluation (scattered vs grid)
- Cross-sections, integration, derivative splines

### D4. `doc/tutorial_special.md` — Polar & spherical domains

- Polar coordinates, boundary functions
- Spherical coordinates, pole handling
- Grid vs scattered variants

### D5. `doc/book_reference.md` — Book equation index

Quick-reference table mapping FITPACK routines to book chapters/equations.

---

## Phase E: Verification

1. `fpm build && fpm test` — all 59 tests pass (no code changes, only comments)
2. `cd project/doxygen && doxygen` — builds without errors/warnings
3. Open `project/doxygen/html/index.html` — verify:
   - Dark mode toggle works
   - MathJax equations render
   - Type hierarchy visible
   - Cross-references link correctly
   - Tutorial pages accessible
4. Push to branch, verify GitHub Actions deploys to gh-pages

---

## Files Created

| File | Description |
|------|-------------|
| `project/doxygen/Doxyfile` | Doxygen configuration |
| `project/doxygen/header.html` | Custom HTML header (dark mode) |
| `project/doxygen/doxygen-awesome.css` | Theme CSS |
| `project/doxygen/doxygen-awesome-darkmode-toggle.js` | Dark mode JS |
| `.github/workflows/deploy-docs.yml` | GitHub Pages deployment |
| `doc/mainpage.md` | Landing page |
| `doc/tutorial_curves.md` | Curve fitting guide |
| `doc/tutorial_surfaces.md` | Surface fitting guide |
| `doc/tutorial_special.md` | Polar/sphere guide |
| `doc/book_reference.md` | Book equation index |

## Files Modified

| File | Changes |
|------|---------|
| `src/fitpack_core.F90` | Add Doxygen blocks to ~112 routines |
| `src/fitpack.f90` | Add module-level Doxygen comment |
| `src/fitpack_fitters.f90` | Add type/method Doxygen comments |
| `src/fitpack_curves.f90` | Add type/method Doxygen comments |
| `src/fitpack_convex_curves.f90` | Add type/method Doxygen comments |
| `src/fitpack_parametric_curves.f90` | Add type/method Doxygen comments |
| `src/fitpack_surfaces.f90` | Add type/method Doxygen comments |
| `src/fitpack_grid_surfaces.f90` | Add type/method Doxygen comments |
| `src/fitpack_parametric_surfaces.f90` | Add type/method Doxygen comments |
| `src/fitpack_polar.f90` | Add type/method Doxygen comments |
| `src/fitpack_gridded_polar.f90` | Add type/method Doxygen comments |
| `src/fitpack_spheres.f90` | Add type/method Doxygen comments |
| `src/fitpack_gridded_sphere.f90` | Add type/method Doxygen comments |
| `.gitignore` | Add `project/doxygen/html/` |

## Execution Order

Phases A → B → C → D → E. Within Phase B, steps B1-B9 are sequential (each
establishes patterns the next follows). Phase C can partially overlap with late
Phase B. Phase D can be written in parallel with Phase C.

Given the size (~112 core routines + 13 OOP types + 5 guide pages), this will
be split into multiple sub-PRs or done as one large documentation PR.
