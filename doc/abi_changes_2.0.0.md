# C ABI changes in 2.0.0

2.0.0 replaces the five hand-written curve binding modules with output from the
[fortran-arrays](https://github.com/perazz/fortran-arrays) binding generator, and extends the
C/C++ interface from the five curve types to all fifteen fitter types. **The C ABI breaks.**
There is no known binary-only consumer; this file is the migration list.

Regenerate it at any time with:

```bash
./scripts/abi_symbols.sh --diff main      # what the generated bindings changed
./scripts/abi_symbols.sh --diff 1.0.0     # everything that changed since the last release
```

## Summary

Against `main` — the state just before the generated bindings landed:

| class | count |
|-------|-------|
| symbols before | 118 |
| symbols in 2.0.0 | 499 |
| identical (name and signature) | 38 |
| same name, changed signature | 39 |
| renamed (old name gone, documented replacement) | 41 |
| new | 422 |

`main` already carries one further ABI change since the 1.0.0 tag, made before this work: the
25 core procedural bindings gained an `fp_` prefix (`curfit_c` became `fp_curfit_c`, and so on
for every `fp_*_c` entry point in `fitpack_core_c.h`). Counting from the 1.0.0 tag instead —
`./scripts/abi_symbols.sh --diff 1.0.0` — folds those in: 13 identical, 39 changed, 66 renamed,
447 new.

The 38 identical symbols are the 25 `fp_*_c` core procedures plus `fitpack_message_c` and
`FITPACK_SUCCESS_c` — all still hand-written in `src/capi/fitpack_core_c.f90`, deliberately not
regenerated — and 11 curve methods whose signature happened to survive unchanged.

## Breaking change: the opaque handle grew from 8 to 24 bytes

```c
/* before */                    /* 2.0.0 */
typedef struct {                typedef struct {
    void* cptr;                     void* cptr;
} fitpack_curve_c;                  bool  is_pointer;      /* non-owning view? */
                                    const char* name_cptr;
                                } fitpack_curve_c;
```

`is_pointer` distinguishes an owning handle from a borrowed view, which is what makes the new
derived-type returns (`cross_section`, `derivative_spline`) and the class hierarchy safe.
Every handle is still initialised from `fitpack_<type>_c_null`; a C caller that declared its
own 8-byte struct must be recompiled.

All the handle structs across the family are layout-identical, which is what lets a derived
type's handle be passed to a parent type's accessor (see the rename table below).

## Breaking change: status and deep-copy arguments

The lifecycle entry points keep their names and gain trailing arguments:

- `allocate`, `destroy`, `move_alloc`, `associate` gain a trailing `fp_status* status` — pass
  `NULL` to keep the old behaviour of stopping on failure.
- `copy` gains both a by-value `bool deep_copy` and the trailing `fp_status*`.

```c
typedef struct fp_status { bool ok; int code; char message[248]; } fp_status;
```

It is declared inline in the generated C headers behind an include-once guard, so no extra
header is needed, and the declaration is C99-clean even with all thirty headers in one
translation unit. Where [fortran-arrays](https://github.com/perazz/fortran-arrays) is on the
include path its own `fx_status` is used instead and `fp_status` becomes an alias of it; the
struct layout is the same either way, so both spellings link against the same binary. The
hand-written `fitpack_core_c.h` does not use it and is unaffected.

## Breaking change: `FITPACK_CAPI_STATIC` on Windows

The generated headers declare the C entry points with `FITPACK_CAPI_EXPORT`, which defaults to
`__declspec(dllimport)` on Windows. This repository ships a static archive, so a Windows
consumer must define `FITPACK_CAPI_STATIC` before including any fitpack header, or the link
fails on `__imp_fitpack_*`. No effect on Linux or macOS.

## Source-level changes after the first bindings drop

These landed on top of the generated-bindings merge, inside the unreleased 2.0.0 window. None
of them touches a symbol name or an argument type, so every table below is unchanged by them:
`const` is not part of the C ABI, and a C++ class name never reaches the archive.

**The status struct is spelled `fp_status`.** 70 C declarations changed `fx_status*` to
`fp_status*`. Same struct, same layout, same symbols — see the section above.

**Constness follows the Fortran receiver's intent.** A method whose Fortran receiver is
`intent(in)` is declared `const <type>_c* self` in C and `const` in C++; anything else is
declared mutable. The previous rule — every function const, every subroutine non-const —
described the return shape, not the effect on the object, and let a caller mutate a fitter
through a `const` handle.

- 66 declarations **lost** `const` on `self`, across all thirteen concrete fitter types:
  `fit` (13), `least_squares` (13), `interpolate` (12), `new_fit` (11), and 17 entry points
  that update the fitter's cached evaluation state — `fitpack_curve_c_curve_eval_one`,
  `fitpack_curve_c_curve_derivative`, `fitpack_curve_c_integral`,
  `fitpack_surface_c_surface_eval_one`, `fitpack_surface_c_surface_derivatives_one`,
  `fitpack_grid_surface_c_gridded_eval_one`, `fitpack_polar_c_polr_eval_one`,
  `fitpack_sphere_c_sphere_eval_one`, `fitpack_convex_curve_c_set_convexity` and their
  siblings on the remaining types.
- 40 declarations **gained** `const`: the `comm_pack` / `core_comm_pack` serializers (27), the
  `_pure` evaluators (4), `zeros` (2), and the derived-type returns `cross_section` (3),
  `derivative_spline` (3) and `grid_eval_many` (1).

The C++ methods carry the same qualifier, so C++ code that called a mutator through a `const`
object or reference now fails to compile — as it should. Two ergonomic wrappers in
`extra_methods/` were `const` for that reason and are not any more: `fpCurve::ddx(x, order,
ierr)` and `fpGridSpline::dfdx(const fpPoint<dim>&, nu, ierr)`.

**`fxFitpackFitter` is `fpFitter`.** The abstract base class now matches the `fp*` family;
`include/fxFitpackFitter.hpp` and `include/fxFitpackFitter_subtypes.hpp` became
`include/fpFitter.hpp` and `include/fpFitter_subtypes.hpp`. The C handle is unchanged
(`fitpack_fitter_c`, `fitpack_fitter_c.h`).

**The Fortran C-API modules moved to `src/capi/`.** Every generated `*_c.f90` / `*_c_types.f90`,
the generated `fitpack_fx_status.f90`, and the hand-written `fitpack_core_c.f90` now live
there. Module names, symbol names and the `include/` layout are untouched, and fpm globs `src/`
recursively, so no consumer sees this.

## Renamed symbols

The generator names a method after the specific procedure that implements it, so generic
resolutions pick up a prefix (`_eval_one` becomes `_curve_eval_one`), and a scalar component is
exposed as a reference accessor (`_degree` becomes `_ref_order`, returning `int32_t*`) rather
than a by-value getter.

A derived type does not re-emit an accessor for a component it inherits: use the parent type's
accessor on the same handle, which is safe because the handle structs are layout-identical and
the Fortran parent occupies offset 0 of the extension. The C++ wrappers do this for you.

| before | 2.0.0 |
|--------|-------|
| `fitpack_closed_curve_c_degree` | `fitpack_parametric_curve_c_ref_order` |
| `fitpack_closed_curve_c_derivative` | `fitpack_closed_curve_c_curve_derivative` |
| `fitpack_closed_curve_c_idim` | `fitpack_parametric_curve_c_ref_idim` |
| `fitpack_closed_curve_c_interpolating` | `fitpack_closed_curve_c_interpolate` |
| `fitpack_closed_curve_c_smoothing` | `fitpack_fitter_c_ref_smoothing` |
| `fitpack_closed_curve_c_ubegin` | `fitpack_parametric_curve_c_ref_ubegin` |
| `fitpack_closed_curve_c_uend` | `fitpack_parametric_curve_c_ref_uend` |
| `fitpack_constrained_curve_c_degree` | `fitpack_parametric_curve_c_ref_order` |
| `fitpack_constrained_curve_c_derivative` | `fitpack_constrained_curve_c_curve_derivative` |
| `fitpack_constrained_curve_c_idim` | `fitpack_parametric_curve_c_ref_idim` |
| `fitpack_constrained_curve_c_interpolating` | `fitpack_constrained_curve_c_interpolate` |
| `fitpack_constrained_curve_c_smoothing` | `fitpack_fitter_c_ref_smoothing` |
| `fitpack_constrained_curve_c_ubegin` | `fitpack_parametric_curve_c_ref_ubegin` |
| `fitpack_constrained_curve_c_uend` | `fitpack_parametric_curve_c_ref_uend` |
| `fitpack_curve_c_all_derivatives` | `fitpack_curve_c_curve_all_derivatives` |
| `fitpack_curve_c_degree` | `fitpack_curve_c_ref_order` |
| `fitpack_curve_c_derivative` | `fitpack_curve_c_curve_derivative` |
| `fitpack_curve_c_eval_many` | `fitpack_curve_c_curve_eval_many` |
| `fitpack_curve_c_eval_one` | `fitpack_curve_c_curve_eval_one` |
| `fitpack_curve_c_fourier` | `fitpack_curve_c_fourier_coefficients` |
| `fitpack_curve_c_get_bc` | `fitpack_curve_c_ref_bc` |
| `fitpack_curve_c_interpolating` | `fitpack_curve_c_interpolate` |
| `fitpack_curve_c_set_bc` | `fitpack_curve_c_ref_bc` |
| `fitpack_curve_c_smoothing` | `fitpack_fitter_c_ref_smoothing` |
| `fitpack_parametric_curve_c_degree` | `fitpack_parametric_curve_c_ref_order` |
| `fitpack_parametric_curve_c_derivative` | `fitpack_parametric_curve_c_curve_derivative` |
| `fitpack_parametric_curve_c_idim` | `fitpack_parametric_curve_c_ref_idim` |
| `fitpack_parametric_curve_c_interpolating` | `fitpack_parametric_curve_c_interpolate` |
| `fitpack_parametric_curve_c_smoothing` | `fitpack_fitter_c_ref_smoothing` |
| `fitpack_parametric_curve_c_ubegin` | `fitpack_parametric_curve_c_ref_ubegin` |
| `fitpack_parametric_curve_c_ubegin_ref` | `fitpack_parametric_curve_c_ref_ubegin` |
| `fitpack_parametric_curve_c_uend` | `fitpack_parametric_curve_c_ref_uend` |
| `fitpack_parametric_curve_c_uend_ref` | `fitpack_parametric_curve_c_ref_uend` |
| `fitpack_periodic_curve_c_all_derivatives` | `fitpack_periodic_curve_c_curve_all_derivatives` |
| `fitpack_periodic_curve_c_degree` | `fitpack_curve_c_ref_order` |
| `fitpack_periodic_curve_c_derivative` | `fitpack_periodic_curve_c_curve_derivative` |
| `fitpack_periodic_curve_c_eval_many` | `fitpack_periodic_curve_c_curve_eval_many` |
| `fitpack_periodic_curve_c_eval_one` | `fitpack_periodic_curve_c_curve_eval_one` |
| `fitpack_periodic_curve_c_fourier` | `fitpack_periodic_curve_c_fourier_coefficients` |
| `fitpack_periodic_curve_c_interpolating` | `fitpack_periodic_curve_c_interpolate` |
| `fitpack_periodic_curve_c_smoothing` | `fitpack_fitter_c_ref_smoothing` |

`_ref_*` accessors return a pointer to the live component: read through it, or assign through
it to replace `set_bc`. The pointer is invalidated by any refit or destroy.

## Same name, changed signature

```
fitpack_closed_curve_c_allocate
    old: this)
    new: this,status)
fitpack_closed_curve_c_copy
    old: this,that)
    new: this,that,deep_copy,status)
fitpack_closed_curve_c_destroy
    old: this)
    new: this,status)
fitpack_closed_curve_c_eval_one
    old: this,u,y)
    new: this,u,ierr,result,n_result)
fitpack_closed_curve_c_fit
    old: this,smoothing,order)
    new: this,smoothing,order,keep_knots)
fitpack_closed_curve_c_move_alloc
    old: to,from)
    new: to,from,status)
fitpack_closed_curve_c_new_fit
    old: this,ndim,npts,x,u,w,smoothing,order)
    new: this,x_n1,x_n2,x,u,w,smoothing,order)
fitpack_closed_curve_c_new_points
    old: this,ndim,npts,x,y,w)
    new: this,x_n1,x_n2,x,u,w)
fitpack_constrained_curve_c_allocate
    old: this)
    new: this,status)
fitpack_constrained_curve_c_copy
    old: this,that)
    new: this,that,deep_copy,status)
fitpack_constrained_curve_c_destroy
    old: this)
    new: this,status)
fitpack_constrained_curve_c_eval_one
    old: this,u,y)
    new: this,u,ierr,result,n_result)
fitpack_constrained_curve_c_fit
    old: this,smoothing,order)
    new: this,smoothing,order,keep_knots)
fitpack_constrained_curve_c_move_alloc
    old: to,from)
    new: to,from,status)
fitpack_constrained_curve_c_new_fit
    old: this,ndim,npts,x,u,w,smoothing,order)
    new: this,x_n1,x_n2,x,u,w,smoothing,order)
fitpack_constrained_curve_c_new_points
    old: this,ndim,npts,x,y,w)
    new: this,x_n1,x_n2,x,u,w)
fitpack_constrained_curve_c_set_constraints
    old: this,nbegin,nend,ddx_begin,ddx_end)
    new: this,ddx_begin_n1,ddx_begin_n2,ddx_begin,ddx_end_n1,ddx_end_n2,ddx_end,ierr)
fitpack_curve_c_allocate
    old: this)
    new: this,status)
fitpack_curve_c_copy
    old: this,that)
    new: this,that,deep_copy,status)
fitpack_curve_c_destroy
    old: this)
    new: this,status)
fitpack_curve_c_fit
    old: this,smoothing,order)
    new: this,smoothing,order,keep_knots)
fitpack_curve_c_move_alloc
    old: to,from)
    new: to,from,status)
fitpack_curve_c_new_fit
    old: this,npts,x,y,w,smoothing)
    new: this,n,x,y,w,smoothing,order)
fitpack_curve_c_new_points
    old: this,npts,x,y,w)
    new: this,n,x,y,w)
fitpack_parametric_curve_c_allocate
    old: this)
    new: this,status)
fitpack_parametric_curve_c_copy
    old: this,that)
    new: this,that,deep_copy,status)
fitpack_parametric_curve_c_destroy
    old: this)
    new: this,status)
fitpack_parametric_curve_c_eval_one
    old: this,u,y)
    new: this,u,ierr,result,n_result)
fitpack_parametric_curve_c_fit
    old: this,smoothing,order)
    new: this,smoothing,order,keep_knots)
fitpack_parametric_curve_c_move_alloc
    old: to,from)
    new: to,from,status)
fitpack_parametric_curve_c_new_fit
    old: this,ndim,npts,x,u,w,smoothing,order)
    new: this,x_n1,x_n2,x,u,w,smoothing,order)
fitpack_parametric_curve_c_new_points
    old: this,ndim,npts,x,y,w)
    new: this,x_n1,x_n2,x,u,w)
fitpack_periodic_curve_c_allocate
    old: this)
    new: this,status)
fitpack_periodic_curve_c_copy
    old: this,that)
    new: this,that,deep_copy,status)
fitpack_periodic_curve_c_destroy
    old: this)
    new: this,status)
fitpack_periodic_curve_c_fit
    old: this,smoothing)
    new: this,smoothing,order,keep_knots)
fitpack_periodic_curve_c_move_alloc
    old: to,from)
    new: to,from,status)
fitpack_periodic_curve_c_new_fit
    old: this,npts,x,y,w,smoothing)
    new: this,n,x,y,w,smoothing,order)
fitpack_periodic_curve_c_new_points
    old: this,npts,x,y,w)
    new: this,n,x,y,w)
```

`ndim, npts` became `x_n1, x_n2` — same meaning (leading and trailing extent of the Fortran
`x(idim,m)`), generated names. `_eval_one` gained the explicit result-buffer convention
`(..., ierr, result, n_result)` shared by every array-returning method.

## New symbols

422 new entry points: ten newly bound fitter types, the abstract base fitter, per-component
`_getcomp_*` and `_ref_*` accessors on every type, and the `comm_*` serialization surface.

| owning type | new symbols |
|-------------|-------------|
| `fitpack_grid_polar_c` | 38 |
| `fitpack_surface_c` | 38 |
| `fitpack_curve_c` | 36 |
| `fitpack_grid_surface_c` | 34 |
| `fitpack_polar_c` | 34 |
| `fitpack_gridded_spline_c` | 32 |
| `fitpack_parametric_curve_c` | 30 |
| `fitpack_sphere_c` | 30 |
| `fitpack_grid_sphere_c` | 29 |
| `fitpack_parametric_surface_c` | 25 |
| `fitpack_periodic_curve_c` | 23 |
| `fitpack_convex_curve_c` | 22 |
| `fitpack_constrained_curve_c` | 21 |
| `fitpack_closed_curve_c` | 15 |
| `fitpack_fitter_c` | 14 |
| `(other)` | 1 |

## Known gap

`fitpack_parametric_surface%new_fit` and `%new_points` take a rank-3 `z(:,:,:)` argument, which
the generator's raw path (capped at rank 2) skips. The previous interface did not bind parametric
surfaces at all, so nothing regresses; the rest of `fpParametricSurface` is bound.

## C++ source-level changes

The C++ wrappers are regenerated too, and the pre-2.0.0 point type is gone:

- `typedef std::vector<FP_REAL> fpPoint` and `flatten_2d_vector()` are **removed** from
  `fitpack_core_c.h`, which is now pure C. The replacement is `fpPoint<dims>` in
  `include/fpPoint.hpp`: a fixed-size, standard-layout point, so a `std::vector<fpPoint<dims>>`
  is bit-identical to the Fortran `x(dims,m)` it feeds. Every `vector<fpPoint>` overload is now
  a `template <FP_SIZE dims>` overload taking `std::vector<fpPoint<dims>>`.
- Methods whose *return type* carries the dimension need it as an explicit template argument:
  `curve.eval<2>(u, &ierr)`, `curve.ddu_all<2>(u, &ierr)`. Where the dimension appears in an
  argument it is deduced: `curve.new_fit(points, u, w, s)`, `curve.constrain_both(a, b)`.
- `fpParametricCurve::interpolate(FP_SIZE order)` is **retired**: it silently ignored its
  argument. `interpolate()` is unchanged; to interpolate at a chosen order, call the generated
  `interpolate(&order)`. `fpCurve::interpolate(order)`, which worked, is kept.
- `fpCurve::eval(xmin, xmax, npts)` no longer defaults `npts` to 100. With the default, an
  `eval(x, ierr)` call with an integer `ierr` is ambiguous against the scalar `eval` under ISO
  C++ (gcc accepts it, clang does not).
- Per-type headers replace the per-module ones: `fitpack_curve_c.h` and
  `fitpack_periodic_curve_c.h` instead of `fitpack_curves_c.h`, and so on. `fitpack.hpp` still
  pulls in the whole C++ surface.
