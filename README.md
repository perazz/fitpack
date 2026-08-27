fitpack
===

This is a Modern Fortran translation of the FITPACK package for curve and surface fitting.
The functions are modernized and translated from the original Fortran77 code [FITPACK](http://www.netlib.org/dierckx) by Paul Dierckx.
The starting code used the double precision version of FITPACK distributed with [scipy](http://www.scipy.org).

An object-oriented interface wrapper was also built. A C/C++ interface is also being built. 

### 1D Spline interpolators:

Class      | Description | Degree
---        | ---         | ---
`fitpack_curve` | 1D spline interpolation of scattered data, $y = s(x)$ | up to 5
`fitpack_periodic_curve` | 1D spline interpolation of scattered data on a periodic domain, $y = s(x), s(0) = s(x_{per})$ | up to 5
`fitpack_parametric_curve` | Parametric 1D curves in N dimensions, $x_i = s_i(u)$, $i=1,\ldots,n$ | up to 5
`fitpack_closed_curve` | Closed parametric 1D curves in N dimensions, $x_i = s_i(u)$, $i=1,\ldots,n$, $x_i(0)=x_i(1)$ | up to 5
`fitpack_constrained_curve` | Parametric 1D curves in N dimensions with value/derivative constraints at the endpoints $x_i = s_i(u)$, $i=1,\ldots,n$, $x_{i}^{(j)}(0)=u_{L,i}^{(j)}$ , $x_{i}^{(j)}(1)=u_{R,i}^{(j)}$, $0\le j \le 2$| up to 5

Here is an example from the `mncurf` 1D generic curve fitting test:

![General Curve Fitting](media/mncurf.gif)

This is an example from the parametric curve test `mnpara`: 

![Parametric Curve Fitting](media/mnpara.gif)

### 2D Spline interpolators:

Class      | Description | Degree
---        | ---         | ---
`fitpack_surface` | 2D spline interpolation of scattered data, $z = s(x,y)$ | up to 5
`fitpack_polar` | 2D spline interpolation of scattered data in a user-defined polar domain $z = s(u,v)$, $u\in[0,1]$, $v\in[-\pi,\pi]$, user-defined domain radius as a function of polar angle $r=r(v)$ | 3
`fitpack_sphere` | 2D spline interpolation of scattered data on a sphere domain $z = s(u,v)$ with latitude $u \in [0,\pi]$, longitude $v \in [-\pi,\pi]$ | 3
`fitpack_grid_surface` | 2D spline interpolation of rectangular 2D data $z = s(x,y)$ with gridded fitting coordinates $x_i, i=1,\ldots,n_x$,  $y_j, j=1,\ldots,n_y$  | up to 5
`fitpack_grid_polar` | 2D spline interpolation of polar data $z = s(u,v)$ in the fixed-radius circular polar domain $u\in[0,r]$, $v\in[-\pi,\pi]$, with user-control of function and derivatives at the origin and the boundaries | 3

### `C`, `C++` interfaces

The C and C++ header-only interfaces are found in the `include` folder. Every fitter class is bound:

Fortran      | C | C++
---        | ---         | ---
`fitpack_fitter` (abstract) | `fitpack_fitter_c` | `fpFitter`
`fitpack_curve` | `fitpack_curve_c` | `fpCurve`
`fitpack_periodic_curve` | `fitpack_periodic_curve_c` | `fpPeriodicCurve`
`fitpack_parametric_curve` | `fitpack_parametric_curve_c` | `fpParametricCurve`
`fitpack_closed_curve` | `fitpack_closed_curve_c` | `fpClosedCurve`
`fitpack_constrained_curve` | `fitpack_constrained_curve_c` | `fpConstrainedCurve`
`fitpack_convex_curve` | `fitpack_convex_curve_c` | `fpConvexCurve`
`fitpack_surface` | `fitpack_surface_c` | `fpSurface`
`fitpack_grid_surface` | `fitpack_grid_surface_c` | `fpGridSurface`
`fitpack_gridded_spline` | `fitpack_gridded_spline_c` | `fpGridSpline`
`fitpack_parametric_surface` | `fitpack_parametric_surface_c` | `fpParametricSurface`
`fitpack_polar` | `fitpack_polar_c` | `fpPolar`
`fitpack_grid_polar` | `fitpack_grid_polar_c` | `fpGridPolar`
`fitpack_sphere` | `fitpack_sphere_c` | `fpSphere`
`fitpack_grid_sphere` | `fitpack_grid_sphere_c` | `fpGridSphere`

The choice to provide a header-only `C++` implementation is motivated by the need to keep the library C-ABI compatible whatever compiler is being used to build it. For example, on macOS, one may build the library with g++/gfortran, that is not ABI-compatible with clang++. So, it is important that no C++ code is compiled together with the Fortran code in the library. The headers require C++17.

#### Generated bindings

Everything under `include/`, and every `src/capi/*_c.f90` except `src/capi/fitpack_core_c.f90`, is produced by the [fortran-arrays](https://github.com/perazz/fortran-arrays) binding generator from `fitpack_bindings.yaml`. **Never hand-edit a generated file**: change the Fortran source, the config, a banner template under `bindings_templates/`, or an ergonomic snippet under `extra_methods/`, then re-run

```bash
FITPACK_BINDINGS_GENERATOR=/path/to/fortran-arrays/tools/bindings/generate-bindings ./scripts/generate_bindings.sh
```

The hand-written parts of the interface are `src/capi/fitpack_core_c.f90` and `include/fitpack_core_c.h` (the 25 `fp_*_c` core procedures; that header is pure C), `include/fpPoint.hpp`, and the snippets in `extra_methods/`.

Every C entry point reports through `fp_status`, FITPACK's own name for the status struct. Where the headers are compiled with [fortran-arrays](https://github.com/perazz/fortran-arrays) on the include path it is the library's `fx_status`; where they are not, the headers define a layout-identical struct of their own, so the two spellings share one ABI. That same build gate (`HAVE_FXARRAY`, auto-detected in `fitpack_config.h`) adds a zero-copy `fxArray<T>` view of every array component next to the `std::vector` copy; the C ABI is identical either way. The published documentation is built with the gate on, so those view methods are listed.

A method's constness follows its Fortran receiver: a routine declared `intent(in)` binds to a `const` C++ method and a `const <type>_c*` handle, and anything else — `fit`, `new_fit`, `interpolate`, `least_squares`, and the evaluators that update the fitter's internal state — binds to the mutable spelling.

On Windows, define `FITPACK_CAPI_STATIC` before including any fitpack header when linking the static archive: the export macro otherwise defaults to `__declspec(dllimport)`.

#### Points: `fpPoint<dims>`

A point of a parametric curve, or a scattered evaluation site of an N-D gridded spline, is an `fpPoint<dims>` — a fixed-size, standard-layout wrapper over `std::array<FP_REAL,dims>`:

```cpp
std::vector<fpPoint<2>> xy = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}};

fpParametricCurve curve;
FP_FLAG ierr = curve.new_fit(xy, 0.05);      // dims deduced from the argument
fpPoint<2> p = curve.eval<2>(0.5, &ierr);    // dims explicit where it is the return type
```

A `std::vector<fpPoint<dims>>` is bit-identical to the Fortran `x(dims,m)` it feeds, so there is no per-point allocation and no scatter loop. This replaces the pre-2.0.0 `typedef std::vector<FP_REAL> fpPoint`, which cost one heap allocation per point and leaked C++ into the nominally-C `fitpack_core_c.h`. See [doc/abi_changes_2.0.0.md](doc/abi_changes_2.0.0.md) for the full migration list; `./scripts/abi_symbols.sh --diff <ref>` regenerates it.

Building, using
===============

An automated build is available via the Fortran Package Manager. To use FITPACK within your FPM project, add the following to your fpm.toml file:

```
[dependencies]
fitpack = { git="https://github.com/perazz/fitpack.git" }
```

Otherwise, a simple command line build script that builds all modules in the src/ folder is possible. 

Several test programs are available through the Fortran Package manager. To run them, just type
```
fpm test
```
 
References
----------
Fitpack contains very robust algorithms for curve interpolation and fitting, based on algorithms described by Paul Dierckx in Ref [1-4]:<br>

[1] P. Dierckx, "An algorithm for smoothing, differentiation and integration of experimental data using spline functions", J.Comp.Appl.Maths 1 (1975) 165-184.

[2] P. Dierckx, "A fast algorithm for smoothing data on a rectangular grid while using spline functions", SIAM J.Numer.Anal. 19 (1982) 1286-1304.

[3] P. Dierckx, "An improved algorithm for curve fitting with spline functions", report tw54, Dept. Computer Science,K.U. Leuven, 1981.

[4] P. Dierckx, "Curve and surface fitting with splines", Monographs on Numerical Analysis, Oxford University Press, 1993.

[5] P. Dierckx, R. Piessens, "Calculation of Fourier coefficients of discrete functions using cubic splines", Journal of Computational and Applied Mathematics  3(3), 207-209, 1977.

See also
--------
- [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl), a Julia interface to FITPACK
- [scipy.interpolate](https://docs.scipy.org/doc/scipy/reference/interpolate.html), which includes interfaces to FITPACK
