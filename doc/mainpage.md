commit a2cd11b23af0988e0e5813ab0525db6a6ae94a3e
Author: Federico Perini <federico.perini@gmail.com>
Date:   Wed Feb 18 22:16:44 2026 +0100

    Fix surfit y-bounds check, percur zero-pivot guard, and NaN-safe tests
    
    Three bug fixes:
    
    1. surfit: y-values were checked against x-bounds instead of y-bounds,
       causing spurious "invalid input" errors when x/y ranges differed.
    
    2. fp_rotate_row_2mat / fp_rotate_row_2mat_vec: add zero-pivot guards
       matching original Netlib fpperi.f. Without these, fpgivs computes
       0/0 = NaN when both pivot and target are zero.
    
    3. fitpack_curves new_points: use percur workspace formula
       nest*(8+5*k) instead of curfit's nest*(7+3*k).
    
    4. Test comparisons: use `.not.(expr <= tol)` instead of `expr > tol`
       so NaN values are caught (IEEE 754: NaN > x is always false).

diff --git a/doc/mainpage.md b/doc/mainpage.md
index 311b4d7..af03e87 100644
--- a/doc/mainpage.md
+++ b/doc/mainpage.md
@@ -76,15 +76,51 @@ fpm build          # build the library
 fpm test           # run all tests
 ```
 
+## Spline Theory
+
+- @ref theory_bsplines &mdash; B-spline basis functions, de Boor evaluation, derivatives, integration, knot insertion
+- @ref theory_curve_fitting &mdash; Smoothing criterion, adaptive knot placement, periodic and parametric curves
+- @ref theory_surface_fitting &mdash; Tensor product splines, scattered vs. gridded fitting, Kronecker product efficiency
+- @ref theory_special_domains &mdash; Polar coordinate transform, spherical coordinates, pole boundary conditions
+
 ## Tutorials
 
-- @ref tutorial_curves — Smoothing, interpolation, derivatives, roots, periodic and parametric curves
-- @ref tutorial_surfaces — Scattered and gridded surface fitting, cross-sections, integration
-- @ref tutorial_special — Polar and spherical domains, pole boundary conditions
-- @ref book_reference — Quick-reference table mapping routines to book chapters and equations
+### Curves
+
+- @ref tutorial_curve &mdash; Univariate smoothing and interpolation with `fitpack_curve`
+- @ref tutorial_periodic_curve &mdash; Periodic splines with `fitpack_periodic_curve`
+- @ref tutorial_parametric_curves &mdash; Parametric, closed, and constrained curves
+- @ref tutorial_convex_curve &mdash; Shape-preserving fitting with convexity constraints
+
+### Surfaces
+
+- @ref tutorial_surface &mdash; Scattered bivariate data with `fitpack_surface`
+- @ref tutorial_grid_surface &mdash; Gridded data and parametric surfaces
+
+### Special Domains
+
+- @ref tutorial_polar &mdash; Disc-shaped domains with `fitpack_polar` and `fitpack_grid_polar`
+- @ref tutorial_sphere &mdash; Spherical data with `fitpack_sphere` and `fitpack_grid_sphere`
+
+## Examples
+
+Compilable example programs are in the `examples/` directory. Build and run with:
+
+```bash
+fpm run --example example_curve
+fpm run --example example_periodic_curve
+fpm run --example example_parametric_curves
+fpm run --example example_convex_curve
+fpm run --example example_surface
+fpm run --example example_grid_surface
+fpm run --example example_polar
+fpm run --example example_sphere
+```
 
 ## Reference
 
+- @ref book_reference &mdash; Quick-reference table mapping routines to book chapters and equations
+
 The algorithms are described in:
 
 > P. Dierckx, *Curve and Surface Fitting with Splines*,
