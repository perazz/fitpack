/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "fitpack", "index.html", [
    [ "fitpack — Modern Fortran Spline Fitting", "index.html", "index" ],
    [ "Book Reference Index", "book_reference.html", [
      [ "Core Algorithms", "book_reference.html#autotoc_md15", null ],
      [ "Curve Fitting (1D)", "book_reference.html#autotoc_md16", null ],
      [ "Curve Evaluation", "book_reference.html#autotoc_md17", null ],
      [ "Surface Fitting (2D)", "book_reference.html#autotoc_md18", null ],
      [ "Surface Evaluation", "book_reference.html#autotoc_md19", null ],
      [ "Polar and Spherical Domains", "book_reference.html#autotoc_md20", null ]
    ] ],
    [ "B-Spline Foundations", "theory_bsplines.html", [
      [ "Piecewise Polynomials and Knot Vectors", "theory_bsplines.html#autotoc_md21", null ],
      [ "B-Spline Basis Functions", "theory_bsplines.html#autotoc_md22", [
        [ "Degree zero", "theory_bsplines.html#autotoc_md23", null ],
        [ "Recurrence for degree", "theory_bsplines.html#autotoc_md24", null ],
        [ "Key properties", "theory_bsplines.html#autotoc_md25", null ]
      ] ],
      [ "Spline as a Linear Combination", "theory_bsplines.html#autotoc_md26", null ],
      [ "Evaluation: the de Boor Algorithm", "theory_bsplines.html#autotoc_md27", null ],
      [ "Derivatives", "theory_bsplines.html#autotoc_md28", null ],
      [ "Integration", "theory_bsplines.html#autotoc_md29", null ],
      [ "Knot Insertion (Oslo Algorithm)", "theory_bsplines.html#autotoc_md30", null ],
      [ "Tensor Products for Bivariate Splines", "theory_bsplines.html#autotoc_md31", null ],
      [ "Further Reading", "theory_bsplines.html#autotoc_md32", null ]
    ] ],
    [ "Curve Fitting Theory", "theory_curve_fitting.html", [
      [ "The Approximation Problem", "theory_curve_fitting.html#autotoc_md33", null ],
      [ "Least-Squares with Fixed Knots", "theory_curve_fitting.html#autotoc_md34", [
        [ "Normal Equations and Banded Structure", "theory_curve_fitting.html#autotoc_md35", null ],
        [ "QR Factorization via Givens Rotations", "theory_curve_fitting.html#autotoc_md36", null ]
      ] ],
      [ "The Smoothing Problem", "theory_curve_fitting.html#autotoc_md37", null ],
      [ "Adaptive Knot Placement", "theory_curve_fitting.html#autotoc_md38", null ],
      [ "Choosing the Smoothing Factor", "theory_curve_fitting.html#autotoc_md39", null ],
      [ "Periodic Splines", "theory_curve_fitting.html#autotoc_md40", null ],
      [ "Parametric Curves", "theory_curve_fitting.html#autotoc_md41", [
        [ "Chord-Length Parameterization", "theory_curve_fitting.html#autotoc_md42", null ],
        [ "Shared Knot Vector", "theory_curve_fitting.html#autotoc_md43", null ],
        [ "Closed Curves", "theory_curve_fitting.html#autotoc_md44", null ],
        [ "Endpoint Constraints", "theory_curve_fitting.html#autotoc_md45", null ]
      ] ],
      [ "Convexity-Constrained Fitting", "theory_curve_fitting.html#autotoc_md46", [
        [ "B-Spline Coefficient Conditions", "theory_curve_fitting.html#autotoc_md47", null ],
        [ "Per-Interval Convexity Flags", "theory_curve_fitting.html#autotoc_md48", null ],
        [ "Quadratic Programming Formulation", "theory_curve_fitting.html#autotoc_md49", null ]
      ] ],
      [ "Summary", "theory_curve_fitting.html#autotoc_md50", null ]
    ] ],
    [ "Surface Fitting Theory", "theory_surface_fitting.html", [
      [ "Tensor Product Splines", "theory_surface_fitting.html#autotoc_md51", null ],
      [ "Scattered Data Fitting (surfit)", "theory_surface_fitting.html#autotoc_md52", [
        [ "The Observation Matrix", "theory_surface_fitting.html#autotoc_md53", null ],
        [ "Least-Squares Formulation", "theory_surface_fitting.html#autotoc_md54", null ],
        [ "Smoothing Formulation", "theory_surface_fitting.html#autotoc_md55", null ],
        [ "Rank Deficiency", "theory_surface_fitting.html#autotoc_md56", null ],
        [ "Adaptive Knot Placement", "theory_surface_fitting.html#autotoc_md57", null ]
      ] ],
      [ "Gridded Data Fitting (regrid)", "theory_surface_fitting.html#autotoc_md58", [
        [ "Kronecker Product Structure", "theory_surface_fitting.html#autotoc_md59", null ],
        [ "Computational Advantage", "theory_surface_fitting.html#autotoc_md60", null ]
      ] ],
      [ "N-Dimensional Gridded Splines", "theory_surface_fitting.html#nd_gridded_splines", [
        [ "Dimension-by-Dimension Solution", "theory_surface_fitting.html#autotoc_md61", null ],
        [ "Alternating-Direction Knot Placement", "theory_surface_fitting.html#autotoc_md62", null ],
        [ "Canonical Data Layout", "theory_surface_fitting.html#autotoc_md63", null ],
        [ "Smoothing Target and Vocabulary", "theory_surface_fitting.html#autotoc_md64", null ],
        [ "Practical Ceiling", "theory_surface_fitting.html#autotoc_md65", null ]
      ] ],
      [ "The Smoothing Norm", "theory_surface_fitting.html#autotoc_md66", null ],
      [ "Choosing the Smoothing Factor", "theory_surface_fitting.html#autotoc_md67", [
        [ "Practical Guidelines", "theory_surface_fitting.html#autotoc_md68", null ]
      ] ],
      [ "Parametric Surfaces", "theory_surface_fitting.html#autotoc_md69", [
        [ "Periodicity", "theory_surface_fitting.html#autotoc_md70", null ]
      ] ],
      [ "Summary", "theory_surface_fitting.html#autotoc_md71", null ]
    ] ],
    [ "Polar and Spherical Domains", "theory_special_domains.html", [
      [ "Polar Coordinates", "theory_special_domains.html#autotoc_md72", [
        [ "The polar domain", "theory_special_domains.html#autotoc_md73", null ],
        [ "Normalized polar transform", "theory_special_domains.html#autotoc_md74", null ],
        [ "FITPACK routines", "theory_special_domains.html#autotoc_md75", null ]
      ] ],
      [ "Continuity at the Origin", "theory_special_domains.html#autotoc_md76", [
        [ "The problem", "theory_special_domains.html#autotoc_md77", null ],
        [ "Continuity orders", "theory_special_domains.html#autotoc_md78", null ],
        [ "Implementation in FITPACK", "theory_special_domains.html#autotoc_md79", null ]
      ] ],
      [ "Boundary Conditions on the Disc Edge", "theory_special_domains.html#autotoc_md80", null ],
      [ "Spherical Coordinates", "theory_special_domains.html#autotoc_md81", [
        [ "The spherical domain", "theory_special_domains.html#autotoc_md82", null ],
        [ "Periodicity in longitude", "theory_special_domains.html#autotoc_md83", null ],
        [ "FITPACK routines", "theory_special_domains.html#autotoc_md84", null ]
      ] ],
      [ "The Pole Problem", "theory_special_domains.html#autotoc_md85", [
        [ "Analogy with the origin problem", "theory_special_domains.html#autotoc_md86", null ],
        [ "Continuity conditions", "theory_special_domains.html#autotoc_md87", null ],
        [ "Pole boundary conditions for gridded spheres", "theory_special_domains.html#autotoc_md88", null ]
      ] ],
      [ "Practical Considerations", "theory_special_domains.html#autotoc_md89", null ]
    ] ],
    [ "Univariate Curve Fitting with fitpack_curve", "tutorial_curve.html", [
      [ "Mathematical Background", "tutorial_curve.html#autotoc_md90", [
        [ "Choosing S", "tutorial_curve.html#autotoc_md91", null ]
      ] ],
      [ "Basic Example", "tutorial_curve.html#autotoc_md92", [
        [ "Step-by-step", "tutorial_curve.html#autotoc_md93", null ]
      ] ],
      [ "Advanced Features", "tutorial_curve.html#autotoc_md94", [
        [ "Weights", "tutorial_curve.html#autotoc_md95", null ],
        [ "Derivatives", "tutorial_curve.html#autotoc_md96", null ],
        [ "Integration", "tutorial_curve.html#autotoc_md97", null ],
        [ "Zero Finding", "tutorial_curve.html#autotoc_md98", null ],
        [ "Interpolation", "tutorial_curve.html#autotoc_md99", null ],
        [ "Smoothing Sweep", "tutorial_curve.html#autotoc_md100", null ]
      ] ],
      [ "API Summary", "tutorial_curve.html#autotoc_md101", null ],
      [ "Complete Example", "tutorial_curve.html#autotoc_md102", null ],
      [ "Error Handling", "tutorial_curve.html#autotoc_md103", null ]
    ] ],
    [ "Periodic Curve Fitting", "tutorial_periodic_curve.html", [
      [ "Overview", "tutorial_periodic_curve.html#autotoc_md104", null ],
      [ "Mathematical Background", "tutorial_periodic_curve.html#autotoc_md105", null ],
      [ "Basic Example", "tutorial_periodic_curve.html#autotoc_md106", null ],
      [ "Periodicity Verification", "tutorial_periodic_curve.html#autotoc_md107", null ],
      [ "Integration Across the Boundary", "tutorial_periodic_curve.html#autotoc_md108", null ],
      [ "Smoothing and Interpolation", "tutorial_periodic_curve.html#autotoc_md109", null ],
      [ "Complete Example", "tutorial_periodic_curve.html#autotoc_md110", null ]
    ] ],
    [ "Parametric, Closed, and Constrained Curves", "tutorial_parametric_curves.html", [
      [ "Overview", "tutorial_parametric_curves.html#autotoc_md111", null ],
      [ "Parametric Curves", "tutorial_parametric_curves.html#autotoc_md112", [
        [ "Chord-Length Parameterization", "tutorial_parametric_curves.html#autotoc_md113", null ],
        [ "API", "tutorial_parametric_curves.html#autotoc_md114", null ],
        [ "Example: Lissajous Figure", "tutorial_parametric_curves.html#autotoc_md115", null ]
      ] ],
      [ "Closed Curves", "tutorial_parametric_curves.html#autotoc_md116", [
        [ "Example: Ellipse", "tutorial_parametric_curves.html#autotoc_md117", null ]
      ] ],
      [ "Constrained Curves", "tutorial_parametric_curves.html#autotoc_md118", [
        [ "Setting Constraints", "tutorial_parametric_curves.html#autotoc_md119", null ],
        [ "Example: Spiral with Prescribed Tangents", "tutorial_parametric_curves.html#autotoc_md120", null ]
      ] ],
      [ "Complete Example", "tutorial_parametric_curves.html#autotoc_md121", null ],
      [ "Error Handling", "tutorial_parametric_curves.html#autotoc_md122", null ],
      [ "Summary", "tutorial_parametric_curves.html#autotoc_md123", null ]
    ] ],
    [ "Convexity-Constrained Curve Fitting", "tutorial_convex_curve.html", [
      [ "Mathematical Background", "tutorial_convex_curve.html#autotoc_md124", null ],
      [ "Basic Example", "tutorial_convex_curve.html#autotoc_md125", [
        [ "Step 1 — Load Data and Set Convexity Flags", "tutorial_convex_curve.html#autotoc_md126", null ],
        [ "Step 2 — Fit with Constraints", "tutorial_convex_curve.html#autotoc_md127", null ],
        [ "Step 3 — Evaluate and Inspect", "tutorial_convex_curve.html#autotoc_md128", null ]
      ] ],
      [ "Comparing Constrained vs Unconstrained", "tutorial_convex_curve.html#autotoc_md129", null ],
      [ "Mixed Constraints", "tutorial_convex_curve.html#autotoc_md130", null ],
      [ "Constructor Shortcut", "tutorial_convex_curve.html#autotoc_md131", null ],
      [ "API Summary", "tutorial_convex_curve.html#autotoc_md132", null ],
      [ "Complete Example", "tutorial_convex_curve.html#autotoc_md133", null ]
    ] ],
    [ "Scattered Surface Fitting Tutorial", "tutorial_surface.html", [
      [ "Mathematical Background", "tutorial_surface.html#autotoc_md134", null ],
      [ "Basic Example", "tutorial_surface.html#autotoc_md135", null ],
      [ "Features", "tutorial_surface.html#autotoc_md136", [
        [ "Point evaluation", "tutorial_surface.html#autotoc_md137", null ],
        [ "Grid evaluation", "tutorial_surface.html#autotoc_md138", null ],
        [ "Partial derivatives", "tutorial_surface.html#autotoc_md139", null ],
        [ "Integration", "tutorial_surface.html#autotoc_md140", null ],
        [ "Cross-sections", "tutorial_surface.html#autotoc_md141", null ],
        [ "Derivative spline", "tutorial_surface.html#autotoc_md142", null ]
      ] ],
      [ "Smoothing Sweep", "tutorial_surface.html#autotoc_md143", null ],
      [ "Fitting Modes", "tutorial_surface.html#autotoc_md144", null ],
      [ "Error Handling", "tutorial_surface.html#autotoc_md145", null ],
      [ "Complete Example", "tutorial_surface.html#autotoc_md146", null ]
    ] ],
    [ "Grid and Parametric Surface Tutorial", "tutorial_grid_surface.html", [
      [ "fitpack_grid_surface", "tutorial_grid_surface.html#autotoc_md147", [
        [ "Creating a Fit", "tutorial_grid_surface.html#autotoc_md148", null ],
        [ "Evaluation", "tutorial_grid_surface.html#autotoc_md149", null ],
        [ "Partial Derivatives", "tutorial_grid_surface.html#autotoc_md150", null ],
        [ "Integration", "tutorial_grid_surface.html#autotoc_md151", null ],
        [ "Cross-Sections and Derivative Splines", "tutorial_grid_surface.html#autotoc_md152", null ],
        [ "Example: Peaks Function on a Grid", "tutorial_grid_surface.html#autotoc_md153", null ]
      ] ],
      [ "fitpack_parametric_surface", "tutorial_grid_surface.html#autotoc_md154", [
        [ "Multi-Component Surfaces", "tutorial_grid_surface.html#autotoc_md155", null ],
        [ "Creating a Fit", "tutorial_grid_surface.html#autotoc_md156", null ],
        [ "Evaluation", "tutorial_grid_surface.html#autotoc_md157", null ],
        [ "Example: Fitting a Torus", "tutorial_grid_surface.html#autotoc_md158", null ]
      ] ],
      [ "Grid vs Scattered Fitting", "tutorial_grid_surface.html#autotoc_md159", null ],
      [ "Complete Example", "tutorial_grid_surface.html#autotoc_md160", null ]
    ] ],
    [ "Polar Domain Fitting Tutorial", "tutorial_polar.html", [
      [ "Overview", "tutorial_polar.html#autotoc_md161", null ],
      [ "Mathematical Background", "tutorial_polar.html#autotoc_md162", [
        [ "Normalized polar coordinates", "tutorial_polar.html#autotoc_md163", null ],
        [ "Continuity at the origin", "tutorial_polar.html#autotoc_md164", null ]
      ] ],
      [ "Scattered Data: fitpack_polar", "tutorial_polar.html#autotoc_md165", [
        [ "Example: smooth function on a unit disc", "tutorial_polar.html#autotoc_md166", null ]
      ] ],
      [ "Gridded Data: fitpack_grid_polar", "tutorial_polar.html#autotoc_md167", [
        [ "Origin boundary conditions", "tutorial_polar.html#autotoc_md168", null ],
        [ "Example: data on a polar grid", "tutorial_polar.html#autotoc_md169", null ]
      ] ],
      [ "Boundary Function", "tutorial_polar.html#autotoc_md170", null ],
      [ "Complete Example", "tutorial_polar.html#autotoc_md171", null ],
      [ "Error Handling", "tutorial_polar.html#autotoc_md172", null ]
    ] ],
    [ "Spherical Spline Fitting Tutorial", "tutorial_sphere.html", [
      [ "Overview", "tutorial_sphere.html#autotoc_md173", null ],
      [ "Mathematical Background", "tutorial_sphere.html#autotoc_md174", [
        [ "The Pole Problem", "tutorial_sphere.html#autotoc_md175", null ],
        [ "Spherical Harmonics Connection", "tutorial_sphere.html#autotoc_md176", null ]
      ] ],
      [ "Scattered Sphere Data — fitpack_sphere", "tutorial_sphere.html#autotoc_md177", [
        [ "Construction and Fitting", "tutorial_sphere.html#autotoc_md178", null ],
        [ "Evaluation", "tutorial_sphere.html#autotoc_md179", null ],
        [ "Example: Spherical Harmonic Test Function", "tutorial_sphere.html#autotoc_md180", null ]
      ] ],
      [ "Gridded Sphere Data — fitpack_grid_sphere", "tutorial_sphere.html#autotoc_md181", [
        [ "Construction and Fitting", "tutorial_sphere.html#autotoc_md182", null ],
        [ "Pole Boundary Conditions", "tutorial_sphere.html#autotoc_md183", null ],
        [ "Example: Function on a Lat-Lon Grid", "tutorial_sphere.html#autotoc_md184", null ]
      ] ],
      [ "Complete Example", "tutorial_sphere.html#autotoc_md185", null ]
    ] ],
    [ "Modules", "namespaces.html", [
      [ "Modules List", "namespaces.html", "namespaces_dup" ],
      [ "Module Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", "namespacemembers_dup" ],
        [ "Functions/Subroutines", "namespacemembers_func.html", "namespacemembers_func" ],
        [ "Variables", "namespacemembers_vars.html", null ]
      ] ]
    ] ],
    [ "Data Types", "annotated.html", [
      [ "Data Types List", "annotated.html", "annotated_dup" ],
      [ "Data Type Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Data Fields", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions/Subroutines", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions/Subroutines", "globals_func.html", null ],
        [ "Typedefs", "globals_type.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"classfp_fitter.html#a23ebf17c80932f5cff8d8696a33e819c",
"classfp_grid_surface.html#af0cd89d402708d1bc7d6c0cefba76dd3",
"classfp_sphere.html#a223ff9f677d6c020f078da5899a185bf",
"functions_vars.html",
"interfacefitpack__parametric__surfaces__c_1_1f__associated.html#a6b1be8dc60687d92061c796652c2c7d4",
"namespacefitpack__core.html#aa452051613d61831a99b2f9ef8e21bc8",
"namespacefitpack__fitters__c.html#add4846e45fcd69fd4a1185f33e06cea6",
"namespacefitpack__gridded__surfaces__c.html#a9f63a96ded709321f19b23a3623e6014",
"namespacefitpack__polar__domains__c.html#a8f689c85309ed55f37f2dc18ffc9d3b4",
"structfitpack__curves_1_1fitpack__curve.html#a6fb1f3284d1718bb94d4da99f65992d9",
"structfitpack__gridded__sphere_1_1fitpack__grid__sphere.html#aaed9c4373c13b313a8c27c145ef81c06",
"structfitpack__parametric__curves_1_1fitpack__parametric__curve.html#a8aef881a8fcc0621022b76622fae7cea",
"structfitpack__surfaces_1_1fitpack__surface.html#a3f42bb481154f8448c67e3543673d753",
"tutorial_polar.html#autotoc_md168"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';