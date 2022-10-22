FRESCO-fitpack
===

This is a Modern Fortran translation of the FITPACK package for curve and surface fitting.
The functions are modernized and translated from the original Fortran77 code [FITPACK](http://www.netlib.org/dierckx) by Paul Dierckx.
An object-oriented interface wrapper is also being built. 

Building, using
===============

An automated build is not available yet. 
- `fitpack_core.f90` contains the refactored package
- `fitpack.f90` contains the object-oriented interface
- `fitpack_tests.f90` contains the original test programs, refactored as subroutines. 

A simple command line build script is: 

```
gfortran src/fitpack_core.f90 src/fitpack_tests.f90 src/fitpack.f90 test/test.f90 -o fitpack_test.exe
```

 
References
----------
Fitpack contains very robust algorithms for curve interpolation and fitting, based on algorithms described by Paul Dierckx in Ref [1-4]:<br>

[1] P. Dierckx, "An algorithm for smoothing, differentiation and integration of experimental data using spline functions", J.Comp.Appl.Maths 1 (1975) 165-184.

[2] P. Dierckx, "A fast algorithm for smoothing data on a rectangular grid while using spline functions", SIAM J.Numer.Anal. 19 (1982) 1286-1304.

[3] P. Dierckx, "An improved algorithm for curve fitting with spline functions", report tw54, Dept. Computer Science,K.U. Leuven, 1981.

[4] P. Dierckx, "Curve and surface fitting with splines", Monographs on Numerical Analysis, Oxford University Press, 1993.
