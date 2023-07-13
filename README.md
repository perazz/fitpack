fortran-fitpack
===

This is a Modern Fortran translation of the FITPACK package for curve and surface fitting.
The functions are modernized and translated from the original Fortran77 code [FITPACK](http://www.netlib.org/dierckx) by Paul Dierckx.
The starting code used the double precision version of FITPACK distributed with [scipy](http://www.scipy.org).
An object-oriented interface wrapper is also being built. 

Building, using
===============

An automated build is available via the Fortran Package Manager.

Otherwise, a simple command line build script is: 
```
gfortran src/fitpack_core.f90 test/fitpack_test_data.f90 test/fitpack_tests.f90 src/fitpack.f90 test/test.f90 -o fitpack_test.exe
```

- `src/fitpack_core.f90` contains the refactored package
- `src/fitpack.f90` contains the object-oriented interface
- `test/fitpack_tests.f90` contains the original test programs, refactored as subroutines. 
- `test/fitpack_test_data.f90` contains support data to the original test programs
- `test/test.f90` is the test driver program 

A simple makefile for the GNU compiler suite is provided in folder `project`; to run it: 

```
cd project/
make -f makefile.gcc
```

The testing executable is just a wrapper to the 29 original test programs. 
These test programs aren't double checked yet, and their output is supposed to be the original F77 program output
 
References
----------
Fitpack contains very robust algorithms for curve interpolation and fitting, based on algorithms described by Paul Dierckx in Ref [1-4]:<br>

[1] P. Dierckx, "An algorithm for smoothing, differentiation and integration of experimental data using spline functions", J.Comp.Appl.Maths 1 (1975) 165-184.

[2] P. Dierckx, "A fast algorithm for smoothing data on a rectangular grid while using spline functions", SIAM J.Numer.Anal. 19 (1982) 1286-1304.

[3] P. Dierckx, "An improved algorithm for curve fitting with spline functions", report tw54, Dept. Computer Science,K.U. Leuven, 1981.

[4] P. Dierckx, "Curve and surface fitting with splines", Monographs on Numerical Analysis, Oxford University Press, 1993.
