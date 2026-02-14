# FITPACK Coding Style Guide

Modern Fortran (2008+) spline fitting library. Build with `fpm build`, test with `fpm test`.

## Types & Kinds

All types are C-compatible via `iso_c_binding`:
```fortran
integer, parameter :: FP_REAL = c_double      ! All reals
integer, parameter :: FP_SIZE = c_int32_t     ! Array sizes
integer, parameter :: FP_FLAG = c_int32_t     ! Error flags
integer, parameter :: FP_BOOL = c_bool        ! Logicals
```

## Error Handling

Return `integer(FP_FLAG)` error codes, never throw/stop in library code:
```fortran
integer(FP_FLAG) :: ierr
call routine(..., ierr)
if (.not.FITPACK_SUCCESS(ierr)) error stop FITPACK_MESSAGE(ierr)
```

Key flags: `FITPACK_OK=0`, negative = success variants, positive = errors.

Use `fitpack_error_handling(ierr, ierr_out, whereAt)` pattern:
- `ierr_out` present: return error code to caller
- `ierr_out` absent: halt with message

## Pure Routines

Core computational routines must be `pure`. Provide dual interfaces:
```fortran
function eval(this,x,ierr) result(y)     ! Non-pure, returns ierr
pure function eval_pure(this,x) result(y) ! Pure, returns NaN on error
```

Use `elemental` for scalar operations: `destroy`, `FITPACK_SUCCESS`, simple getters.

## Arrays

- **Public API**: assumed-shape `x(:)`, size relationships via `y(size(x))`
- **Core routines**: explicit-shape `x(m)` with size parameters for performance
- **Derived types**: `allocatable` arrays, always provide `destroy` method

## C Interface

Opaque pointer pattern for types:
```fortran
type, bind(C) :: fitpack_curve_c
    type(c_ptr) :: cptr = c_null_ptr
end type
```

C-callable routines:
```fortran
pure subroutine splev_c(...) bind(C, name='splev_c')
    real(FP_REAL), intent(in), value :: x  ! Scalars by value
    real(FP_REAL), intent(out) :: y(*)     ! Arrays as pointers
```

## Naming

- **snake_case** everywhere
- Types: `fitpack_curve`, `fitpack_surface`
- Procedures: `curve_fit_automatic_knots`, `new_fit`, `destroy`
- Constants: `FITPACK_OK`, `OUTSIDE_EXTRAPOLATE`, `MAX_K`
- Named literals: `zero`, `one`, `half`, `pi` (not magic numbers)

## Module Structure

```
fitpack (public API)
├── fitpack_core (types, constants, core algorithms)
├── fitpack_curves, fitpack_surfaces, ... (domain modules)
└── fitpack_*_c (C bindings)
```

Default `private`, explicit `public` declarations. One type per domain module.

## Type-Bound Procedures

Use generics for overloading:
```fortran
type :: fitpack_curve
contains
    procedure :: eval_one, eval_many
    generic :: eval => eval_one, eval_many
end type
```

## Testing

Tests return `logical` success flag:
```fortran
logical function test_feature(iunit) result(success)
    success = FITPACK_SUCCESS(ierr)
end function
```

## Git

- NEVER mention Claude in any commits or PRs

## Don'ts

- No preprocessor directives
- No `iso_fortran_env` (use `iso_c_binding` kinds)
- No `stop`/`error stop` in library routines (return error flags)
- No magic numbers (use named constants)
