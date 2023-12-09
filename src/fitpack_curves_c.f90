! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   (C) Federico Perini, 12/09/2023
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************
module fitpack_curves_c
    use fitpack_core
    use fitpack_curves
    use iso_c_binding
    implicit none
    private

    !> Opaque pointer C type for `fitpack_curve` objects
    type, public, bind(C) :: fp_curve_c
        type(c_ptr) :: cptr = c_null_ptr
    end type fp_curve_c

    public :: fp_curve_c_allocate
    public :: fp_curve_c_destroy
    public :: fp_curve_c_pointer
    public :: fp_curve_c_new_points
    public :: fp_curve_c_new_fit
    public :: fp_curve_c_fit
    public :: fp_curve_c_interpolating
    public :: fp_curve_c_eval_one
    public :: fp_curve_c_eval_many
    public :: fp_curve_c_integral
    public :: fp_curve_c_fourier
    public :: fp_curve_c_derivative
    public :: fp_curve_c_all_derivatives
    public :: fp_curve_c_smoothing
    public :: fp_curve_c_mse
    public :: fp_curve_c_degree

    contains

     !> Helper function: return Fortran pointer
     subroutine fp_curve_c_get_pointer(this,fptr)
        type(fp_curve_c), intent(in) :: this
        type(fitpack_curve), pointer, intent(out) :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr,fptr)
        else
            nullify(fptr)
        end if

     end subroutine fp_curve_c_get_pointer

     !> C wrapper: deallocate fit
     subroutine fp_curve_c_destroy(this) bind(C,name='fp_curve_c_destroy')
         type(fp_curve_c), intent(inout) :: this

         type(fitpack_curve), pointer :: fthis
         integer :: ierr

         call fp_curve_c_get_pointer(this, fthis)

         ! Clear memory from Fortran
         if (associated(fthis)) then
              call fthis%destroy()
              deallocate(fthis,stat=ierr)

              ! Should never happen
              if (ierr/=0) stop 'fp_curve_c_destroy deallocation error'
         end if

         ! Nullify C pointer
         this%cptr = c_null_ptr

     end subroutine fp_curve_c_destroy

     !> C wrapper: initialize fit
     subroutine fp_curve_c_allocate(this) bind(C,name='fp_curve_c_allocate')
         type(fp_curve_c), intent(inout) :: this

         type(fitpack_curve), pointer :: fcurve

         call fp_curve_c_pointer(this,fcurve)

     end subroutine fp_curve_c_allocate

     !> C wrapper: initialize fit
     subroutine fp_curve_c_pointer(this,fcurve)
         type(fp_curve_c), intent(inout) :: this
         type(fitpack_curve), pointer, intent(out) :: fcurve
         integer :: ierr

         ! Check if fcurve is already associated
         call fp_curve_c_get_pointer(this,fcurve)

         ! Already associated: do nothing
         if (associated(fcurve)) return

         ! Allocate if necessary
         allocate(fcurve,stat=ierr)
         if (ierr/=0) stop 'fp_curve_c_allocate allocation error'

         ! Locate new fortran object
         call fcurve%destroy()
         this%cptr = c_loc(fcurve)

     end subroutine fp_curve_c_pointer

     !> Wrapper to new_points
     subroutine fp_curve_c_new_points(this,npts,x,y,w) bind(c,name='fp_curve_c_new_points')
         type(fp_curve_c), intent(inout) :: this
         integer(FP_SIZE), intent(in), value :: npts
         real(FP_REAL), intent(in) :: x(npts),y(npts)
         real(FP_REAL), optional, intent(in) :: w(npts)

         !> Local variables
         type(fitpack_curve), pointer :: fcurve

         !> Get object; allocate it in case
         call fp_curve_c_pointer(this, fcurve)

         call fcurve%new_points(x,y,w)

     end subroutine fp_curve_c_new_points

     !> Wrapper to new_fit
     integer(FP_SIZE) function fp_curve_c_new_fit(this,npts,x,y,w,smoothing) result(ierr) bind(c,name='fp_curve_c_new_fit')
         type(fp_curve_c), intent(inout) :: this
         integer(FP_SIZE), intent(in), value :: npts
         real(FP_REAL), intent(in) :: x(npts),y(npts)
         real(FP_REAL), optional, intent(in) :: w(npts)
         real(FP_REAL), optional, intent(in) :: smoothing

         !> Local variables
         type(fitpack_curve), pointer :: fcurve

         !> Get object; allocate it in case
         call fp_curve_c_pointer(this, fcurve)

         ierr = fcurve%new_fit(x,y,w,smoothing)

     end function fp_curve_c_new_fit

     !> Wrapper to curve_fit_automatic_knots
     integer(FP_SIZE) function fp_curve_c_fit(this,smoothing) result(ierr) bind(c,name='fp_curve_c_fit')
        type(fp_curve_c), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this, fcurve)

        ierr = fcurve%fit(smoothing)

     end function fp_curve_c_fit

     !> Wrapper to interpolating_curve
     integer(FP_SIZE) function fp_curve_c_interpolating(this) result(ierr) bind(c,name='fp_curve_FP_SIZEeprolating')
        type(fp_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this, fcurve)

        ierr = fcurve%interpolate()

     end function fp_curve_c_interpolating

     !> Wrapper to curve_eval_one
     real(FP_REAL) function fp_curve_c_eval_one(this,x,ierr) result(y) bind(c,name='fp_curve_c_eval_one')
        type(fp_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: x
        integer(FP_SIZE), optional, intent(out) :: ierr


        type(fitpack_curve), pointer :: fcurve
        integer :: ierr0

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this, fcurve)

        y = fcurve%eval(x,ierr0)

        if (present(ierr)) ierr = ierr0

     end function fp_curve_c_eval_one

     !> Wrapper to curve_eval_many
     subroutine fp_curve_c_eval_many(this,npts,x,y,ierr) bind(c,name='fp_curve_c_eval_many')
        type(fp_curve_c), intent(inout) :: this
        integer(FP_SIZE), intent(in), value :: npts
        real(FP_REAL), intent(in) :: x(npts)
        real(FP_REAL), intent(out) :: y(npts)
        integer(FP_SIZE), optional, intent(out) :: ierr

        type(fitpack_curve), pointer :: fcurve
        integer :: ierr0

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this, fcurve)

        y = fcurve%eval(x,ierr0)

        if (present(ierr)) ierr = ierr0

     end subroutine fp_curve_c_eval_many

     !> Wrapper to integral
     real(FP_REAL) function fp_curve_c_integral(this,from,to) result(y) bind(c,name='fp_curve_c_integral')
        type(fp_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: from,to

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this, fcurve)

        y = fcurve%integral(from,to)

     end function fp_curve_c_integral

     !> Wrapper to fourier_coefficients
     subroutine fp_curve_c_fourier(this,nparm,alpha,A,B,ierr) bind(c,name='fp_curve_c_fourier')
        type(fp_curve_c), intent(inout) :: this
        integer(FP_SIZE), intent(in), value :: nparm
        real(FP_REAL), intent(in) :: alpha(nparm)
        real(FP_REAL), intent(out) :: a(nparm),b(nparm)
        integer(FP_SIZE), optional, intent(out) :: ierr

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this, fcurve)

        call fcurve%fourier_coefficients(alpha,a,b,ierr)

     end subroutine fp_curve_c_fourier

     !> Wrapper to curve_derivative
     real(FP_REAL) function fp_curve_c_derivative(this,x,order,ierr) result(ddx) bind(c,name='fp_curve_c_derivative')
        type(fp_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: x
        integer(FP_SIZE), intent(in), value :: order
        integer(FP_SIZE), optional, intent(out) :: ierr

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this,fcurve)

        ddx = fcurve%dfdx(x,order,ierr)

     end function fp_curve_c_derivative

     !> Wrapper to curve_all_derivatives
     integer(FP_FLAG) function fp_curve_c_all_derivatives(this,x,ddx) result(ierr) bind(c,name='fp_curve_c_all_derivatives')
        type(fp_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: x
        real(FP_REAL), intent(out) :: ddx(*)

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this,fcurve)

        ddx(1:fcurve%order+1) = fcurve%dfdx_all(x,ierr)

     end function fp_curve_c_all_derivatives

     !> Get smoothing
     real(FP_REAL) function fp_curve_c_smoothing(this) bind(c,name='fp_curve_c_smoothing')
        type(fp_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this,fcurve)

        fp_curve_c_smoothing = fcurve%smoothing

     end function fp_curve_c_smoothing

     !> Get MSE
     real(FP_REAL) function fp_curve_c_mse(this) bind(c,name='fp_curve_c_mse')
        type(fp_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this,fcurve)

        fp_curve_c_mse = fcurve%fp

     end function fp_curve_c_mse

     !> Get spline degree
     integer(FP_SIZE) function fp_curve_c_degree(this) bind(c,name='fp_curve_c_degree')
        type(fp_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fp_curve_c_pointer(this,fcurve)

        fp_curve_c_degree = fcurve%order

     end function fp_curve_c_degree

end module fitpack_curves_c
