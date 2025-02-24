! **************************************************************************************************
!
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   Refactored by Federico Perini, 10/6/2022
!   Based on the netlib library by Paul Dierckx
!
! **************************************************************************************************
!
!  fitpack_curves_c
!> @brief C interface to fitpack_curves
!
! **************************************************************************************************
!
!> @author Federico Perini <federico.perini@gmail.com>
!> @since     12/09/2023
!
! **************************************************************************************************
module fitpack_curves_c
     use fitpack_curves, only: fitpack_curve
     use fitpack_core, only: FP_SIZE, FP_REAL, FP_FLAG, FP_BOOL
     use, intrinsic :: iso_c_binding
     implicit none
     private

     !> Public Fortran interface
     public :: fitpack_curve_c_get_pointer
     public :: fitpack_curve_c_pointer
     public :: fitpack_curve_c_allocate
     public :: fitpack_curve_c_destroy
     public :: fitpack_curve_c_copy
     public :: fitpack_curve_c_new_points
     public :: fitpack_curve_c_new_fit
     public :: fitpack_curve_c_fit
     public :: fitpack_curve_c_interpolating
     public :: fitpack_curve_c_eval_one
     public :: fitpack_curve_c_eval_many
     public :: fitpack_curve_c_integral
     public :: fitpack_curve_c_fourier
     public :: fitpack_curve_c_derivative
     public :: fitpack_curve_c_all_derivatives
     public :: fitpack_curve_c_smoothing
     public :: fitpack_curve_c_mse
     public :: fitpack_curve_c_degree

     !> Opaque-pointer C derived type
     type, public, bind(C) :: fitpack_curve_c
        type(c_ptr) :: cptr = c_null_ptr
     end type fitpack_curve_c


     !> Null-initialized object
     type(fitpack_curve_c), parameter, public :: fitpack_curve_c_null = fitpack_curve_c(cptr=c_null_ptr)

     contains

     !> Get pointer to embedded Fortran object
     subroutine fitpack_curve_c_get_pointer(this,fptr)
        type(fitpack_curve_c), intent(in) :: this
        type(fitpack_curve), pointer, intent(inout) :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr,fptr)
        else
            nullify(fptr)
        end if
     end subroutine fitpack_curve_c_get_pointer

     !> Get embedded Fortran object; allocate it if not done yet
     subroutine fitpack_curve_c_pointer(this,fthis)
        type(fitpack_curve_c), intent(inout) :: this
        type(fitpack_curve), pointer, intent(inout) :: fthis
        integer :: ierr

        ! Check if fthis is already associated
        call fitpack_curve_c_get_pointer(this,fthis)

        ! Already associated: do nothing
        if (associated(fthis)) return

        ! Allocate if necessary
        allocate(fthis,stat=ierr)
        if (ierr/=0) stop 'fitpack_curve_c_pointer allocation error'

        ! Locate new fortran object
        call fthis%destroy()
        this%cptr = c_loc(fthis)

     end subroutine fitpack_curve_c_pointer

     !> Instantiate Fortran object
     subroutine fitpack_curve_c_allocate(this) bind(C,name='fitpack_curve_c_allocate')
        type(fitpack_curve_c), intent(inout) :: this
        type(fitpack_curve), pointer :: fthis
        call fitpack_curve_c_pointer(this,fthis)
     end subroutine fitpack_curve_c_allocate

     !> Deallocate Fortran object
     subroutine fitpack_curve_c_destroy(this) bind(C,name='fitpack_curve_c_destroy')
        type(fitpack_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fthis
        integer :: ierr

        call fitpack_curve_c_get_pointer(this, fthis)
        if (associated(fthis)) then
            call fthis%destroy()
            deallocate(fthis,stat=ierr)
            if (ierr/=0) stop 'fitpack_curve_c_destroy deallocation error'
        end if

        ! Nullify C pointers
        this = fitpack_curve_c_null

     end subroutine fitpack_curve_c_destroy

     !> Copy contents from another object
     subroutine fitpack_curve_c_copy(this,that) bind(C,name='fitpack_curve_c_copy')
         type(fitpack_curve_c), intent(inout) :: this,that

         type(fitpack_curve), pointer :: fthis,fthat

         call fitpack_curve_c_get_pointer(that,fthat)

         if (associated(fthat)) then
             ! Rhs has an object: copy it
             call fitpack_curve_c_pointer(this,fthis)
             fthis = fthat
         else
            ! Rhs not allocated: return empty object
            call fitpack_curve_c_destroy(this)
         end if

     end subroutine fitpack_curve_c_copy

      !> Move allocation from one object to another
      subroutine fitpack_curve_c_move_alloc(to,from) bind(C,name='fitpack_curve_c_move_alloc')
         type(fitpack_curve_c), intent(inout) :: to,from

         ! Ensure the receiver starts deallocated
         call fitpack_curve_c_destroy(to)

         ! Whether rhs has an object or not, move ownership to 'to'
         to  %cptr = from%cptr
         from%cptr = c_null_ptr

      end subroutine fitpack_curve_c_move_alloc

     !> Wrapper to new_points
     subroutine fitpack_curve_c_new_points(this,npts,x,y,w) bind(c,name='fitpack_curve_c_new_points')
         type(fitpack_curve_c), intent(inout) :: this
         integer(FP_SIZE), intent(in), value :: npts
         real(FP_REAL), intent(in) :: x(npts),y(npts)
         real(FP_REAL), optional, intent(in) :: w(npts)

         !> Local variables
         type(fitpack_curve), pointer :: fcurve

         !> Get object; allocate it in case
         call fitpack_curve_c_pointer(this, fcurve)

         call fcurve%new_points(x,y,w)

     end subroutine fitpack_curve_c_new_points

     !> Wrapper to new_fit
     integer(FP_FLAG) function fitpack_curve_c_new_fit(this,npts,x,y,w,smoothing) result(ierr) &
                                bind(c,name='fitpack_curve_c_new_fit')
         type(fitpack_curve_c), intent(inout) :: this
         integer(FP_SIZE), intent(in), value :: npts
         real(FP_REAL), intent(in) :: x(npts),y(npts)
         real(FP_REAL), optional, intent(in) :: w(npts)
         real(FP_REAL), optional, intent(in) :: smoothing

         !> Local variables
         type(fitpack_curve), pointer :: fcurve

         !> Get object; allocate it in case
         call fitpack_curve_c_pointer(this, fcurve)

         ierr = fcurve%new_fit(x,y,w,smoothing)

     end function fitpack_curve_c_new_fit

     !> Wrapper to curve_fit_automatic_knots
     integer(FP_FLAG) function fitpack_curve_c_fit(this,smoothing,order) result(ierr) &
                               bind(c,name='fitpack_curve_c_fit')
        type(fitpack_curve_c), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        integer(FP_SIZE), optional, intent(in) :: order

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this, fcurve)

        ierr = fcurve%fit(smoothing,order)

     end function fitpack_curve_c_fit

     !> Wrapper to interpolating_curve
     integer(FP_FLAG) function fitpack_curve_c_interpolating(this) result(ierr) bind(c,name='fitpack_curve_c_interpolating')
        type(fitpack_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this, fcurve)

        ierr = fcurve%interpolate()

     end function fitpack_curve_c_interpolating

     !> Wrapper to curve_eval_one
     real(FP_REAL) function fitpack_curve_c_eval_one(this,x,ierr) result(y) bind(c,name='fitpack_curve_c_eval_one')
        type(fitpack_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: x
        integer(FP_FLAG), optional, intent(out) :: ierr

        type(fitpack_curve), pointer :: fcurve
        integer :: ierr0

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this, fcurve)

        y = fcurve%eval(x,ierr0)

        if (present(ierr)) ierr = ierr0

     end function fitpack_curve_c_eval_one

     !> Wrapper to curve_eval_many
     subroutine fitpack_curve_c_eval_many(this,npts,x,y,ierr) bind(c,name='fitpack_curve_c_eval_many')
        type(fitpack_curve_c), intent(inout) :: this
        integer(FP_SIZE), intent(in), value :: npts
        real(FP_REAL), intent(in) :: x(npts)
        real(FP_REAL), intent(out) :: y(npts)
        integer(FP_FLAG), optional, intent(out) :: ierr

        type(fitpack_curve), pointer :: fcurve
        integer(FP_FLAG) :: ierr0

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this, fcurve)

        y = fcurve%eval(x,ierr0)

        if (present(ierr)) ierr = ierr0

     end subroutine fitpack_curve_c_eval_many

     !> Wrapper to integral
     real(FP_REAL) function fitpack_curve_c_integral(this,from,to) result(y) bind(c,name='fitpack_curve_c_integral')
        type(fitpack_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: from,to

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this, fcurve)

        y = fcurve%integral(from,to)

     end function fitpack_curve_c_integral

     !> Wrapper to fourier_coefficients
     subroutine fitpack_curve_c_fourier(this,nparm,alpha,A,B,ierr) bind(c,name='fitpack_curve_c_fourier')
        type(fitpack_curve_c), intent(inout) :: this
        integer(FP_SIZE), intent(in), value :: nparm
        real(FP_REAL), intent(in) :: alpha(nparm)
        real(FP_REAL), intent(out) :: a(nparm),b(nparm)
        integer(FP_FLAG), optional, intent(out) :: ierr

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this, fcurve)

        call fcurve%fourier_coefficients(alpha,a,b,ierr)

     end subroutine fitpack_curve_c_fourier

     !> Wrapper to curve_derivative
     real(FP_REAL) function fitpack_curve_c_derivative(this,x,order,ierr) result(ddx) bind(c,name='fitpack_curve_c_derivative')
        type(fitpack_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: x
        integer(FP_SIZE), intent(in), value :: order
        integer(FP_SIZE), optional, intent(out) :: ierr

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this,fcurve)

        ddx = fcurve%dfdx(x,order,ierr)

     end function fitpack_curve_c_derivative

     !> Wrapper to curve_all_derivatives
     integer(FP_FLAG) function fitpack_curve_c_all_derivatives(this,x,ddx) result(ierr) &
                               bind(c,name='fitpack_curve_c_all_derivatives')
        type(fitpack_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: x
        real(FP_REAL), intent(out) :: ddx(*)

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this,fcurve)

        ddx(1:fcurve%order+1) = fcurve%dfdx_all(x,ierr)

     end function fitpack_curve_c_all_derivatives

     !> Get smoothing
     real(FP_REAL) function fitpack_curve_c_smoothing(this) bind(c,name='fitpack_curve_c_smoothing')
        type(fitpack_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this,fcurve)

        fitpack_curve_c_smoothing = fcurve%smoothing

     end function fitpack_curve_c_smoothing

     !> Get MSE
     real(FP_REAL) function fitpack_curve_c_mse(this) bind(c,name='fitpack_curve_c_mse')
        type(fitpack_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this,fcurve)

        fitpack_curve_c_mse = fcurve%fp

     end function fitpack_curve_c_mse

     !> Get spline degree
     integer(FP_SIZE) function fitpack_curve_c_degree(this) bind(c,name='fitpack_curve_c_degree')
        type(fitpack_curve_c), intent(inout) :: this

        type(fitpack_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_curve_c_pointer(this,fcurve)

        fitpack_curve_c_degree = fcurve%order

     end function fitpack_curve_c_degree


end module fitpack_curves_c
