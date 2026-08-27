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
!  fitpack_constrained_curves_c
!> @brief C interface to fitpack_constrained_curves
!
! **************************************************************************************************
!
!> @author Federico Perini <federico.perini@gmail.com>
!> @since     01/07/2024
!
! **************************************************************************************************
module fitpack_constrained_curves_c
     use fitpack_parametric_curves, only: fitpack_constrained_curve
     use fitpack_core, only: FP_SIZE, FP_REAL, FP_FLAG, FP_BOOL, FITPACK_INPUT_ERROR
     use, intrinsic :: iso_c_binding
     implicit none
     private

     !> Public Fortran interface
     public :: fitpack_constrained_curve_c_get_pointer
     public :: fitpack_constrained_curve_c_pointer
     public :: fitpack_constrained_curve_c_allocate
     public :: fitpack_constrained_curve_c_destroy
     public :: fitpack_constrained_curve_c_copy
     public :: fitpack_constrained_curve_c_move_alloc
     public :: fitpack_constrained_curve_c_new_points
     public :: fitpack_constrained_curve_c_set_default_parameters
     public :: fitpack_constrained_curve_c_new_fit
     public :: fitpack_constrained_curve_c_fit
     public :: fitpack_constrained_curve_c_interpolating
     public :: fitpack_constrained_curve_c_eval_one
     public :: fitpack_constrained_curve_c_derivative
     public :: fitpack_constrained_curve_c_smoothing
     public :: fitpack_constrained_curve_c_mse
     public :: fitpack_constrained_curve_c_degree
     public :: fitpack_constrained_curve_c_idim
     public :: fitpack_constrained_curve_c_ubegin
     public :: fitpack_constrained_curve_c_uend

     !> Opaque-pointer C derived type
     type, public, bind(C) :: fitpack_constrained_curve_c
        type(c_ptr) :: cptr = c_null_ptr
     end type fitpack_constrained_curve_c


     !> Null-initialized object
     type(fitpack_constrained_curve_c), parameter, public :: &
        fitpack_constrained_curve_c_null = fitpack_constrained_curve_c(cptr=c_null_ptr)

     contains

     !> Get pointer to embedded Fortran object
     subroutine fitpack_constrained_curve_c_get_pointer(this,fptr)
        type(fitpack_constrained_curve_c), intent(in) :: this
        type(fitpack_constrained_curve), pointer, intent(inout) :: fptr

        if (c_associated(this%cptr)) then
            call c_f_pointer(this%cptr,fptr)
        else
            nullify(fptr)
        end if
     end subroutine fitpack_constrained_curve_c_get_pointer

     !> Get embedded Fortran object; allocate it if not done yet
     subroutine fitpack_constrained_curve_c_pointer(this,fthis)
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve), pointer, intent(inout) :: fthis
        integer :: ierr

        ! Check if fthis is already associated
        call fitpack_constrained_curve_c_get_pointer(this,fthis)

        ! Already associated: do nothing
        if (associated(fthis)) return

        ! Allocate if necessary
        allocate(fthis,stat=ierr)
        if (ierr/=0) stop 'fitpack_constrained_curve_c_pointer allocation error'

        ! Locate new fortran object
        call fthis%destroy()
        this%cptr = c_loc(fthis)

     end subroutine fitpack_constrained_curve_c_pointer

     !> Instantiate Fortran object
     subroutine fitpack_constrained_curve_c_allocate(this) &
                bind(C,name='fitpack_constrained_curve_c_allocate')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        type(fitpack_constrained_curve), pointer :: fthis
        call fitpack_constrained_curve_c_pointer(this,fthis)
     end subroutine fitpack_constrained_curve_c_allocate

     !> Deallocate Fortran object
     subroutine fitpack_constrained_curve_c_destroy(this) &
                bind(C,name='fitpack_constrained_curve_c_destroy')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fthis
        integer :: ierr

        call fitpack_constrained_curve_c_get_pointer(this, fthis)
        if (associated(fthis)) then
            call fthis%destroy()
            deallocate(fthis,stat=ierr)
            if (ierr/=0) stop 'fitpack_constrained_curve_c_destroy deallocation error'
        end if

        ! Nullify C pointers
        this = fitpack_constrained_curve_c_null

     end subroutine fitpack_constrained_curve_c_destroy

     !> Copy contents from another object
     subroutine fitpack_constrained_curve_c_copy(this,that) &
                bind(C,name='fitpack_constrained_curve_c_copy')
         type(fitpack_constrained_curve_c), intent(inout) :: this,that

         type(fitpack_constrained_curve), pointer :: fthis,fthat

         call fitpack_constrained_curve_c_get_pointer(that,fthat)

         if (associated(fthat)) then
             ! Rhs has an object: copy it
             call fitpack_constrained_curve_c_pointer(this,fthis)
             fthis = fthat
         else
            ! Rhs not allocated: return empty object
            call fitpack_constrained_curve_c_destroy(this)
         end if

     end subroutine fitpack_constrained_curve_c_copy

      !> Move allocation from one object to another
      subroutine fitpack_constrained_curve_c_move_alloc(to,from) &
                 bind(C,name='fitpack_constrained_curve_c_move_alloc')
         type(fitpack_constrained_curve_c), intent(inout) :: to,from

         ! Ensure the receiver starts deallocated
         call fitpack_constrained_curve_c_destroy(to)

         ! Whether rhs has an object or not, move ownership to 'to'
         to  %cptr = from%cptr
         from%cptr = c_null_ptr

      end subroutine fitpack_constrained_curve_c_move_alloc

     !> Wrapper to new_points
     subroutine fitpack_constrained_curve_c_new_points(this,ndim,npts,x,y,w) &
                bind(c,name='fitpack_constrained_curve_c_new_points')
         type(fitpack_constrained_curve_c), intent(inout) :: this
         integer(FP_SIZE), intent(in), value :: ndim,npts
         real(FP_REAL), intent(in) :: x(ndim,npts),y(ndim,npts)
         real(FP_REAL), optional, intent(in) :: w(ndim,npts)

         !> Local variables
         type(fitpack_constrained_curve), pointer :: fcurve

         !> Get object; allocate it in case
         call fitpack_constrained_curve_c_pointer(this, fcurve)

         call fcurve%new_points(x,y,w)

     end subroutine fitpack_constrained_curve_c_new_points

     subroutine fitpack_constrained_curve_c_set_default_parameters(this) &
         bind(c,name='fitpack_constrained_curve_c_set_default_parameters')
         type(fitpack_constrained_curve_c), intent(inout) :: this
         type(fitpack_constrained_curve), pointer :: fcurve
         call fitpack_constrained_curve_c_get_pointer(this,fcurve)
         if (associated(fcurve)) call fcurve%set_default_parameters()
     end subroutine fitpack_constrained_curve_c_set_default_parameters

     !> Wrapper to new_fit
     integer(FP_FLAG) function fitpack_constrained_curve_c_new_fit(this,ndim,npts,x,u,w,smoothing,order) result(ierr) &
                                bind(c,name='fitpack_constrained_curve_c_new_fit')
         type(fitpack_constrained_curve_c), intent(inout) :: this
         integer(FP_SIZE), value,    intent(in) :: ndim,npts
         real   (FP_REAL),           intent(in) :: x(ndim,npts)
         real   (FP_REAL), optional, intent(in) :: u(npts)
         real   (FP_REAL), optional, intent(in) :: w(npts)
         real   (FP_REAL), optional, intent(in) :: smoothing
         integer(FP_SIZE), optional, intent(in) :: order

         !> Local variables
         type(fitpack_constrained_curve), pointer :: fcurve

         !> Get object; allocate it in case
         call fitpack_constrained_curve_c_pointer(this, fcurve)

         ierr = fcurve%new_fit(x,u,w,smoothing,order)

     end function fitpack_constrained_curve_c_new_fit

     !> Wrapper to curve_fit_automatic_knots
     integer(FP_FLAG) function fitpack_constrained_curve_c_fit(this,smoothing,order) result(ierr) &
                               bind(c,name='fitpack_constrained_curve_c_fit')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        real   (FP_REAL), optional, intent(in) :: smoothing
        integer(FP_SIZE), optional, intent(in) :: order

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_pointer(this, fcurve)

        ierr = fcurve%fit(smoothing,order)

     end function fitpack_constrained_curve_c_fit

     !> Wrapper to interpolating_curve
     integer(FP_FLAG) function fitpack_constrained_curve_c_interpolating(this) result(ierr) &
                                bind(c,name='fitpack_constrained_curve_c_interpolating')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_pointer(this, fcurve)

        ierr = fcurve%interpolate()

     end function fitpack_constrained_curve_c_interpolating

     !> Wrapper to curve_eval_one
     integer(FP_FLAG) function fitpack_constrained_curve_c_eval_one(this,u,y) result(ierr) &
                            bind(c,name='fitpack_constrained_curve_c_eval_one')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        real(FP_REAL), intent(in), value :: u
        real(FP_REAL), intent(out) :: y(*)

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_pointer(this, fcurve)

        y(1:fcurve%idim) = fcurve%eval(u,ierr)

     end function fitpack_constrained_curve_c_eval_one

     !> Evaluate k-th derivative of the curve at point u
     integer(FP_FLAG) function fitpack_constrained_curve_c_derivative(this,u,order,ddx) result(ierr) &
                            bind(c,name='fitpack_constrained_curve_c_derivative')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        real   (FP_REAL), intent(in), value :: u
        integer(FP_SIZE), intent(in), value :: order
        real   (FP_REAL), intent(out)       :: ddx(*)

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) ddx(1:fcurve%idim) = fcurve%dfdx(u,order,ierr)

     end function fitpack_constrained_curve_c_derivative

     !> Get smoothing
     real(FP_REAL) function fitpack_constrained_curve_c_smoothing(this) &
                            bind(c,name='fitpack_constrained_curve_c_smoothing')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) then
           fitpack_constrained_curve_c_smoothing = fcurve%smoothing
        else
           fitpack_constrained_curve_c_smoothing = -1.0_FP_REAL
        endif

     end function fitpack_constrained_curve_c_smoothing


     !> Get dimensions
     integer(FP_SIZE) function fitpack_constrained_curve_c_idim(this) &
                            bind(c,name='fitpack_constrained_curve_c_idim')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) then
           fitpack_constrained_curve_c_idim = fcurve%idim
        else
           fitpack_constrained_curve_c_idim = 0
        end if

     end function fitpack_constrained_curve_c_idim

     !> Get begin endpoint
     real(FP_REAL) function fitpack_constrained_curve_c_ubegin(this) &
                            bind(c,name='fitpack_constrained_curve_c_ubegin')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) then
           fitpack_constrained_curve_c_ubegin = fcurve%ubegin
        else
           fitpack_constrained_curve_c_ubegin = -huge(0.0_FP_REAL)
        end if

     end function fitpack_constrained_curve_c_ubegin

     !> Get end endpoint
     real(FP_REAL) function fitpack_constrained_curve_c_uend(this) &
                            bind(c,name='fitpack_constrained_curve_c_uend')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) then
           fitpack_constrained_curve_c_uend = fcurve%uend
        else
           fitpack_constrained_curve_c_uend = -huge(0.0_FP_REAL)
        end if

     end function fitpack_constrained_curve_c_uend

     !> Get MSE
     real(FP_REAL) function fitpack_constrained_curve_c_mse(this) &
                            bind(c,name='fitpack_constrained_curve_c_mse')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object; allocate it in case
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) then
           fitpack_constrained_curve_c_mse = fcurve%fp
        else
           fitpack_constrained_curve_c_mse = 0.0_FP_REAL
        endif

     end function fitpack_constrained_curve_c_mse

     !> Get spline degree
     integer(FP_SIZE) function fitpack_constrained_curve_c_degree(this) &
                            bind(c,name='fitpack_constrained_curve_c_degree')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) then
           fitpack_constrained_curve_c_degree = fcurve%order
        else
           fitpack_constrained_curve_c_degree = 0
        end if

     end function fitpack_constrained_curve_c_degree

     !> Reset constraints
     subroutine fitpack_constrained_curve_c_clean_constraints(this) &
                bind(C,name='fitpack_constrained_curve_c_clean_constraints')
        type(fitpack_constrained_curve_c), intent(inout) :: this

        type(fitpack_constrained_curve), pointer :: fcurve

        !> Get object
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)

        if (associated(fcurve)) call fcurve%clean_constraints()

     end subroutine fitpack_constrained_curve_c_clean_constraints

     !> Set constraints _and reset all previous ones_: a missing "ddx_end" means: no endpoint constraints
     !> Begin/end point constraints: (1:idim,0)=function; (1:idim,i)=i-th derivative
     integer(FP_FLAG) function fitpack_constrained_curve_c_set_constraints(this,nbegin,nend,ddx_begin,ddx_end) &
                               result(ierr) bind(C,name='fitpack_constrained_curve_c_set_constraints')
        type(fitpack_constrained_curve_c), intent(inout) :: this
        integer(FP_SIZE), intent(in), value :: nbegin,nend ! Number of constraints on begin, end points
        real(FP_REAL), optional, intent(in) :: ddx_begin(*)
        real(FP_REAL), optional, intent(in) :: ddx_end  (*)

        type(fitpack_constrained_curve), pointer :: fcurve
        integer(FP_SIZE) :: idim
        real(FP_REAL), allocatable, dimension(:,:) :: con_begin,con_end

        ierr = FITPACK_INPUT_ERROR

        !> Get dimensions
        idim      = fitpack_constrained_curve_c_idim(this); if (idim<=0) return

        !> Get object
        call fitpack_constrained_curve_c_get_pointer(this,fcurve)
        if (.not.associated(fcurve)) return

        !> Use allocatable copies for optional arguments
        if (present(ddx_begin)) allocate(con_begin(idim,0:nbegin-1),source=reshape(ddx_begin(:idim*nbegin),[idim,nbegin]))
        if (present(ddx_end  )) allocate(con_end  (idim,0:nend-1  ),source=reshape(ddx_end  (:idim*nend),[idim,nend]))

        call fcurve%set_constraints(con_begin,con_end,ierr)

     end function fitpack_constrained_curve_c_set_constraints





end module fitpack_constrained_curves_c
