! **************************************************************************************************
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
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************
module fitpack_curves
    use fitpack_core
    implicit none
    private

    public :: fitpack_curve
    public :: fitpack_periodic_curve

    integer, parameter :: MAX_K = 5

    !> A public type describing a curve fitter y = c(x)
    type :: fitpack_curve

        !> The data points
        integer :: m = 0
        real(RKIND), allocatable :: x(:),y(:)

        !> Spline degree
        integer :: order = 3

        !> Interval boundaries
        real(RKIND) :: xleft,xright

        ! Node weights
        real(RKIND), allocatable :: sp(:),w(:)

        ! Estimated and actual number of knots and their allocations
        integer                  :: nest  = 0
        integer                  :: lwrk  = 0
        integer, allocatable     :: iwrk(:)
        real(RKIND), allocatable :: wrk(:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(RKIND) :: smoothing = 1000.0_RKIND

        ! Actual curve MSE
        real(RKIND) :: fp = zero

        ! Curve extrapolation behavior
        integer     :: bc = OUTSIDE_NEAREST_BND

        ! Knots
        integer     :: knots = 0
        real(RKIND), allocatable :: t(:)  ! Knot location

        ! Spline coefficients [knots-order-1]
        real(RKIND), allocatable :: c(:)

        ! Runtime flag
        integer :: iopt = 0

        contains

           !> Clean memory
           procedure :: destroy

           !> Set new points
           procedure :: new_points

           !> Generate new fit
           procedure :: new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit         => curve_fit_automatic_knots
           procedure :: interpolate => interpolating_curve

           !> Evaluate curve at given coordinates
           procedure, private :: curve_eval_one
           procedure, private :: curve_eval_many
           generic :: eval => curve_eval_one,curve_eval_many

           !> Evaluate derivative at given coordinates
           procedure, private :: curve_derivative
           procedure, private :: curve_derivatives
           generic   :: dfdx => curve_derivative,curve_derivatives

           !> Properties: MSE
           procedure, non_overridable :: mse => curve_error

    end type fitpack_curve

    !> Derived type describing a periodic curve. No changes are made to the storage,
    !> but the appropriate package functions will be called depending on the type
    type, extends(fitpack_curve) :: fitpack_periodic_curve

    end type fitpack_periodic_curve

    ! Default constructor
    interface fitpack_curve
       module procedure new_from_points
    end interface fitpack_curve

    contains

    ! A default constructor
    type(fitpack_curve) function new_from_points(x,y,w,ierr) result(this)
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(x,y,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new curve fit')

    end function new_from_points

    ! Fit a new curve
    integer function new_fit(this,x,y,w,smoothing)
        class(fitpack_curve), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights
        real(RKIND), optional, intent(in) :: smoothing

        call this%new_points(x,y,w)

        new_fit = this%fit(smoothing)

    end function new_fit

    elemental subroutine destroy(this)
       class(fitpack_curve), intent(inout) :: this
       integer :: ierr
       this%m = 0
       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%sp,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this% wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)
       this%xleft = zero
       this%xright = zero

       this%smoothing = 1000.0_RKIND
       this%order     = 3
       this%iopt      = 0
       this%nest      = 0
       this%lwrk      = 0
       this%knots     = 0
       this%fp        = 0.0_RKIND
       this%bc        = OUTSIDE_NEAREST_BND

    end subroutine destroy

    subroutine new_points(this,x,y,w)
        class(fitpack_curve), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights

        integer, allocatable :: isort(:)
        integer, parameter   :: SAFE = 10

        associate(m=>this%m,nest=>this%nest,lwrk=>this%lwrk)

        call this%destroy()

        m = size(x)

        ! Ensure x are sorted
        isort = fitpack_argsort(x)
        allocate(this%x,source=x(isort))
        allocate(this%y,source=y(isort))

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w(isort))
        else
           allocate(this%w(m),source=1.0_RKIND)
        endif
        allocate(this%sp(m),source=0.0_RKIND)

        ! Setup boundaries
        this%xleft  = minval(x)
        this%xright = maxval(x)

        ! Reset run flag
        this%iopt = 0

        ! Reset estimated knots
        nest = max(SAFE*2*MAX_K+2,ceiling(1.4*this%m))
        allocate(this%iwrk(nest),this%t(nest),this%c(nest))

        ! Setup working space.
        lwrk = (m*(MAX_K+1)+nest*(7+3*MAX_K))
        allocate(this%wrk(lwrk),source=0.0_RKIND)

        endassociate

    end subroutine new_points

    real(RKIND) function curve_eval_one(this,x,ierr) result(y)
        class(fitpack_curve), intent(inout)  :: this
        real(RKIND),          intent(in)     :: x      ! Evaluation point
        integer, optional,    intent(out)    :: ierr   ! Optional error flag

        real(RKIND) :: y1(1)

        y1 = curve_eval_many(this,[x],ierr)
        y  = y1(1)

    end function curve_eval_one


    ! Curve evaluation driver
    function curve_eval_many(this,x,ierr) result(y)
        class(fitpack_curve), intent(inout)  :: this
        real(RKIND),          intent(in)     :: x(:)   ! Evaluation points
        integer, optional,    intent(out)    :: ierr   ! Optional error flag
        real(RKIND) :: y(size(x))

        integer :: npts,ier

        npts = size(x)

        !  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
        !  a spline s(x) of degree k, given in its b-spline representation.
        !
        !  calling sequence:
        !     call splev(t,n,c,k,x,y,m,e,ier)
        !
        !  input parameters:
        !    t    : array,length n, which contains the position of the knots.
        !    n    : integer, giving the total number of knots of s(x).
        !    c    : array,length n, which contains the b-spline coefficients.
        !    k    : integer, giving the degree of s(x).
        !    x    : array,length m, which contains the points where s(x) must
        !           be evaluated.
        !    m    : integer, giving the number of points where s(x) must be
        !           evaluated.
        !    e    : integer, if 0 the spline is extrapolated from the end
        !           spans for points not in the support, if 1 the spline
        !           evaluates to zero for those points, if 2 ier is set to
        !           1 and the subroutine returns, and if 3 the spline evaluates
        !           to the value of the nearest boundary point.

        call splev(t=this%t,&                      ! the position of the knots
                   n=this%knots,&                  ! total number of knots of s(x)
                   c=this%c,&                      ! the b-spline coefficients
                   k=this%order,&                  ! the degree of s(x)
                   x=x,m=npts,&                    ! the points where s(x) must be evaluated.
                   y=y, &                          ! the predictions
                   e=this%bc,           &          ! What to do outside mapped knot range
                   ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate 1d spline')

    end function curve_eval_many

    ! Interpolating curve
    integer function interpolating_curve(this) result(ierr)
        class(fitpack_curve), intent(inout) :: this

        ! Set zero smoothing
        ierr = curve_fit_automatic_knots(this,zero)

    end function interpolating_curve

    ! Curve fitting driver: automatic number of knots
    integer function curve_fit_automatic_knots(this,smoothing) result(ierr)
        class(fitpack_curve), intent(inout) :: this
        real(RKIND), optional, intent(in) :: smoothing

        integer :: loop,nit

        real(RKIND), parameter :: smoothing_trajectory(*) = [1000.d0,60.d0,30.d0]
        real(RKIND), dimension(size(smoothing_trajectory)) :: smooth_now

        if (present(smoothing)) then
            smooth_now = smoothing
            nit        = 1
        else
            smooth_now = smoothing_trajectory
            nit        = size(smoothing_trajectory)
        end if

        ! First iteration lets solver decide knots
        this%iopt = 0

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            select type (curve => this)

               class is (fitpack_periodic_curve)

                  ! Call fitting function
                  call percur(curve%iopt,                      &  ! option
                              curve%m,curve%x,curve%y,curve%w, &  ! points
                              curve%order,curve%smoothing,     &  ! spline accuracy
                              curve%nest,curve%knots,curve%t,  &  ! spline output
                              curve%c,curve%fp,                &  ! spline output
                              curve%wrk,curve%lwrk,curve%iwrk, &  ! memory
                              ierr)                           ! Error flag

               class default

                  ! Call curvfit
                  call curfit(curve%iopt,                      &  ! option
                              curve%m,curve%x,curve%y,curve%w, &  ! points
                              curve%xleft,curve%xright,        &  ! x range
                              curve%order,curve%smoothing,     &  ! spline accuracy
                              curve%nest,curve%knots,curve%t,  &  ! spline output
                              curve%c,curve%fp,                &  ! spline output
                              curve%wrk,curve%lwrk,curve%iwrk, &  ! memory
                              ierr)                           ! Error flag

            end select

        end do

    end function curve_fit_automatic_knots

    ! Return fitting MSE
    elemental real(RKIND) function curve_error(this)
       class(fitpack_curve), intent(in) :: this
       curve_error = this%fp
    end function curve_error

    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    function curve_derivatives(this, x, order, ierr) result(ddx)
       class(fitpack_curve), intent(inout) :: this
       real(RKIND),          intent(in)    :: x(:)   ! Evaluation point (scalar)
       integer, optional,    intent(in)    :: order  ! Derivative order. Default 1
       integer, optional,    intent(out)   :: ierr  ! Optional error flag
       real(RKIND), dimension(size(x))     :: ddx

       integer :: ddx_order,m,ierr0

       if (present(order)) then
          ddx_order = max(0,order)
       else
          ddx_order = 1
       end if

       ierr0 = FITPACK_OK


       m = size(x); if (m<=0) goto 1

       !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m the derivative of
       !  order nu of a spline s(x) of degree k, given in its b-spline representation.
       call splder(this%t,     & ! Position of the knots
                   this%knots, & ! Number of knots
                   this%c,     & ! spline coefficients
                   this%order, & ! spline degree
                   ddx_order,  & ! derivative order (0<=order<=this%order)  x,y,m,e,wrk,ier)
                   x,          & ! Array of points where this should be evaluated
                   ddx,        & ! Evaluated derivatives
                   m,          & ! Number of input points
                   this%bc,    & ! Extrapolation behavior
                   this%wrk,   & ! Temporary working space
                   ierr0)        ! Output flag

       1 call fitpack_error_handling(ierr0,ierr,'evaluate derivative')

    end function curve_derivatives

    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    real(RKIND) function curve_derivative(this, x, order, ierr) result(ddx)
       class(fitpack_curve), intent(inout)  :: this
       real(RKIND),          intent(in)     :: x      ! Evaluation point (scalar)
       integer, optional,    intent(in)     :: order  ! Derivative order. Default 1
       integer, optional,    intent(out)    :: ierr   ! Optional error flag

       real(RKIND) :: ddxa(1)

       ! Use the array-based wrapper
       ddxa = curve_derivatives(this,[x],order,ierr)
       ddx  = ddxa(1)

    end function curve_derivative


end module fitpack_curves
