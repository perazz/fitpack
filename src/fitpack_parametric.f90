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
module fitpack_parametric_curves
    use fitpack_core
    implicit none
    private

    public :: fitpack_parametric_curve

    integer, parameter :: MAX_K = 5

    !> A public type describing a parametric curve fitter defined by points x(:,i) in the idim-dimensional
    !> space, attached to a set of strictly increasing parameter values u(i)
    type :: fitpack_parametric_curve

        !> Number of points
        integer :: m = 0

        !> Number of dimensions
        integer :: idim = 0

        !> The data points
        real(RKIND), allocatable :: x(:,:) ! [idim x m]

        !> Parameter values: one for each point. They may be optional, in which case, they will
        !> be internally calculated by fitpack
        logical :: has_params = .false.
        real(RKIND), allocatable :: u(:)

        !> Spline degree
        integer :: order = 3

        !> Interval boundaries
        real(RKIND) :: ubegin = zero
        real(RKIND) :: uend = zero

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

        ! Knots
        integer     :: knots = 0
        real(RKIND), allocatable :: t(:)  ! Knot location

        ! Spline coefficients [knots-order-1]
        real(RKIND), allocatable :: c(:)

        ! Runtime flag
        integer :: iopt = IOPT_NEW_SMOOTHING

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
!
!           !> Evaluate derivative at given coordinates
!           procedure, private :: curve_derivative
!           procedure, private :: curve_derivatives
!           generic   :: dfdx => curve_derivative,curve_derivatives
!
           !> Properties: MSE
           procedure, non_overridable :: mse => curve_error

    end type fitpack_parametric_curve


    ! Default constructor
    interface fitpack_parametric_curve
       module procedure new_from_points
    end interface fitpack_parametric_curve

    contains

    ! A default constructor
    type(fitpack_parametric_curve) function new_from_points(x,u,w,ierr) result(this)
        real(RKIND), intent(in) :: x(:,:)
        real(RKIND), optional, intent(in) :: u(size(x,2)) ! parameter values
        real(RKIND), optional, intent(in) :: w(size(x,2)) ! node weights
        integer, optional, intent(out) :: ierr

        integer :: ierr0

   !     ierr0 = this%new_fit(x,y,w)
   stop 'not implemented'

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new curve fit')

    end function new_from_points

    ! Fit a new curve
    integer function new_fit(this,x,u,w,smoothing,order)
        class(fitpack_parametric_curve), intent(inout) :: this
        real(RKIND), intent(in) :: x(:,:)
        real(RKIND), optional, intent(in) :: u(size(x,2)) ! parameter values
        real(RKIND), optional, intent(in) :: w(size(x,2)) ! node weights
        real(RKIND), optional, intent(in) :: smoothing
        integer    , optional, intent(in) :: order

        call this%new_points(x,u,w)

        new_fit = this%fit(smoothing,order)

    end function new_fit

    elemental subroutine destroy(this)
       class(fitpack_parametric_curve), intent(inout) :: this
       integer :: ierr
       this%m    = 0
       this%idim = 0
       this%has_params = .false.
       deallocate(this%x,stat=ierr)
       deallocate(this%u,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%sp,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this% wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)
       this%ubegin = zero
       this%uend = zero

       this%smoothing = 1000.0_RKIND
       this%order     = 3
       this%iopt      = IOPT_NEW_SMOOTHING
       this%nest      = 0
       this%lwrk      = 0
       this%knots     = 0
       this%fp        = 0.0_RKIND

    end subroutine destroy

    subroutine new_points(this,x,u,w)
        class(fitpack_parametric_curve), intent(inout) :: this
        real(RKIND), intent(in) :: x(:,:)
        real(RKIND), optional, intent(in) :: u(size(x,2)) ! parameter values
        real(RKIND), optional, intent(in) :: w(size(x,2)) ! node weights

        integer, allocatable :: isort(:)
        integer, parameter   :: SAFE = 2

        associate(m=>this%m,idim=>this%idim,nest=>this%nest,lwrk=>this%lwrk)

        call this%destroy()

        ! Curve size
        m    = size(x,2)
        idim = size(x,1)

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w)
        else
           allocate(this%w(m),source=one)
        endif
        allocate(this%sp(m),source=zero)

        ! Check if parameters were provided, and if so, ensure x's are sorted
        this%has_params = present(u)
        if (this%has_params) then
            isort = fitpack_argsort(u)
            allocate(this%u,source=u(isort))
            allocate(this%x,source=x(:,isort))
            this%w(:) = this%w(isort)

            this%ubegin  = minval(u)
            this%uend    = maxval(u)


        else
            allocate(this%u(m),source=zero)
            this%x = x

            this%ubegin  = zero
            this%uend    = one

        end if

        ! Reset run flag
        this%iopt = IOPT_NEW_SMOOTHING

        ! Reset estimated knots.
        ! In most practical situation nest=m/2 will be sufficient.
        ! always large enough is nest=m+k+1, the number of knots needed for interpolation (s=0).
        nest = SAFE*(m+MAX_ORDER+1)
        allocate(this%iwrk(nest),this%t(nest),this%c(nest*idim))

        ! Setup working space.
        lwrk = (m*(MAX_K+1)+nest*(6+idim+3*MAX_K))
        allocate(this%wrk(lwrk),source=zero)

        endassociate

    end subroutine new_points

    function curve_eval_one(this,u,ierr) result(y)
        class(fitpack_parametric_curve), intent(inout)  :: this
        real(RKIND),          intent(in)     :: u      ! Evaluation point
        integer, optional,    intent(out)    :: ierr   ! Optional error flag
        real(RKIND) :: y(this%idim)

        real(RKIND) :: y1(this%idim,1)

        y1 = curve_eval_many(this,[u],ierr)
        y  = y1(:,1)

    end function curve_eval_one

    ! Curve evaluation driver
    function curve_eval_many(this,u,ierr) result(x)
        class(fitpack_parametric_curve), intent(inout)  :: this
        real(RKIND),          intent(in)     :: u(:)   ! Evaluation points (parameter value)
        integer, optional,    intent(out)    :: ierr   ! Optional error flag
        real(RKIND) :: x(this%idim,size(u))

        integer :: npts,ier

        npts = size(u)

        !  subroutine curev evaluates in a number of points u(i),i=1,2,...,m a spline curve s(u)
        !  of degree k and dimension idim, given in its b-spline representation.
        !
        !  input parameters:
        !    idim : integer, giving the dimension of the spline curve.
        !    t    : array,length n, which contains the position of the knots.
        !    n    : integer, giving the total number of knots of s(u).
        !    c    : array,length nc, which contains the b-spline coefficients.
        !    nc   : integer, giving the total number of coefficients of s(u).
        !    k    : integer, giving the degree of s(u).
        !    u    : array,length m, which contains the points where s(u) must be evaluated.
        !    m    : integer, giving the number of points where s(u) must be evaluated.
        !    mx   : integer, giving the dimension of the array x. mx >= m*idim

        call curev(idim=this%idim,&                ! dimension of the spline curve
                   t=this%t,&                      ! the position of the knots
                   n=this%knots,&                  ! total number of knots of s(u(x))
                   c=this%c,&                      ! the b-spline coefficients
                   nc=size(this%c),&               ! number of coefficients of s(u(x))
                   k=this%order,&                  ! the degree of s(u(x))
                   u=u,m=npts,&                    ! the points where s(u) must be evaluated.
                   x=x,mx=size(x), &               ! the predictions
                   ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate 1d spline')

    end function curve_eval_many

    ! Interpolating curve
    integer function interpolating_curve(this) result(ierr)
        class(fitpack_parametric_curve), intent(inout) :: this

        ! Set zero smoothing
        ierr = curve_fit_automatic_knots(this,zero)

    end function interpolating_curve

    ! Curve fitting driver: automatic number of knots
    integer function curve_fit_automatic_knots(this,smoothing,order) result(ierr)
        class(fitpack_parametric_curve), intent(inout) :: this
        real(RKIND), optional, intent(in) :: smoothing
        integer    , optional, intent(in) :: order

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

        ! Set order
        if (present(order)) this%order = order

        ! First iteration lets solver decide knots
        this%iopt = 0

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! If this is the first call

            ! Call fitting function
            call parcur(this%iopt,                    &  ! option
                        merge(1,0,this%has_params),   &  ! have input parameter values?
                        this%idim,this%m,             &  ! Number of dimensions
                        this%u,                       &  ! Parameter array
                        size(this%x),                 &  ! Unrolled size of x array
                        this%x,this%w,                &  ! points and weights
                        this%ubegin,this%uend,        &  ! Parameter range
                        this%order,this%smoothing,    &  ! spline accuracy
                        this%nest,this%knots,this%t,  &  ! spline output
                        size(this%c),this%c,this%fp,  &  ! spline output
                        this%wrk,this%lwrk,this%iwrk, &  ! memory
                        ierr)                            ! Error flag

            ! After any successful call, parameters have surely been computed.
            if (FITPACK_SUCCESS(ierr)) this%has_params = .true.

        end do

    end function curve_fit_automatic_knots

    ! Return fitting MSE
    elemental real(RKIND) function curve_error(this)
       class(fitpack_parametric_curve), intent(in) :: this
       curve_error = this%fp
    end function curve_error
!
!    !> Evaluate k-th derivative of the curve at points x
!    !> Use 1st derivative if order not present
!    function curve_derivatives(this, x, order, ierr) result(ddx)
!       class(fitpack_parametric_curve), intent(inout) :: this
!       real(RKIND),          intent(in)    :: x(:)   ! Evaluation point (scalar)
!       integer, optional,    intent(in)    :: order  ! Derivative order. Default 1
!       integer, optional,    intent(out)   :: ierr  ! Optional error flag
!       real(RKIND), dimension(size(x))     :: ddx
!
!       integer :: ddx_order,m,ierr0
!
!       if (present(order)) then
!          ddx_order = max(0,order)
!       else
!          ddx_order = 1
!       end if
!
!       ierr0 = FITPACK_OK
!
!
!       m = size(x); if (m<=0) goto 1
!
!       !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m the derivative of
!       !  order nu of a spline s(x) of degree k, given in its b-spline representation.
!       call splder(this%t,     & ! Position of the knots
!                   this%knots, & ! Number of knots
!                   this%c,     & ! spline coefficients
!                   this%order, & ! spline degree
!                   ddx_order,  & ! derivative order (0<=order<=this%order)  x,y,m,e,wrk,ier)
!                   x,          & ! Array of points where this should be evaluated
!                   ddx,        & ! Evaluated derivatives
!                   m,          & ! Number of input points
!                   this%bc,    & ! Extrapolation behavior
!                   this%wrk,   & ! Temporary working space
!                   ierr0)        ! Output flag
!
!       1 call fitpack_error_handling(ierr0,ierr,'evaluate derivative')
!
!    end function curve_derivatives
!
!    !> Evaluate k-th derivative of the curve at points x
!    !> Use 1st derivative if order not present
!    real(RKIND) function curve_derivative(this, x, order, ierr) result(ddx)
!       class(fitpack_parametric_curve), intent(inout)  :: this
!       real(RKIND),          intent(in)     :: x      ! Evaluation point (scalar)
!       integer, optional,    intent(in)     :: order  ! Derivative order. Default 1
!       integer, optional,    intent(out)    :: ierr   ! Optional error flag
!
!       real(RKIND) :: ddxa(1)
!
!       ! Use the array-based wrapper
!       ddxa = curve_derivatives(this,[x],order,ierr)
!       ddx  = ddxa(1)
!
!    end function curve_derivative


end module fitpack_parametric_curves
