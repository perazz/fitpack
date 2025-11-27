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
    use ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    private

    public :: fitpack_curve
    public :: fitpack_periodic_curve

    integer, parameter :: MAX_K = 5

    !> A public type describing a curve fitter y = c(x)
    type :: fitpack_curve

        !> The data points
        integer(FP_SIZE) :: m = 0
        real(FP_REAL), allocatable :: x(:),y(:)

        !> Spline degree
        integer(FP_SIZE) :: order = 3

        !> Interval boundaries
        real(FP_REAL) :: xleft,xright

        ! Node weights
        real(FP_REAL), allocatable :: sp(:),w(:)

        ! Estimated and actual number of knots and their allocations
        integer(FP_SIZE)               :: nest  = 0
        integer(FP_SIZE)               :: lwrk  = 0
        integer(FP_SIZE), allocatable  :: iwrk(:)
        real(FP_REAL), allocatable     :: wrk(:),wrk_fou(:,:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(FP_REAL) :: smoothing = 1000.0_FP_REAL

        ! Actual curve MSE
        real(FP_REAL) :: fp = zero

        ! Curve extrapolation behavior
        integer(FP_FLAG) :: bc = OUTSIDE_NEAREST_BND

        ! Knots
        integer(FP_SIZE) :: knots = 0
        real(FP_REAL), allocatable :: t(:)  ! Knot location

        ! Spline coefficients [knots-order-1]
        real(FP_REAL), allocatable :: c(:)

        ! Runtime flag
        integer(FP_SIZE) :: iopt = 0

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
           procedure, private :: curve_eval_one_noerr
           procedure, private :: curve_eval_many
           procedure, private :: curve_eval_many_pure
           generic :: eval => curve_eval_one,curve_eval_one_noerr, &
                              curve_eval_many,curve_eval_many_pure

           !> Integrate
           procedure :: integral

           !> Fourier coefficients
           procedure :: fourier_coefficients

           !> Find the zeros of a spline s(x), only if it is cubic
           procedure :: zeros

           !> Evaluate derivative at given coordinates
           procedure, private :: curve_derivative
           procedure, private :: curve_derivatives
           procedure, private :: curve_all_derivatives
           procedure, private :: curve_all_derivatives_pure
           generic   :: dfdx     => curve_derivative,curve_derivatives
           generic   :: dfdx_all => curve_all_derivatives,curve_all_derivatives_pure

           !> Properties: MSE
           procedure, non_overridable :: mse => curve_error

           !> Parallel communication interface (size/pack/expand)
           procedure :: comm_size   => curve_comm_size
           procedure :: comm_pack   => curve_comm_pack
           procedure :: comm_expand => curve_comm_expand

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
        real(FP_REAL), intent(in) :: x(:),y(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ierr0

        ierr0 = this%new_fit(x,y,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new curve fit')

    end function new_from_points

    ! Fit a new curve
    integer(FP_FLAG) function new_fit(this,x,y,w,smoothing,order)
        class(fitpack_curve), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights
        real(FP_REAL), optional, intent(in) :: smoothing
        integer(FP_SIZE), optional, intent(in) :: order

        call this%new_points(x,y,w)

        new_fit = this%fit(smoothing,order)

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
       deallocate(this%wrk_fou,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)
       this%xleft = zero
       this%xright = zero

       this%smoothing = 1000.0_FP_REAL
       this%order     = 3
       this%iopt      = 0
       this%nest      = 0
       this%lwrk      = 0
       this%knots     = 0
       this%fp        = 0.0_FP_REAL
       this%bc        = OUTSIDE_NEAREST_BND

    end subroutine destroy

    subroutine new_points(this,x,y,w)
        class(fitpack_curve), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights

        integer, allocatable :: isort(:)
        integer, parameter   :: SAFE = 10

        associate(m=>this%m,nest=>this%nest,lwrk=>this%lwrk)

        call this%destroy()

        m = size(x)

        ! Ensure x are sorted
        isort = fitpack_argsort(x)
        allocate(this%x(m),source=x(isort))
        allocate(this%y(m),source=y(isort))

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w(isort))
        else
           allocate(this%w(m),source=1.0_FP_REAL)
        endif
        allocate(this%sp(m),source=0.0_FP_REAL)

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
        allocate(this%wrk(lwrk),source=0.0_FP_REAL)
        allocate(this%wrk_fou(nest,2))

        endassociate

    end subroutine new_points

    real(FP_REAL) function curve_eval_one(this,x,ierr) result(y)
        class(fitpack_curve), intent(inout) :: this
        real(FP_REAL),        intent(in)    :: x      ! Evaluation point
        integer(FP_FLAG),     intent(out)   :: ierr   ! Optional error flag

        real(FP_REAL) :: y1(1)

        y1 = curve_eval_many(this,[x],ierr)
        y  = y1(1)

    end function curve_eval_one

    pure real(FP_REAL) function curve_eval_one_noerr(this,x) result(y)
        class(fitpack_curve), intent(in) :: this
        real(FP_REAL),        intent(in) :: x      ! Evaluation point

        real(FP_REAL) :: y1(1)

        y1 = curve_eval_many_pure(this,[x])
        y  = y1(1)

    end function curve_eval_one_noerr

    ! Curve evaluation driver
    function curve_eval_many(this,x,ierr) result(y)
        class(fitpack_curve), intent(inout) :: this
        real(FP_REAL),        intent(in)    :: x(:)   ! Evaluation points
        integer(FP_FLAG),     intent(out)   :: ierr   ! Optional error flag
        real(FP_REAL) :: y(size(x))

        integer(FP_SIZE) :: npts
        integer(FP_FLAG) :: ier

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

    ! Curve evaluation driver
    pure function curve_eval_many_pure(this,x) result(y)
        class(fitpack_curve), intent(in)  :: this
        real(FP_REAL), intent(in) :: x(:)   ! Evaluation points
        real(FP_REAL) :: y(size(x))

        integer(FP_SIZE) :: npts
        integer(FP_FLAG) :: ier

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

        ! For the pure version, return NaNs on error
        if (.not.FITPACK_SUCCESS(ier)) y = ieee_value(0.0_FP_REAL,ieee_quiet_nan)

    end function curve_eval_many_pure

    ! Interpolating curve
    integer(FP_FLAG) function interpolating_curve(this,order) result(ierr)
        class(fitpack_curve), intent(inout) :: this
        integer(FP_SIZE), optional, intent(in) :: order

        ! Set zero smoothing
        this%iopt = IOPT_NEW_SMOOTHING
        
        ! Update order if necessary
        ierr = curve_fit_automatic_knots(this,smoothing=zero,order=order)

    end function interpolating_curve

    ! Curve fitting driver: automatic number of knots
    integer(FP_FLAG) function curve_fit_automatic_knots(this,smoothing,order) result(ierr)
        class(fitpack_curve), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        integer(FP_SIZE), optional, intent(in) :: order

        integer(FP_SIZE) :: loop,nit
        real(FP_REAL) :: smooth_now(3)

        !> Get smoothing trajectory
        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots
        if (this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        !> Set/update order
        if (present(order)) this%order = order

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            select type (curve => this)

               type is (fitpack_periodic_curve)

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

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function curve_fit_automatic_knots

    ! Return fitting MSE
    elemental real(FP_REAL) function curve_error(this)
       class(fitpack_curve), intent(in) :: this
       curve_error = this%fp
    end function curve_error

    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    function curve_derivatives(this, x, order, ierr) result(ddx)
       class(fitpack_curve), intent(inout) :: this
       real(FP_REAL),        intent(in)    :: x(:)   ! Evaluation point (scalar)
       integer,              intent(in)    :: order  ! Derivative order. Default 1
       integer(FP_FLAG), optional, intent(out)   :: ierr  ! Optional error flag
       real(FP_REAL), dimension(size(x))   :: ddx

       integer(FP_SIZE) :: ddx_order,m
       integer(FP_FLAG) :: ierr0

       ! Order 0 = spline value
       ddx_order = max(0,order)

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


    !> Evaluate ALL derivatives of the curve at points x
    !>              (j-1)
    !>      d(j) = s     (x) , j=1,2,...,k1
    !>  of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline representation.
    function curve_all_derivatives(this, x, ierr) result(ddx)
       class(fitpack_curve), intent(inout) :: this
       real(FP_REAL),        intent(in)    :: x   ! Evaluation point (scalar)
       integer(FP_FLAG),     intent(out)   :: ierr  ! Optional error flag
       real(FP_REAL), dimension(this%order+1)  :: ddx

       integer(FP_FLAG) :: ierr0

       ierr0 = FITPACK_OK

       !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m the derivative of
       !  order nu of a spline s(x) of degree k, given in its b-spline representation.
       call spalde(this%t,       & ! Position of the knots
                   this%knots,   & ! Number of knots
                   this%c,       & ! spline coefficients
                   this%order+1, & ! spline order (=degree+1)
                   x,            & ! Point where this should be evaluated
                   ddx,          & ! Evaluated derivatives
                   ierr0)        ! Output flag

       call fitpack_error_handling(ierr0,ierr,'evaluate all derivatives')

    end function curve_all_derivatives

    !> Evaluate ALL derivatives of the curve at points x
    !>              (j-1)
    !>      d(j) = s     (x) , j=1,2,...,k1
    !>  of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline representation.
    function curve_all_derivatives_pure(this, x) result(ddx)
       class(fitpack_curve), intent(in)        :: this
       real(FP_REAL),        intent(in)        :: x   ! Evaluation point (scalar)
       real(FP_REAL), dimension(this%order+1)  :: ddx

       integer(FP_FLAG) :: ierr0

       ierr0 = FITPACK_OK

       !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m the derivative of
       !  order nu of a spline s(x) of degree k, given in its b-spline representation.
       call spalde(this%t,       & ! Position of the knots
                   this%knots,   & ! Number of knots
                   this%c,       & ! spline coefficients
                   this%order+1, & ! spline order (=degree+1)
                   x,            & ! Point where this should be evaluated
                   ddx,          & ! Evaluated derivatives
                   ierr0)        ! Output flag

       if (.not.FITPACK_SUCCESS(ierr0)) ddx = ieee_value(0.0_FP_REAL,ieee_quiet_nan)

    end function curve_all_derivatives_pure

    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    real(FP_REAL) function curve_derivative(this, x, order, ierr) result(ddx)
       class(fitpack_curve), intent(inout) :: this
       real(FP_REAL),        intent(in)    :: x      ! Evaluation point (scalar)
       integer,              intent(in)    :: order  ! Derivative order. Default 1
       integer(FP_FLAG), optional, intent(out)   :: ierr   ! Optional error flag

       real(FP_REAL) :: ddxa(1)

       ! Use the array-based wrapper
       ddxa = curve_derivatives(this,[x],order,ierr)
       ddx  = ddxa(1)

    end function curve_derivative

    ! Calculates the integral of the spline function in interval [from,to]
    real(FP_REAL) function integral(this,from,to)
       class(fitpack_curve), intent(inout) :: this
       real(FP_REAL), intent(in) :: from,to

       integral = splint(this%t, &  ! array of knots
                         this%knots, & ! number of knots
                         this%c    , & ! array of spline coefficients
                         this%order, & ! degree of the spline
                         from,to,    & ! endpoints of the integration interval
                         this%wrk)     ! working space

    end function integral

    !> Fourier coefficients: compute fourier coefficients from the interior of the spline
    !> represenation (must be >=10 knots) in the form:
    !>                 /t(n-3)
    !>    A(i) =      |        s(x)*sin(alfa(i)*x) dx    and
    !>           t(4)/
    !>                 /t(n-3)
    !>    B(i) =      |        s(x)*cos(alfa(i)*x) dx, i=1,...,size(alfa),
    !>           t(4)/
    !> for user defined alpha(:))
    subroutine fourier_coefficients(this,alpha,A,B,ierr)
        class(fitpack_curve), intent(inout) :: this
        real(FP_REAL),    intent(in) :: alpha(:)
        real(FP_REAL),    intent(out), dimension(size(alpha)) :: A,B
        integer(FP_FLAG), intent(out), optional :: ierr

        integer(FP_FLAG) :: ier

        call fourco(this%t,this%knots,  & ! Spline knots)
                    this%c,             & ! spline coefficients
                    alpha,size(alpha),  & ! Parameters alpha(i)
                    A,B,                & ! Computed Fourier coefficients
                    this%wrk_fou(:,1),  & ! Working space
                    this%wrk_fou(:,2),  & ! Working space
                    ier)

        call fitpack_error_handling(ier,ierr,'fourier_coefficients')

    end subroutine fourier_coefficients

    !> Find the zeros of a cubic spline s(x)
    function zeros(this,ierr)
        class(fitpack_curve), intent(in) :: this
        real(FP_REAL), allocatable :: zeros(:)
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_SIZE), parameter :: maxit = 10
        integer(FP_SIZE) :: mest,m,nit,i
        integer(FP_FLAG) :: ier
        real(FP_REAL), allocatable :: tmp_zeros(:)

        ! Estimate a relatively safe number of zeros
        m    = this%knots
        ier  = FITPACK_INSUFFICIENT_STORAGE
        nit  = 0

        ! Iterate until all zeros fit the array
        storage_attempts: do while (nit<maxit .and. ier==FITPACK_INSUFFICIENT_STORAGE)

            nit = nit+1

            ! Allocate storage with the new estimate
            mest = 2*m
            deallocate(tmp_zeros,stat=i)
            allocate(tmp_zeros(mest))

            call sproot(this%t,this%knots,   & ! Spline knots
                        this%c,              & ! the spline coefficients
                        tmp_zeros,mest,      & ! Temporary storage
                        m,                   & ! Actual number of zeros
                        ier)

        end do storage_attempts

        !> Always return an allocated array
        if (m>0) then
            allocate(zeros(m),source=tmp_zeros(1:m))
        else
            allocate(zeros(0))
        endif

        call fitpack_error_handling(ier,ierr,'compute zeros')

    end function zeros

    ! =================================================================================================
    ! PARALLEL COMMUNICATION (size/pack/expand)
    ! =================================================================================================

    !> Return communication buffer size (number of FP_REAL elements)
    !> This counts storage for all curve data needed to reconstruct the spline
    pure integer(FP_SIZE) function curve_comm_size(this)
        class(fitpack_curve), intent(in) :: this

        ! Scalar integers: m, order, knots, bc, iopt, nest, lwrk (7 values)
        ! Scalar reals: xleft, xright, smoothing, fp (4 values)
        ! Use 11 FP_REAL slots for 11 scalar values
        curve_comm_size = 11 &
                        + FP_COMM_SIZE(this%x) &
                        + FP_COMM_SIZE(this%y) &
                        + FP_COMM_SIZE(this%sp) &
                        + FP_COMM_SIZE(this%w) &
                        + FP_COMM_SIZE(this%t) &
                        + FP_COMM_SIZE(this%c) &
                        + FP_COMM_SIZE(this%wrk) &
                        + FP_COMM_SIZE(this%iwrk) &
                        + FP_COMM_SIZE(this%wrk_fou)

    end function curve_comm_size

    !> Pack curve data into communication buffer
    pure subroutine curve_comm_pack(this, buffer)
        class(fitpack_curve), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)

        integer(FP_SIZE) :: pos, n

        ! Pack scalar integers as reals
        buffer(1)  = real(this%m, FP_COMM)
        buffer(2)  = real(this%order, FP_COMM)
        buffer(3)  = real(this%knots, FP_COMM)
        buffer(4)  = real(this%bc, FP_COMM)
        buffer(5)  = real(this%iopt, FP_COMM)
        buffer(6)  = real(this%nest, FP_COMM)
        buffer(7)  = real(this%lwrk, FP_COMM)
        buffer(8)  = this%xleft
        buffer(9)  = this%xright
        buffer(10) = this%smoothing
        buffer(11) = this%fp
        pos = 12

        ! Pack 1D real arrays using FP_COMM_PACK
        call FP_COMM_PACK(this%x, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_PACK(this%y, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_PACK(this%sp, buffer(pos:));  pos = pos + FP_COMM_SIZE(this%sp)
        call FP_COMM_PACK(this%w, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_PACK(this%t, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_PACK(this%c, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%c)
        call FP_COMM_PACK(this%wrk, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk)
        call FP_COMM_PACK(this%iwrk, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%iwrk)
        call FP_COMM_PACK(this%wrk_fou, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk_fou)

    end subroutine curve_comm_pack

    !> Expand curve data from communication buffer
    pure subroutine curve_comm_expand(this, buffer)
        class(fitpack_curve), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)

        integer(FP_SIZE) :: pos

        ! Expand scalar integers
        this%m     = nint(buffer(1), FP_SIZE)
        this%order = nint(buffer(2), FP_SIZE)
        this%knots = nint(buffer(3), FP_SIZE)
        this%bc    = nint(buffer(4), FP_FLAG)
        this%iopt  = nint(buffer(5), FP_SIZE)
        this%nest  = nint(buffer(6), FP_SIZE)
        this%lwrk  = nint(buffer(7), FP_SIZE)
        this%xleft     = buffer(8)
        this%xright    = buffer(9)
        this%smoothing = buffer(10)
        this%fp        = buffer(11)
        pos = 12

        ! Expand arrays using FP_COMM_EXPAND
        ! After expand, FP_COMM_SIZE returns size based on new allocation
        call FP_COMM_EXPAND(this%x, buffer(pos:));       pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_EXPAND(this%y, buffer(pos:));       pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_EXPAND(this%sp, buffer(pos:));      pos = pos + FP_COMM_SIZE(this%sp)
        call FP_COMM_EXPAND(this%w, buffer(pos:));       pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_EXPAND(this%t, buffer(pos:));       pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_EXPAND(this%c, buffer(pos:));       pos = pos + FP_COMM_SIZE(this%c)
        call FP_COMM_EXPAND(this%wrk, buffer(pos:));     pos = pos + FP_COMM_SIZE(this%wrk)
        call FP_COMM_EXPAND(this%iwrk, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%iwrk)
        call FP_COMM_EXPAND(this%wrk_fou, buffer(pos:))

    end subroutine curve_comm_expand

end module fitpack_curves
