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
    public :: fitpack_closed_curve
    public :: fitpack_constrained_curve

    integer, parameter :: MAX_K = 5

    !> A public type describing a parametric curve fitter defined by points x(:,i) in the idim-dimensional
    !> space, attached to a set of strictly increasing parameter values u(i)
    type :: fitpack_parametric_curve

        !> Number of points
        integer :: m = 0

        !> Number of dimensions
        integer :: idim = 0

        !> The data points
        real(FP_REAL), allocatable :: x(:,:) ! [idim x m]

        !> Parameter values: one for each point. They may be optional, in which case, they will
        !> be internally calculated by fitpack
        logical :: has_params = .false.
        real(FP_REAL), allocatable :: u(:)

        !> Spline degree
        integer :: order = 3

        !> Interval boundaries
        real(FP_REAL) :: ubegin = zero
        real(FP_REAL) :: uend = zero

        ! Node weights
        real(FP_REAL), allocatable :: sp(:),w(:)

        ! Estimated and actual number of knots and their allocations
        integer                  :: nest  = 0
        integer                  :: lwrk  = 0
        integer, allocatable     :: iwrk(:)
        real(FP_REAL), allocatable :: wrk(:)

        ! Space for derivative evaluation
        real(FP_REAL), allocatable :: dd(:,:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(FP_REAL) :: smoothing = 1000.0_FP_REAL

        ! Actual curve MSE
        real(FP_REAL) :: fp = zero

        ! Knots
        integer     :: knots = 0
        real(FP_REAL), allocatable :: t(:)  ! Knot location

        ! Spline coefficients [knots-order-1]
        real(FP_REAL), allocatable :: c(:)

        ! Runtime flag
        integer :: iopt = IOPT_NEW_SMOOTHING

        contains

           !> Clean memory
           procedure :: destroy

           !> Set new points
           procedure :: new_points
           procedure :: set_default_parameters

           !> Generate new fit
           procedure :: new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit         => curve_fit_automatic_knots
           procedure :: interpolate => interpolating_curve

           !> Evaluate curve at given coordinates
           procedure :: curve_eval_one
           procedure :: curve_eval_many
           generic :: eval => curve_eval_one,curve_eval_many

           !> Evaluate derivative at given coordinates
           procedure :: curve_derivative
           procedure :: curve_derivatives
           procedure :: curve_all_derivatives
           generic   :: dfdx => curve_derivative,curve_derivatives
           generic   :: dfdx_all => curve_all_derivatives

           !> Properties: MSE
           procedure, non_overridable :: mse => curve_error

    end type fitpack_parametric_curve

    !> Derived type describing a closed parametric curve. No changes are made to the storage,
    !> but the appropriate package functions will be called depending on the type
    type, extends(fitpack_parametric_curve) :: fitpack_closed_curve

    end type fitpack_closed_curve

    !> Derived type describing a parametric curve with constraints at the boundaries.
    !> The dimensional splines will satisfy the following boundary constraints
    !>                     (l)
    !>       if ib >= 0 :  sj   (u(1)) = db(1:idim,0:ib-1), ib = boundary constraint on the (i-1)-th derivative
    !>   and                (l)
    !>       if ie >= 0 :  sj   (u(m)) = de(1:idim,0:ie-1), ie = boundary constraint on the (i-1)-th derivative
    type, extends(fitpack_parametric_curve) :: fitpack_constrained_curve

        !> Left boundary derivatives
        integer                    :: ib = 0 ! Number of boundary constraints: up to (ib-1)-th derivative
        real(FP_REAL), allocatable :: deriv_begin(:,:) ! 0:ib-1
        !> Right boundary derivatives
        integer                    :: ie = 0 ! Number of boundary constraints: up to (ib-1)-th derivative
        real(FP_REAL), allocatable :: deriv_end(:,:)

        !> On exit xx contains the coordinates of the data points to which a spline curve with zero
        !> derivative constraints has been determined. if the computation mode iopt =1 is used xx should
        !> be left unchanged between calls.
        real(FP_REAL), allocatable :: xx(:,:)

        !> On exit cp contains the b-spline coefficients of a polynomial curve which satisfies the
        !> boundary constraints. if the computation mode iopt =1 is used cp should be left unchanged
        !> between calls.
        real(FP_REAL), allocatable :: cp(:,:)

        contains

           !> Clear memory
           procedure :: destroy => con_destroy

           !> Curve constraints
           procedure :: clean_constraints
           procedure ::   set_constraints

    end type fitpack_constrained_curve

    ! Default constructor
    interface fitpack_parametric_curve
       module procedure new_from_points
    end interface fitpack_parametric_curve

    contains

    ! A default constructor
    type(fitpack_parametric_curve) function new_from_points(x,u,w,ierr) result(this)
        real(FP_REAL), intent(in) :: x(:,:)
        real(FP_REAL), optional, intent(in) :: u(size(x,2)) ! parameter values
        real(FP_REAL), optional, intent(in) :: w(size(x,2)) ! node weights
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(x,u,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new curve fit')

    end function new_from_points

    ! Fit a new curve
    integer function new_fit(this,x,u,w,smoothing,order)
        class(fitpack_parametric_curve), intent(inout) :: this
        real   (FP_REAL), intent(in) :: x(:,:)
        real   (FP_REAL), optional, intent(in) :: u(size(x,2)) ! parameter values
        real   (FP_REAL), optional, intent(in) :: w(size(x,2)) ! node weights
        real   (FP_REAL), optional, intent(in) :: smoothing
        integer(FP_SIZE), optional, intent(in) :: order

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
       deallocate(this%dd,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)
       this%ubegin = zero
       this%uend = zero

       this%smoothing = 1000.0_FP_REAL
       this%order     = 3
       this%iopt      = IOPT_NEW_SMOOTHING
       this%nest      = 0
       this%lwrk      = 0
       this%knots     = 0
       this%fp        = 0.0_FP_REAL

    end subroutine destroy

    elemental subroutine con_destroy(this)
       class(fitpack_constrained_curve), intent(inout) :: this

       integer :: ierr
       call destroy(this)
       call clean_constraints(this)
       deallocate(this%xx,stat=ierr)
       deallocate(this%cp,stat=ierr)

    end subroutine con_destroy

    subroutine new_points(this,x,u,w)
        class(fitpack_parametric_curve), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:,:)
        real(FP_REAL), optional, intent(in) :: u(size(x,2))  ! parameter values
        real(FP_REAL), optional, intent(in) :: w(size(x,2))  ! node weights

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

            call set_default_parameters(this)

        end if

        ! Reset run flag
        this%iopt = IOPT_NEW_SMOOTHING

        ! Reset estimated knots.
        ! In most practical situation nest=m/2 will be sufficient.
        ! always large enough is nest=m+k+1, the number of knots needed for interpolation (s=0).
        select type (curv => this)

           type is (fitpack_constrained_curve)

              ! on exit cp will contain the b-spline coefficients of a polynomial curve which satisfies
              ! the boundary constraints. if the computation mode iopt =1 is used cp should be left
              ! unchanged between calls.
              allocate(curv%xx(idim,m),&
                       curv%cp(idim,2*(MAX_ORDER+1)))

              nest = SAFE*(m+MAX_ORDER+1+max(0,MAX_ORDER)+max(0,MAX_ORDER))
              lwrk = (m*(MAX_ORDER+1)+nest*(6+idim+3*MAX_ORDER))

              ! Ensure constraint arrays are allocated and no constraints set
              call set_constraints(curv)

           type is (fitpack_closed_curve)

              nest = SAFE*(m+2*MAX_ORDER)
              lwrk = (m*(MAX_K+1)+nest*(7+idim+5*MAX_K))

           class default

              nest = SAFE*(m+MAX_ORDER+1)
              lwrk = (m*(MAX_K+1)+nest*(6+idim+3*MAX_K))

        end select

        ! Setup working space.
        allocate(this%iwrk(nest),this%t(nest),this%c(nest*idim))
        allocate(this%wrk(lwrk),source=zero)

        ! Setup space for derivative evaluatiuon
        allocate(this%dd(idim,0:MAX_ORDER),source=zero)

        endassociate

    end subroutine new_points

    elemental subroutine clean_constraints(this)
       class(fitpack_constrained_curve), intent(inout) :: this
       integer :: ierr
       deallocate(this%deriv_begin,stat=ierr)
       deallocate(this%deriv_end,stat=ierr)
       this%ib = 0
       this%ie = 0
    end subroutine clean_constraints

    ! A call to set_constraints will RESET ALL contraints: a missing "ddx_end" means: no constraints
    ! at the endpoint
    subroutine set_constraints(this,ddx_begin,ddx_end,ierr)
        class(fitpack_constrained_curve), intent(inout) :: this

        !> Begin point constraints: (:,0)=function; (:,i)=i-th derivative
        real(FP_REAL), optional, intent(in) :: ddx_begin(:,0:)
        real(FP_REAL), optional, intent(in) :: ddx_end  (:,0:)

        integer, optional, intent(out) :: ierr

        integer :: ier

        ier = FITPACK_OK

        ! Prepare deallocated space
        call clean_constraints(this)

        if (present(ddx_begin)) then
            if (size(ddx_begin,1)/=this%idim .or. .not.size(ddx_begin,2)>0) then
                ier = FITPACK_INVALID_CONSTRAINT
                call fitpack_error_handling(ier,ierr,'constrained_curve: begin point constraint')
                return
            else
                this%ib = size(ddx_begin,2)
                allocate(this%deriv_begin(this%idim,0:ubound(ddx_begin,2)),source=ddx_begin)
            end if
        else
            ! Do not leave array unallocated
            this%ib = 0
            allocate(this%deriv_begin(this%idim,1))
        end if

        if (present(ddx_end)) then
            if (size(ddx_end,1)/=this%idim .or. .not.size(ddx_end,2)>0) then
                ier = FITPACK_INVALID_CONSTRAINT
                call fitpack_error_handling(ier,ierr,'constrained_curve: endpoint constraint')
                return
            else
                this%ie = size(ddx_end,2)
                allocate(this%deriv_end(this%idim,0:ubound(ddx_end,2)),source=ddx_end)
            end if
        else
            ! Do not leave array unallocated
            this%ie = 0
            allocate(this%deriv_end(this%idim,1))
        end if

        ! Add point constraints
        call fitpack_error_handling(ier,ierr,'constrained_curve: set_constraints')

    end subroutine set_constraints

    function curve_eval_one(this,u,ierr) result(y)
        class(fitpack_parametric_curve), intent(inout)  :: this
        real(FP_REAL),          intent(in)     :: u      ! Evaluation point
        integer, optional,    intent(out)    :: ierr   ! Optional error flag
        real(FP_REAL) :: y(this%idim)

        real(FP_REAL) :: y1(this%idim,1)

        y1 = curve_eval_many(this,[u],ierr)
        y  = y1(:,1)

    end function curve_eval_one

    ! Curve evaluation driver
    function curve_eval_many(this,u,ierr) result(x)
        class(fitpack_parametric_curve), intent(inout)  :: this
        real(FP_REAL),          intent(in)     :: u(:)   ! Evaluation points (parameter value)
        integer, optional,    intent(out)    :: ierr   ! Optional error flag
        real(FP_REAL) :: x(this%idim,size(u))

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
        real(FP_REAL), optional, intent(in) :: smoothing
        integer(FP_SIZE), optional, intent(in) :: order

        integer :: loop,nit
        real(FP_REAL) :: smooth_now(3)

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots
        if (this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        ! Set order
        if (present(order)) this%order = order

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! Call fitting function
            select type (curv => this)

               type is (fitpack_closed_curve)

                  call clocur(curv%iopt,                    &  ! option
                              merge(1,0,curv%has_params),   &  ! have input parameter values?
                              curv%idim,curv%m,             &  ! Number of dimensions
                              curv%u,                       &  ! Parameter array
                              size(curv%x),                 &  ! Unrolled size of x array
                              curv%x,curv%w,                &  ! points and weights
                              curv%order,curv%smoothing,    &  ! spline accuracy
                              curv%nest,curv%knots,curv%t,  &  ! spline output
                              size(curv%c),curv%c,curv%fp,  &  ! spline output
                              curv%wrk,curv%lwrk,curv%iwrk, &  ! memory
                              ierr)

               type is (fitpack_parametric_curve)

                  call parcur(curv%iopt,                    &  ! option
                              merge(1,0,curv%has_params),   &  ! have input parameter values?
                              curv%idim,curv%m,             &  ! Number of dimensions
                              curv%u,                       &  ! Parameter array
                              size(curv%x),                 &  ! Unrolled size of x array
                              curv%x,curv%w,                &  ! points and weights
                              curv%ubegin,curv%uend,        &  ! Parameter range
                              curv%order,curv%smoothing,    &  ! spline accuracy
                              curv%nest,curv%knots,curv%t,  &  ! spline output
                              size(curv%c),curv%c,curv%fp,  &  ! spline output
                              curv%wrk,curv%lwrk,curv%iwrk, &  ! memory
                              ierr)

               type is (fitpack_constrained_curve)

                  call concur(curv%iopt,                    &  ! option
                              curv%idim,curv%m,             &  ! Number of dimensions
                              curv%u,                       &  ! Parameter array
                              size(curv%x),                 &  ! Unrolled size of x array
                              curv%x,curv%xx,curv%w,        &  ! points and weights
                              curv%ib,curv%deriv_begin,size(curv%deriv_begin), & ! BEGIN point constraints
                              curv%ie,curv%deriv_end  ,size(curv%deriv_end),   & ! ENDpoint constraints
                              curv%order,curv%smoothing,    &  ! spline accuracy
                              curv%nest,curv%knots,curv%t,  &  ! spline output
                              size(curv%c),curv%c,          &  ! spline output
                              size(curv%cp),curv%cp,        &  ! spline output
                              curv%fp,                      &  ! spline output
                              curv%wrk,curv%lwrk,curv%iwrk, &  ! memory
                              ierr)

            end select

            ! After any successful call, parameters have surely been computed.
            if (FITPACK_SUCCESS(ierr)) then
                this%has_params = .true.
                this%iopt       = IOPT_OLD_FIT
            endif

        end do

    end function curve_fit_automatic_knots

    ! Return fitting MSE
    elemental real(FP_REAL) function curve_error(this)
       class(fitpack_parametric_curve), intent(in) :: this
       curve_error = this%fp
    end function curve_error

    !> Evaluate k-th derivative of the curve at point u
    function curve_derivative(this, u, order, ierr) result(ddx)
       class(fitpack_parametric_curve), intent(inout) :: this
       real(FP_REAL),          intent(in)    :: u      ! Evaluation points (parameter)
       integer,              intent(in)    :: order  ! Derivative order. 0=function; 1:k=i-th derivative
       integer, optional,    intent(out)   :: ierr   ! Optional error flag
       real(FP_REAL), dimension(this%idim)   :: ddx

       integer :: ddx_order,ierr0

       ! Choose order
       ddx_order = max(0,order)

       ierr0 = FITPACK_OK

       !  subroutine cualde evaluates at the point u all the derivatives
       !                (l)
       !     d(j,l) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim
       !  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim.
       call cualde(this%idim,    & ! Number of dimensions
                   this%t,       & ! Position of the knots
                   this%knots,   & ! Number of knots
                   this%c,       & ! Spline coefficients
                   size(this%c), & ! Number of coefficients
                   this%order+1, & ! k1 = order of s(u) (order = degree+1)
                   u,            & ! Where the derivatives must be evaluated
                   this%dd,      & ! Space for derivative evaluation
                   size(this%dd),& ! Its size
                   ierr0)          ! Output flag

       ! Derivative order is 0:k <- extract derivative
       ddx = this%dd(:,ddx_order)

       call fitpack_error_handling(ierr0,ierr,'evaluate derivative')

    end function curve_derivative

    !> Evaluate all derivatives (0:k) of the curve at point u
    function curve_all_derivatives(this, u, ierr) result(ddx)
       class(fitpack_parametric_curve), intent(inout) :: this
       real(FP_REAL),          intent(in)    :: u      ! Evaluation points (parameter)
       integer, optional,    intent(out)   :: ierr   ! Optional error flag
       real(FP_REAL), dimension(this%idim,0:this%order) :: ddx

       integer :: ierr0

       ierr0 = FITPACK_OK

       !  subroutine cualde evaluates at the point u all the derivatives
       !                (l)
       !     d(j,l) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim
       !  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim.
       call cualde(this%idim,    & ! Number of dimensions
                   this%t,       & ! Position of the knots
                   this%knots,   & ! Number of knots
                   this%c,       & ! Spline coefficients
                   size(this%c), & ! Number of coefficients
                   this%order+1, & ! k1 = order of s(u) (order = degree+1)
                   u,            & ! Where the derivatives must be evaluated
                   this%dd,      & ! Space for derivative evaluation
                   size(this%dd),& ! Its size
                   ierr0)          ! Output flag

       ! Derivative order is 0:k <- extract derivative
       ddx = this%dd(:,0:this%order)

       call fitpack_error_handling(ierr0,ierr,'evaluate all derivatives')

    end function curve_all_derivatives

    ! Set normalized coordinates in [0,1] when not provided by the user
    subroutine set_default_parameters(this)
        class(fitpack_parametric_curve), intent(inout) :: this

        integer :: i,m

        associate(u=>this%u,x=>this%x)

        ! Number of points
        m = size(this%x,dim=2)

        ! Point coordinates are stored in x(:), offset by idim values
        u(1) = zero
        do i=2,m
           u(i) = u(i-1) + norm2(x(:,i)-x(:,i-1))
        end do
        if (u(m)>zero) u(2:) = u(2:)/u(m)
        u(m)  = one

        this%ubegin = u(1)
        this%uend   = u(m)

        this%has_params = .true.

        endassociate

    end subroutine set_default_parameters

    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    function curve_derivatives(this, u, order, ierr) result(ddx)
       class(fitpack_parametric_curve), intent(inout)  :: this
       real(FP_REAL),          intent(in)     :: u(:)   ! Evaluation points (parameter)
       integer,              intent(in)     :: order  ! Derivative order. Default 1
       integer, optional,    intent(out)    :: ierr   ! Optional error flag
       real(FP_REAL), dimension(this%idim,size(u)) :: ddx

       integer :: i,ierr0

       do i=1,size(u)
          ddx(:,i) = curve_derivative(this,u(i),order,ierr0)
          if (.not.FITPACK_SUCCESS(ierr0)) exit
       end do

       if (present(ierr)) ierr = ierr0

    end function curve_derivatives


end module fitpack_parametric_curves
