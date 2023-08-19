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
module fitpack_parametric_surfaces
    use fitpack_core
    implicit none
    private

    public :: fitpack_parametric_surface

    !> A public type describing a bicubic parametric surface fitter defined by points z(j,i,:) in the
    !> idim-dimensional space, organized on a grid of strictly increasing parameter values u(i), v(j)
    type :: fitpack_parametric_surface

        !> Number of dimensions
        integer :: idim = 0

        !> Parameter grid values
        real(RKIND), allocatable :: u(:),v(:)

        !> idim-dimensional function values [size(v),size(u),idim]
        real(RKIND), allocatable :: z(:,:,:)

        !> Flags to determine whether dimensions u,v are periodic
        logical :: periodic_dim(2) = .false.

        !> Estimated and actual number of knots and their allocations
        integer :: nest(2)  = 0
        integer :: nmax     = 0
        integer                  :: lwrk = 0, liwrk = 0
        integer, allocatable     :: iwrk(:)
        real(RKIND), allocatable :: wrk (:)

!        ! Space for derivative evaluation
!        real(RKIND), allocatable :: dd(:,:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(RKIND) :: smoothing = 1000.0_RKIND

        ! Actual curve MSE
        real(RKIND) :: fp = zero

        ! Knots
        integer     :: knots(2) = 0
        real(RKIND), allocatable :: t(:,:)  ! Knot locations (:,1)=u; (:,2)=v

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

           !> Evaluate derivative at given coordinates
           procedure, private :: curve_derivative
           procedure, private :: curve_derivatives
           procedure, private :: curve_all_derivatives
           generic   :: dfdx => curve_derivative,curve_derivatives
           generic   :: dfdx_all => curve_all_derivatives

           !> Properties: MSE
           procedure, non_overridable :: mse => curve_error

    end type fitpack_parametric_surface

    ! Default constructor
    interface fitpack_parametric_surface
       module procedure new_from_points
    end interface fitpack_parametric_surface

    contains

    ! A default constructor
    type(fitpack_parametric_surface) function new_from_points(u,v,z,periodic_BC,ierr) result(this)
        real(RKIND), intent(in) :: u(:),v(:),z(:,:,:)
        logical    , optional, intent(in) :: periodic_BC(2)
        integer    , optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(u,v,z,periodic_BC=periodic_BC)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new curve fit')

    end function new_from_points

    ! Fit a new curve
    integer function new_fit(this,u,v,z,smoothing,periodic_BC)
        class(fitpack_parametric_surface), intent(inout) :: this
        real(RKIND), intent(in) :: u(:),v(:),z(:,:,:)
        real(RKIND), optional, intent(in) :: smoothing
        logical    , optional, intent(in) :: periodic_BC(2)

        call this%new_points(u,v,z,periodic_BC)

        new_fit = this%fit(smoothing)

    end function new_fit

    elemental subroutine destroy(this)
       class(fitpack_parametric_surface), intent(inout) :: this
       integer :: ierr
       this%idim = 0
       this%periodic_dim = .false.
       deallocate(this%u,stat=ierr)
       deallocate(this%v,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this% wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)

       this%smoothing = 1000.0_RKIND
       this%iopt      = IOPT_NEW_SMOOTHING
       this%nest      = 0
       this%lwrk      = 0
       this%liwrk     = 0
       this%knots     = 0
       this%fp        = 0.0_RKIND

    end subroutine destroy

    subroutine new_points(this,u,v,z,periodic_BC)
        class(fitpack_parametric_surface), intent(inout) :: this
        real(RKIND), intent(in) :: u(:),v(:),z(:,:,:)
        logical, optional, intent(in) :: periodic_BC(2)

        integer, parameter   :: SAFE = 2
        integer :: m(2),q,clen


        associate(idim=>this%idim,nest=>this%nest,nmax=>this%nmax,lwrk=>this%lwrk)

        call this%destroy()

        ! Surface size
        m    = [size(u),size(v)]
        idim = size(z,3)

        ! Optional periodic BC
        if (present(periodic_BC)) this%periodic_dim = periodic_BC

        ! Should we do bounds checking here?
        !if (size(z,1)/=m(2) .or. size(z,2)/=m(1)) then
        !end if

        ! Load data
        allocate(this%u,source=u)
        allocate(this%v,source=v)
        allocate(this%z,source=z)

        ! Reset run flag
        this%iopt = IOPT_NEW_SMOOTHING

        ! Knot space: overestimate (mv+4+2*ipar(2))
        nest = SAFE*(m+4+2)
        nmax = maxval(nest)
        allocate(this%t(nmax,2),source=zero)

        ! Spline coefficients
        clen = product(nest-4)*idim
        allocate(this%c(clen),source=zero)

        this%fp = zero

        ! Working space
        this%liwrk = 3+sum(m+nest)
        allocate(this%iwrk(this%liwrk),source=0)

        ! wrk
        ! lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+4*(mu+mv)+q*idim
        ! where q is the larger of mv and nuest.
        q = max(m(2),nest(1))
        this%lwrk = 4+nest(1)*(m(2)*idim+11+4)+nest(2)*(11+4)+4*sum(m)+q*idim
        allocate(this%wrk(this%lwrk),source=zero)

        endassociate

    end subroutine new_points

    function curve_eval_one(this,u,ierr) result(y)
        class(fitpack_parametric_surface), intent(inout)  :: this
        real(RKIND),          intent(in)     :: u      ! Evaluation point
        integer, optional,    intent(out)    :: ierr   ! Optional error flag
        real(RKIND) :: y(this%idim)

        real(RKIND) :: y1(this%idim,1)

        y1 = curve_eval_many(this,[u],ierr)
        y  = y1(:,1)

    end function curve_eval_one

    ! Curve evaluation driver
    function curve_eval_many(this,u,ierr) result(x)
        class(fitpack_parametric_surface), intent(inout)  :: this
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

!        call curev(idim=this%idim,&                ! dimension of the spline curve
!                   t=this%t,&                      ! the position of the knots
!                   n=this%knots,&                  ! total number of knots of s(u(x))
!                   c=this%c,&                      ! the b-spline coefficients
!                   nc=size(this%c),&               ! number of coefficients of s(u(x))
!                   k=this%order,&                  ! the degree of s(u(x))
!                   u=u,m=npts,&                    ! the points where s(u) must be evaluated.
!                   x=x,mx=size(x), &               ! the predictions
!                   ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate 1d spline')

    end function curve_eval_many

    ! Interpolating curve
    integer function interpolating_curve(this) result(ierr)
        class(fitpack_parametric_surface), intent(inout) :: this

        ! Set zero smoothing
        ierr = curve_fit_automatic_knots(this,zero)

    end function interpolating_curve

    ! Curve fitting driver: automatic number of knots
    integer function curve_fit_automatic_knots(this,smoothing) result(ierr)
        class(fitpack_parametric_surface), intent(inout) :: this
        real(RKIND), optional, intent(in) :: smoothing

        integer :: loop,nit,ipar(2)

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

            ! Call fitting function
            select type (surf => this)

               class is (fitpack_parametric_surface)

                  ipar = merge(1,0,surf%periodic_dim)

                  call parsur(surf%iopt,                    &  ! option
                              ipar,                         &  ! periodic dimension [u,v] flags
                              surf%idim,                    &  ! Number of dimensions
                              size(surf%u),surf%u,          &  ! U-parameter coordinates
                              size(Surf%v),surf%v,          &  ! V-parameter coordinates
                              surf%z,                       &  ! Values on the U-V grid [size(v),size(u),idim]
                              surf%smoothing,               &  ! Spline accuracy
                              surf%nest(1),surf%nest(2),    &  ! spline output
                              surf%knots(1),surf%t(:,1),    &  ! u knots (out)
                              surf%knots(2),surf%t(:,2),    &  ! v knots (out)
                              surf%c,surf%fp,               &  ! spline output. size(c)>=(nxest-kx-1)*(nyest-ky-1)
                              surf%wrk,surf%lwrk,           &  ! memory
                              surf%iwrk,surf%liwrk,         &  ! memory
                              ierr)

            end select

            ! After any successful call, parameters have surely been computed.
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function curve_fit_automatic_knots

    ! Return fitting MSE
    elemental real(RKIND) function curve_error(this)
       class(fitpack_parametric_surface), intent(in) :: this
       curve_error = this%fp
    end function curve_error

    !> Evaluate k-th derivative of the curve at point u
    function curve_derivative(this, u, order, ierr) result(ddx)
       class(fitpack_parametric_surface), intent(inout) :: this
       real(RKIND),          intent(in)    :: u      ! Evaluation points (parameter)
       integer,              intent(in)    :: order  ! Derivative order. 0=function; 1:k=i-th derivative
       integer, optional,    intent(out)   :: ierr   ! Optional error flag
       real(RKIND), dimension(this%idim)   :: ddx

       integer :: ddx_order,ierr0

       ! Choose order
       ddx_order = max(0,order)

       ierr0 = FITPACK_OK

       !  subroutine cualde evaluates at the point u all the derivatives
       !                (l)
       !     d(j,l) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim
       !  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim.
!       call cualde(this%idim,    & ! Number of dimensions
!                   this%t,       & ! Position of the knots
!                   this%knots,   & ! Number of knots
!                   this%c,       & ! Spline coefficients
!                   size(this%c), & ! Number of coefficients
!                   this%order+1, & ! k1 = order of s(u) (order = degree+1)
!                   u,            & ! Where the derivatives must be evaluated
!                   this%dd,      & ! Space for derivative evaluation
!                   size(this%dd),& ! Its size
!                   ierr0)          ! Output flag

       ! Derivative order is 0:k <- extract derivative
       !ddx = this%dd(:,1+ddx_order)

       call fitpack_error_handling(ierr0,ierr,'evaluate derivative')

    end function curve_derivative

    !> Evaluate all derivatives (0:k) of the curve at point u
    function curve_all_derivatives(this, u, ierr) result(ddx)
       class(fitpack_parametric_surface), intent(inout) :: this
       real(RKIND),          intent(in)    :: u      ! Evaluation points (parameter)
       integer, optional,    intent(out)   :: ierr   ! Optional error flag
       real(RKIND), dimension(this%idim,0:3) :: ddx

       integer :: ierr0

       ierr0 = FITPACK_OK
!
!       !  subroutine cualde evaluates at the point u all the derivatives
!       !                (l)
!       !     d(j,l) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim
!       !  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim.
!       call cualde(this%idim,    & ! Number of dimensions
!                   this%t,       & ! Position of the knots
!                   this%knots,   & ! Number of knots
!                   this%c,       & ! Spline coefficients
!                   size(this%c), & ! Number of coefficients
!                   this%order+1, & ! k1 = order of s(u) (order = degree+1)
!                   u,            & ! Where the derivatives must be evaluated
!                   this%dd,      & ! Space for derivative evaluation
!                   size(this%dd),& ! Its size
!                   ierr0)          ! Output flag

       ! Derivative order is 0:k <- extract derivative
       !ddx = this%dd(:,0:this%order)

       call fitpack_error_handling(ierr0,ierr,'evaluate all derivatives')

    end function curve_all_derivatives


    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    function curve_derivatives(this, u, order, ierr) result(ddx)
       class(fitpack_parametric_surface), intent(inout)  :: this
       real(RKIND),          intent(in)     :: u(:)   ! Evaluation points (parameter)
       integer,              intent(in)     :: order  ! Derivative order. Default 1
       integer, optional,    intent(out)    :: ierr   ! Optional error flag
       real(RKIND), dimension(this%idim,size(u)) :: ddx

       integer :: i,ierr0

       do i=1,size(u)
          ddx(:,i) = curve_derivative(this,u(i),order,ierr0)
          if (.not.FITPACK_SUCCESS(ierr0)) exit
       end do

       if (present(ierr)) ierr = ierr0

    end function curve_derivatives


end module fitpack_parametric_surfaces
