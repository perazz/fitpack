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
module fitpack_sphere_domains
    use fitpack_core
    implicit none
    private

    public :: fitpack_sphere

    !> A public type describing a smooth bicubic surface fitter r = s(teta,phi),
    !  0 <= teta <= pi   "Latitude"
    !> 0 <= phi <= 2*pi  "Longitude"
    !> which is arbitrarily scattered over the sphere domain
    type :: fitpack_sphere

        !> Scattered data points
        integer :: m = 0
        real(RKIND), allocatable :: theta(:),phi(:)

        !> Function values at the points
        real(RKIND), allocatable :: r(:)

        ! Node weights
        real(RKIND), allocatable :: w(:)

        ! Internal Storage
        integer                  :: lwrk1 = 0, lwrk2 = 0, liwrk = 0
        integer, allocatable     :: iwrk(:)
        real(RKIND), allocatable :: wrk1(:),wrk2(:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(RKIND) :: smoothing = 1000.d0

        ! Actual curve MSE
        real(RKIND) :: fp = zero

        ! Knots: extimated max number
        integer :: nest(2)  = 0
        integer :: nmax     = 0
        integer :: knots(2) = 0
        real(RKIND), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

        ! Coefficients of the spline approximation
        real(RKIND), allocatable :: c(:)

        ! Runtime flag
        integer :: iopt = IOPT_NEW_SMOOTHING ! -> iopt(1)

        contains

           !> Clean memory
           procedure :: destroy    => sphere_destroy

           !> Set new points
           procedure :: new_points => sphere_new_points

           !> Generate new fit
           procedure :: new_fit    => sphere_new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit           => surface_fit_automatic_knots
           procedure :: least_squares => surface_fit_least_squares
           procedure :: interpolate   => surface_fit_interpolating

           !> Evaluate polar domain at given x,y coordinates
           procedure, private :: sphere_eval_one
           procedure, private :: sphere_eval_many
           generic :: eval => sphere_eval_one,sphere_eval_many

    end type fitpack_sphere

    interface fitpack_sphere
       module procedure sphere_new_from_points
    end interface fitpack_sphere

    contains

    ! Fit a surface to least squares of the current knots
    integer function surface_fit_least_squares(this) result(ierr)
       class(fitpack_sphere), intent(inout) :: this

       this%iopt = IOPT_NEW_LEASTSQUARES
       ierr = this%fit()

    end function surface_fit_least_squares

    ! Find interpolating surface
    integer function surface_fit_interpolating(this) result(ierr)
        class(fitpack_sphere), intent(inout) :: this

        ! Set zero smoothing
        ierr = surface_fit_automatic_knots(this,smoothing=zero)

    end function surface_fit_interpolating

    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer function surface_fit_automatic_knots(this,smoothing) result(ierr)
        class(fitpack_sphere), intent(inout) :: this
        real(RKIND), optional, intent(in) :: smoothing

        integer :: loop,nit
        real(RKIND) :: smooth_now(3)

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! Call fitting routine
            call sphere(this%iopt,                   &  ! Continuation and BCs
                       this%m,this%theta,this%phi,   &  ! scattered points
                       this%r,                       &  ! spline values at the points
                       this%w,                       &  ! Set of weights
                       this%smoothing,               &  ! Smoothing parameter
                       this%nest(1),this%nest(2),    &  ! Knot space
                       epsilon(zero),                &  ! Threshold for the rank of the linear system
                       this%knots(1),this%t(:,1),    &  ! theta  knots (out)
                       this%knots(2),this%t(:,2),    &  ! phi knots (out)
                       this%c,this%fp,               &  ! Spline representation and MSE
                       this%wrk1,this%lwrk1,         &  ! memory
                       this%wrk2,this%lwrk2,         &  ! memory
                       this%iwrk,this%liwrk,         &  ! memory
                       ierr)                            ! Error flag

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function surface_fit_automatic_knots


    elemental subroutine sphere_destroy(this)
       class(fitpack_sphere), intent(inout) :: this
       integer :: ierr
       this%m = 0
       deallocate(this%theta,stat=ierr)
       deallocate(this%phi,stat=ierr)
       deallocate(this%r,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this%wrk1,stat=ierr)
       deallocate(this%wrk2,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)

       this%smoothing = 1000.0_RKIND
       this%nest      = 0
       this%lwrk1     = 0
       this%lwrk2     = 0
       this%liwrk     = 0
       this%knots     = 0
       this%fp        = zero

    end subroutine sphere_destroy

    subroutine sphere_new_points(this,theta,phi,r,w)
        class(fitpack_sphere), intent(inout) :: this
        real(RKIND), intent(in) :: theta(:),phi(size(theta)),r(size(theta))
        real(RKIND), optional, intent(in) :: w(size(theta)) ! node weights

        integer :: clen,u,v
        integer, parameter :: SAFE = 2

        associate(m=>this%m,nest=>this%nest,nmax=>this%nmax)

        call this%destroy()

        m = size(theta)

        ! Copy scattered points
        this%theta = theta
        this%phi   = phi
        this%r     = r

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w)
        else
           allocate(this%w(m),source=one)
        endif

        ! Reset run flag
        this%iopt = 0

        ! Knot space: overestimate upper bound
        ! ntest, npest >=8. In most situations, ntest = npest = 8+sqrt(m/2) will be sufficient.
        nest = SAFE*(8+ceiling(sqrt(m*half)))
        nmax = maxval(nest)
        allocate(this%t(nmax,2),source=zero)

        ! Spline coefficients
        clen = product(nest-4)
        allocate(this%c(clen),source=zero)

        this%fp = zero

        ! Working space
        this%liwrk = m+product(nest-7)
        allocate(this%iwrk(this%liwrk),source=0)

        ! wrk1
        u = nest(1)-7
        v = nest(2)-7
        this%lwrk1 = 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
        allocate(this%wrk1(this%lwrk1),source=zero)

        ! wrk2
        this%lwrk2 = 48+21*v+7*u*v+4*(u-1)*v**2
        allocate(this%wrk2(this%lwrk2),source=zero)

        endassociate

    end subroutine sphere_new_points

    ! A default constructor
    type(fitpack_sphere) function sphere_new_from_points(theta,phi,r,w,ierr) result(this)
        real(RKIND), intent(in) :: theta(:),phi(size(theta)),r(size(theta))
        real(RKIND), optional, intent(in) :: w(size(theta)) ! node weights
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(theta,phi,r,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new surface fit')

    end function sphere_new_from_points

    function sphere_eval_many(this,theta,phi,ierr) result(r)
        class(fitpack_sphere), intent(inout)  :: this
        real(RKIND), intent(in) :: theta(:),phi(:)  ! Evaluation points
        real(RKIND) :: r(size(phi),size(theta))
        integer, optional, intent(out) :: ierr ! Optional error flag

        integer :: ier

        !  evaluation of the spline approximation.
        !  Assume cubic spline in both directions

        ! On successful exit r(j,i) contains the value of s(x,y) at point
        ! (x(i),y(j)),i=1,...,mx; j=1,...,my.
        call bispev(tx=this%t(:,1),nx=this%knots(1), &
                    ty=this%t(:,2),ny=this%knots(2), &
                    c=this%c, &
                    kx=3,ky=3, &
                    x=theta,mx=size(theta), &
                    y=phi,my=size(phi), &
                    z=r, &
                    wrk=this%wrk2,lwrk=this%lwrk2, &
                    iwrk=this%iwrk,kwrk=this%liwrk,ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate sphere spline')

    end function sphere_eval_many

    ! Curve evaluation driver
    real(RKIND) function sphere_eval_one(this,theta,phi,ierr) result(r)
        class(fitpack_sphere), intent(inout)  :: this
        real(RKIND),          intent(in)     :: theta,phi ! Evaluation point
        integer, optional,    intent(out)    :: ierr      ! Optional error flag
        real(RKIND) :: r1(1,1)

        r1 = sphere_eval_many(this,[theta],[phi],ierr)
        r = r1(1,1)

    end function sphere_eval_one

    ! Fit a new curve
    integer function sphere_new_fit(this,theta,phi,r,w,smoothing)
        class(fitpack_sphere), intent(inout) :: this
        real(RKIND), intent(in) :: theta(:),phi(size(theta)),r(size(theta))
        real(RKIND), optional, intent(in) :: w(size(theta)) ! node weights
        real(RKIND), optional, intent(in) :: smoothing

        call this%new_points(theta,phi,r,w)

        sphere_new_fit = this%fit(smoothing)

    end function sphere_new_fit

end module fitpack_sphere_domains
