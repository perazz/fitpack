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
module fitpack_polar_domains
    use fitpack_core
    use fitpack_fitters
    implicit none
    private

    public :: fitpack_polar

    !> A public type describing a polar fitter z = f(x,y) to scattered polar data,
    !> which is arbitrarily scattered over the polar domain x**2+y**2 <= rad(atan(y/x))**2
    !> through the transform:  x = u*rad(v)*cos(v),
    !>                         y = u*rad(v)*sin(v)
    !> the approximation problem is reduced to the determination of a bi-cubic spline s(u,v) fitting a
    !> corresponding set of data points (u(i),v(i),z(i)) on the rectangle 0<=u<=1,-pi<=v<=pi.

    !> rad is an external real function defining the boundary of the approximation domain, i.e
    !>          x = rad(v)*cos(v) , y = rad(v)*sin(v), -pi <= v <= pi.
    !> It can be a sphere or anything else.
    type, extends(fitpack_fitter) :: fitpack_polar

        !> Scattered data points
        integer :: m = 0
        real(FP_REAL), allocatable :: x(:),y(:)

        !> Function values at the points
        real(FP_REAL), allocatable :: z(:)

        !> Coordinates of the data points in rectangular coordinates (u,v)
        real(FP_REAL), allocatable :: u(:),v(:)

        !> Function that describes shape of the polar boundary: radius as a function of
        !> angle theta:
        !> x = rad(v)*cos(v) , y = rad(v)*sin(v), -pi <= v <= pi.
        procedure(fitpack_polar_boundary), nopass, pointer :: rad => null()

        ! Node weights
        real(FP_REAL), allocatable :: w(:)

        ! Internal Storage
        integer                  :: lwrk2 = 0
        real(FP_REAL), allocatable :: wrk2(:)

        ! Curve behavior
        integer :: bc_continuity_origin = 2 ! Continuity at origin (C0, C1, C2)
        integer :: bc_boundary = OUTSIDE_EXTRAPOLATE ! Extrapolate (=0) or go to zero (=1)

        ! Knots: extimated max number
        integer :: nest(2)  = 0
        integer :: nmax     = 0
        integer :: knots(2) = 0
        real(FP_REAL), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

        contains

           !> Clean memory
           procedure :: destroy    => polar_destroy

           !> Set new points
           procedure :: new_points => polar_new_points

           !> Generate new fit
           procedure :: new_fit    => polr_new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit           => surface_fit_automatic_knots
           procedure :: least_squares => surface_fit_least_squares
           procedure :: interpolate   => surface_fit_interpolating


           !> Evaluate polar domain at given x,y coordinates
           procedure, private :: polr_eval_one
           procedure, private :: polr_eval_many
           generic :: eval => polr_eval_one,polr_eval_many

           !> Communication interface
           procedure :: comm_size   => polar_comm_size
           procedure :: comm_pack   => polar_comm_pack
           procedure :: comm_expand => polar_comm_expand

    end type fitpack_polar

    interface fitpack_polar
       module procedure polr_new_from_points
    end interface fitpack_polar

    contains

    ! Fit a surface to least squares of the current knots
    integer function surface_fit_least_squares(this,smoothing,reset_knots) result(ierr)
       class(fitpack_polar), intent(inout) :: this
       real(FP_REAL), optional, intent(in) :: smoothing
       logical, optional, intent(in) :: reset_knots

       logical :: do_reset

       ! Optionally recompute knots via a smoothing fit first
       do_reset = .false.; if (present(reset_knots)) do_reset = reset_knots
       if (do_reset) then
           this%iopt = IOPT_NEW_SMOOTHING
           ierr = this%fit(smoothing)
           if (.not.FITPACK_SUCCESS(ierr)) return
       end if

       this%iopt = IOPT_NEW_LEASTSQUARES
       ierr = this%fit()

    end function surface_fit_least_squares

    ! Find interpolating surface
    integer function surface_fit_interpolating(this,reset_knots) result(ierr)
        class(fitpack_polar), intent(inout) :: this
        logical, optional, intent(in) :: reset_knots

        logical :: do_reset

        do_reset = .true.; if (present(reset_knots)) do_reset = reset_knots
        if (do_reset) this%iopt = IOPT_NEW_SMOOTHING
        ierr = surface_fit_automatic_knots(this,smoothing=zero,keep_knots=.not.do_reset)

    end function surface_fit_interpolating

    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer function surface_fit_automatic_knots(this,smoothing,keep_knots) result(ierr)
        class(fitpack_polar), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        logical, optional, intent(in) :: keep_knots

        integer :: loop,nit,iopt(3)
        real(FP_REAL) :: smooth_now(3)
        logical :: do_guard

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots (unless caller wants to keep them)
        do_guard = .true.; if (present(keep_knots)) do_guard = .not.keep_knots
        if (do_guard .and. this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! [Continuation, Origin BC, Boundary BC]
            iopt = [this%iopt,this%bc_continuity_origin,this%bc_boundary]

            ! Call fitting routine
            call polar(iopt,                         &  ! Continuation and BCs
                       this%m,this%x,this%y,this%z,  &  ! scattered points, their coordinates and values
                       this%w,                       &  ! Set of weights
                       this%rad,                     &  ! Polar domain boundary
                       this%smoothing,               &  ! Smoothing parameter
                       this%nest(1),this%nest(2),    &  ! Knot space
                       epsilon(zero),                &  ! Threshold for the rank of the linear system
                       this%knots(1),this%t(:,1),    &  ! u (0:1) knots (out)
                       this%knots(2),this%t(:,2),    &  ! v (-pi:pi) knots (out)
                       this%u,this%v,                &  ! co-ordinates of the i-th data point w.r.t the rectangular domain
                       this%c,this%fp,               &  ! Spline representation and MSE
                       this%wrk,this%lwrk,         &  ! memory
                       this%wrk2,this%lwrk2,         &  ! memory
                       this%iwrk,this%liwrk,         &  ! memory
                       ierr)                            ! Error flag

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function surface_fit_automatic_knots


    elemental subroutine polar_destroy(this)
       class(fitpack_polar), intent(inout) :: this
       integer :: ierr
       call this%destroy_base()
       this%m = 0
       nullify(this%rad)
       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%u,stat=ierr)
       deallocate(this%v,stat=ierr)
       deallocate(this%wrk2,stat=ierr)
       deallocate(this%t,stat=ierr)

       this%nest      = 0
       this%lwrk2     = 0
       this%knots     = 0
       this%bc_continuity_origin = 2
       this%bc_boundary = OUTSIDE_EXTRAPOLATE

    end subroutine polar_destroy

    subroutine polar_new_points(this,x,y,z,boundary,w,boundary_bc)
        class(fitpack_polar), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x)),z(size(x))
        procedure(fitpack_polar_boundary) :: boundary
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights
        integer    , optional, intent(in) :: boundary_BC

        integer :: clen,l,k,p,q,iopt2_min,iopt2_max,iopt3_min,iopt3_max
        integer, parameter :: SAFE = 2

        associate(m=>this%m,nest=>this%nest,nmax=>this%nmax)

        call this%destroy()

        m = size(x)

        ! Copy scattered points
        allocate(this%x(m),source=x)
        allocate(this%y(m),source=y)
        allocate(this%z(m),source=z)

        ! Allocate rectangular coordinates
        allocate(this%u(m),this%v(m))

        ! Assign boundary function
        this%rad => boundary

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w)
        else
           allocate(this%w(m),source=one)
        endif

        ! Reset run flag
        this%iopt = 0

        ! Reset boundary BC
        if (present(boundary_BC)) this%bc_boundary = boundary_BC

        ! Knot space: overestimate upper bound
        ! nuest, nvest >=8. In most situations, nuest = nvest = 8+sqrt(m/2) will be sufficient.
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

        ! wrk
        iopt2_max = 2
        iopt3_max = 1
        iopt2_min = 0
        iopt3_min = 0
        k = nest(1)-7
        l = nest(2)-7
        p = 1+iopt2_max*(iopt2_max+3)/2
        q = k+2-iopt2_min-iopt3_min
        this%lwrk = 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
        allocate(this%wrk(this%lwrk),source=zero)

        ! wrk2
        this%lwrk2 = (p+l*q+1)*(4*l+p)+p+l*q
        allocate(this%wrk2(this%lwrk2),source=zero)

        endassociate

    end subroutine polar_new_points

    ! A default constructor
    type(fitpack_polar) function polr_new_from_points(x,y,z,boundary,w,boundary_bc,ierr) result(this)
        real(FP_REAL), intent(in) :: x(:),y(size(x)),z(size(x))
        procedure(fitpack_polar_boundary) :: boundary
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights
        integer, optional, intent(in) :: boundary_bc
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(x,y,z,boundary,w,boundary_bc)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new surface fit')

    end function polr_new_from_points

    real(FP_REAL) function polr_eval_one(this,x,y,ierr) result(z)
        class(fitpack_polar), intent(inout)  :: this
        real(FP_REAL),          intent(in)     :: x,y    ! Evaluation point
        integer, optional,    intent(out)    :: ierr   ! Optional error flag

        z = evapol(this%t(:,1),this%knots(1), & ! u-direction knots
                   this%t(:,2),this%knots(2), & ! v-direction knots
                   this%c,                    & ! b-spline coefficients
                   this%rad,                  & ! Boundary function
                   x,y)                         ! Coordinates of the evaluation point

        call fitpack_error_handling(FITPACK_OK,ierr,'evaluate polar spline')

    end function polr_eval_one

    ! Curve evaluation driver
    function polr_eval_many(this,x,y,ierr) result(z)
        class(fitpack_polar), intent(inout)  :: this
        real(FP_REAL),          intent(in)     :: x(:),y(size(x)) ! Evaluation points
        integer, optional,    intent(out)    :: ierr            ! Optional error flag
        real(FP_REAL) :: z(size(x))

        integer :: i

        do i=1,size(x)
            z(i) = polr_eval_one(this,x(i),y(i),ierr)
        end do

    end function polr_eval_many

    ! Fit a new curve
    integer function polr_new_fit(this,x,y,z,boundary,w,boundary_bc,smoothing)
        class(fitpack_polar), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x)),z(size(x))
        procedure(fitpack_polar_boundary) :: boundary
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights
        integer    , optional, intent(in) :: boundary_bc
        real(FP_REAL), optional, intent(in) :: smoothing

        call this%new_points(x,y,z,boundary,w,boundary_bc)

        polr_new_fit = this%fit(smoothing)

    end function polr_new_fit

    ! Communication: buffer size
    elemental integer(FP_SIZE) function polar_comm_size(this)
        class(fitpack_polar), intent(in) :: this
        polar_comm_size = this%core_comm_size() &
                        + 9 &
                        + FP_COMM_SIZE(this%x) &
                        + FP_COMM_SIZE(this%y) &
                        + FP_COMM_SIZE(this%z) &
                        + FP_COMM_SIZE(this%u) &
                        + FP_COMM_SIZE(this%v) &
                        + FP_COMM_SIZE(this%w) &
                        + FP_COMM_SIZE(this%wrk2) &
                        + FP_COMM_SIZE(this%t)
    end function polar_comm_size

    ! Communication: pack into buffer
    pure subroutine polar_comm_pack(this, buffer)
        class(fitpack_polar), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_pack(buffer)
        pos = this%core_comm_size() + 1

        buffer(pos) = real(this%m, FP_COMM);                       pos = pos + 1
        buffer(pos) = real(this%lwrk2, FP_COMM);                   pos = pos + 1
        buffer(pos) = real(this%bc_continuity_origin, FP_COMM);    pos = pos + 1
        buffer(pos) = real(this%bc_boundary, FP_COMM);             pos = pos + 1
        buffer(pos) = real(this%nest(1), FP_COMM);                 pos = pos + 1
        buffer(pos) = real(this%nest(2), FP_COMM);                 pos = pos + 1
        buffer(pos) = real(this%nmax, FP_COMM);                    pos = pos + 1
        buffer(pos) = real(this%knots(1), FP_COMM);                pos = pos + 1
        buffer(pos) = real(this%knots(2), FP_COMM);                pos = pos + 1

        call FP_COMM_PACK(this%x, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_PACK(this%y, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_PACK(this%z, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_PACK(this%u, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%u)
        call FP_COMM_PACK(this%v, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%v)
        call FP_COMM_PACK(this%w, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_PACK(this%wrk2, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk2)
        call FP_COMM_PACK(this%t, buffer(pos:))
    end subroutine polar_comm_pack

    ! Communication: expand from buffer
    pure subroutine polar_comm_expand(this, buffer)
        class(fitpack_polar), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_expand(buffer)
        pos = this%core_comm_size() + 1

        this%m                       = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%lwrk2                   = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%bc_continuity_origin    = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%bc_boundary             = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nest(1)                 = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nest(2)                 = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nmax                    = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(1)                = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(2)                = nint(buffer(pos), FP_SIZE);  pos = pos + 1

        call FP_COMM_EXPAND(this%x, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_EXPAND(this%y, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_EXPAND(this%z, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_EXPAND(this%u, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%u)
        call FP_COMM_EXPAND(this%v, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%v)
        call FP_COMM_EXPAND(this%w, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_EXPAND(this%wrk2, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk2)
        call FP_COMM_EXPAND(this%t, buffer(pos:))
    end subroutine polar_comm_expand

end module fitpack_polar_domains
