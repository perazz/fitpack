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
module fitpack_gridded_polar
    use fitpack_core
    implicit none
    private

    public :: fitpack_grid_polar

    !> A public type describing a polar fitter z = f(x,y) to GRIDDED polar data,
    !> which is distributed across the domain x**2+y**2 <= rad(atan(y/x))**2
    !> through the transform:  x = u*rad*cos(v),
    !>                         y = u*rad*sin(v)
    !> Gridded data is provided in terms of (u,v) coordinates and the *constant* boundary radius, r
    type :: fitpack_grid_polar

        !> Coordinates of the data points in grid coordinates (u,v) and domain size
        real(FP_REAL), allocatable :: u(:),v(:)
        real(FP_REAL) :: r = zero

        !> Function values at the points
        real(FP_REAL), allocatable :: z(:,:)    ! Function values z(iu,iv)

        ! Origin of the polar domain
        !> Function value at the origin
        real(FP_REAL) :: z0 = zero              ! Function value at the origin

        !> Fit spline to z0 exactly at the origin
        logical :: z0_present = .false.
        logical :: z0_exact   = .false.

        !> Decide if the origin should have zero-gradient BCs in both u and v
        logical :: z0_zero_gradient = .false.

        !> Node weights are not allowed

        ! Estimated and actual number of knots and their allocations
        integer :: nest(2)  = 0
        integer :: nmax     = 0
        integer                  :: lwrk = 0, liwrk = 0
        integer, allocatable     :: iwrk(:)
        real(FP_REAL), allocatable :: wrk (:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(FP_REAL) :: smoothing = 1000.d0

        ! Actual curve MSE
        real(FP_REAL) :: fp = zero

        ! Knots
        integer     :: knots(2) = 0
        real(FP_REAL), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

        ! Spline coefficients [knots-order-1]
        real(FP_REAL), allocatable :: c(:)

        ! Curve behavior
        integer :: bc_continuity_origin = 1 ! Continuity at origin (C0, C1)
        integer :: bc_boundary = OUTSIDE_EXTRAPOLATE ! Extrapolate (=0) or go to zero (=1)

        ! Runtime flag
        integer :: iopt = IOPT_NEW_SMOOTHING ! -> iopt(1)

        contains

           !> Clean memory
           procedure :: destroy       => surf_destroy

           !> Set new points
           procedure :: new_points    => surf_new_points

           !> Generate new fit
           procedure :: new_fit       => surf_new_fit

           !> Set origin BC
           procedure :: set_origin_BC

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit           => polr_fit_automatic_knots
           procedure :: least_squares => polr_fit_least_squares
           procedure :: interpolate   => polr_fit_interpolating

           !> Evaluate gridded domain at given x,y coordinates
           procedure, private :: gridded_eval_one
           procedure, private :: gridded_eval_many
           generic :: eval => gridded_eval_one,gridded_eval_many

           !> Write to disk
           procedure :: write => gridded_to_disk

    end type fitpack_grid_polar

    interface fitpack_grid_polar
       module procedure surf_new_from_points
    end interface fitpack_grid_polar

    contains

    ! Fit a surface to least squares of the current knots
    integer function polr_fit_least_squares(this) result(ierr)
       class(fitpack_grid_polar), intent(inout) :: this

       this%iopt = IOPT_NEW_LEASTSQUARES
       ierr = this%fit()

    end function polr_fit_least_squares

    ! Find interpolating surface
    integer function polr_fit_interpolating(this) result(ierr)
        class(fitpack_grid_polar), intent(inout) :: this

        ! Set zero smoothing
        ierr = polr_fit_automatic_knots(this,smoothing=zero)

    end function polr_fit_interpolating


    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer function polr_fit_automatic_knots(this,smoothing) result(ierr)
        class(fitpack_grid_polar), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing

        integer :: loop,nit,iopt(3),ider(2)
        real(FP_REAL) :: smooth_now(3)

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots
        if (this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! [Continuation, Origin continuity, Boundary BC]
            iopt = [this%iopt,this%bc_continuity_origin,this%bc_boundary]

            ! [Fit z0 exactly, z0 zero gradient BC]
            ider = [merge(merge(1,0,this%z0_exact),-1,this%z0_present), &
                    merge(1,0,this%z0_zero_gradient .and. this%bc_continuity_origin>0)]

            call pogrid(iopt,                         &  ! Continuation and BCs
                        ider,                         &  ! Origin point behavior
                        size(this%u),this%u,          &  ! U grid
                        size(this%v),this%v,          &  ! V grid
                        this%z,this%z0,               &  ! Gridded function values
                        this%r,                       &  ! Polar domain boundary
                        this%smoothing,               &  ! Smoothing parameter
                        this%nest(1),this%nest(2),    &  ! Knot space
                        this%knots(1),this%t(:,1),    &  ! u (0:1) knots (out)
                        this%knots(2),this%t(:,2),    &  ! v (-pi:pi) knots (out)
                        this%c,this%fp,               &  ! Spline representation and MSE
                        this%wrk,this%lwrk,           &  ! memory
                        this%iwrk,this%liwrk,         &  ! memory
                        ierr)                            ! Error flag

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function polr_fit_automatic_knots


    elemental subroutine surf_destroy(this)
       class(fitpack_grid_polar), intent(inout) :: this
       integer :: ierr
       this%r = zero
       this%z0 = zero
       this%z0_exact = .false.
       this%z0_present = .false.
       this%z0_zero_gradient = .false.
       deallocate(this%u,stat=ierr)
       deallocate(this%v,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this%wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)

       this%smoothing = 1000.0_FP_REAL
       this%iopt      = 0
       this%bc_boundary = OUTSIDE_EXTRAPOLATE
       this%bc_continuity_origin = 1
       this%nest      = 0
       this%nmax      = 0
       this%lwrk      = 0
       this%liwrk     = 0
       this%knots     = 0
       this%fp        = zero

    end subroutine surf_destroy

    subroutine surf_new_points(this,u,v,r,z,z0)
        class(fitpack_grid_polar), intent(inout) :: this
        real(FP_REAL), intent(in) :: u(:),v(:),r ! polar domain
        real(FP_REAL), intent(in) :: z(size(v),size(u)) ! Gridded values
        real(FP_REAL), optional, intent(in) :: z0 ! Origin value (optional)

        integer :: clen,m(2),q
        integer, parameter :: SAFE = 2

        associate(nest=>this%nest,nmax=>this%nmax)

        call this%destroy()

        m = [size(u),size(v)] ! /= shape(z), == shape(transpose(z))

        ! Set domain radius
        this%r = r

        ! Copy grid and data
        allocate(this%u,source=u)
        allocate(this%v,source=v)
        allocate(this%z,source=z)

        call set_origin_BC(this,z0)

        ! Reset run flag
        this%iopt = IOPT_NEW_SMOOTHING

        ! Knot space: overestimate (mu+5+iopt(2)+iopt(3),mv+7)
        nest = SAFE*(m + 7)
        nmax = maxval(nest)
        allocate(this%t(nmax,2),source=zero)

        ! Spline coefficients
        clen = product(nest-4) ! nest-order-1, fixed order==3
        allocate(this%c(clen),source=zero)

        this%fp = zero

        ! Working space
        this%liwrk = 4+sum(m+nest)
        allocate(this%iwrk(this%liwrk),source=0)

        ! wrk
        ! lwrk >= 8+nuest*(mv+nvest+3)+nvest*21+4*mu+6*mv+q
        q = max(m(2)+nest(2),nest(1))
        this%lwrk = 8+nest(1)*(m(2)+nest(2)+3)+21*nest(2)+4*m(1)+6*m(2)+q
        allocate(this%wrk(this%lwrk),source=zero)

        endassociate

    end subroutine surf_new_points

    ! A default constructor
    type(fitpack_grid_polar) function surf_new_from_points(u,v,r,z,z0,ierr) result(this)
        real(FP_REAL), intent(in) :: u(:),v(:),r ! polar domain
        real(FP_REAL), intent(in) :: z(size(v),size(u)) ! Gridded values
        real(FP_REAL), optional, intent(in) :: z0 ! Origin value (optional)
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(u,v,r,z,z0)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new gridded surface fit')

    end function surf_new_from_points

    ! Fit a new curve
    integer function surf_new_fit(this,u,v,r,z,z0,smoothing)
        class(fitpack_grid_polar), intent(inout) :: this
        real(FP_REAL), intent(in) :: u(:),v(:),r ! polar domain
        real(FP_REAL), intent(in) :: z(size(v),size(u)) ! Gridded values
        real(FP_REAL), optional, intent(in) :: z0 ! Origin value (optional)
        real(FP_REAL), optional, intent(in) :: smoothing

        call this%new_points(u,v,r,z,z0)

        surf_new_fit = this%fit(smoothing)

    end function surf_new_fit

    function gridded_eval_many(this,u,v,ierr) result(f)
        class(fitpack_grid_polar), intent(inout)  :: this
        real(FP_REAL), intent(in) :: u(:),v(:)  ! Evaluation grid points (polar coordinates)
        real(FP_REAL) :: f(size(v),size(u))
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
                    x=u,mx=size(u), &
                    y=v,my=size(v), &
                    z=f, & ! output in format (j,i)
                    wrk=this%wrk,lwrk=this%lwrk, &
                    iwrk=this%iwrk,kwrk=this%liwrk,ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate gridded surface')

    end function gridded_eval_many

    ! Curve evaluation driver
    real(FP_REAL) function gridded_eval_one(this,u,v,ierr) result(f)
        class(fitpack_grid_polar), intent(inout)  :: this
        real(FP_REAL),          intent(in)      :: u,v ! Evaluation point (grid polar coordinates)
        integer, optional,    intent(out)     :: ierr      ! Optional error flag
        real(FP_REAL) :: f1(1,1)

        f1 = gridded_eval_many(this,[u],[v],ierr)
        f  = f1(1,1)

    end function gridded_eval_one

    subroutine set_origin_BC(this,z0,exact,differentiable)
        class(fitpack_grid_polar), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: z0 ! Function value at origin
        logical, optional, intent(in) :: exact,differentiable

        this%z0_present = present(z0)
        if (present(exact)) then
            this%z0_exact = exact
        else
            this%z0_exact = .false.
        endif

        if (present(differentiable)) then
            this%bc_continuity_origin = merge(1,0,differentiable)
        else
            this%bc_continuity_origin = 0
        end if

    end subroutine set_origin_BC

    !> Print gridded polar data to disk
    subroutine gridded_to_disk(this,fileName)
        class(fitpack_grid_polar), intent(inout) :: this
        character(*), intent(in) :: fileName

        integer :: iunit,ierr,j,i

        open(newunit=iunit,file=fileName,action='write',form='formatted',iostat=ierr)

        ! Write header
        write(iunit,1) 'u',(named('v',j)       ,j=1,size(this%v))
        write(iunit,1) ' ',(numbered(this%v(j)),j=1,size(this%v))

        do i=1,size(this%u)
            write(iunit,2) this%u(i),this%z(:,i)
        end do

        close(iunit)

        1 format('!',a12,*(1x,a12))
        2 format(*(1x,1pe12.5))

        contains

           character(12) function named(name,i)
              character(*), intent(in) :: name
              integer, intent(in) :: i
              write(named,"(a,'_',i0)") name,i
              named = adjustr(named)
           end function named
           character(12) function numbered(v)
              real(FP_REAL), intent(in) :: v
              write(numbered,"(1pe12.5)") v
              numbered = adjustr(numbered)
           end function numbered

    end subroutine gridded_to_disk

end module fitpack_gridded_polar
