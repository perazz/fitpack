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
module fitpack_grid_surfaces
    use fitpack_core, only: FITPACK_SUCCESS,FP_REAL,FP_SIZE,FP_FLAG,zero,IOPT_NEW_SMOOTHING,IOPT_OLD_FIT, &
                            IOPT_NEW_LEASTSQUARES,bispev,fitpack_error_handling,get_smoothing,regrid, pardtc
    implicit none
    private

    public :: fitpack_grid_surface

    !> A public type describing a surface fitter z = s(x,y) to gridded x,y data
    type :: fitpack_grid_surface

        !> The data points
        integer :: m = 0
        real(FP_REAL), allocatable :: x(:),y(:) ! Grid values in x, y dimensions
        real(FP_REAL), allocatable :: z(:,:)    ! Function values z(iy,ix)

        !> Spline degree
        integer(FP_SIZE) :: order(2) = 3

        !> Interval boundaries
        real(FP_REAL) :: left(2),right(2)

        !> Node weights are not allowed

        ! Estimated and actual number of knots and their allocations
        integer(FP_SIZE) :: nest(2)  = 0
        integer(FP_SIZE) :: nmax     = 0
        integer(FP_SIZE)                  :: lwrk = 0, liwrk = 0
        integer(FP_SIZE), allocatable     :: iwrk(:)
        real(FP_REAL), allocatable :: wrk (:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(FP_REAL) :: smoothing = 1000.d0

        ! Actual curve MSE
        real(FP_REAL) :: fp = zero

        ! Knots
        integer(FP_SIZE) :: knots(2) = 0
        real(FP_REAL), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

        ! Spline coefficients [knots-order-1]
        real(FP_REAL), allocatable :: c(:)

        ! Runtime flag
        integer(FP_FLAG) :: iopt = IOPT_NEW_SMOOTHING

        contains

           !> Clean memory
           procedure :: destroy       => surf_destroy

           !> Set new points
           procedure :: new_points    => surf_new_points

           !> Generate new fit
           procedure :: new_fit       => surf_new_fit

           !> Generate/update fitting curve, with optional /Users/federico/code/fitpack/test/fitpack_curve_tests.f90smoothing
           procedure :: fit           => surface_fit_automatic_knots
           procedure :: least_squares => surface_fit_least_squares
           procedure :: interpolate   => surface_fit_interpolating

           !> Evaluate gridded domain at given x,y coordinates
           procedure, private :: gridded_eval_one
           procedure, private :: gridded_eval_many
           generic :: eval => gridded_eval_one,gridded_eval_many

    end type fitpack_grid_surface

    interface fitpack_grid_surface
       module procedure surf_new_from_points
    end interface fitpack_grid_surface

    contains

    ! Fit a surface to least squares of the current knots
    integer function surface_fit_least_squares(this) result(ierr)
       class(fitpack_grid_surface), intent(inout) :: this

       this%iopt = IOPT_NEW_LEASTSQUARES
       ierr = this%fit()

    end function surface_fit_least_squares

    ! Find interpolating surface
    integer function surface_fit_interpolating(this) result(ierr)
        class(fitpack_grid_surface), intent(inout) :: this

        ! Set zero smoothing
        ierr = surface_fit_automatic_knots(this,smoothing=zero)

    end function surface_fit_interpolating


    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer(FP_FLAG) function surface_fit_automatic_knots(this,smoothing,order) result(ierr)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        integer, optional, intent(in) :: order

        integer(FP_SIZE) :: loop,nit
        real(FP_REAL) :: smooth_now(3)

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots
        if (this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        ! User may want to change the order for both x and y
        if (present(order)) this%order = order
        
        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)
            
            call regrid(this%iopt,                   &  ! [-1]=lsq on given knots; [0,1]=smoothing spline
                        size(this%x),this%x,         &  ! x coordinate of the grid points
                        size(this%y),this%y,         &  ! y coordinate of the grid points
                        this%z,                      &  ! z(ix,jy) gridded data points
                        this%left(1),this%right(1),  &  ! x range
                        this%left(2),this%right(2),  &  ! y range
                        this%order(1),this%order(2), &  ! [1:5] x,y spline order. Recommended: bicubic (x=y=3)
                        this%smoothing,              &  ! spline accuracy (iopt>=0)
                        this%nest(1),this%nest(2),   &  ! estimated number of knots and storage nxest >= 2*(kx+1), nyest >= 2*(ky+1)
                        this%knots(1),this%t(:,1),   &  ! x knots (out)
                        this%knots(2),this%t(:,2),   &  ! y knots (out)
                        this%c,this%fp,              &  ! spline output. size(c)>=(nxest-kx-1)*(nyest-ky-1)
                        this%wrk,this%lwrk,          &  ! memory
                        this%iwrk,this%liwrk,        &  ! memory
                        ierr)                           ! Error flag

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function surface_fit_automatic_knots


    elemental subroutine surf_destroy(this)
       class(fitpack_grid_surface), intent(inout) :: this
       integer :: ierr
       this%m = 0
       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this%wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)
       this%left  = zero
       this%right = zero

       this%smoothing = 1000.0_FP_REAL
       this%order     = 3
       this%iopt      = 0
       this%nest      = 0
       this%nmax      = 0
       this%lwrk      = 0
       this%liwrk     = 0
       this%knots     = 0
       this%fp        = zero

    end subroutine surf_destroy

    subroutine surf_new_points(this,x,y,z)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(:),z(size(y),size(x))

        integer(FP_SIZE) :: clen,u,m(2)
        integer(FP_SIZE), parameter :: SAFE = 2

        associate(nest=>this%nest,nmax=>this%nmax,order=>this%order)

        call this%destroy()

        m = [size(x),size(y)]

        ! Ensure x are sorted
        allocate(this%x,source=x)
        allocate(this%y,source=y)
        allocate(this%z,source=z)

        ! Setup boundaries
        this%left(1)  = minval(x,1)
        this%left(2)  = minval(y,1)
        this%right(1) = maxval(x,1)
        this%right(2) = maxval(y,1)

        ! Reset run flag
        this%iopt = IOPT_NEW_SMOOTHING

        ! Knot space: overestimate (2*order+1 => order+m+1)
        nest = SAFE*(order + m + 1)
        nmax = maxval(nest)
        nest = nmax
        allocate(this%t(nmax,2),source=zero)

        ! Spline coefficients
        clen = product(nest-order-1)
        allocate(this%c(clen),source=zero)

        this%fp = zero

        ! Working space
        this%liwrk = 3+sum(m+nest)
        allocate(this%iwrk(this%liwrk),source=0)

        ! wrk
        ! lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+ my*(ky+1) +u
        u = max(m(2),nest(1))
        ! Do not use sum() or it wil segfault gfortran 13
        this%lwrk = 4+u+ (nest(1)*(m(1)+2*order(1)+5)) + (m(1)*(order(1)+1)) &
                        + (nest(2)*(m(2)+2*order(2)+5)) + (m(2)*(order(2)+1))
        allocate(this%wrk(this%lwrk),source=zero)

        endassociate

    end subroutine surf_new_points

    ! A default constructor
    type(fitpack_grid_surface) function surf_new_from_points(x,y,z,ierr) result(this)
        real(FP_REAL), intent(in) :: x(:),y(:),z(size(y),size(x))
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ierr0

        ierr0 = this%new_fit(x,y,z)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new gridded surface fit')

    end function surf_new_from_points

    ! Fit a new curve
    integer function surf_new_fit(this,x,y,z,smoothing,order)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(:),z(size(y),size(x))
        real(FP_REAL), optional, intent(in) :: smoothing
        integer    , optional, intent(in) :: order

        call this%new_points(x,y,z)

        surf_new_fit = this%fit(smoothing,order)

    end function surf_new_fit

    type(fitpack_grid_surface) function surf_derivative(this, nux, nuy, ierr) result(f)
        class(fitpack_grid_surface), intent(in) :: this
        integer, intent(in) :: nux, nuy
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ierr0
        ! everything is copied to the output 
        f = this

        !call pardtc(tx,nx,ty,ny,c,kx,ky,nux,nuy,newc,ierr0)

        ! Error handling
        call fitpack_error_handling(ierr0, ierr, 'derivative of the spline')

    end function

    function gridded_eval_many(this,x,y,ierr) result(f)
        class(fitpack_grid_surface), intent(inout)  :: this
        real(FP_REAL), intent(in) :: x(:),y(:)  ! Evaluation points
        real(FP_REAL) :: f(size(y),size(x))
        integer(FP_FLAG), optional, intent(out) :: ierr ! Optional error flag

        integer(FP_FLAG) :: ier
        integer(FP_SIZE) :: mx, my, kx, ky
        integer(FP_SIZE) :: lwrk, liwrk
        
        !  evaluation of the spline approximation.
        !  Assume cubic spline in both directions

        ! On successful exit r(j,i) contains the value of s(x,y) at point
        ! (x(i),y(j)),i=1,...,mx; j=1,...,my.

        mx = size(x)
        my = size(y)
        kx = this%order(1)
        ky = this%order(2)
        ! lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
        
        if (this%lwrk < mx*(kx+1)+my*(ky+1)) then
            if (allocated(this%wrk)) deallocate(this%wrk)
            this%lwrk = mx*(kx+1)+my*(ky+1)
            allocate(this%wrk(this%lwrk))
        endif
        if (this%liwrk < mx+my) then
            if (allocated(this%iwrk)) deallocate(this%iwrk)
            this%liwrk = mx+my
            allocate(this%iwrk(this%liwrk))
        endif

        call bispev(tx=this%t(:,1),nx=this%knots(1), &
                    ty=this%t(:,2),ny=this%knots(2), &
                    c=this%c, &
                    kx=kx,ky=ky, &
                    x=x,mx=size(x), &
                    y=y,my=size(y), &
                    z=f, &
                    wrk=this%wrk,lwrk=this%lwrk, &
                    iwrk=this%iwrk,kwrk=this%liwrk,ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate gridded surface')

    end function gridded_eval_many

    ! Curve evaluation driver
    real(FP_REAL) function gridded_eval_one(this,x,y,ierr) result(f)
        class(fitpack_grid_surface), intent(inout)  :: this
        real(FP_REAL),               intent(in)      :: x,y ! Evaluation point
        integer(FP_FLAG), optional, intent(out)     :: ierr      ! Optional error flag
        real(FP_REAL) :: f1(1,1)

        f1 = gridded_eval_many(this,[x],[y],ierr)
        f  = f1(1,1)

    end function gridded_eval_one

end module fitpack_grid_surfaces
