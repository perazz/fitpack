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
    use fitpack_core
    implicit none
    private

    public :: fitpack_grid_surface

    !> A public type describing a surface fitter z = s(x,y) to gridded x,y data
    type :: fitpack_grid_surface

        !> The data points
        integer :: m = 0
        real(RKIND), allocatable :: x(:),y(:) ! Grid values in x, y dimensions
        real(RKIND), allocatable :: z(:,:)    ! Function values z(ix,iy)

        !> Spline degree
        integer :: order(2) = 3

        !> Interval boundaries
        real(RKIND) :: left(2),right(2)

        !> Node weights are not allowed

        ! Estimated and actual number of knots and their allocations
        integer :: nest(2)  = 0
        integer :: nmax     = 0
        integer                  :: lwrk = 0, liwrk = 0
        integer, allocatable     :: iwrk(:)
        real(RKIND), allocatable :: wrk (:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(RKIND) :: smoothing = 1000.d0

        ! Actual curve MSE
        real(RKIND) :: fp = zero

        ! Knots
        integer     :: knots(2) = 0
        real(RKIND), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

        ! Spline coefficients [knots-order-1]
        real(RKIND), allocatable :: c(:)

        ! Runtime flag
        integer :: iopt = 0

        contains

           !> Clean memory
           procedure :: destroy    => surf_destroy

           !> Set new points
           procedure :: new_points => surf_new_points

           !> Generate new fit
           procedure :: new_fit    => surf_new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit        => surface_fit_automatic_knots

    end type fitpack_grid_surface

    interface fitpack_grid_surface
       module procedure surf_new_from_points
    end interface fitpack_grid_surface

    contains

    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer function surface_fit_automatic_knots(this,smoothing,order) result(ierr)
        class(fitpack_grid_surface), intent(inout) :: this
        real(RKIND), optional, intent(in) :: smoothing
        integer, optional, intent(in) :: order

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

        ! User may want to change the order for both x and y
        if (present(order)) this%order = order

        ! First iteration lets solver decide knots
        this%iopt = 0

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

       this%smoothing = 1000.0_RKIND
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
        real(RKIND), intent(in) :: x(:),y(:),z(size(x),size(y))

        integer :: clen,u,m(2)
        integer, parameter :: SAFE = 2


        associate(nest=>this%nest,nmax=>this%nmax,order=>this%order)

        call this%destroy()

        m = shape(z)

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
        this%iopt = 0

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
        this%lwrk = 4+sum(nest*(m+2*order+5))+sum(m*(order+1)) +u
        allocate(this%wrk(this%lwrk),source=zero)

        endassociate

    end subroutine surf_new_points

    ! A default constructor
    type(fitpack_grid_surface) function surf_new_from_points(x,y,z,ierr) result(this)
        real(RKIND), intent(in) :: x(:),y(:),z(size(x),size(y))
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(x,y,z)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new gridded surface fit')

    end function surf_new_from_points

    ! Fit a new curve
    integer function surf_new_fit(this,x,y,z,smoothing,order)
        class(fitpack_grid_surface), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(:),z(size(x),size(y))
        real(RKIND), optional, intent(in) :: smoothing
        integer    , optional, intent(in) :: order

        call this%new_points(x,y,z)

        surf_new_fit = this%fit(smoothing,order)

    end function surf_new_fit

end module fitpack_grid_surfaces
