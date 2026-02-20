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
!> @brief OOP wrapper for bivariate surface fitting to data on a rectangular grid.
!!
!! Provides fitpack_grid_surface, a derived type for fitting tensor-product B-spline
!! surfaces \f$ z = s(x, y) \f$ to data values given on a rectangular grid
!! \f$ (x_i, y_j) \f$. The underlying core routine is regrid, which exploits the
!! grid structure for faster fitting than surfit. Supports evaluation, partial
!! derivatives, integration, cross-section extraction, and derivative-spline computation.
!!
!! @see Dierckx, Ch. 5, §5.4 (pp. 98–103); regrid, bispev, parder, pardeu, dblint, profil
module fitpack_grid_surfaces
    use fitpack_core, only: FITPACK_SUCCESS,FP_REAL,FP_SIZE,FP_FLAG,FP_COMM,zero,IOPT_NEW_SMOOTHING,IOPT_OLD_FIT, &
                            IOPT_NEW_LEASTSQUARES,bispev,bispeu,fitpack_error_handling,get_smoothing,regrid, &
                            parder,pardeu,FITPACK_INPUT_ERROR, &
                            dblint,profil,pardtc, &
                            FP_COMM_SIZE,FP_COMM_PACK,FP_COMM_EXPAND
    use fitpack_fitters
    use fitpack_curves, only: fitpack_curve
    implicit none
    private

    public :: fitpack_grid_surface

    !> @brief Bivariate surface fitter \f$ z = s(x, y) \f$ for gridded data.
    !!
    !! Stores grid vectors \f$ x_i \f$ (\f$ i = 1, \ldots, m_x \f$) and
    !! \f$ y_j \f$ (\f$ j = 1, \ldots, m_y \f$) together with gridded function values
    !! \f$ z(j, i) \f$, the fitted tensor-product B-spline representation, and the fitting
    !! state. Uses regrid, which is significantly more efficient than surfit for gridded input.
    !!
    !! @see Dierckx, Ch. 5, §5.4 (pp. 98–103)
    type, extends(fitpack_fitter) :: fitpack_grid_surface

        !> The data points
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
        ! (lwrk/wrk inherited from fitpack_fitter)

        ! Knots
        integer(FP_SIZE) :: knots(2) = 0
        real(FP_REAL), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

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

           !> Evaluate at scattered (x,y) points (bispeu). For gridded output, use eval_ongrid.
           procedure, private :: gridsurf_eval_one
           procedure, private :: gridsurf_eval_many
           generic :: eval => gridsurf_eval_one,gridsurf_eval_many

           !> Evaluate on a rectangular grid x(:) x y(:) (bispev). Returns f(ny,nx).
           procedure, private :: gridded_eval_one
           procedure, private :: gridded_eval_many
           generic :: eval_ongrid => gridded_eval_one,gridded_eval_many
           
           !> Evaluate derivatives at given coordinates
           procedure, private :: gridded_derivatives_gridded
           procedure, private :: gridded_derivatives_many
           procedure, private :: gridded_derivatives_one
           generic   :: dfdx => gridded_derivatives_one,gridded_derivatives_many
           generic   :: dfdx_ongrid => gridded_derivatives_gridded

           !> Double integration over a rectangular domain
           procedure :: integral => gridsurf_integral

           !> Extract a 1D cross-section curve from the surface
           procedure :: cross_section => gridsurf_cross_section

           !> Compute the derivative spline coefficients
           procedure :: derivative_spline => gridsurf_derivative_spline

           !> Parallel communication
           procedure :: comm_size   => gridsurf_comm_size
           procedure :: comm_pack   => gridsurf_comm_pack
           procedure :: comm_expand => gridsurf_comm_expand

    end type fitpack_grid_surface

    interface fitpack_grid_surface
       module procedure surf_new_from_points
    end interface fitpack_grid_surface

    contains

    !> @brief Fit a least-squares gridded surface with fixed knots.
    !!
    !! @see regrid
    integer function surface_fit_least_squares(this,smoothing,reset_knots) result(ierr)
       class(fitpack_grid_surface), intent(inout) :: this
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

    !> @brief Fit an interpolating gridded surface (\f$ s = 0 \f$).
    !!
    !! @see regrid
    integer function surface_fit_interpolating(this,reset_knots) result(ierr)
        class(fitpack_grid_surface), intent(inout) :: this
        logical, optional, intent(in) :: reset_knots

        logical :: do_reset

        do_reset = .true.; if (present(reset_knots)) do_reset = reset_knots
        if (do_reset) this%iopt = IOPT_NEW_SMOOTHING
        ierr = surface_fit_automatic_knots(this,smoothing=zero,keep_knots=.not.do_reset)

    end function surface_fit_interpolating


    !> @brief Fit a smoothing gridded surface with automatic knot placement.
    !!
    !! Uses the regrid core routine, which exploits the rectangular grid structure for
    !! efficiency.
    !!
    !! @see regrid
    integer(FP_FLAG) function surface_fit_automatic_knots(this,smoothing,order,keep_knots) result(ierr)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        integer, optional, intent(in) :: order
        logical, optional, intent(in) :: keep_knots

        integer(FP_SIZE) :: loop,nit
        real(FP_REAL) :: smooth_now(3)
        logical :: do_guard

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots (unless caller wants to keep them)
        do_guard = .true.; if (present(keep_knots)) do_guard = .not.keep_knots
        if (do_guard .and. this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

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


    !> @brief Destroy a grid surface object and release all allocated memory.
    elemental subroutine surf_destroy(this)
       class(fitpack_grid_surface), intent(inout) :: this
       integer :: ierr

       call this%destroy_base()

       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%t,stat=ierr)
       this%left  = zero
       this%right = zero

       this%order     = 3
       this%nest      = 0
       this%nmax      = 0
       this%knots     = 0

    end subroutine surf_destroy

    !> @brief Load new gridded data and allocate workspace.
    !!
    !! @param[in,out] this  The grid surface (destroyed and reinitialized).
    !! @param[in]     x     Grid coordinates in the x direction.
    !! @param[in]     y     Grid coordinates in the y direction.
    !! @param[in]     z     Function values `z(j,i)` = \f$ f(x_i, y_j) \f$.
    !!
    !! @see regrid
    subroutine surf_new_points(this,x,y,z)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(:),z(size(y),size(x))

        integer(FP_SIZE) :: clen,u,m(2)
        integer(FP_SIZE), parameter :: SAFE = 2

        associate(nest=>this%nest,nmax=>this%nmax,order=>this%order)

        call this%destroy()

        m = [size(x),size(y)]

        ! Ensure x are sorted
        allocate(this%x(m(1)),source=x)
        allocate(this%y(m(2)),source=y)
        allocate(this%z(m(2),m(1)),source=z)

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

    !> @brief Construct a grid surface from gridded data and perform a default fit.
    !!
    !! @see regrid
    type(fitpack_grid_surface) function surf_new_from_points(x,y,z,ierr) result(this)
        real(FP_REAL), intent(in) :: x(:),y(:),z(size(y),size(x))
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ierr0

        ierr0 = this%new_fit(x,y,z)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new gridded surface fit')

    end function surf_new_from_points

    !> @brief Load new gridded data and perform a fresh fit.
    !!
    !! @see regrid
    integer function surf_new_fit(this,x,y,z,smoothing,order)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(:),z(size(y),size(x))
        real(FP_REAL), optional, intent(in) :: smoothing
        integer    , optional, intent(in) :: order

        call this%new_points(x,y,z)

        surf_new_fit = this%fit(smoothing,order)

    end function surf_new_fit

    !> @brief Evaluate the grid surface at scattered \f$ (x_i, y_i) \f$ points.
    !!
    !! @see bispeu
    function gridsurf_eval_many(this,x,y,ierr) result(f)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x))
        real(FP_REAL) :: f(size(x))
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ier
        integer(FP_SIZE) :: min_lwrk
        real(FP_REAL), allocatable :: wrk(:)

        ! bispeu workspace: lwrk >= kx+ky+2
        min_lwrk = sum(this%order) + 2
        allocate(wrk(min_lwrk))

        call bispeu(tx=this%t(:,1),nx=this%knots(1),   &
                    ty=this%t(:,2),ny=this%knots(2),   &
                    c=this%c,                          &
                    kx=this%order(1),ky=this%order(2), &
                    x=x,y=y,z=f,m=size(x),            &
                    wrk=wrk,lwrk=min_lwrk,             &
                    ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate grid surface at scattered points')

    end function gridsurf_eval_many

    !> @brief Evaluate the grid surface at a single \f$ (x, y) \f$ point.
    !!
    !! @see bispeu
    real(FP_REAL) function gridsurf_eval_one(this,x,y,ierr) result(f)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x,y
        integer(FP_FLAG), optional, intent(out) :: ierr

        real(FP_REAL) :: z1(1)

        z1 = gridsurf_eval_many(this,[x],[y],ierr)
        f  = z1(1)

    end function gridsurf_eval_one

    !> @brief Evaluate the grid surface on a rectangular evaluation grid.
    !!
    !! Returns `f(j,i)` = \f$ s(x_i, y_j) \f$.
    !!
    !! @see bispev
    function gridded_eval_many(this,x,y,ierr) result(f)
        class(fitpack_grid_surface), intent(inout)  :: this
        real(FP_REAL), intent(in) :: x(:),y(:)  ! Evaluation points
        real(FP_REAL) :: f(size(y),size(x))
        integer(FP_FLAG), optional, intent(out) :: ierr ! Optional error flag

        integer(FP_FLAG) :: ier

        !  evaluation of the spline approximation.
        !  Assume cubic spline in both directions

        ! On successful exit f(j,i) contains the value of s(x,y) at point
        ! (x(i),y(j)),i=1,...,mx; j=1,...,my.
        call bispev(tx=this%t(:,1),nx=this%knots(1), &
                    ty=this%t(:,2),ny=this%knots(2), &
                    c=this%c, &
                    kx=this%order(1),ky=this%order(2), &
                    x=x,mx=size(x), &
                    y=y,my=size(y), &
                    z=f, &
                    wrk=this%wrk,lwrk=this%lwrk, &
                    iwrk=this%iwrk,kwrk=this%liwrk,ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate gridded surface')

    end function gridded_eval_many

    !> @brief Evaluate the grid surface at a single grid point via bispev.
    real(FP_REAL) function gridded_eval_one(this,x,y,ierr) result(f)
        class(fitpack_grid_surface), intent(inout)  :: this
        real(FP_REAL),               intent(in)      :: x,y ! Evaluation point
        integer(FP_FLAG), optional, intent(out)     :: ierr      ! Optional error flag
        real(FP_REAL) :: f1(1,1)

        f1 = gridded_eval_many(this,[x],[y],ierr)
        f  = f1(1,1)

    end function gridded_eval_one

    !> @brief Evaluate partial derivatives on a rectangular grid.
    !!
    !! @see parder
    function gridded_derivatives_gridded(this,x,y,dx,dy,ierr) result(f)
        class(fitpack_grid_surface), intent(inout)  :: this
        
        ! Grid evaluation ranges
        real(FP_REAL),    intent(in) :: x(:),y(:)  
        
        ! Order of the partial derivatives w.r.t. x and y
        integer(FP_SIZE), intent(in) :: dx,dy      
        
        ! f(j,i) contains the value of the specified partial derivative of s(x,y) at (x(i),y(j))
        real(FP_REAL) :: f(size(y),size(x))
        
        ! Optional error flag
        integer(FP_FLAG), optional, intent(out) :: ierr 

        integer(FP_FLAG) :: ier        
        integer(FP_SIZE) :: min_lwrk,min_kwrk,m(2)
        real(FP_REAL), allocatable :: min_wrk(:)
        integer(FP_SIZE), allocatable :: min_iwrk(:)
        
        ! Check if the internal working storage is large enough: reallocate it if necessary
        m = [size(x),size(y)]
        
        ! Assert real working storage
        min_lwrk = sum(m*(this%order+1)) + product(this%knots-this%order-1)   
        if (min_lwrk>this%lwrk) then 
            allocate(min_wrk(min_lwrk),source=0.0_FP_REAL)
            call move_alloc(from=min_wrk,to=this%wrk)
            this%lwrk = min_lwrk
        end if
        
        ! Assert integer working storage
        min_kwrk = sum(m)
        if (min_kwrk>this%liwrk) then 
            allocate(min_iwrk(min_kwrk),source=0_FP_SIZE)
            call move_alloc(from=min_iwrk,to=this%iwrk)
            this%liwrk = min_kwrk
        end if
        
        ! On successful exit f(j,i) contains the value of s(x,y) at point
        ! (x(i),y(j)),i=1,...,mx; j=1,...,my.
        call parder(tx=this%t(:,1),nx=this%knots(1),   &    ! position of the knots in the x-direction
                    ty=this%t(:,2),ny=this%knots(2),   &    ! position of the knots in the y-direction
                    c=this%c,                          &    ! the b-spline coefficients                  
                    kx=this%order(1),ky=this%order(2), &    ! the degrees of the spline
                    nux=dx,nuy=dy,                     &    ! order of the derivatives 0<=nux<kx, 0<=nuy<ky 
                    x=x,mx=size(x),                    &    ! x grid points: tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2:mx.
                    y=y,my=size(y),                    &    ! y grid points: ty(ky+1)<=y(i-1)<=y(i)<=ty(ny-ky), i=2:mx.
                    z=f,                               &    ! Value of the partial derivative
                    wrk=this%wrk,lwrk=this%lwrk,       &    ! memory
                    iwrk=this%iwrk,kwrk=this%liwrk,    &    ! memory
                    ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate gridded derivatives')

    end function gridded_derivatives_gridded

    !> @brief Evaluate partial derivatives at scattered \f$ (x_i, y_i) \f$ points.
    !!
    !! @see pardeu
    function gridded_derivatives_many(this,x,y,dx,dy,ierr) result(f)
        class(fitpack_grid_surface), intent(inout)  :: this
        
        ! Grid evaluation ranges
        real(FP_REAL),    intent(in) :: x(:),y(:)  
        
        ! Order of the partial derivatives w.r.t. x and y
        integer(FP_SIZE), intent(in) :: dx,dy      
        
        ! f(j,i) contains the value of the specified partial derivative of s(x,y) at (x(i),y(j))
        real(FP_REAL) :: f(size(x))
        
        ! Optional error flag
        integer(FP_FLAG), optional, intent(out) :: ierr 

        integer(FP_FLAG) :: ier        
        integer(FP_SIZE) :: min_lwrk,min_kwrk,m(2)
        real(FP_REAL), allocatable :: min_wrk(:)
        integer(FP_SIZE), allocatable :: min_iwrk(:)
        
        ! Check if the input arrays have the same size
        m = [size(x),size(y)]
        if (m(1)/=m(2)) then 
            
            ier = FITPACK_INPUT_ERROR
            
        else
        
            ! Assert real working storage
            min_lwrk = sum(m*(this%order+1)) + product(this%knots-this%order-1)   
            if (min_lwrk>this%lwrk) then 
                allocate(min_wrk(min_lwrk),source=0.0_FP_REAL)
                call move_alloc(from=min_wrk,to=this%wrk)
                this%lwrk = min_lwrk
            end if
            
            ! Assert integer working storage
            min_kwrk = sum(m)
            if (min_kwrk>this%liwrk) then 
                allocate(min_iwrk(min_kwrk),source=0_FP_SIZE)
                call move_alloc(from=min_iwrk,to=this%iwrk)
                this%liwrk = min_kwrk
            end if
            
            ! On successful exit f(j,i) contains the value of s(x,y) at point
            ! (x(i),y(j)),i=1,...,mx; j=1,...,my.
            call pardeu(tx=this%t(:,1),nx=this%knots(1),   &    ! position of the knots in the x-direction
                        ty=this%t(:,2),ny=this%knots(2),   &    ! position of the knots in the y-direction
                        c=this%c,                          &    ! the b-spline coefficients                  
                        kx=this%order(1),ky=this%order(2), &    ! the degrees of the spline
                        nux=dx,nuy=dy,                     &    ! order of the derivatives 0<=nux<kx, 0<=nuy<ky 
                        x=x,y=y,                           &    ! list of (x,y) points
                        z=f,                               &    ! Value of the partial derivative
                        m=m(1),                            &    ! Number of input points
                        wrk=this%wrk,lwrk=this%lwrk,       &    ! memory
                        iwrk=this%iwrk,kwrk=this%liwrk,    &    ! memory
                        ier=ier)

            
        end if
        
        call fitpack_error_handling(ier,ierr,'evaluate bivariate derivatives')

    end function gridded_derivatives_many

    !> @brief Evaluate a partial derivative at a single \f$ (x, y) \f$ point.
    !!
    !! @see pardeu
    real(FP_REAL) function gridded_derivatives_one(this,x,y,dx,dy,ierr) result(f)
        class(fitpack_grid_surface), intent(inout) :: this
        
        ! Evaluation point
        real(FP_REAL),          intent(in)      :: x,y    

        ! Order of the partial derivatives w.r.t. x and y
        integer(FP_SIZE), intent(in) :: dx,dy              
        
        ! Optional error flag
        integer(FP_FLAG), optional, intent(out) :: ierr   

        real(FP_REAL) :: z1(1),x1(1),y1(1)
        
        x1 = x
        y1 = y

        z1 = gridded_derivatives_many(this,x1,y1,dx,dy,ierr)
        f  = z1(1)

    end function gridded_derivatives_one

    !> @brief Compute the double integral of the grid surface over a rectangular domain.
    !!
    !! @see dblint
    real(FP_REAL) function gridsurf_integral(this, lower, upper)
        class(fitpack_grid_surface), intent(in) :: this
        real(FP_REAL), intent(in) :: lower(2), upper(2)

        real(FP_REAL) :: wrk(this%knots(1)+this%knots(2)-this%order(1)-this%order(2)-2)

        gridsurf_integral = dblint(this%t(:,1), this%knots(1), &
                                   this%t(:,2), this%knots(2), &
                                   this%c, &
                                   this%order(1), this%order(2), &
                                   lower(1), upper(1), lower(2), upper(2), wrk)

    end function gridsurf_integral

    !> @brief Extract a 1D cross-section curve from the grid surface.
    !!
    !! If `along_y=.true.`, returns \f$ f(y) = s(u, y) \f$; otherwise \f$ g(x) = s(x, u) \f$.
    !!
    !! @see profil
    function gridsurf_cross_section(this, u, along_y, ierr) result(curve)
        class(fitpack_grid_surface), intent(in) :: this
        real(FP_REAL), intent(in) :: u
        logical, intent(in) :: along_y
        integer(FP_FLAG), optional, intent(out) :: ierr
        type(fitpack_curve) :: curve

        integer(FP_FLAG) :: ier
        integer(FP_SIZE) :: iopt, nu, nc
        real(FP_REAL), allocatable :: cu(:)

        if (along_y) then
            iopt = 0
            nu   = this%knots(2)
            nc   = this%knots(2) - this%order(2) - 1
        else
            iopt = 1
            nu   = this%knots(1)
            nc   = this%knots(1) - this%order(1) - 1
        end if

        allocate(cu(nu))

        call profil(iopt, this%t(:,1), this%knots(1), &
                          this%t(:,2), this%knots(2), &
                          this%c, &
                          this%order(1), this%order(2), &
                          u, nu, cu, ier)

        if (FITPACK_SUCCESS(ier)) then
            curve%order = merge(this%order(2), this%order(1), along_y)
            curve%knots = nu
            curve%nest  = nu
            if (along_y) then
                allocate(curve%t(nu), source=this%t(1:nu, 2))
            else
                allocate(curve%t(nu), source=this%t(1:nu, 1))
            end if
            allocate(curve%c(nu), source=zero)
            curve%c(:nc) = cu(:nc)
            curve%xleft  = curve%t(curve%order + 1)
            curve%xright = curve%t(curve%knots - curve%order)
        end if

        call fitpack_error_handling(ier, ierr, 'grid surface cross-section')

    end function gridsurf_cross_section

    !> @brief Compute the B-spline representation of a partial derivative surface.
    !!
    !! Returns a new grid surface with reduced degrees and trimmed knots.
    !!
    !! @see pardtc
    function gridsurf_derivative_spline(this, nux, nuy, ierr) result(dsurf)
        class(fitpack_grid_surface), intent(in) :: this
        integer(FP_SIZE), intent(in) :: nux, nuy
        integer(FP_FLAG), optional, intent(out) :: ierr
        type(fitpack_grid_surface) :: dsurf

        integer(FP_FLAG) :: ier
        integer(FP_SIZE) :: nx, ny, newkx, newky, newnx, newny, nc_old, nc_new
        real(FP_REAL), allocatable :: newc(:)

        nx = this%knots(1)
        ny = this%knots(2)
        nc_old = (nx - this%order(1) - 1) * (ny - this%order(2) - 1)

        allocate(newc(nc_old))

        call pardtc(this%t(:,1), nx, this%t(:,2), ny, &
                    this%c, this%order(1), this%order(2), &
                    nux, nuy, newc, ier)

        if (FITPACK_SUCCESS(ier)) then
            newkx  = this%order(1) - nux
            newky  = this%order(2) - nuy
            newnx  = nx - 2*nux
            newny  = ny - 2*nuy
            nc_new = (newnx - newkx - 1) * (newny - newky - 1)

            dsurf%order  = [newkx, newky]
            dsurf%knots  = [newnx, newny]
            dsurf%nmax   = max(newnx, newny)
            dsurf%nest   = dsurf%nmax

            allocate(dsurf%t(dsurf%nmax, 2), source=zero)
            dsurf%t(1:newnx, 1) = this%t(1+nux:nx-nux, 1)
            dsurf%t(1:newny, 2) = this%t(1+nuy:ny-nuy, 2)

            allocate(dsurf%c(nc_new), source=newc(1:nc_new))

            dsurf%left  = [dsurf%t(newkx+1, 1), dsurf%t(newky+1, 2)]
            dsurf%right = [dsurf%t(newnx-newkx, 1), dsurf%t(newny-newky, 2)]

            ! Allocate evaluation workspace (bispev needs lwrk >= mx*(kx+1)+my*(ky+1), kwrk >= mx+my)
            dsurf%lwrk  = newnx*(newkx+1) + newny*(newky+1) + nc_new
            dsurf%liwrk = newnx + newny
            allocate(dsurf%wrk(dsurf%lwrk), source=zero)
            allocate(dsurf%iwrk(dsurf%liwrk), source=0_FP_SIZE)
        end if

        call fitpack_error_handling(ier, ierr, 'grid surface derivative spline')

    end function gridsurf_derivative_spline


    ! =================================================================================================
    ! PARALLEL COMMUNICATION
    ! =================================================================================================

    !> @brief Return the communication buffer size for the grid surface.
    elemental integer(FP_SIZE) function gridsurf_comm_size(this)
        class(fitpack_grid_surface), intent(in) :: this
        ! Base fields + grid-surface-specific scalars:
        ! order(2), left(2), right(2), nest(2), nmax, knots(2) = 11
        gridsurf_comm_size = this%core_comm_size() &
                           + 11 &
                           + FP_COMM_SIZE(this%x) &
                           + FP_COMM_SIZE(this%y) &
                           + FP_COMM_SIZE(this%z) &
                           + FP_COMM_SIZE(this%t)
    end function gridsurf_comm_size

    !> @brief Pack grid surface data into a communication buffer.
    pure subroutine gridsurf_comm_pack(this, buffer)
        class(fitpack_grid_surface), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_pack(buffer)
        pos = this%core_comm_size() + 1

        buffer(pos) = real(this%order(1), FP_COMM);  pos = pos + 1
        buffer(pos) = real(this%order(2), FP_COMM);  pos = pos + 1
        buffer(pos) = this%left(1);                   pos = pos + 1
        buffer(pos) = this%left(2);                   pos = pos + 1
        buffer(pos) = this%right(1);                  pos = pos + 1
        buffer(pos) = this%right(2);                  pos = pos + 1
        buffer(pos) = real(this%nest(1), FP_COMM);   pos = pos + 1
        buffer(pos) = real(this%nest(2), FP_COMM);   pos = pos + 1
        buffer(pos) = real(this%nmax, FP_COMM);      pos = pos + 1
        buffer(pos) = real(this%knots(1), FP_COMM);  pos = pos + 1
        buffer(pos) = real(this%knots(2), FP_COMM);  pos = pos + 1

        call FP_COMM_PACK(this%x, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_PACK(this%y, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_PACK(this%z, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_PACK(this%t, buffer(pos:))
    end subroutine gridsurf_comm_pack

    !> @brief Expand grid surface data from a communication buffer.
    pure subroutine gridsurf_comm_expand(this, buffer)
        class(fitpack_grid_surface), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_expand(buffer)
        pos = this%core_comm_size() + 1

        this%order(1) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%order(2) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%left(1)  = buffer(pos);                  pos = pos + 1
        this%left(2)  = buffer(pos);                  pos = pos + 1
        this%right(1) = buffer(pos);                  pos = pos + 1
        this%right(2) = buffer(pos);                  pos = pos + 1
        this%nest(1)  = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nest(2)  = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nmax     = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(1) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(2) = nint(buffer(pos), FP_SIZE);  pos = pos + 1

        call FP_COMM_EXPAND(this%x, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_EXPAND(this%y, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_EXPAND(this%z, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_EXPAND(this%t, buffer(pos:))
    end subroutine gridsurf_comm_expand

end module fitpack_grid_surfaces
