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
    use fitpack_core, only: FITPACK_SUCCESS,FP_REAL,FP_SIZE,FP_FLAG,FP_COMM,zero,IOPT_NEW_SMOOTHING,IOPT_OLD_FIT, &
                            IOPT_NEW_LEASTSQUARES,bispev,fitpack_error_handling,get_smoothing,regrid, &
                            parder,pardeu,FITPACK_INPUT_ERROR, &
                            FP_COMM_SIZE,FP_COMM_PACK,FP_COMM_EXPAND
    use fitpack_fitters
    implicit none
    private

    public :: fitpack_grid_surface

    !> A public type describing a surface fitter z = s(x,y) to gridded x,y data
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
        integer(FP_SIZE)                  :: lwrk = 0, liwrk = 0
        real(FP_REAL), allocatable :: wrk(:)

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

           !> Evaluate gridded domain at given x,y coordinates
           procedure, private :: gridded_eval_one
           procedure, private :: gridded_eval_many
           generic :: eval => gridded_eval_one,gridded_eval_many
           
           !> Evaluate derivatives at given coordinates
           procedure, private :: gridded_derivatives_gridded
           procedure, private :: gridded_derivatives_many
           procedure, private :: gridded_derivatives_one
           generic   :: dfdx => gridded_derivatives_one,gridded_derivatives_many
           generic   :: dfdx_ongrid => gridded_derivatives_gridded

           !> Parallel communication
           procedure :: comm_size   => gridsurf_comm_size
           procedure :: comm_pack   => gridsurf_comm_pack
           procedure :: comm_expand => gridsurf_comm_expand

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

       call this%destroy_base()

       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       this%left  = zero
       this%right = zero

       this%order     = 3
       this%nest      = 0
       this%nmax      = 0
       this%lwrk      = 0
       this%liwrk     = 0
       this%knots     = 0

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

    ! Curve evaluation driver
    real(FP_REAL) function gridded_eval_one(this,x,y,ierr) result(f)
        class(fitpack_grid_surface), intent(inout)  :: this
        real(FP_REAL),               intent(in)      :: x,y ! Evaluation point
        integer(FP_FLAG), optional, intent(out)     :: ierr      ! Optional error flag
        real(FP_REAL) :: f1(1,1)

        f1 = gridded_eval_many(this,[x],[y],ierr)
        f  = f1(1,1)

    end function gridded_eval_one

    !> Evaluate derivatives on a grid domain
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

    !> Evaluate derivatives on a list of (x(i),y(i)) points
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

    ! =================================================================================================
    ! PARALLEL COMMUNICATION
    ! =================================================================================================

    elemental integer(FP_SIZE) function gridsurf_comm_size(this)
        class(fitpack_grid_surface), intent(in) :: this
        ! Base fields + grid-surface-specific scalars:
        ! order(2), left(2), right(2), nest(2), nmax, lwrk, liwrk, knots(2) = 13
        gridsurf_comm_size = this%core_comm_size() &
                           + 13 &
                           + FP_COMM_SIZE(this%x) &
                           + FP_COMM_SIZE(this%y) &
                           + FP_COMM_SIZE(this%z) &
                           + FP_COMM_SIZE(this%t) &
                           + FP_COMM_SIZE(this%wrk)
    end function gridsurf_comm_size

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
        buffer(pos) = real(this%lwrk, FP_COMM);      pos = pos + 1
        buffer(pos) = real(this%liwrk, FP_COMM);     pos = pos + 1
        buffer(pos) = real(this%knots(1), FP_COMM);  pos = pos + 1
        buffer(pos) = real(this%knots(2), FP_COMM);  pos = pos + 1

        call FP_COMM_PACK(this%x, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_PACK(this%y, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_PACK(this%z, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_PACK(this%t, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_PACK(this%wrk, buffer(pos:))
    end subroutine gridsurf_comm_pack

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
        this%lwrk     = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%liwrk    = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(1) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(2) = nint(buffer(pos), FP_SIZE);  pos = pos + 1

        call FP_COMM_EXPAND(this%x, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_EXPAND(this%y, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_EXPAND(this%z, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_EXPAND(this%t, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_EXPAND(this%wrk, buffer(pos:))
    end subroutine gridsurf_comm_expand

end module fitpack_grid_surfaces
