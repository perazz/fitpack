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
module fitpack_surfaces
    use fitpack_core
    use fitpack_fitters
    implicit none
    private

    public :: fitpack_surface

    !> A public type describing a bivariate surface fitter z = s(x,y) to scattered x,y data
    type, extends(fitpack_fitter) :: fitpack_surface

        !> The input data points
        integer :: m = 0
        real(FP_REAL), allocatable :: x(:),y(:),z(:)

        !> Spline degree
        integer :: order(2) = 3

        !> Interval boundaries
        real(FP_REAL) :: left(2),right(2)

        ! Node weights
        real(FP_REAL), allocatable :: w(:)

        ! Estimated and actual number of knots and their allocations
        integer :: nest(2)  = 0
        integer :: nmax = 0
        integer                  :: lwrk1 = 0, lwrk2 = 0
        real(FP_REAL), allocatable :: wrk1(:),wrk2(:)

        ! Curve extrapolation behavior
        integer     :: bc = OUTSIDE_NEAREST_BND

        ! Knots
        integer     :: knots(2) = 0
        real(FP_REAL), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

        contains

           !> Clean memory
           procedure :: destroy    => surf_destroy

           !> Set new points
           procedure :: new_points => surf_new_points

           !> Generate new fit
           procedure :: new_fit    => surf_new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit           => surface_fit_automatic_knots
           procedure :: interpolate   => surface_fit_interpolating
           procedure :: least_squares => surface_fit_least_squares

           !> Evaluate surface at given coordinates
           procedure, private :: surface_eval_one
           procedure, private :: surface_eval_many
           generic   :: eval => surface_eval_one,surface_eval_many

           !> Evaluate surface on a grid
           procedure, private :: surface_eval_gridded
           generic   :: eval_ongrid => surface_eval_gridded

           !> Evaluate derivatives at given coordinates
           procedure, private :: surface_derivatives_gridded
           procedure, private :: surface_derivatives_many
           procedure, private :: surface_derivatives_one
           generic   :: dfdx => surface_derivatives_one,surface_derivatives_many
           generic   :: dfdx_ongrid => surface_derivatives_gridded

           !> Parallel communication interface
           procedure :: comm_size   => surf_comm_size
           procedure :: comm_pack   => surf_comm_pack
           procedure :: comm_expand => surf_comm_expand

    end type fitpack_surface

    interface fitpack_surface
       module procedure surf_new_from_points
    end interface fitpack_surface

    contains

    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer(FP_FLAG) function surface_fit_automatic_knots(this,smoothing,order) result(ierr)
        class(fitpack_surface), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        integer, optional, intent(in) :: order

        integer :: loop,nit
        real(FP_REAL) :: smooth_now(3)

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots
        if (this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        ! User may want to change the order for both x and y
        if (present(order)) this%order = order

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! Call curvfit
            call surfit(this%iopt,                   &  ! [-1]=lsq on given knots; [0,1]=smoothing spline
                        this%m,this%x,this%y,this%z, &  ! points and their coordinates
                        this%w,                      &  ! weights
                        this%left(1),this%right(1),  &  ! x range
                        this%left(2),this%right(2),  &  ! y range
                        this%order(1),this%order(2), &  ! [1:5] x,y spline order. Recommended: bicubic (x=y=3)
                        this%smoothing,              &  ! spline accuracy (iopt>=0)
                        this%nest(1),this%nest(2),   &  ! estimated number of knots and storage nxest >= 2*(kx+1), nyest >= 2*(ky+1)
                        this%nmax,                   &  ! nmax>=max(nxest,nyest) size of tx,ty
                        epsilon(zero),               &  ! numeric limits of the current kind
                        this%knots(1),this%t(:,1),   &  ! x knots (out)
                        this%knots(2),this%t(:,2),   &  ! y knots (out)
                        this%c,this%fp,              &  ! spline output. size(c)>=(nxest-kx-1)*(nyest-ky-1)
                        this%wrk1,this%lwrk1,        &  ! memory
                        this%wrk2,this%lwrk2,        &  ! memory
                        this%iwrk,this%liwrk,        &  ! memory
                        ierr)                           ! Error flag

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function surface_fit_automatic_knots


    ! Find interpolating surface
    integer(FP_FLAG) function surface_fit_interpolating(this) result(ierr)
        class(fitpack_surface), intent(inout) :: this
        ierr = surface_fit_automatic_knots(this,smoothing=zero)
    end function surface_fit_interpolating

    ! Fit a surface to least squares of the current knots
    integer(FP_FLAG) function surface_fit_least_squares(this) result(ierr)
       class(fitpack_surface), intent(inout) :: this
       this%iopt = IOPT_NEW_LEASTSQUARES
       ierr = this%fit()
    end function surface_fit_least_squares

    !> Evaluate surface on a list of (x(i),y(i)) scattered points using bispeu
    function surface_eval_many(this,x,y,ierr) result(f)
        class(fitpack_surface), intent(inout) :: this
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

        call fitpack_error_handling(ier,ierr,'evaluate bivariate surface')

    end function surface_eval_many

    !> Evaluate surface at a single (x,y) point
    real(FP_REAL) function surface_eval_one(this,x,y,ierr) result(f)
        class(fitpack_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x,y
        integer(FP_FLAG), optional, intent(out) :: ierr

        real(FP_REAL) :: z1(1)

        z1 = surface_eval_many(this,[x],[y],ierr)
        f  = z1(1)

    end function surface_eval_one

    !> Evaluate surface on a grid domain using bispev
    function surface_eval_gridded(this,x,y,ierr) result(f)
        class(fitpack_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(:)
        real(FP_REAL) :: f(size(y),size(x))
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ier
        integer(FP_SIZE) :: min_lwrk,min_kwrk,m(2)
        real(FP_REAL), allocatable :: min_wrk(:)
        integer(FP_SIZE), allocatable :: min_iwrk(:)

        m = [size(x),size(y)]

        ! Assert real working storage
        min_lwrk = sum(m*(this%order+1)) + product(this%knots-this%order-1)
        if (min_lwrk>this%lwrk1) then
            allocate(min_wrk(min_lwrk),source=0.0_FP_REAL)
            call move_alloc(from=min_wrk,to=this%wrk1)
            this%lwrk1 = min_lwrk
        end if

        ! Assert integer working storage
        min_kwrk = sum(m)
        if (min_kwrk>this%liwrk) then
            allocate(min_iwrk(min_kwrk),source=0_FP_SIZE)
            call move_alloc(from=min_iwrk,to=this%iwrk)
            this%liwrk = min_kwrk
        end if

        call bispev(tx=this%t(:,1),nx=this%knots(1),   &
                    ty=this%t(:,2),ny=this%knots(2),   &
                    c=this%c,                          &
                    kx=this%order(1),ky=this%order(2), &
                    x=x,mx=size(x),                    &
                    y=y,my=size(y),                    &
                    z=f,                               &
                    wrk=this%wrk1,lwrk=this%lwrk1,     &
                    iwrk=this%iwrk,kwrk=this%liwrk,    &
                    ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate gridded surface')

    end function surface_eval_gridded

    elemental subroutine surf_destroy(this)
       class(fitpack_surface), intent(inout) :: this
       integer :: ierr
       call this%destroy_base()
       this%m = 0
       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%wrk1,stat=ierr)
       deallocate(this%wrk2,stat=ierr)
       deallocate(this%t,stat=ierr)
       this%left  = zero
       this%right = zero

       this%order     = 3
       this%nest      = 0
       this%nmax      = 0
       this%lwrk1     = 0
       this%lwrk2     = 0
       this%knots     = 0
       this%bc        = OUTSIDE_NEAREST_BND

    end subroutine surf_destroy

    subroutine surf_new_points(this,x,y,z,w)
        class(fitpack_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x)),z(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights

        integer :: clen,uv(2),km,bxy(2),b1,b2
        integer, parameter :: SAFE = 2


        associate(m=>this%m,nest=>this%nest,nmax=>this%nmax,order=>this%order)

        call this%destroy()

        m = size(x)

        ! Ensure x are sorted
        allocate(this%x(m),source=x)
        allocate(this%y(m),source=y)
        allocate(this%z(m),source=z)

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w)
        else
           allocate(this%w(m),source=one)
        endif

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
        this%liwrk = m+product(nest-2*order-1)
        allocate(this%iwrk(this%liwrk),source=0)

        ! wrk1
        uv  = nest-order-1
        km  = maxval(order,1)+1
        bxy(1) = order(1)*uv(2)+order(2)+1
        bxy(2) = order(2)*uv(1)+order(1)+1
        if (bxy(1)<=bxy(2)) then
            b1 = bxy(1); b2 = b1 + uv(2) - order(2)
        else
            b1 = bxy(2); b2 = b1 + uv(1) - order(1)
        end if
        this%lwrk1 = product(uv)*(2+b1+b2)+2*(sum(uv)+km*(m+nmax)+nmax-sum(order))+b2+1
        allocate(this%wrk1(this%lwrk1),source=zero)

        ! wrk2
        this%lwrk2 = product(uv)*(b2+1)+b2
        allocate(this%wrk2(this%lwrk2),source=zero)

        endassociate

    end subroutine surf_new_points

    ! A default constructor
    type(fitpack_surface) function surf_new_from_points(x,y,z,w,ierr) result(this)
        real(FP_REAL), intent(in) :: x(:),y(size(x)),z(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ierr0

        ierr0 = this%new_fit(x,y,z,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new surface fit')

    end function surf_new_from_points

    ! Fit a new curve
    integer(FP_FLAG) function surf_new_fit(this,x,y,z,w,smoothing,order)
        class(fitpack_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x)),z(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights
        real(FP_REAL), optional, intent(in) :: smoothing
        integer    , optional, intent(in) :: order

        call this%new_points(x,y,z,w)

        surf_new_fit = this%fit(smoothing,order)

    end function surf_new_fit

    !> Evaluate derivatives on a grid domain
    function surface_derivatives_gridded(this,x,y,dx,dy,ierr) result(f)
        class(fitpack_surface), intent(inout)  :: this
        
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
        if (min_lwrk>this%lwrk1) then 
            allocate(min_wrk(min_lwrk),source=0.0_FP_REAL)
            call move_alloc(from=min_wrk,to=this%wrk1)
            this%lwrk1 = min_lwrk
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
                    wrk=this%wrk1,lwrk=this%lwrk1,     &    ! memory
                    iwrk=this%iwrk,kwrk=this%liwrk,    &    ! memory
                    ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate gridded derivatives')

    end function surface_derivatives_gridded

    !> Evaluate derivatives on a list of (x(i),y(i)) points
    function surface_derivatives_many(this,x,y,dx,dy,ierr) result(f)
        class(fitpack_surface), intent(inout)  :: this
        
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
            if (min_lwrk>this%lwrk1) then 
                allocate(min_wrk(min_lwrk),source=0.0_FP_REAL)
                call move_alloc(from=min_wrk,to=this%wrk1)
                this%lwrk1 = min_lwrk
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
                        wrk=this%wrk1,lwrk=this%lwrk1,     &    ! memory
                        iwrk=this%iwrk,kwrk=this%liwrk,    &    ! memory
                        ier=ier)

            
        end if
        
        call fitpack_error_handling(ier,ierr,'evaluate bivariate derivatives')

    end function surface_derivatives_many

    real(FP_REAL) function surface_derivatives_one(this,x,y,dx,dy,ierr) result(f)
        class(fitpack_surface), intent(inout) :: this
        
        ! Evaluation point
        real(FP_REAL),          intent(in)      :: x,y    

        ! Order of the partial derivatives w.r.t. x and y
        integer(FP_SIZE), intent(in) :: dx,dy              
        
        ! Optional error flag
        integer(FP_FLAG), optional, intent(out) :: ierr   

        real(FP_REAL) :: z1(1),x1(1),y1(1)
        
        x1 = x
        y1 = y

        z1 = surface_derivatives_many(this,x1,y1,dx,dy,ierr)
        f  = z1(1)

    end function surface_derivatives_one


    ! =================================================================================================
    ! PARALLEL COMMUNICATION
    ! =================================================================================================

    elemental integer(FP_SIZE) function surf_comm_size(this)
        class(fitpack_surface), intent(in) :: this
        ! Base fields + surface-specific scalars:
        ! m, order(2), left(2), right(2), nest(2), nmax, lwrk1, lwrk2, bc, knots(2) = 15
        surf_comm_size = this%core_comm_size() &
                       + 15 &
                       + FP_COMM_SIZE(this%x) &
                       + FP_COMM_SIZE(this%y) &
                       + FP_COMM_SIZE(this%z) &
                       + FP_COMM_SIZE(this%w) &
                       + FP_COMM_SIZE(this%t) &
                       + FP_COMM_SIZE(this%wrk1) &
                       + FP_COMM_SIZE(this%wrk2)
    end function surf_comm_size

    pure subroutine surf_comm_pack(this, buffer)
        class(fitpack_surface), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_pack(buffer)
        pos = this%core_comm_size() + 1

        buffer(pos) = real(this%m, FP_COMM);         pos = pos + 1
        buffer(pos) = real(this%order(1), FP_COMM);  pos = pos + 1
        buffer(pos) = real(this%order(2), FP_COMM);  pos = pos + 1
        buffer(pos) = this%left(1);                   pos = pos + 1
        buffer(pos) = this%left(2);                   pos = pos + 1
        buffer(pos) = this%right(1);                  pos = pos + 1
        buffer(pos) = this%right(2);                  pos = pos + 1
        buffer(pos) = real(this%nest(1), FP_COMM);   pos = pos + 1
        buffer(pos) = real(this%nest(2), FP_COMM);   pos = pos + 1
        buffer(pos) = real(this%nmax, FP_COMM);      pos = pos + 1
        buffer(pos) = real(this%lwrk1, FP_COMM);     pos = pos + 1
        buffer(pos) = real(this%lwrk2, FP_COMM);     pos = pos + 1
        buffer(pos) = real(this%bc, FP_COMM);        pos = pos + 1
        buffer(pos) = real(this%knots(1), FP_COMM);  pos = pos + 1
        buffer(pos) = real(this%knots(2), FP_COMM);  pos = pos + 1

        call FP_COMM_PACK(this%x, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_PACK(this%y, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_PACK(this%z, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_PACK(this%w, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_PACK(this%t, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_PACK(this%wrk1, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk1)
        call FP_COMM_PACK(this%wrk2, buffer(pos:))
    end subroutine surf_comm_pack

    pure subroutine surf_comm_expand(this, buffer)
        class(fitpack_surface), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_expand(buffer)
        pos = this%core_comm_size() + 1

        this%m        = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%order(1) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%order(2) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%left(1)  = buffer(pos);                  pos = pos + 1
        this%left(2)  = buffer(pos);                  pos = pos + 1
        this%right(1) = buffer(pos);                  pos = pos + 1
        this%right(2) = buffer(pos);                  pos = pos + 1
        this%nest(1)  = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nest(2)  = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nmax     = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%lwrk1    = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%lwrk2    = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%bc       = nint(buffer(pos), FP_FLAG);  pos = pos + 1
        this%knots(1) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(2) = nint(buffer(pos), FP_SIZE);  pos = pos + 1

        call FP_COMM_EXPAND(this%x, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%x)
        call FP_COMM_EXPAND(this%y, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%y)
        call FP_COMM_EXPAND(this%z, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_EXPAND(this%w, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_EXPAND(this%t, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_EXPAND(this%wrk1, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk1)
        call FP_COMM_EXPAND(this%wrk2, buffer(pos:))
    end subroutine surf_comm_expand

end module fitpack_surfaces
