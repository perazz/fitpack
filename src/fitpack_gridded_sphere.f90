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
module fitpack_gridded_sphere
    use fitpack_core
    implicit none
    private

    public :: fitpack_grid_sphere

    !> A public type describing a sphere fitter z = s(u,v) to GRIDDED sphere data,
    !> on the latitude-longitude grid (u(i),v(j)), i=1,...,mu ; j=1,...,mv , spgrid determines a smooth
    !> u = latitude, 0<=u<=pi,
    !> v = 2*pi-periodic longitude, vb<=v<=ve (vb = v(1), ve=vb+2*pi).
    type :: fitpack_grid_sphere

        !> Coordinates of the data points in grid coordinates (u,v) and domain size
        real(RKIND), allocatable :: u(:),v(:)

        !> Function values at the data points z(iv,iu) note (v,u) indices
        real(RKIND), allocatable :: z(:,:)

        !> North pole BC (u=0, index=1), South pole BC (u=pi, index=2)
        real(RKIND) :: pole_z0(2)      = zero      ! Function value at north pole
        logical     :: pole_present(2) = .false.   ! Is pole data provided?
        logical     :: pole_exct(2)    = .false.   ! Should that be treated as exact
        integer     :: pole_continuity(2) = 1      ! Continuity at pole BC (C0 or C1)
        logical     :: pole_zero_grad(2) = .true.  ! Should gradients vanish

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
        integer :: iopt = IOPT_NEW_SMOOTHING ! -> iopt(1)



        contains

           !> Clean memory
           procedure :: destroy       => spgrid_destroy

           !> Set new points
           procedure :: new_points    => spgrid_new_points

           !> Generate new fit
           procedure :: new_fit       => spgrid_new_fit

           !> Set pole BCs
           procedure :: BC_north_pole
           procedure :: BC_south_pole

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit           => spgrid_fit_automatic_knots
           procedure :: least_squares => spgrid_fit_least_squares
           procedure :: interpolate   => spgrid_fit_interpolating

           !> Evaluate gridded domain at given x,y coordinates
           procedure, private :: gridded_eval_one
           procedure, private :: gridded_eval_many
           generic :: eval => gridded_eval_one,gridded_eval_many

           !> Write to disk
           procedure :: write => gridded_to_disk

    end type fitpack_grid_sphere

    interface fitpack_grid_sphere
       module procedure spgrid_new_from_points
    end interface fitpack_grid_sphere

    contains

    ! Fit a surface to least squares of the current knots
    integer function spgrid_fit_least_squares(this) result(ierr)
       class(fitpack_grid_sphere), intent(inout) :: this

       this%iopt = IOPT_NEW_LEASTSQUARES
       ierr = this%fit()

    end function spgrid_fit_least_squares

    ! Find interpolating surface
    integer function spgrid_fit_interpolating(this) result(ierr)
        class(fitpack_grid_sphere), intent(inout) :: this

        ! Set zero smoothing
        ierr = spgrid_fit_automatic_knots(this,smoothing=zero)

    end function spgrid_fit_interpolating


    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer function spgrid_fit_automatic_knots(this,smoothing) result(ierr)
        class(fitpack_grid_sphere), intent(inout) :: this
        real(RKIND), optional, intent(in) :: smoothing

        integer :: loop,nit,iopt(3),ider(4)

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

            ! [Continuation, North pole continuity, South pole continuity]
            iopt = [this%iopt,this%pole_continuity]

            ! [Fit z0 exactly, z0 zero gradient BC]
            ider([1,3]) = merge(merge(1,0,this%pole_exct),-1,this%pole_present)
            ider([2,4]) = merge(1,0,this%pole_zero_grad .and. this%pole_continuity>0)

            call spgrid(iopt,                            &  ! Continuation and BCs
                        ider,                            &  ! Origin point behavior
                        size(this%u),this%u,             &  ! U grid
                        size(this%v),this%v,             &  ! V grid
                        this%z,                          &  ! Sphere radius
                        this%pole_z0(1),this%pole_z0(2), & ! Data value at the poles
                        this%smoothing,                  &  ! Smoothing parameter
                        this%nest(1),this%nest(2),       &  ! Knot space
                        this%knots(1),this%t(:,1),       &  ! u (0:pi) knots (out)
                        this%knots(2),this%t(:,2),       &  ! v (-pi:pi) knots (out)
                        this%c,this%fp,                  &  ! Spline representation and MSE
                        this%wrk,this%lwrk,              &  ! memory
                        this%iwrk,this%liwrk,            &  ! memory
                        ierr)                            ! Error flag

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function spgrid_fit_automatic_knots

    elemental subroutine spgrid_destroy(this)
       class(fitpack_grid_sphere), intent(inout) :: this
       integer :: ierr

       this%pole_z0         = zero
       this%pole_present    = .false.
       this%pole_exct       = .false.
       this%pole_continuity = 1
       this%pole_zero_grad  = .false.

       deallocate(this%u,stat=ierr)
       deallocate(this%v,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this%wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)

       this%smoothing = 1000.0_RKIND
       this%iopt      = 0
       this%nest      = 0
       this%nmax      = 0
       this%lwrk      = 0
       this%liwrk     = 0
       this%knots     = 0
       this%fp        = zero

    end subroutine spgrid_destroy

    subroutine spgrid_new_points(this,u,v,z)
        class(fitpack_grid_sphere), intent(inout) :: this
        real(RKIND), intent(in) :: u(:),v(:) ! sphere domain
        real(RKIND), intent(in) :: z(size(v),size(u)) ! Gridded values

        integer :: clen,m(2),q
        integer, parameter :: SAFE = 2

        associate(nest=>this%nest,nmax=>this%nmax)

        call this%destroy()

        m = [size(u),size(v)] ! /= shape(z), == shape(transpose(z))

        ! Copy grid and data
        allocate(this%u,source=u)
        allocate(this%v,source=v)
        allocate(this%z,source=z)

        ! Reset run flag
        this%iopt = IOPT_NEW_SMOOTHING

        ! Knot space
        nest = SAFE*(m + 8)
        nmax = maxval(nest)
        allocate(this%t(nmax,2),source=zero)

        ! Spline coefficients
        clen = product(nest-4) ! nest-order-1, fixed order==3
        allocate(this%c(clen),source=zero)

        this%fp = zero

        ! Working space
        this%liwrk = 5+sum(m+nest)
        allocate(this%iwrk(this%liwrk),source=0)

        ! wrk
        q = max(m(2)+nest(2),nest(1))
        this%lwrk = 12+nest(1)*(m(2)+nest(2)+3)+nest(2)*24+4*m(1)+8*m(2)+q
        allocate(this%wrk(this%lwrk),source=zero)

        endassociate

    end subroutine spgrid_new_points

    ! A default constructor
    type(fitpack_grid_sphere) function spgrid_new_from_points(u,v,z,ierr) result(this)
        real(RKIND), intent(in) :: u(:),v(:) ! sphere domain
        real(RKIND), intent(in) :: z(size(v),size(u)) ! Gridded values
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(u,v,z)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new gridded surface fit')

    end function spgrid_new_from_points

    ! Fit a new curve
    integer function spgrid_new_fit(this,u,v,z,smoothing)
        class(fitpack_grid_sphere), intent(inout) :: this
        real(RKIND), intent(in) :: u(:),v(:) ! sphere domain
        real(RKIND), intent(in) :: z(size(v),size(u)) ! Gridded values
        real(RKIND), optional, intent(in) :: smoothing

        call this%new_points(u,v,z)

        spgrid_new_fit = this%fit(smoothing)

    end function spgrid_new_fit

    function gridded_eval_many(this,u,v,ierr) result(f)
        class(fitpack_grid_sphere), intent(inout)  :: this
        real(RKIND), intent(in) :: u(:),v(:)  ! Evaluation grid points (polar coordinates)
        real(RKIND) :: f(size(v),size(u))
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
    real(RKIND) function gridded_eval_one(this,u,v,ierr) result(f)
        class(fitpack_grid_sphere), intent(inout)  :: this
        real(RKIND),          intent(in)      :: u,v ! Evaluation point (grid polar coordinates)
        integer, optional,    intent(out)     :: ierr      ! Optional error flag
        real(RKIND) :: f1(1,1)

        f1 = gridded_eval_many(this,[u],[v],ierr)
        f  = f1(1,1)

    end function gridded_eval_one

    !> Set pole BC
    subroutine pole_BC(this,pole,z0,exact,differentiable,zero_grad)
        class(fitpack_grid_sphere), intent(inout) :: this
        integer, intent(in) :: pole
        real(RKIND), optional, intent(in) :: z0 ! Function value at origin
        logical, optional, intent(in) :: exact,differentiable,zero_grad

        this%pole_present(pole) = present(z0)
        if (present(exact)) then
            this%pole_exct(pole) = exact
        else
            this%pole_exct(pole) = .false.
        endif

        if (present(differentiable)) then
            this%pole_continuity(pole) = merge(1,0,differentiable)
        else
            this%pole_continuity(pole) = 0
        end if

        if (present(zero_grad)) then
            this%pole_zero_grad(pole) = zero_grad
        else
            this%pole_zero_grad = .true.
        end if

    end subroutine pole_BC

    !> North pole BC
    subroutine BC_north_pole(this,z0,exact,differentiable,zero_grad)
        class(fitpack_grid_sphere), intent(inout) :: this
        real(RKIND), optional, intent(in) :: z0 ! Function value at origin
        logical, optional, intent(in) :: exact,differentiable,zero_grad
        call pole_BC(this,1,z0,exact,differentiable,zero_grad)
    end subroutine BC_north_pole

    !> South pole BC
    subroutine BC_south_pole(this,z0,exact,differentiable,zero_grad)
        class(fitpack_grid_sphere), intent(inout) :: this
        real(RKIND), optional, intent(in) :: z0 ! Function value at origin
        logical, optional, intent(in) :: exact,differentiable,zero_grad
        call pole_BC(this,1,z0,exact,differentiable,zero_grad)
    end subroutine BC_south_pole

    !> Print gridded polar data to disk
    subroutine gridded_to_disk(this,fileName)
        class(fitpack_grid_sphere), intent(inout) :: this
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
              real(RKIND), intent(in) :: v
              write(numbered,"(1pe12.5)") v
              numbered = adjustr(numbered)
           end function numbered

    end subroutine gridded_to_disk

end module fitpack_gridded_sphere
