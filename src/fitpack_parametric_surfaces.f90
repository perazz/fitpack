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
    use fitpack_fitters
    implicit none
    private

    public :: fitpack_parametric_surface

    !> A public type describing a bicubic parametric surface fitter defined by points z(j,i,:) in the
    !> idim-dimensional space, organized on a grid of strictly increasing parameter values u(i), v(j)
    type, extends(fitpack_fitter) :: fitpack_parametric_surface

        !> Number of dimensions
        integer :: idim = 0

        !> Parameter grid values
        real(FP_REAL), allocatable :: u(:),v(:)

        !> idim-dimensional function values [size(v),size(u),idim]
        real(FP_REAL), allocatable :: z(:,:,:)

        !> Flags to determine whether dimensions u,v are periodic
        logical :: periodic_dim(2) = .false.

        !> Estimated and actual number of knots and their allocations
        integer :: nest(2)  = 0
        integer :: nmax     = 0
        ! (lwrk/wrk inherited from fitpack_fitter)

!        ! Space for derivative evaluation
!        real(FP_REAL), allocatable :: dd(:,:)

        ! Knots
        integer     :: knots(2) = 0
        real(FP_REAL), allocatable :: t(:,:)  ! Knot locations (:,1)=u; (:,2)=v

        contains

           !> Clean memory
           procedure :: destroy

           !> Set new points
           procedure :: new_points

           !> Generate new fit
           procedure :: new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit           => surf_fit_automatic_knots
           procedure :: interpolate   => interpolating_curve
           procedure :: least_squares => surface_fit_least_squares

           !> Evaluate curve at given coordinates
           procedure, private :: surf_eval_one
           procedure, private :: surf_eval_grid
           generic :: eval => surf_eval_one,surf_eval_grid

           !> Parallel communication interface
           procedure :: comm_size   => parsurf_comm_size
           procedure :: comm_pack   => parsurf_comm_pack
           procedure :: comm_expand => parsurf_comm_expand

    end type fitpack_parametric_surface

    ! Default constructor
    interface fitpack_parametric_surface
       module procedure new_from_points
    end interface fitpack_parametric_surface

    contains

    ! A default constructor
    type(fitpack_parametric_surface) function new_from_points(u,v,z,periodic_BC,ierr) result(this)
        real(FP_REAL), intent(in) :: u(:),v(:),z(:,:,:)
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
        real(FP_REAL), intent(in) :: u(:),v(:),z(:,:,:)
        real(FP_REAL), optional, intent(in) :: smoothing
        logical    , optional, intent(in) :: periodic_BC(2)

        call this%new_points(u,v,z,periodic_BC)

        new_fit = this%fit(smoothing)

    end function new_fit

    elemental subroutine destroy(this)
       class(fitpack_parametric_surface), intent(inout) :: this
       integer :: ierr
       call this%destroy_base()
       this%idim = 0
       this%periodic_dim = .false.
       deallocate(this%u,stat=ierr)
       deallocate(this%v,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%t,stat=ierr)

       this%nest      = 0
       this%knots     = 0

    end subroutine destroy

    subroutine new_points(this,u,v,z,periodic_BC)
        class(fitpack_parametric_surface), intent(inout) :: this
        real(FP_REAL), intent(in) :: u(:),v(:),z(:,:,:)
        logical, optional, intent(in) :: periodic_BC(2)

        integer, parameter   :: SAFE = 2
        integer :: m(2)

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
        allocate(this%u(m(1)),source=u)
        allocate(this%v(m(2)),source=v)
        allocate(this%z(size(z,1),size(z,2),size(z,3)),source=z)

        ! Reset run flag
        this%iopt = IOPT_NEW_SMOOTHING

        ! Knot space: overestimate (mv+4+2*ipar(2))
        nest = SAFE*(m+4+2)
        call allocate_knot_storage(this,nest)

        this%fp = zero

        endassociate

    end subroutine new_points

    pure subroutine allocate_knot_storage(this,knots)
       class(fitpack_parametric_surface), intent(inout) :: this
       integer, intent(in) :: knots(2)

       integer :: clen,q,m(2)
       real(FP_REAL), allocatable :: t(:,:),c(:),wrk(:)
       integer, allocatable :: iwrk(:)

       associate(nmax=>this%nmax,idim=>this%idim)

       m = [size(this%u),size(this%v)]

       nmax = maxval(knots)
       allocate(t(nmax,2),source=zero)
       call move_alloc(from=t,to=this%t)

       ! Spline coefficients
       clen = product(knots-4)*idim
       allocate(c(clen),source=zero)
       call move_alloc(from=c,to=this%c)

       ! Working space
       this%liwrk = 3+sum(m+knots)
       allocate(iwrk(this%liwrk),source=0)
       call move_alloc(from=iwrk,to=this%iwrk)

       ! wrk
       ! lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+4*(mu+mv)+q*idim
       ! where q is the larger of mv and nuest.
       q = max(m(2),knots(1))
       this%lwrk = 4+knots(1)*(m(2)*idim+11+4)+knots(2)*(11+4)+4*sum(m)+q*idim
       allocate(wrk(this%lwrk),source=zero)
       call move_alloc(from=wrk,to=this%wrk)

       endassociate

    end subroutine allocate_knot_storage

    function surf_eval_one(this,u,v,ierr) result(y)
        class(fitpack_parametric_surface), intent(inout)  :: this
        real(FP_REAL),          intent(in)     :: u,v      ! Evaluation point
        integer, optional,    intent(out)    :: ierr     ! Optional error flag
        real(FP_REAL) :: y(this%idim)

        real(FP_REAL) :: y1(1,1,this%idim)

        y1 = surf_eval_grid(this,[u],[v],ierr)
        y  = y1(1,1,:)

    end function surf_eval_one

    ! subroutine surev evaluate parameter values on a grid (u(i),v(j)),i=1,...,mu; j=1,...,mv a bicubic spline
    ! surface of dimension idim, given in the b-spline representation.

    ! on successful exit f(mu*mv*(l-1)+mv*(i-1)+j) contains the l-th co-ordinate of the bicubic
    ! spline surface at the point (u(i),v(j)),l=1,...,idim,i=1,...,mu;j=1,...,mv.

    ! f(1:mv,1:mu,1:idim)


    function surf_eval_grid(this,u,v,ierr) result(f)
        class(fitpack_parametric_surface), intent(inout)  :: this
        real(FP_REAL),          intent(in)     :: u(:),v(:) ! Evaluation grid (parameter range)
        integer, optional,    intent(out)    :: ierr      ! Optional error flag
        real(FP_REAL) :: f(size(v),size(u),this%idim)

        integer :: ier

        call surev(idim=this%idim,                  &  ! dimension of the spline surface
                   tu=this%t(:,1),nu=this%knots(1), &  ! knots in the u-direction
                   tv=this%t(:,2),nv=this%knots(2), &  ! knots in the v-direction
                   c=this%c,                        &  ! the b-spline coefficients
                   u=u,mu=size(u),                  &  ! u co-ordinates of the grid points along the u-axis.
                   v=v,mv=size(v),                  &  ! v co-ordinates of the grid points along the v-axis.
                   f=f,mf=size(f),                  &  ! Array of the results
                   wrk=this%wrk,lwrk=this%lwrk,     &  ! workspace
                   iwrk=this%iwrk,kwrk=this%liwrk,  &  ! workspace
                   ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate parametric surface')

    end function surf_eval_grid

    ! Interpolating curve
    integer function interpolating_curve(this,reset_knots) result(ierr)
        class(fitpack_parametric_surface), intent(inout) :: this
        logical, optional, intent(in) :: reset_knots

        logical :: do_reset

        do_reset = .true.; if (present(reset_knots)) do_reset = reset_knots
        if (do_reset) this%iopt = IOPT_NEW_SMOOTHING
        ierr = surf_fit_automatic_knots(this,zero,keep_knots=.not.do_reset)

    end function interpolating_curve

    ! Least-squares surface on current or given knots
    integer function surface_fit_least_squares(this,u_knots,v_knots,smoothing,reset_knots) result(ierr)
       class(fitpack_parametric_surface), intent(inout) :: this
       real(FP_REAL), optional, intent(in) :: u_knots(:)
       real(FP_REAL), optional, intent(in) :: v_knots(:)
       real(FP_REAL), optional, intent(in) :: smoothing
       logical, optional, intent(in) :: reset_knots

       integer :: nuuser,nu,nvuser,nv,new_knots(2)
       logical :: do_reset

       ! Optionally recompute knots via a smoothing fit first
       do_reset = .false.; if (present(reset_knots)) do_reset = reset_knots
       if (do_reset) then
           this%iopt = IOPT_NEW_SMOOTHING
           ierr = this%fit(smoothing)
           if (.not.FITPACK_SUCCESS(ierr)) return
       end if

       this%iopt = IOPT_NEW_LEASTSQUARES

       if (present(u_knots)) then
           nuuser = size(u_knots)
           nu     = nuuser+8
       else
           nuuser = 0
           nu     = this%knots(1)
       endif

       if (present(v_knots)) then
           nvuser = size(v_knots)
           nv     = nvuser+8
       else
           nvuser = 0
           nv     = this%knots(2)
       endif

       ! Reallocate knots if necessary
       new_knots = [nu,nv]
       if (any(new_knots>this%nmax)) call allocate_knot_storage(this,new_knots)
       this%knots = new_knots

       ! Copy knots
       if (present(u_knots)) this%t(4+1:4+nuuser,1) = u_knots
       if (present(v_knots)) this%t(4+1:4+nvuser,2) = v_knots

       ! Fit
       this%iopt = IOPT_NEW_LEASTSQUARES
       ierr = this%fit()

    end function surface_fit_least_squares

    ! Curve fitting driver: automatic number of knots
    integer function surf_fit_automatic_knots(this,smoothing,periodic,keep_knots) result(ierr)
        class(fitpack_parametric_surface), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        logical, optional, intent(in) :: periodic(2)
        logical, optional, intent(in) :: keep_knots

        integer :: loop,nit,ipar(2)
        real(FP_REAL) :: smooth_now(3)
        logical :: do_guard

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots (unless caller wants to keep them)
        do_guard = .true.; if (present(keep_knots)) do_guard = .not.keep_knots
        if (do_guard .and. this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        ! Optionally set periodicity
        if (present(periodic)) this%periodic_dim = periodic

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

    end function surf_fit_automatic_knots

    ! =================================================================================================
    ! PARALLEL COMMUNICATION
    ! =================================================================================================

    elemental integer(FP_SIZE) function parsurf_comm_size(this)
        class(fitpack_parametric_surface), intent(in) :: this
        ! Base fields + parametric-surface-specific scalars:
        ! idim, periodic_dim(2), nest(2), nmax, knots(2) = 8
        parsurf_comm_size = this%core_comm_size() &
                          + 8 &
                          + FP_COMM_SIZE(this%u) &
                          + FP_COMM_SIZE(this%v) &
                          + FP_COMM_SIZE(this%z) &
                          + FP_COMM_SIZE(this%t)
    end function parsurf_comm_size

    pure subroutine parsurf_comm_pack(this, buffer)
        class(fitpack_parametric_surface), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_pack(buffer)
        pos = this%core_comm_size() + 1

        buffer(pos) = real(this%idim, FP_COMM);                              pos = pos + 1
        buffer(pos) = real(merge(1, 0, this%periodic_dim(1)), FP_COMM);      pos = pos + 1
        buffer(pos) = real(merge(1, 0, this%periodic_dim(2)), FP_COMM);      pos = pos + 1
        buffer(pos) = real(this%nest(1), FP_COMM);                           pos = pos + 1
        buffer(pos) = real(this%nest(2), FP_COMM);                           pos = pos + 1
        buffer(pos) = real(this%nmax, FP_COMM);                              pos = pos + 1
        buffer(pos) = real(this%knots(1), FP_COMM);                          pos = pos + 1
        buffer(pos) = real(this%knots(2), FP_COMM);                          pos = pos + 1

        call FP_COMM_PACK(this%u, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%u)
        call FP_COMM_PACK(this%v, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%v)
        call FP_COMM_PACK(this%z, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_PACK(this%t, buffer(pos:))
    end subroutine parsurf_comm_pack

    pure subroutine parsurf_comm_expand(this, buffer)
        class(fitpack_parametric_surface), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_expand(buffer)
        pos = this%core_comm_size() + 1

        this%idim            = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%periodic_dim(1) = nint(buffer(pos)) /= 0;     pos = pos + 1
        this%periodic_dim(2) = nint(buffer(pos)) /= 0;     pos = pos + 1
        this%nest(1)         = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nest(2)         = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nmax            = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(1)        = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(2)        = nint(buffer(pos), FP_SIZE);  pos = pos + 1

        call FP_COMM_EXPAND(this%u, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%u)
        call FP_COMM_EXPAND(this%v, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%v)
        call FP_COMM_EXPAND(this%z, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%z)
        call FP_COMM_EXPAND(this%t, buffer(pos:))
    end subroutine parsurf_comm_expand

end module fitpack_parametric_surfaces
