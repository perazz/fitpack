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
!> @brief OOP wrapper for constrained-convexity curve fitting.
!!
!! Provides fitpack_convex_curve, an extension of fitpack_curve that enforces local
!! convexity or concavity constraints during spline fitting. Always uses cubic splines
!! (\f$ k = 3 \f$). The underlying core routines are concon (automatic knot placement)
!! and cocosp (least-squares with user-supplied knots).
!!
!! @see Dierckx, Ch. 8, §8.3–8.4 (pp. 173–196); concon, cocosp
module fitpack_convex_curves
    use fitpack_core
    use fitpack_fitters
    use fitpack_curves
    implicit none
    private

    public :: fitpack_convex_curve

    !> @brief Cubic spline fitter with pointwise convexity/concavity constraints.
    !!
    !! Extends fitpack_curve with per-data-point convexity flags: each point can be
    !! constrained to lie on a locally concave (\f$ v_i = 1 \f$), convex
    !! (\f$ v_i = -1 \f$), or unconstrained (\f$ v_i = 0 \f$) portion of the spline.
    !! The fitting routines solve a constrained quadratic programming problem to determine
    !! knot positions and coefficients that satisfy these shape constraints while minimizing
    !! the smoothing functional.
    !!
    !! @see Dierckx, Ch. 8, §8.3–8.4 (pp. 173–196)
    type, extends(fitpack_curve) :: fitpack_convex_curve

        !> Convexity constraints at data points: 1=concave, -1=convex, 0=unconstrained
        real(FP_REAL), allocatable    :: v(:)

        !> Spline values at data points after fit
        real(FP_REAL), allocatable    :: sx(:)

        !> Active constraint flags at knots
        logical(FP_BOOL), allocatable :: bind(:)

        !> Tree storage estimate for QP solver
        integer(FP_SIZE) :: maxtr  = 100

        !> Max number of zero-curvature knots
        integer(FP_SIZE) :: maxbin = 10

    contains

        !> Clean memory
        procedure :: destroy       => convex_destroy

        !> Set new points
        procedure :: new_points    => convex_new_points

        !> Generate/update fitting curve with automatic knots (concon)
        procedure :: fit           => convex_fit_automatic_knots

        !> Least-squares fit with given knots (cocosp)
        procedure :: least_squares => convex_fit_least_squares

        !> Set convexity constraint values
        procedure :: set_convexity

        !> Parallel communication interface
        procedure :: comm_size     => convex_comm_size
        procedure :: comm_pack     => convex_comm_pack
        procedure :: comm_expand   => convex_comm_expand

    end type fitpack_convex_curve

    ! Default constructor
    interface fitpack_convex_curve
       module procedure convex_new_from_points
    end interface fitpack_convex_curve

    contains

    !> @brief Construct a convex curve from data points with convexity constraints and fit it.
    !!
    !! @param[in]  x          Independent variable values.
    !! @param[in]  y          Dependent variable values, same size as `x`.
    !! @param[in]  v          Convexity constraints: 1=concave, -1=convex, 0=unconstrained.
    !! @param[in]  w          Optional weights (positive), same size as `x`.
    !! @param[in]  smoothing  Optional smoothing factor.
    !! @param[out] ierr       Optional error flag.
    !! @return Fitted convex curve object.
    !!
    !! @see concon
    type(fitpack_convex_curve) function convex_new_from_points(x,y,v,w,smoothing,ierr) result(this)
        real(FP_REAL), intent(in) :: x(:),y(size(x))
        real(FP_REAL), intent(in) :: v(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x))
        real(FP_REAL), optional, intent(in) :: smoothing
        integer(FP_FLAG), optional, intent(out) :: ierr

        integer(FP_FLAG) :: ierr0

        call this%new_points(x,y,w)
        ierr0 = this%set_convexity(v)
        if (FITPACK_SUCCESS(ierr0)) ierr0 = this%fit(smoothing)

        call fitpack_error_handling(ierr0,ierr,'new convex curve fit')

    end function convex_new_from_points

    !> @brief Destroy the convex curve object and release all allocated memory.
    elemental subroutine convex_destroy(this)
       class(fitpack_convex_curve), intent(inout) :: this
       integer :: ierr

       ! Destroy parent
       call this%fitpack_curve%destroy()

       ! Destroy convex-specific fields
       deallocate(this%v,stat=ierr)
       deallocate(this%sx,stat=ierr)
       deallocate(this%bind,stat=ierr)
       this%maxtr  = 100
       this%maxbin = 10

    end subroutine convex_destroy

    !> @brief Load new data points and allocate workspace for convex fitting.
    !!
    !! @param[in,out] this  The convex curve object (destroyed and reinitialized).
    !! @param[in]     x     Independent variable values.
    !! @param[in]     y     Dependent variable values, same size as `x`.
    !! @param[in]     w     Optional positive weights, same size as `x`.
    !!
    !! @see concon, cocosp
    subroutine convex_new_points(this,x,y,w)
        class(fitpack_convex_curve), intent(inout) :: this
        real(FP_REAL), intent(in) :: x(:),y(size(x))
        real(FP_REAL), optional, intent(in) :: w(size(x)) ! node weights

        integer, allocatable :: isort(:)

        associate(m=>this%m,nest=>this%nest,lwrk=>this%lwrk,maxbin=>this%maxbin,maxtr=>this%maxtr)

        call this%destroy()

        m = size(x)

        ! Ensure x are sorted
        isort = fitpack_argsort(x)
        allocate(this%x(m),source=x(isort))
        allocate(this%y(m),source=y(isort))

        ! Set up weights
        if (present(w)) then
           allocate(this%w(m),source=w(isort))
        else
           allocate(this%w(m),source=1.0_FP_REAL)
        endif
        allocate(this%sp(m),source=0.0_FP_REAL)

        ! Setup boundaries
        this%xleft  = minval(x)
        this%xright = maxval(x)

        ! Force cubic (concon/cocosp always cubic)
        this%order = 3

        ! Reset run flag
        this%iopt = 0

        ! Estimated knots: m+4 is always sufficient for concon
        nest = max(8, m + 4)
        this%liwrk = maxtr*4 + 2*(maxbin+1)
        allocate(this%iwrk(this%liwrk))
        allocate(this%t(nest),this%c(nest))

        ! concon workspace: m*4 + nest*8 + maxbin*(maxbin+nest+1)
        lwrk = m*4 + nest*8 + maxbin*(maxbin+nest+1)
        allocate(this%wrk(lwrk),source=0.0_FP_REAL)
        allocate(this%wrk_fou(nest,2))

        ! Allocate convex-specific arrays
        allocate(this%v(m),source=zero)
        allocate(this%sx(m),source=zero)
        allocate(this%bind(nest),source=FP_FALSE)

        endassociate

    end subroutine convex_new_points

    !> @brief Set pointwise convexity constraints.
    !!
    !! @param[in,out] this  The convex curve object.
    !! @param[in]     v     Constraint values: 1=concave, -1=convex, 0=unconstrained. Must have size `m`.
    !! @return Error flag (FITPACK_INPUT_ERROR if size mismatch).
    integer(FP_FLAG) function set_convexity(this, v) result(ierr)
        class(fitpack_convex_curve), intent(inout) :: this
        real(FP_REAL), intent(in) :: v(:)

        if (size(v) /= this%m) then
            ierr = FITPACK_INPUT_ERROR
            return
        end if

        this%v = v
        ierr = FITPACK_OK

    end function set_convexity

    ! Remap concon/cocosp error codes 1-3 to distinct FITPACK constants
    elemental integer(FP_FLAG) function remap_concon_error(ier) result(remapped)
        integer(FP_FLAG), intent(in) :: ier

        select case (ier)
            case (1);  remapped = CONCON_MAXBIN
            case (2);  remapped = CONCON_MAXTR
            case (3);  remapped = CONCON_QP_FAIL
            case default; remapped = ier
        end select

    end function remap_concon_error

    !> @brief Fit a smoothing spline with convexity constraints and automatic knot placement.
    !!
    !! Uses the concon core routine which solves a constrained quadratic programming problem
    !! to place knots while respecting local convexity/concavity constraints.
    !!
    !! @param[in,out] this       The convex curve object.
    !! @param[in]     smoothing  Optional smoothing factor.
    !! @param[in]     order      Ignored (always cubic).
    !! @param[in]     keep_knots If `.true.`, continue from current knot set.
    !! @return Error flag.
    !!
    !! @see concon
    integer(FP_FLAG) function convex_fit_automatic_knots(this,smoothing,order,keep_knots) result(ierr)
        class(fitpack_convex_curve), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        integer(FP_SIZE), optional, intent(in) :: order
        logical, optional, intent(in) :: keep_knots

        logical :: do_keep

        ! Set smoothing if provided
        if (present(smoothing)) this%smoothing = smoothing

        ! Map keep_knots -> iopt (0=new, 1=continue)
        do_keep = .false.; if (present(keep_knots)) do_keep = keep_knots
        if (.not.do_keep) this%iopt = IOPT_NEW_SMOOTHING

        ! Ignore order argument (always cubic)

        call concon(this%iopt, this%m, this%x, this%y, this%w, this%v, &
                    this%smoothing, this%nest, this%maxtr, this%maxbin, &
                    this%knots, this%t, this%c, this%fp, this%sx, this%bind, &
                    this%wrk, this%lwrk, this%iwrk, this%liwrk, ierr)

        ! Remap error codes
        ierr = remap_concon_error(ierr)

        ! On success, set iopt to continue mode
        if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

    end function convex_fit_automatic_knots

    !> @brief Fit a least-squares spline with convexity constraints on the current knots.
    !!
    !! Uses the cocosp core routine. Constraint values at interior knots are derived from
    !! the nearest data-point convexity flags.
    !!
    !! @param[in,out] this        The convex curve object.
    !! @param[in]     smoothing   Optional smoothing for initial knot placement (if `reset_knots`).
    !! @param[in]     reset_knots If `.true.`, recompute knots via concon first.
    !! @return Error flag.
    !!
    !! @see cocosp
    integer(FP_FLAG) function convex_fit_least_squares(this,smoothing,reset_knots) result(ierr)
        class(fitpack_convex_curve), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        logical, optional, intent(in) :: reset_knots

        logical :: do_reset
        integer(FP_SIZE) :: j, n6
        real(FP_REAL) :: e(this%nest)

        ! Optionally recompute knots via concon first
        do_reset = .false.; if (present(reset_knots)) do_reset = reset_knots
        if (do_reset) then
            ierr = this%fit(smoothing)
            if (.not.FITPACK_SUCCESS(ierr)) return
        end if

        ! Derive constraint values e(j) at interior knots from v(i) at data points.
        ! For each interior knot t(j+3), find the nearest data point and use its v value.
        n6 = this%knots - 6
        e = zero
        do j = 1, n6
            e(j) = nearest_constraint(this%t(j+3), this%x, this%v, this%m)
        end do

        ! cocosp workspace (m*4+n*7+maxbin*(maxbin+n+1)) fits within concon's
        ! (m*4+nest*8+maxbin*(maxbin+nest+1)) since n<=nest and 7<8.
        ! Reuse this%wrk, this%iwrk, this%bind.
        this%bind = FP_FALSE

        call cocosp(this%m, this%x, this%y, this%w, &
                    this%knots, this%t, e, this%maxtr, this%maxbin, &
                    this%c, this%fp, this%sx, this%bind, &
                    this%wrk, this%lwrk, this%iwrk, this%liwrk, ierr)

        ! Remap error codes
        ierr = remap_concon_error(ierr)

    end function convex_fit_least_squares

    ! Find the constraint value for the nearest data point to a knot location
    pure real(FP_REAL) function nearest_constraint(t_knot, x, v, m) result(vc)
        real(FP_REAL), intent(in) :: t_knot
        integer(FP_SIZE), intent(in) :: m
        real(FP_REAL), intent(in) :: x(m), v(m)

        integer(FP_SIZE) :: i, imin
        real(FP_REAL) :: dmin, d

        dmin = huge(one)
        imin = 1
        do i = 1, m
            d = abs(x(i) - t_knot)
            if (d < dmin) then
                dmin = d
                imin = i
            end if
        end do
        vc = v(imin)

    end function nearest_constraint

    ! =================================================================================================
    ! PARALLEL COMMUNICATION (size/pack/expand)
    ! =================================================================================================

    !> @brief Return the communication buffer size for the convex curve.
    elemental integer(FP_SIZE) function convex_comm_size(this)
        class(fitpack_convex_curve), intent(in) :: this

        ! Parent curve size + convex-specific scalars (maxtr, maxbin = 2)
        ! + convex-specific arrays (v, sx, bind)
        convex_comm_size = this%fitpack_curve%comm_size() &
                         + 2 &
                         + FP_COMM_SIZE(this%v) &
                         + FP_COMM_SIZE(this%sx) &
                         + FP_COMM_SIZE(this%bind)

    end function convex_comm_size

    !> @brief Pack convex curve data into a communication buffer.
    pure subroutine convex_comm_pack(this, buffer)
        class(fitpack_convex_curve), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)

        integer(FP_SIZE) :: pos

        ! Pack parent fields first
        call this%fitpack_curve%comm_pack(buffer)
        pos = this%fitpack_curve%comm_size() + 1

        ! Pack convex-specific scalars
        buffer(pos)   = real(this%maxtr, FP_COMM);   pos = pos + 1
        buffer(pos)   = real(this%maxbin, FP_COMM);  pos = pos + 1

        ! Pack convex-specific arrays
        call FP_COMM_PACK(this%v, buffer(pos:));      pos = pos + FP_COMM_SIZE(this%v)
        call FP_COMM_PACK(this%sx, buffer(pos:));     pos = pos + FP_COMM_SIZE(this%sx)
        call FP_COMM_PACK(this%bind, buffer(pos:))

    end subroutine convex_comm_pack

    !> @brief Expand convex curve data from a communication buffer.
    pure subroutine convex_comm_expand(this, buffer)
        class(fitpack_convex_curve), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)

        integer(FP_SIZE) :: pos

        ! Expand parent fields first
        call this%fitpack_curve%comm_expand(buffer)
        pos = this%fitpack_curve%comm_size() + 1

        ! Expand convex-specific scalars
        this%maxtr  = nint(buffer(pos), FP_SIZE);   pos = pos + 1
        this%maxbin = nint(buffer(pos), FP_SIZE);   pos = pos + 1

        ! Expand convex-specific arrays
        call FP_COMM_EXPAND(this%v, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%v)
        call FP_COMM_EXPAND(this%sx, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%sx)
        call FP_COMM_EXPAND(this%bind, buffer(pos:))

    end subroutine convex_comm_expand

end module fitpack_convex_curves
