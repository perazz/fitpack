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
!> @brief OOP wrapper for bivariate spline fitting on the sphere to scattered data.
!!
!! Provides fitpack_sphere, a derived type for fitting bicubic splines to data scattered
!! over the unit sphere, parameterized by colatitude \f$ \theta \in [0, \pi] \f$ and
!! longitude \f$ \phi \in [0, 2\pi] \f$. The fitted surface
!! \f$ r = s(\theta, \phi) \f$ satisfies appropriate pole constraints to ensure smoothness
!! at \f$ \theta = 0 \f$ and \f$ \theta = \pi \f$. Uses the sphere core routine.
!!
!! @see Dierckx, Ch. 11, §11.2 (pp. 263–269); sphere
module fitpack_sphere_domains
    use fitpack_core
    use fitpack_fitters
    implicit none
    private

    public :: fitpack_sphere

    !> @brief Bicubic spline fitter for data scattered on the sphere.
    !!
    !! Stores the scattered spherical data \f$ (\theta_i, \phi_i, r_i) \f$, optional
    !! weights, and the fitted B-spline representation. The smoothing parameter controls
    !! closeness of fit versus smoothness. Pole constraints ensure regularity at
    !! \f$ \theta = 0 \f$ (north pole) and \f$ \theta = \pi \f$ (south pole).
    !!
    !! @see Dierckx, Ch. 11, §11.2 (pp. 263–269)
    type, extends(fitpack_fitter) :: fitpack_sphere

        !> Scattered data points
        integer :: m = 0
        real(FP_REAL), allocatable :: theta(:),phi(:)

        !> Function values at the points
        real(FP_REAL), allocatable :: r(:)

        ! Node weights
        real(FP_REAL), allocatable :: w(:)

        ! Internal Storage
        integer                  :: lwrk2 = 0
        real(FP_REAL), allocatable :: wrk2(:)

        ! Knots: extimated max number
        integer :: nest(2)  = 0
        integer :: nmax     = 0
        integer :: knots(2) = 0
        real(FP_REAL), allocatable :: t(:,:) ! Knot locations (:,1)=x; (:,2)=y

        contains

           !> Clean memory
           procedure :: destroy    => sphere_destroy

           !> Set new points
           procedure :: new_points => sphere_new_points

           !> Generate new fit
           procedure :: new_fit    => sphere_new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit           => surface_fit_automatic_knots
           procedure :: least_squares => surface_fit_least_squares
           procedure :: interpolate   => surface_fit_interpolating

           !> Evaluate polar domain at given x,y coordinates
           procedure, private :: sphere_eval_one
           procedure, private :: sphere_eval_many
           generic :: eval => sphere_eval_one,sphere_eval_many

           !> Parallel communication interface
           procedure :: comm_size   => sphere_comm_size
           procedure :: comm_pack   => sphere_comm_pack
           procedure :: comm_expand => sphere_comm_expand

    end type fitpack_sphere

    interface fitpack_sphere
       module procedure sphere_new_from_points
    end interface fitpack_sphere

    contains

    !> @brief Fit a least-squares spherical surface with fixed knots.
    !!
    !! @see sphere
    integer function surface_fit_least_squares(this,smoothing,reset_knots) result(ierr)
       class(fitpack_sphere), intent(inout) :: this
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

    !> @brief Fit an interpolating spherical surface (\f$ s = 0 \f$).
    !!
    !! @param[in] reset_knots  If `.true.` (default), start with fresh knot placement.
    !! @return Error flag.
    !!
    !! @see sphere
    integer function surface_fit_interpolating(this,reset_knots) result(ierr)
        class(fitpack_sphere), intent(inout) :: this
        logical, optional, intent(in) :: reset_knots

        logical :: do_reset

        do_reset = .true.; if (present(reset_knots)) do_reset = reset_knots
        if (do_reset) this%iopt = IOPT_NEW_SMOOTHING
        ierr = surface_fit_automatic_knots(this,smoothing=zero,keep_knots=.not.do_reset)

    end function surface_fit_interpolating

    !> @brief Fit a smoothing spherical surface with automatic knot placement.
    !!
    !! Iterates over the smoothing schedule, calling the sphere core routine to fit a
    !! bicubic spline \f$ r = s(\theta, \phi) \f$ to scattered data on the unit sphere.
    !!
    !! @param[in] smoothing   Smoothing factor (\f$ s \ge 0 \f$); default uses stored value.
    !! @param[in] keep_knots  If `.true.`, reuse the current knot set.
    !! @return Error flag.
    !!
    !! @see sphere
    integer function surface_fit_automatic_knots(this,smoothing,keep_knots) result(ierr)
        class(fitpack_sphere), intent(inout) :: this
        real(FP_REAL), optional, intent(in) :: smoothing
        logical, optional, intent(in) :: keep_knots

        integer :: loop,nit
        real(FP_REAL) :: smooth_now(3)
        logical :: do_guard

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        !> Ensure we start with new knots (unless caller wants to keep them)
        do_guard = .true.; if (present(keep_knots)) do_guard = .not.keep_knots
        if (do_guard .and. this%iopt==IOPT_OLD_FIT) this%iopt = IOPT_NEW_SMOOTHING

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! Call fitting routine
            call sphere(this%iopt,                   &  ! Continuation and BCs
                       this%m,this%theta,this%phi,   &  ! scattered points
                       this%r,                       &  ! spline values at the points
                       this%w,                       &  ! Set of weights
                       this%smoothing,               &  ! Smoothing parameter
                       this%nest(1),this%nest(2),    &  ! Knot space
                       epsilon(zero),                &  ! Threshold for the rank of the linear system
                       this%knots(1),this%t(:,1),    &  ! theta  knots (out)
                       this%knots(2),this%t(:,2),    &  ! phi knots (out)
                       this%c,this%fp,               &  ! Spline representation and MSE
                       this%wrk,this%lwrk,         &  ! memory
                       this%wrk2,this%lwrk2,         &  ! memory
                       this%iwrk,this%liwrk,         &  ! memory
                       ierr)                            ! Error flag

            ! If fit was successful, set iopt to "old"
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT

        end do

    end function surface_fit_automatic_knots


    !> @brief Release all allocated memory and reset the sphere fitter to its default state.
    elemental subroutine sphere_destroy(this)
       class(fitpack_sphere), intent(inout) :: this
       integer :: ierr
       call this%destroy_base()
       this%m = 0
       deallocate(this%theta,stat=ierr)
       deallocate(this%phi,stat=ierr)
       deallocate(this%r,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%wrk2,stat=ierr)
       deallocate(this%t,stat=ierr)

       this%nest      = 0
       this%lwrk2     = 0
       this%knots     = 0

    end subroutine sphere_destroy

    !> @brief Load new scattered spherical data and allocate working storage.
    !!
    !! Replaces any previous data with the given colatitude \f$ \theta_i \in [0, \pi] \f$,
    !! longitude \f$ \phi_i \in [0, 2\pi] \f$, and radial values \f$ r_i \f$.
    !!
    !! @param[in] theta  Colatitude coordinates.
    !! @param[in] phi    Longitude coordinates (same size as theta).
    !! @param[in] r      Function values at each point.
    !! @param[in] w      Optional weights (default: uniform).
    subroutine sphere_new_points(this,theta,phi,r,w)
        class(fitpack_sphere), intent(inout) :: this
        real(FP_REAL), intent(in) :: theta(:),phi(size(theta)),r(size(theta))
        real(FP_REAL), optional, intent(in) :: w(size(theta)) ! node weights

        integer :: clen,u,v
        integer, parameter :: SAFE = 2

        associate(m=>this%m,nest=>this%nest,nmax=>this%nmax)

        call this%destroy()

        m = size(theta)

        ! Copy scattered points
        this%theta = theta
        this%phi   = phi
        this%r     = r

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w)
        else
           allocate(this%w(m),source=one)
        endif

        ! Reset run flag
        this%iopt = 0

        ! Knot space: overestimate upper bound
        ! ntest, npest >=8. In most situations, ntest = npest = 8+sqrt(m/2) will be sufficient.
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
        u = nest(1)-7
        v = nest(2)-7
        this%lwrk = 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
        allocate(this%wrk(this%lwrk),source=zero)

        ! wrk2
        this%lwrk2 = 48+21*v+7*u*v+4*(u-1)*v**2
        allocate(this%wrk2(this%lwrk2),source=zero)

        endassociate

    end subroutine sphere_new_points

    !> @brief Construct a fitpack_sphere from scattered data and perform an initial fit.
    !!
    !! @param[in]  theta  Colatitude coordinates \f$ \theta_i \in [0, \pi] \f$.
    !! @param[in]  phi    Longitude coordinates \f$ \phi_i \in [0, 2\pi] \f$.
    !! @param[in]  r      Function values at each point.
    !! @param[in]  w      Optional weights.
    !! @param[out] ierr   Optional error flag; if absent, halts on error.
    !!
    !! @see sphere
    type(fitpack_sphere) function sphere_new_from_points(theta,phi,r,w,ierr) result(this)
        real(FP_REAL), intent(in) :: theta(:),phi(size(theta)),r(size(theta))
        real(FP_REAL), optional, intent(in) :: w(size(theta)) ! node weights
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(theta,phi,r,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new surface fit')

    end function sphere_new_from_points

    !> @brief Evaluate the spherical spline on a grid of colatitude and longitude values.
    !!
    !! Returns \f$ r(j,i) = s(\theta_i, \phi_j) \f$ for all combinations of the input
    !! theta and phi vectors (tensor-product evaluation).
    !!
    !! @param[in]  theta  Colatitude evaluation points.
    !! @param[in]  phi    Longitude evaluation points.
    !! @param[out] ierr   Optional error flag.
    !! @return Grid of spline values, dimensioned `(size(phi), size(theta))`.
    !!
    !! @see bispev
    function sphere_eval_many(this,theta,phi,ierr) result(r)
        class(fitpack_sphere), intent(inout)  :: this
        real(FP_REAL), intent(in) :: theta(:),phi(:)  ! Evaluation points
        real(FP_REAL) :: r(size(phi),size(theta))
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
                    x=theta,mx=size(theta), &
                    y=phi,my=size(phi), &
                    z=r, &
                    wrk=this%wrk2,lwrk=this%lwrk2, &
                    iwrk=this%iwrk,kwrk=this%liwrk,ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate sphere spline')

    end function sphere_eval_many

    !> @brief Evaluate the spherical spline at a single point.
    !!
    !! @param[in]  theta  Colatitude \f$ \theta \in [0, \pi] \f$.
    !! @param[in]  phi    Longitude \f$ \phi \in [0, 2\pi] \f$.
    !! @param[out] ierr   Optional error flag.
    !! @return Spline value \f$ s(\theta, \phi) \f$.
    !!
    !! @see bispev
    real(FP_REAL) function sphere_eval_one(this,theta,phi,ierr) result(r)
        class(fitpack_sphere), intent(inout)  :: this
        real(FP_REAL),          intent(in)     :: theta,phi ! Evaluation point
        integer, optional,    intent(out)    :: ierr      ! Optional error flag
        real(FP_REAL) :: r1(1,1)

        r1 = sphere_eval_many(this,[theta],[phi],ierr)
        r = r1(1,1)

    end function sphere_eval_one

    !> @brief Load new scattered data and fit a smoothing spherical spline in one call.
    !!
    !! @param[in] theta      Colatitude coordinates.
    !! @param[in] phi        Longitude coordinates.
    !! @param[in] r          Function values.
    !! @param[in] w          Optional weights.
    !! @param[in] smoothing  Optional smoothing factor.
    !! @return Error flag.
    !!
    !! @see sphere
    integer function sphere_new_fit(this,theta,phi,r,w,smoothing)
        class(fitpack_sphere), intent(inout) :: this
        real(FP_REAL), intent(in) :: theta(:),phi(size(theta)),r(size(theta))
        real(FP_REAL), optional, intent(in) :: w(size(theta)) ! node weights
        real(FP_REAL), optional, intent(in) :: smoothing

        call this%new_points(theta,phi,r,w)

        sphere_new_fit = this%fit(smoothing)

    end function sphere_new_fit

    ! =================================================================================================
    ! PARALLEL COMMUNICATION
    ! =================================================================================================

    !> @brief Return the communication buffer size for parallel pack/expand.
    elemental integer(FP_SIZE) function sphere_comm_size(this)
        class(fitpack_sphere), intent(in) :: this
        ! Base fields + sphere-specific scalars:
        ! m, lwrk2, nest(2), nmax, knots(2) = 7
        sphere_comm_size = this%core_comm_size() &
                         + 7 &
                         + FP_COMM_SIZE(this%theta) &
                         + FP_COMM_SIZE(this%phi) &
                         + FP_COMM_SIZE(this%r) &
                         + FP_COMM_SIZE(this%w) &
                         + FP_COMM_SIZE(this%wrk2) &
                         + FP_COMM_SIZE(this%t)
    end function sphere_comm_size

    !> @brief Pack the sphere fitter state into a communication buffer.
    pure subroutine sphere_comm_pack(this, buffer)
        class(fitpack_sphere), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_pack(buffer)
        pos = this%core_comm_size() + 1

        buffer(pos) = real(this%m, FP_COMM);        pos = pos + 1
        buffer(pos) = real(this%lwrk2, FP_COMM);    pos = pos + 1
        buffer(pos) = real(this%nest(1), FP_COMM);  pos = pos + 1
        buffer(pos) = real(this%nest(2), FP_COMM);  pos = pos + 1
        buffer(pos) = real(this%nmax, FP_COMM);     pos = pos + 1
        buffer(pos) = real(this%knots(1), FP_COMM); pos = pos + 1
        buffer(pos) = real(this%knots(2), FP_COMM); pos = pos + 1

        call FP_COMM_PACK(this%theta, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%theta)
        call FP_COMM_PACK(this%phi, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%phi)
        call FP_COMM_PACK(this%r, buffer(pos:));     pos = pos + FP_COMM_SIZE(this%r)
        call FP_COMM_PACK(this%w, buffer(pos:));     pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_PACK(this%wrk2, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk2)
        call FP_COMM_PACK(this%t, buffer(pos:))
    end subroutine sphere_comm_pack

    !> @brief Restore the sphere fitter state from a communication buffer.
    pure subroutine sphere_comm_expand(this, buffer)
        class(fitpack_sphere), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_expand(buffer)
        pos = this%core_comm_size() + 1

        this%m        = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%lwrk2    = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nest(1)  = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nest(2)  = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%nmax     = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(1) = nint(buffer(pos), FP_SIZE);  pos = pos + 1
        this%knots(2) = nint(buffer(pos), FP_SIZE);  pos = pos + 1

        call FP_COMM_EXPAND(this%theta, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%theta)
        call FP_COMM_EXPAND(this%phi, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%phi)
        call FP_COMM_EXPAND(this%r, buffer(pos:));     pos = pos + FP_COMM_SIZE(this%r)
        call FP_COMM_EXPAND(this%w, buffer(pos:));     pos = pos + FP_COMM_SIZE(this%w)
        call FP_COMM_EXPAND(this%wrk2, buffer(pos:)); pos = pos + FP_COMM_SIZE(this%wrk2)
        call FP_COMM_EXPAND(this%t, buffer(pos:))
    end subroutine sphere_comm_expand

end module fitpack_sphere_domains
