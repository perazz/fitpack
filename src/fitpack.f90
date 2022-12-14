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
module fitpack
    use fitpack_core
    implicit none
    private

    ! Public interface
    public :: RKIND
    public :: fitpack_curve
    public :: fitpack_surface

    ! Max fitting degree
    integer, parameter :: MAX_K = 5

    !> A public type describing a surface fitter z = s(x,y)
    type :: fitpack_surface

        !> The data points
        integer :: m = 0
        real(RKIND), allocatable :: x(:),y(:),z(:)

        !> Spline degree
        integer :: order(2) = 3

        !> Interval boundaries
        real(RKIND) :: left(2),right(2)

        ! Node weights
        real(RKIND), allocatable :: w(:)

        ! Estimated and actual number of knots and their allocations
        integer :: nest(2)  = 0
        integer :: nmax = 0
        integer                  :: lwrk1 = 0, lwrk2 = 0, liwrk = 0
        integer, allocatable     :: iwrk(:)
        real(RKIND), allocatable :: wrk1(:),wrk2(:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(RKIND) :: smoothing = 1000.d0

        ! Actual curve MSE
        real(RKIND) :: fp = zero

        ! Curve extrapolation behavior
        integer     :: bc = OUTSIDE_NEAREST_BND

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

    end type fitpack_surface

    !> A public type describing a curve fitter y = c(x)
    type :: fitpack_curve

        !> The data points
        integer :: m = 0
        real(RKIND), allocatable :: x(:),y(:)

        !> Spline degree
        integer :: order = 3

        !> Interval boundaries
        real(RKIND) :: xleft,xright

        ! Node weights
        real(RKIND), allocatable :: sp(:),w(:)

        ! Estimated and actual number of knots and their allocations
        integer                  :: nest  = 0
        integer                  :: lwrk  = 0
        integer, allocatable     :: iwrk(:)
        real(RKIND), allocatable :: wrk(:)

        ! Curve fit smoothing parameter (fit vs. points MSE)
        real(RKIND) :: smoothing = 1000.0_RKIND

        ! Actual curve MSE
        real(RKIND) :: fp = zero

        ! Curve extrapolation behavior
        integer     :: bc = OUTSIDE_NEAREST_BND

        ! Knots
        integer     :: knots = 0
        real(RKIND), allocatable :: t(:)  ! Knot location

        ! Spline coefficients [knots-order-1]
        real(RKIND), allocatable :: c(:)

        ! Runtime flag
        integer :: iopt = 0

        contains

           !> Clean memory
           procedure :: destroy

           !> Set new points
           procedure :: new_points

           !> Generate new fit
           procedure :: new_fit

           !> Generate/update fitting curve, with optional smoothing
           procedure :: fit         => curve_fit_automatic_knots
           procedure :: interpolate => interpolating_curve

           !> Evaluate curve at given coordinates
           procedure, private :: curve_eval_one
           procedure, private :: curve_eval_many
           generic :: eval => curve_eval_one,curve_eval_many

           !> Evaluate derivative at given coordinates
           procedure, private :: curve_derivative
           procedure, private :: curve_derivatives
           generic   :: dfdx => curve_derivative,curve_derivatives

           !> Properties: MSE
           procedure, non_overridable :: mse => curve_error

    end type fitpack_curve

    ! Default constructor
    interface fitpack_curve
       module procedure new_from_points
    end interface fitpack_curve

    interface fitpack_surface
       module procedure surf_new_from_points
    end interface fitpack_surface


    contains

    ! A default constructor
    type(fitpack_curve) function new_from_points(x,y,w,ierr) result(this)
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(x,y,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new curve fit')

    end function new_from_points

    ! Fit a new curve
    integer function new_fit(this,x,y,w,smoothing)
        class(fitpack_curve), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights
        real(RKIND), optional, intent(in) :: smoothing

        call this%new_points(x,y,w)

        new_fit = this%fit(smoothing)

    end function new_fit

    elemental subroutine destroy(this)
       class(fitpack_curve), intent(inout) :: this
       integer :: ierr
       this%m = 0
       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%sp,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this% wrk,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)
       this%xleft = zero
       this%xright = zero

       this%smoothing = 1000.0_RKIND
       this%order     = 3
       this%iopt      = 0
       this%nest      = 0
       this%lwrk      = 0
       this%knots     = 0
       this%fp        = 0.0_RKIND
       this%bc        = OUTSIDE_NEAREST_BND

    end subroutine destroy

    subroutine new_points(this,x,y,w)
        class(fitpack_curve), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights

        integer, allocatable :: isort(:)
        integer, parameter   :: SAFE = 10

        associate(m=>this%m,nest=>this%nest,lwrk=>this%lwrk)

        call this%destroy()

        m = size(x)

        ! Ensure x are sorted
        isort = RKIND_argsort(x)
        allocate(this%x,source=x(isort))
        allocate(this%y,source=y(isort))

        ! set up uniform weights
        if (present(w)) then
           allocate(this%w(m),source=w(isort))
        else
           allocate(this%w(m),source=1.0_RKIND)
        endif
        allocate(this%sp(m),source=0.0_RKIND)

        ! Setup boundaries
        this%xleft  = minval(x)
        this%xright = maxval(x)

        ! Reset run flag
        this%iopt = 0

        ! Reset estimated knots
        nest = max(SAFE*2*MAX_K+2,ceiling(1.4*this%m))
        allocate(this%iwrk(nest),this%t(nest),this%c(nest))

        ! Setup working space.
        lwrk = (m*(MAX_K+1)+nest*(7+3*MAX_K))
        allocate(this%wrk(lwrk),source=0.0_RKIND)

        endassociate

    end subroutine new_points

    real(RKIND) function curve_eval_one(this,x,ierr) result(y)
        class(fitpack_curve), intent(inout)  :: this
        real(RKIND),          intent(in)     :: x      ! Evaluation point
        integer, optional,    intent(out)    :: ierr   ! Optional error flag

        real(RKIND) :: y1(1)

        y1 = curve_eval_many(this,[x],ierr)
        y  = y1(1)

    end function curve_eval_one


    ! Curve evaluation driver
    function curve_eval_many(this,x,ierr) result(y)
        class(fitpack_curve), intent(inout)  :: this
        real(RKIND),          intent(in)     :: x(:)   ! Evaluation points
        integer, optional,    intent(out)    :: ierr   ! Optional error flag
        real(RKIND) :: y(size(x))

        integer :: npts,ier

        npts = size(x)

        !  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
        !  a spline s(x) of degree k, given in its b-spline representation.
        !
        !  calling sequence:
        !     call splev(t,n,c,k,x,y,m,e,ier)
        !
        !  input parameters:
        !    t    : array,length n, which contains the position of the knots.
        !    n    : integer, giving the total number of knots of s(x).
        !    c    : array,length n, which contains the b-spline coefficients.
        !    k    : integer, giving the degree of s(x).
        !    x    : array,length m, which contains the points where s(x) must
        !           be evaluated.
        !    m    : integer, giving the number of points where s(x) must be
        !           evaluated.
        !    e    : integer, if 0 the spline is extrapolated from the end
        !           spans for points not in the support, if 1 the spline
        !           evaluates to zero for those points, if 2 ier is set to
        !           1 and the subroutine returns, and if 3 the spline evaluates
        !           to the value of the nearest boundary point.

        call splev(t=this%t,&                      ! the position of the knots
                   n=this%knots,&                  ! total number of knots of s(x)
                   c=this%c,&                      ! the b-spline coefficients
                   k=this%order,&                  ! the degree of s(x)
                   x=x,m=npts,&                    ! the points where s(x) must be evaluated.
                   y=y, &                          ! the predictions
                   e=this%bc,           &          ! What to do outside mapped knot range
                   ier=ier)

        call fitpack_error_handling(ier,ierr,'evaluate 1d spline')

    end function curve_eval_many

    ! Interpolating curve
    integer function interpolating_curve(this) result(ierr)
        class(fitpack_curve), intent(inout) :: this

        ! Set zero smoothing
        ierr = curve_fit_automatic_knots(this,zero)

    end function interpolating_curve

    ! Curve fitting driver: automatic number of knots
    integer function curve_fit_automatic_knots(this,smoothing) result(ierr)
        class(fitpack_curve), intent(inout) :: this
        real(RKIND), optional, intent(in) :: smoothing

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

        ! First iteration lets solver decide knots
        this%iopt = 0

        do loop=1,nit

            ! Set current smoothing
            this%smoothing = smooth_now(loop)

            ! Call curvfit
            call curfit(this%iopt,                   &  ! option
                        this%m,this%x,this%y,this%w, &  ! points
                        this%xleft,this%xright,      &  ! x range
                        this%order,this%smoothing,   &  ! spline accuracy
                        this%nest,this%knots,this%t, &  ! spline output
                        this%c,this%fp,              &  ! spline output
                        this%wrk,this%lwrk,this%iwrk,&  ! memory
                        ierr)                           ! Error flag
        end do

    end function curve_fit_automatic_knots



    ! Return fitting MSE
    elemental real(RKIND) function curve_error(this)
       class(fitpack_curve), intent(in) :: this
       curve_error = this%fp
    end function curve_error

    ! Utilities: argsort
    ! Return indices of sorted array
    pure function RKIND_argsort(list) result(ilist)
        real(RKIND), dimension(:), intent(in) :: list
        integer(RSIZE), dimension(size(list)) :: ilist

        real(RKIND), allocatable :: copy(:)
        integer(RSIZE) :: i

        ! Prepare data
        allocate(copy(size(list)),source=list)
        forall(i=1:size(list,kind=RSIZE)) ilist(i) = i

        ! Perform sort
        call RKIND_quicksort_andlist(copy,ilist)

        deallocate(copy)

    end function RKIND_argsort

    ! Quicksort
    pure recursive subroutine RKIND_quicksort_andlist(list,ilist,down)

      real(RKIND), dimension(:), intent(inout) :: list
      integer(RSIZE), dimension(size(list)), intent(inout) :: ilist
      logical, optional, intent(in) :: down

      integer(RSIZE)  :: i, j, n
      real(RKIND)  :: chosen
      integer(RSIZE), parameter :: SMALL_SIZE = 8
      logical         :: descending

      descending = .false.; if (present(down)) descending = down
      n = size(list,kind=RSIZE)

      choose_sorting_algorithm: if (n <= SMALL_SIZE) then

         ! Use interchange sort for small lists
         do i = 1, n - 1
            do j = i + 1, n
               if (toBeSwapped(list(i),list(j),.false.)) then
                   call swap_data(list(i),list(j))
                   call swap_size(ilist(i),ilist(j))
               end if
            end do
         end do

      else
         ! Use partition (quick) sort if the list is big
         chosen = list(int(n/2))
         i = 0
         j = n + 1

         scan_lists: do
            ! Scan list from left end
            ! until element >= chosen is found
            scan_from_left: do
              i = i + 1
              if (toBeSwapped(list(i),chosen,.true.) .or. i==n) exit scan_from_left
            end do scan_from_left

            ! Scan list from right end
            ! until element <= chosen is found

            scan_from_right: do
               j = j - 1
               if (toBeSwapped(chosen,list(j),.true.) .or. j==1) exit scan_from_right
            end do scan_from_right

            swap_element: if (i < j) then
                ! Swap two out of place elements
                call swap_data(list(i),list(j))
                call swap_size(ilist(i),ilist(j))
            else if (i == j) then
                i = i + 1
                exit
            else
                exit
            endif swap_element

         end do scan_lists

         if (1 < j) call RKIND_quicksort_andlist(list(:j),ilist(:j),down)
         if (i < n) call RKIND_quicksort_andlist(list(i:),ilist(i:),down)

      end if choose_sorting_algorithm ! test for small array

      contains

         elemental logical function toBeSwapped(a,b,orEqual)
            real(RKIND), intent(in) :: a,b
            logical, intent(in) :: orEqual

            toBeSwapped = merge(is_before(a,b),is_after(a,b),descending)
            if (orEqual .and. a==b) toBeSwapped = .true.

         end function toBeSwapped

    end subroutine RKIND_quicksort_andlist

    elemental subroutine swap_data(a,b)
      real(RKIND), intent(inout) :: a, b
      real(RKIND)                :: tmp
      tmp = a
      a   = b
      b   = tmp
      return
    end subroutine swap_data

    elemental subroutine swap_size(a,b)
      integer(RSIZE), intent(inout) :: a, b
      integer(RSIZE)                :: tmp
      tmp = a
      a   = b
      b   = tmp
      return
    end subroutine swap_size

    elemental logical function is_before(a,b)
       real(RKIND), intent(in) :: a,b
       is_before = a<b
    end function is_before

    elemental logical function is_after(a,b)
       real(RKIND), intent(in) :: a,b
       is_after = a>b
    end function is_after

    elemental logical function is_ge(a,b)
       real(RKIND), intent(in) :: a,b
       is_ge = a>=b
    end function is_ge

    elemental logical function is_le(a,b)
       real(RKIND), intent(in) :: a,b
       is_le = a<=b
    end function is_le

    ! Flow control: on output flag present, return it;
    ! otherwise, halt on error
    subroutine fitpack_error_handling(ierr,ierr_out,whereAt)
        integer, intent(in) :: ierr
        integer, optional, intent(out) :: ierr_out
        character(*), intent(in) :: whereAt


        if (present(ierr_out)) then
            ierr_out = ierr
        elseif (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[fitpack] at '//trim(whereAt)//' failed with error '//FITPACK_MESSAGE(ierr)
            stop ierr
        end if
    end subroutine fitpack_error_handling

    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    function curve_derivatives(this, x, order, ierr) result(ddx)
       class(fitpack_curve), intent(inout) :: this
       real(RKIND),          intent(in)    :: x(:)   ! Evaluation point (scalar)
       integer, optional,    intent(in)    :: order  ! Derivative order. Default 1
       integer, optional,    intent(out)   :: ierr  ! Optional error flag
       real(RKIND), dimension(size(x))     :: ddx

       integer :: ddx_order,m,ierr0

       if (present(order)) then
          ddx_order = max(0,order)
       else
          ddx_order = 1
       end if

       ierr0 = FITPACK_OK


       m = size(x); if (m<=0) goto 1

       !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m the derivative of
       !  order nu of a spline s(x) of degree k, given in its b-spline representation.

       call splder(this%t, & ! Position of the knots
                   this%knots, & ! Number of knots
                   this%c,     & ! spline coefficients
                   this%order, & ! spline degree
                   ddx_order,  & ! derivative order (0<=order<=this%order)  x,y,m,e,wrk,ier)
                   x,          & ! Array of points where this should be evaluated
                   ddx,        & ! Evaluated derivatives
                   m,          & ! Number of input points
                   this%bc,    & ! Extrapolation behavior
                   this%wrk,   & ! Temporary working space
                   ierr0)        ! Output flag

       1 call fitpack_error_handling(ierr0,ierr,'evaluate derivative')

    end function curve_derivatives

    !> Evaluate k-th derivative of the curve at points x
    !> Use 1st derivative if order not present
    real(RKIND) function curve_derivative(this, x, order, ierr) result(ddx)
       class(fitpack_curve), intent(inout)  :: this
       real(RKIND),          intent(in)     :: x      ! Evaluation point (scalar)
       integer, optional,    intent(in)     :: order  ! Derivative order. Default 1
       integer, optional,    intent(out)    :: ierr   ! Optional error flag

       real(RKIND) :: ddxa(1)

       ! Use the array-based wrapper
       ddxa = curve_derivatives(this,[x],order,ierr)
       ddx  = ddxa(1)

    end function curve_derivative

    ! Fit a surface z = s(x,y) defined on a meshgrid: x[1:n], y[1:m]
    integer function surface_fit_automatic_knots(this,smoothing,order) result(ierr)
        class(fitpack_surface), intent(inout) :: this
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

        end do

    end function surface_fit_automatic_knots


    elemental subroutine surf_destroy(this)
       class(fitpack_surface), intent(inout) :: this
       integer :: ierr
       this%m = 0
       deallocate(this%x,stat=ierr)
       deallocate(this%y,stat=ierr)
       deallocate(this%z,stat=ierr)
       deallocate(this%w,stat=ierr)
       deallocate(this%iwrk,stat=ierr)
       deallocate(this%wrk1,stat=ierr)
       deallocate(this%wrk2,stat=ierr)
       deallocate(this%t,stat=ierr)
       deallocate(this%c,stat=ierr)
       this%left  = zero
       this%right = zero

       this%smoothing = 1000.0_RKIND
       this%order     = 3
       this%iopt      = 0
       this%nest      = 0
       this%nmax      = 0
       this%lwrk1     = 0
       this%lwrk2     = 0
       this%liwrk     = 0
       this%knots     = 0
       this%fp        = zero
       this%bc        = OUTSIDE_NEAREST_BND

    end subroutine surf_destroy

    subroutine surf_new_points(this,x,y,z,w)
        class(fitpack_surface), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x)),z(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights

        integer :: clen,uv(2),km,bxy(2),b1,b2
        integer, parameter :: SAFE = 2


        associate(m=>this%m,nest=>this%nest,nmax=>this%nmax,order=>this%order)

        call this%destroy()

        m = size(x)

        ! Ensure x are sorted
        allocate(this%x,source=x)
        allocate(this%y,source=y)
        allocate(this%z,source=z)

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
        real(RKIND), intent(in) :: x(:),y(size(x)),z(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights
        integer, optional, intent(out) :: ierr

        integer :: ierr0

        ierr0 = this%new_fit(x,y,z,w)

        ! Error handling
        call fitpack_error_handling(ierr0,ierr,'new surface fit')

    end function surf_new_from_points

    ! Fit a new curve
    integer function surf_new_fit(this,x,y,z,w,smoothing,order)
        class(fitpack_surface), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x)),z(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights
        real(RKIND), optional, intent(in) :: smoothing
        integer    , optional, intent(in) :: order

        call this%new_points(x,y,z,w)

        surf_new_fit = this%fit(smoothing,order)

    end function surf_new_fit

end module fitpack
