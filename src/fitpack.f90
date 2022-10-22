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

    ! Max fitting degree    
    integer, parameter :: MAX_K = 5

    ! A public type describing a curve fitter
    type, public :: fitpack_fitter

        ! The data points
        integer :: m = 0
        real(RKIND), allocatable :: x(:),y(:)

        ! Spline degree
        integer :: order = 3

        ! Interval boundaries
        real(RKIND) :: xleft,xright

        ! Node weights
        real(RKIND), allocatable :: sp(:),w(:)

        ! Estimated and actual number of knots and their allocations
        integer :: nest  = 0
        integer :: lwrk  = 0
        integer, allocatable :: iwrk(:)
        real(RKIND), allocatable :: wrk(:)

        ! Spline representation
        real(RKIND) :: smoothing = 1000.d0
        real(RKIND) :: fp = 0.d0
        integer      :: knots = 0
        real(RKIND), allocatable :: t(:)  ! Position
        real(RKIND), allocatable :: c(:)  ! Spline coefficients [knots-order-1]

        ! Runtime flag
        integer :: iopt = 0

        contains

           !> Clean memory
           procedure :: destroy

           !> Set new points
           procedure :: new_points

           !> Generate new fit
           procedure :: new_fit

           !> Generate/update fit
           procedure :: fit => curve_fit_automatic_knots

           !> Evaluate curve at given coordinates
           procedure :: eval => curve_eval

    end type fitpack_fitter

    contains

    integer function new_fit(this,x,y,w)
        class(fitpack_fitter), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights

        call this%new_points(x,y,w)
        new_fit = this%fit()

    end function new_fit

    elemental subroutine destroy(this)
       class(fitpack_fitter), intent(inout) :: this
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
       this%xleft = 0.d0
       this%xright = 0.d0

       this%smoothing = 1000.d0
       this%order = 3
       this%iopt = 0
       this%nest = 0
       this%lwrk = 0
       this%knots = 0
       this%fp = 0.d0


    end subroutine destroy

    subroutine new_points(this,x,y,w)
        class(fitpack_fitter), intent(inout) :: this
        real(RKIND), intent(in) :: x(:),y(size(x))
        real(RKIND), optional, intent(in) :: w(size(x)) ! node weights

        integer, allocatable :: isort(:)

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
        nest = max(2*MAX_K+2,ceiling(1.4*this%m))
        allocate(this%iwrk(nest),this%t(nest),this%c(nest))

        ! Setup working space.
        lwrk = (m*(MAX_K+1)+nest*(7+3*MAX_K))
        allocate(this%wrk(lwrk),source=0.0_RKIND)

        endassociate

    end subroutine new_points

    ! Curve evaluation driver
    integer function curve_eval(this,x,y) result(ierr)
        class(fitpack_fitter), intent(inout) :: this
        real(RKIND), intent(in) :: x(:)
        real(RKIND), intent(out) :: y(size(x))

        integer :: npts

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
                   e=SPLINE_NEAREST_BND,&          ! What to do outside mapped knot range
                   ier=ierr)

!        call curev(idim=1,&                        ! dimension of the spline curve.
!                   t=this%t,&                      ! the position of the knots
!                   n=this%knots,&                  ! total number of knots of s(u)
!                   c=this%c,&                      ! the b-spline coefficients
!                   nc=this%knots-this%order+1, &   ! number of coefficients of s(u)
!                   k=this%order,&                  ! the degree of s(u)
!                   u=x,m=npts,&                    ! the points where s(u) must be evaluated.
!                   x=y,&
!                   mx=npts,&                       ! the dimension of array x
!                   ier=ierr)

    end function curve_eval

    ! Curve fitting driver: automatic number of knots
    integer function curve_fit_automatic_knots(this,smoothing) result(ierr)
        class(fitpack_fitter), intent(inout) :: this
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

end module fitpack
