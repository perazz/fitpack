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
module fitpack_core
    use iso_fortran_env, only: real64,int32
    implicit none
    private

    ! Precision and array size
    integer, parameter, public :: RKIND = real64
    integer, parameter, public :: RSIZE = int32

    ! Polar
    public :: evapol ! Evaluate polar spline s(u,v) at (x,y)
    public :: polar,pogrid,percur,cocosp,concon,concur,cualde
    public :: dblint,parcur,parder,pardeu,pardtc

    ! Fourier coefficients
    public :: fourco

    ! Closed-curve
    public :: clocur

    ! Surface
    public :: surfit

    public :: parsur
    public :: surev   ! Evaluate parametric surface on an (u(i),v(j)) grid

    ! Curve
    public :: curev,curfit

    ! B-Spline
    public :: bispev  ! Evaluate 2d spline on a meshgrid
    public :: bispeu  ! Evaluate 2d spline on arbitrary (x,y) points
    public :: insert,sphere,spgrid


    public :: splev   ! Evaluate 1d spline value
    public :: splder  ! Evaluate 1d spline derivative
    public :: spalde  ! Evaluate all 1d derivatives at x
    public :: splint  ! Evaluate integral below spline
    public :: sproot  ! Find zeroes of a 1d spline

    public :: regrid
    public :: profil

    ! Spline behavior for points not in the support
    integer, parameter, public :: OUTSIDE_EXTRAPOLATE = 0 ! extrapolated from the end spans
    integer, parameter, public :: OUTSIDE_ZERO        = 1 ! spline evaluates to zero
    integer, parameter, public :: OUTSIDE_NOT_ALLOWED = 2 ! an error flag is returned
    integer, parameter, public :: OUTSIDE_NEAREST_BND = 3 ! evaluate to value of nearest boundary point

                        public :: FITPACK_MESSAGE
                        public :: FITPACK_SUCCESS
    integer, parameter, public :: FITPACK_OK                   = 0  ! ok for spline, abs(fp-s)/s <= tol=0.001
    integer, parameter, public :: FITPACK_INTERPOLATING_OK     = -1 ! ok for interpolating spline, fp=0
    integer, parameter, public :: FITPACK_LEASTSQUARES_OK      = -2 ! ok for weighted least-squares polynomial of degree k.
    integer, parameter, public :: FITPACK_INSUFFICIENT_STORAGE = 1
    integer, parameter, public :: FITPACK_S_TOO_SMALL          = 2
    integer, parameter, public :: FITPACK_MAXIT                = 3
    integer, parameter, public :: FITPACK_TOO_MANY_KNOTS       = 4
    integer, parameter, public :: FITPACK_OVERLAPPING_KNOTS    = 5
    integer, parameter, public :: FITPACK_INVALID_RANGE        = 6
    integer, parameter, public :: FITPACK_INPUT_ERROR          = 10

    ! Internal Parameters
    integer    , parameter, public :: SIZ_IDIM = 10                ! Max number of dimensions
    integer    , parameter, public :: SIZ_K    = 19                ! Max order (for array allocation)
    real(RKIND), parameter, public :: one    = 1.0_RKIND
    real(RKIND), parameter, public :: zero   = 0.0_RKIND
    real(RKIND), parameter, public :: half   = 0.5_RKIND
    real(RKIND), parameter, public :: onep5  = 1.5_RKIND
    real(RKIND), parameter, public :: fourth = 0.25_RKIND
    real(RKIND), parameter, public :: two    = 2.0_RKIND
    real(RKIND), parameter, public :: three  = 3.0_RKIND
    real(RKIND), parameter, public :: four   = 4.0_RKIND
    real(RKIND), parameter, public :: six    = 6.0_RKIND
    real(RKIND), parameter, public :: ten    = 10.0_RKIND
    real(RKIND), parameter, public :: pi     = atan2(zero,-one)
    real(RKIND), parameter, public :: pi2    = 2*pi
    real(RKIND), parameter, public :: pi4    = 4*pi
    real(RKIND), parameter, public :: smallnum03 = 1.0e-03_RKIND
    real(RKIND), parameter, public :: smallnum10 = 1.0e-10_RKIND

    abstract interface
       ! Function defining the boundary of the curve approximation domain
       pure real(RKIND) function boundary(v) result(rad)
          import RKIND
          real(RKIND), intent(in) :: v
       end function boundary
    end interface

    contains

      ! Wrapper for the error flag
      pure function FITPACK_MESSAGE(ierr) result(msg)
         integer, intent(in) :: ierr
         character(len=:), allocatable :: msg

         select case (ierr)
            case (FITPACK_OK); msg = 'Success!'
            case (FITPACK_INTERPOLATING_OK); msg = 'Success! (interpolation)'
            case (FITPACK_LEASTSQUARES_OK); msg = 'Success! (least-squares)'
            case (FITPACK_INSUFFICIENT_STORAGE); msg = 'Insufficient Storage'
            case (FITPACK_S_TOO_SMALL); msg = 'Smoothing parameter is too small'
            case (FITPACK_MAXIT); msg = 'Infinite loop detected'
            case (FITPACK_TOO_MANY_KNOTS); msg = 'More knots than data points'
            case (FITPACK_OVERLAPPING_KNOTS); msg = 'Overlapping knots found'
            case (FITPACK_INVALID_RANGE); msg = 'Invalid variable range'
            case (FITPACK_INPUT_ERROR); msg = 'Invalid input'
            case default; msg = 'UNKNOWN ERROR'
         end select

      end function FITPACK_MESSAGE

      ! Wrapper for OK
      elemental logical function FITPACK_SUCCESS(ierr)
         integer, intent(in) :: ierr
         FITPACK_SUCCESS = ierr<=FITPACK_OK
      end function FITPACK_SUCCESS

      pure subroutine bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk,ier)

      !  subroutine bispeu evaluates on a set of points (x(i),y(i)),i=1,...,m
      !  a bivariate spline s(x,y) of degrees kx and ky, given in the
      !  b-spline representation.
      !
      !  calling sequence:
      !     call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk,iwrk,kwrk,ier)
      !
      !  input parameters:
      !   tx    : real array, length nx, which contains the position of the knots in the x-direction.
      !   nx    : integer, giving the total number of knots in the x-direction
      !   ty    : real array, length ny, which contains the position of the knots in the y-direction.
      !   ny    : integer, giving the total number of knots in the y-direction
      !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the b-spline coefficients.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   x     : real array of dimension (mx).
      !   y     : real array of dimension (my).
      !   m     : on entry m must specify the number points. m >= 1.
      !   wrk   : real array of dimension lwrk. used as workspace.
      !   lwrk  : integer, specifying the dimension of wrk. lwrk >= kx+ky+2
      !
      !  output parameters:
      !   z     : real array of dimension m.
      !           on successful exit z(i) contains the value of s(x,y)
      !           at the point (x(i),y(i)), i=1,...,m.
      !   ier   : integer error flag
      !
      !  restrictions:
      !   m >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
      !   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
      !   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
      !
      !  other subroutines required:
      !    fpbisp,fpbspl
      !
      !  ..scalar arguments..
      integer, intent(in) :: nx,ny,kx,ky,m,lwrk
      integer, intent(out) :: ier

      !  ..array arguments..
      real(RKIND), intent(in)    :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(m),y(m)
      real(RKIND), intent(inout) :: wrk(lwrk)
      real(RKIND), intent(out)   :: z(m)

      !  ..local scalars..
      integer :: iwrk(2),i,lwest

      !  Check inputs
      lwest = kx+ky+2

      if (lwrk<lwest .or. m<1) then

         ier = FITPACK_INPUT_ERROR
         return

      else

         ier = FITPACK_OK
         do i=1,m
            call fpbisp(tx,nx,ty,ny,c,kx,ky,x(i),1,y(i),1,z(i),wrk(1),wrk(kx+2),iwrk(1),iwrk(2))
         end do

      end if

      end subroutine bispeu


      !  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...,my a bivariate spline
      !  s(x,y) of degrees kx and ky, given in the b-spline representation.
      pure subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk,ier)

      !
      !  input parameters:
      !   tx    : real array, length nx, which contains the position of the knots in the x-direction.
      !   nx    : integer, giving the total number of knots in the x-direction
      !   ty    : real array, length ny, which contains the position of the knots in the y-direction.
      !   ny    : integer, giving the total number of knots in the y-direction
      !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the b-spline coefficients.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   x     : real array of dimension (mx).
      !           before entry x(i) must be set to the x co-ordinate of the i-th grid point along the x-axis.
      !           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
      !   mx    : on entry mx must specify the number of grid points along the x-axis. mx >=1.
      !   y     : real array of dimension (my).
      !           before entry y(j) must be set to the y co-ordinate of the j-th grid point along the y-axis.
      !           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
      !   my    : on entry my must specify the number of grid points along the y-axis. my >=1.
      !   wrk   : real array of dimension lwrk. used as workspace.
      !   lwrk  : integer, specifying the dimension of wrk.
      !           lwrk >= mx*(kx+1)+my*(ky+1)
      !   iwrk  : integer array of dimension kwrk. used as workspace.
      !   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
      !
      !  output parameters:
      !   z     : real array of dimension (mx*my).
      !           on successful exit z(my*(i-1)+j) contains the value of s(x,y)
      !           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
      !   ier   : integer error flag
      !
      !  restrictions:
      !   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
      !   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
      !   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
      !
      !  other subroutines required:
      !    fpbisp,fpbspl
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory
      !                6 (1972) 50-62.
      !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
      !                applics 10 (1972) 134-149.
      !    dierckx p. : curve and surface fitting with splines, monographs on
      !                 numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1987
      !
      !  ..scalar arguments..
      integer, intent(in) :: nx,ny,kx,ky,mx,my,lwrk,kwrk
      integer, intent(out) :: ier
      !  ..array arguments..
      integer, intent(inout) :: iwrk(kwrk)
      real(RKIND), intent(in)    :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my)
      real(RKIND), intent(out)   :: z(mx*my)
      real(RKIND), intent(inout) :: wrk(lwrk)
      !  ..local scalars..
      integer :: iw,lwest
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier   = FITPACK_INPUT_ERROR
      lwest = (kx+1)*mx+(ky+1)*my

      if (lwrk<lwest .or. kwrk<(mx+my) .or. mx<1 .or. my<1) return
      if (mx>1 .and. any(x(2:mx)<x(1:mx-1))) return
      if (my>1 .and. any(y(2:my)<y(1:my-1))) return

      ! Evaluate spline
      ier = FITPACK_OK
      iw  = mx*(kx+1)+1
      call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(1),wrk(iw),iwrk(1),iwrk(mx+1))

      end subroutine bispev


      recursive subroutine clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest, &
       n,t,nc,c,fp,wrk,lwrk,iwrk,ier)

      !  given the ordered set of m points x(i) in the idim-dimensional space
      !  with x(1)=x(m), and given also a corresponding set of strictly in-
      !  creasing values u(i) and the set of positive numbers w(i),i=1,2,...,m
      !  subroutine clocur determines a smooth approximating closed spline
      !  curve s(u), i.e.
      !      x1 = s1(u)
      !      x2 = s2(u)       u(1) <= u <= u(m)
      !      .........
      !      xidim = sidim(u)
      !  with sj(u),j=1,2,...,idim periodic spline functions of degree k with
      !  common knots t(j),j=1,2,...,n.
      !  if ipar=1 the values u(i),i=1,2,...,m must be supplied by the user.
      !  if ipar=0 these values are chosen automatically by clocur as
      !      v(1) = 0
      !      v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
      !      u(i) = v(i)/v(m) ,i=1,2,...,m
      !  if iopt=-1 clocur calculates the weighted least-squares closed spline
      !  curve according to a given set of knots.
      !  if iopt>=0 the number of knots of the splines sj(u) and the position
      !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
      !  ness of s(u) is then achieved by minimalizing the discontinuity
      !  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
      !  n-k-1. the amount of smoothness is determined by the condition that
      !  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
      !  negative constant, called the smoothing factor.
      !  the fit s(u) is given in the b-spline representation and can be
      !  evaluated by means of subroutine curev.
      !
      !  calling sequence:
      !     call clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,
      !    * fp,wrk,lwrk,iwrk,ier)
      !
      !  parameters:
      !   iopt  : integer flag. on entry iopt must specify whether a weighted
      !           least-squares closed spline curve (iopt=-1) or a smoothing
      !           closed spline curve (iopt=0 or 1) must be determined. if
      !           iopt=0 the routine will start with an initial set of knots
      !           t(i)=u(1)+(u(m)-u(1))*(i-k-1),i=1,2,...,2*k+2. if iopt=1 the
      !           routine will continue with the knots found at the last call.
      !           attention: a call with iopt=1 must always be immediately
      !           preceded by another call with iopt=1 or iopt=zero
      !           unchanged on exit.
      !   ipar  : integer flag. on entry ipar must specify whether (ipar=1)
      !           the user will supply the parameter values u(i),or whether
      !           (ipar=0) these values are to be calculated by clocur.
      !           unchanged on exit.
      !   idim  : integer. on entry idim must specify the dimension of the
      !           curve. 0 < idim < 11.
      !           unchanged on exit.
      !   m     : integer. on entry m must specify the number of data points.
      !           m > 1. unchanged on exit.
      !   u     : real array of dimension at least (m). in case ipar=1,before
      !           entry, u(i) must be set to the i-th value of the parameter
      !           variable u for i=1,2,...,m. these values must then be
      !           supplied in strictly ascending order and will be unchanged
      !           on exit. in case ipar=0, on exit,the array will contain the
      !           values u(i) as determined by clocur.
      !   mx    : integer. on entry mx must specify the actual dimension of
      !           the array x as declared in the calling (sub)program. mx must
      !           not be too small (see x). unchanged on exit.
      !   x     : real array of dimension at least idim*m.
      !           before entry, x(idim*(i-1)+j) must contain the j-th coord-
      !           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
      !           idim. since first and last data point must coincide it
      !           means that x(j)=x(idim*(m-1)+j),j=1,2,...,idim.
      !           unchanged on exit.
      !   w     : real array of dimension at least (m). before entry, w(i)
      !           must be set to the i-th value in the set of weights. the
      !           w(i) must be strictly positive. w(m) is not used.
      !           unchanged on exit. see also further comments.
      !   k     : integer. on entry k must specify the degree of the splines.
      !           1<=k<=5. it is recommended to use cubic splines (k=3).
      !           the user is strongly dissuaded from choosing k even,together
      !           with a small s-value. unchanged on exit.
      !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
      !           factor. s >=zero unchanged on exit.
      !           for advice on the choice of s see further comments.
      !   nest  : integer. on entry nest must contain an over-estimate of the
      !           total number of knots of the splines returned, to indicate
      !           the storage space available to the routine. nest >=2*k+2.
      !           in most practical situation nest=m/2 will be sufficient.
      !           always large enough is nest=m+2*k, the number of knots
      !           needed for interpolation (s=0). unchanged on exit.
      !   n     : integer.
      !           unless ier = 10 (in case iopt >=0), n will contain the
      !           total number of knots of the smoothing spline curve returned
      !           if the computation mode iopt=1 is used this value of n
      !           should be left unchanged between subsequent calls.
      !           in case iopt=-1, the value of n must be specified on entry.
      !   t     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the knots of the
      !           spline curve,i.e. the position of the interior knots t(k+2),
      !           t(k+3),..,t(n-k-1) as well as the position of the additional
      !           t(1),t(2),..,t(k+1)=u(1) and u(m)=t(n-k),...,t(n) needed for
      !           the b-spline representation.
      !           if the computation mode iopt=1 is used, the values of t(1),
      !           t(2),...,t(n) should be left unchanged between subsequent
      !           calls. if the computation mode iopt=-1 is used, the values
      !           t(k+2),...,t(n-k-1) must be supplied by the user, before
      !           entry. see also the restrictions (ier=10).
      !   nc    : integer. on entry nc must specify the actual dimension of
      !           the array c as declared in the calling (sub)program. nc
      !           must not be too small (see c). unchanged on exit.
      !   c     : real array of dimension at least (nest*idim).
      !           on successful exit, this array will contain the coefficients
      !           in the b-spline representation of the spline curve s(u),i.e.
      !           the b-spline coefficients of the spline sj(u) will be given
      !           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
      !   fp    : real. unless ier = 10, fp contains the weighted sum of
      !           squared residuals of the spline curve returned.
      !   wrk   : real array of dimension at least m*(k+1)+nest*(7+idim+5*k).
      !           used as working space. if the computation mode iopt=1 is
      !           used, the values wrk(1),...,wrk(n) should be left unchanged
      !           between subsequent calls.
      !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
      !           the array wrk as declared in the calling (sub)program. lwrk
      !           must not be too small (see wrk). unchanged on exit.
      !   iwrk  : integer array of dimension at least (nest).
      !           used as working space. if the computation mode iopt=1 is
      !           used,the values iwrk(1),...,iwrk(n) should be left unchanged
      !           between subsequent calls.
      !   ier   : integer. unless the routine detects an error, ier contains a
      !           non-positive value on exit, i.e.
      !    ier=0  : normal return. the close curve returned has a residual
      !             sum of squares fp such that abs(fp-s)/s <= tol with tol a
      !             relative tolerance set to zero001 by the program.
      !    ier=-1 : normal return. the curve returned is an interpolating
      !             spline curve (fp=0).
      !    ier=-2 : normal return. the curve returned is the weighted least-
      !             squares point,i.e. each spline sj(u) is a constant. in
      !             this extreme case fp gives the upper bound fp0 for the
      !             smoothing factor s.
      !    ier=1  : error. the required storage space exceeds the available
      !             storage space, as specified by the parameter nest.
      !             probably causes : nest too small. if nest is already
      !             large (say nest > m/2), it may also indicate that s is
      !             too small
      !             the approximation returned is the least-squares closed
      !             curve according to the knots t(1),t(2),...,t(n). (n=nest)
      !             the parameter fp gives the corresponding weighted sum of
      !             squared residuals (fp>s).
      !    ier=2  : error. a theoretically impossible result was found during
      !             the iteration process for finding a smoothing curve with
      !             fp = s. probably causes : s too small.
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=3  : error. the maximal number of iterations maxit (set to 20
      !             by the program) allowed for finding a smoothing curve
      !             with fp=s has been reached. probably causes : s too small
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=10 : error. on entry, the input data are controlled on validity
      !             the following restrictions must be satisfied.
      !             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,2,...,m
      !             0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(7+idim+5*k),
      !             nc>=nest*idim, x(j)=x(idim*(m-1)+j), j=1,2,...,idim
      !             if ipar=0: sum j=1,idim (x(i*idim+j)-x((i-1)*idim+j))**2>0
      !                        i=1,2,...,m-1.
      !             if ipar=1: u(1)<u(2)<...<u(m)
      !             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k)
      !                         u(1)<t(k+2)<t(k+3)<...<t(n-k-1)<u(m)
      !                            (u(1)=0 and u(m)=1 in case ipar=0)
      !                       the schoenberg-whitney conditions, i.e. there
      !                       must be a subset of data points uu(j) with
      !                       uu(j) = u(i) or u(i)+(u(m)-u(1)) such that
      !                         t(j) < uu(j) < t(j+k+1), j=k+1,...,n-k-1
      !             if iopt>=0: s>=0
      !                         if s=0 : nest >= m+2*k
      !             if one of these conditions is found to be violated,control
      !             is immediately repassed to the calling program. in that
      !             case there is no approximation returned.
      !
      !  further comments:
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the curve will be too smooth and signal will be
      !   lost ; if s is too small the curve will pick up too much noise. in
      !   the extreme cases the program will return an interpolating curve if
      !   s=0 and the weighted least-squares point if s is very large.
      !   between these extremes, a properly chosen s will result in a good
      !   compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the weights w(i). if these are
      !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
      !   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
      !   sqrt(2*m)). if nothing is known about the statistical error in x(i)
      !   each w(i) can be set equal to one and s determined by trial and
      !   error, taking account of the comments above. the best is then to
      !   start with a very large value of s ( to determine the weighted
      !   least-squares point and the upper bound fp0 for s) and then to
      !   progressively decrease the value of s ( say by a factor 10 in the
      !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
      !   approximating curve shows more detail) to obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=zero
      !   if iopt=1 the program will continue with the set of knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if clocur is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   curve underlying the data. but, if the computation mode iopt=1 is
      !   used, the knots returned may also depend on the s-values at previous
      !   calls (if these were smaller). therefore, if after a number of
      !   trials with different s-values and iopt=1, the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   clocur once more with the selected value for s but now with iopt=zero
      !   indeed, clocur may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !
      !   the form of the approximating curve can strongly be affected  by
      !   the choice of the parameter values u(i). if there is no physical
      !   reason for choosing a particular parameter u, often good results
      !   will be obtained with the choice of clocur(in case ipar=0), i.e.
      !        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
      !   where
      !        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
      !   other possibilities for q(i) are
      !        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
      !        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
      !        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
      !        q(i)= 1
      !
      !
      !  other subroutines required:
      !    fpbacp,fpbspl,fpchep,fpclos,fpdisc,fpgivs,fpknot,fprati,fprota
      !
      !  references:
      !   dierckx p. : algorithms for smoothing data with periodic and
      !                parametric splines, computer graphics and image
      !                processing 20 (1982) 171-184.
      !   dierckx p. : algorithms for smoothing data with periodic and param-
      !                etric splines, report tw55, dept. computer science,
      !                k.u.leuven, 1981.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !  latest update : march 1987
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) s,fp
      integer iopt,ipar,idim,m,mx,k,nest,n,nc,lwrk,ier
      !  ..array arguments..
      real(RKIND) u(m),x(mx),w(m),t(nest),c(nc),wrk(lwrk)
      integer iwrk(nest)
      !  ..local scalars..
      real(RKIND) per,dist
      integer i,ia1,ia2,ib,ifp,ig1,ig2,iq,iz,i1,i2,j1,j2,k1,k2,lwest,m1,nmin,ncc,j
      !  we set up the parameters tol and maxit
      integer, parameter :: maxit = 20
      real(RKIND), parameter :: tol = smallnum03
      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if(iopt<(-1) .or. iopt>1) go to 90
      if(ipar<0 .or. ipar>1) go to 90
      if(idim<=0 .or. idim>10) go to 90
      if(k<=0 .or. k>5) go to 90
      k1 = k+1
      k2 = k1+1
      nmin = 2*k1
      if(m<2 .or. nest<nmin) go to 90
      ncc = nest*idim
      if(mx<m*idim .or. nc<ncc) go to 90
      lwest = m*k1+nest*(7+idim+5*k)
      if(lwrk<lwest) go to 90
      i1 = idim
      i2 = m*idim
      do 5 j=1,idim
         if(x(i1)/=x(i2)) go to 90
         i1 = i1-1
         i2 = i2-1
   5  continue
      if(ipar/=0 .or. iopt>0) go to 40
      i1 = 0
      i2 = idim
      u(1) = zero
      do 20 i=2,m
         dist = zero
         do 10 j1=1,idim
            i1 = i1+1
            i2 = i2+1
            dist = dist+(x(i2)-x(i1))**2
  10     continue
         u(i) = u(i-1)+sqrt(dist)
  20  continue
      if(u(m)<=0.) go to 90
      do 30 i=2,m
         u(i) = u(i)/u(m)
  30  continue
      u(m) = one
  40  if(w(1)<=zero) go to 90
      m1 = m-1
      do 50 i=1,m1
         if(u(i)>=u(i+1) .or. w(i)<=zero) go to 90
  50  continue
      if(iopt>=0) go to 70
      if(n<=nmin .or. n>nest) go to 90
      per = u(m)-u(1)
      j1 = k1
      t(j1) = u(1)
      i1 = n-k
      t(i1) = u(m)
      j2 = j1
      i2 = i1
      do 60 i=1,k
         i1 = i1+1
         i2 = i2-1
         j1 = j1+1
         j2 = j2-1
         t(j2) = t(i2)-per
         t(i1) = t(j1)+per
  60  continue
      ier = fpchep(u,m,t,n,k)
      if (ier==FITPACK_OK) go to 80
      go to 90
  70  if(s<0.) go to 90
      if(s==zero .and. nest<(m+2*k)) go to 90
      ier = 0
      ! we partition the working space and determine the spline approximation.
  80  ifp = 1
      iz = ifp+nest
      ia1 = iz+ncc
      ia2 = ia1+nest*k1
      ib = ia2+nest*k
      ig1 = ib+nest*k2
      ig2 = ig1+nest*k2
      iq = ig2+nest*k1
      call fpclos(iopt,idim,m,u,mx,x,w,k,s,nest,tol,maxit,k1,k2,n,t, &
       ncc,c,fp,wrk(ifp),wrk(iz),wrk(ia1),wrk(ia2),wrk(ib),wrk(ig1), &
       wrk(ig2),wrk(iq),iwrk,ier)
  90  return
      end subroutine clocur


      recursive subroutine cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq, &
       sx,bind,wrk,lwrk,iwrk,kwrk,ier)

      !  given the set of data points (x(i),y(i)) and the set of positive
      !  numbers w(i),i=1,2,...,m, subroutine cocosp determines the weighted
      !  least-squares cubic spline s(x) with given knots t(j),j=1,2,...,n
      !  which satisfies the following concavity/convexity conditions
      !      s''(t(j+3))*e(j) <= 0, j=1,2,...n-6
      !  the fit is given in the b-spline representation( b-spline coef-
      !  ficients c(j),j=1,2,...n-4) and can be evaluated by means of
      !  subroutine splev.
      !
      !  calling sequence:
      !     call cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,
      !    * lwrk,iwrk,kwrk,ier)
      !
      !  parameters:
      !    m   : integer. on entry m must specify the number of data points.
      !          m > 3. unchanged on exit.
      !    x   : real array of dimension at least (m). before entry, x(i)
      !          must be set to the i-th value of the independent variable x,
      !          for i=1,2,...,m. these values must be supplied in strictly
      !          ascending order. unchanged on exit.
      !    y   : real array of dimension at least (m). before entry, y(i)
      !          must be set to the i-th value of the dependent variable y,
      !          for i=1,2,...,m. unchanged on exit.
      !    w   : real array of dimension at least (m). before entry, w(i)
      !          must be set to the i-th value in the set of weights. the
      !          w(i) must be strictly positive. unchanged on exit.
      !    n   : integer. on entry n must contain the total number of knots
      !          of the cubic spline. m+4>=n>=8. unchanged on exit.
      !    t   : real array of dimension at least (n). before entry, this
      !          array must contain the knots of the spline, i.e. the position
      !          of the interior knots t(5),t(6),...,t(n-4) as well as the
      !          position of the boundary knots t(1),t(2),t(3),t(4) and t(n-3)
      !          t(n-2),t(n-1),t(n) needed for the b-spline representation.
      !          unchanged on exit. see also the restrictions (ier=10).
      !    e   : real array of dimension at least (n). before entry, e(j)
      !          must be set to 1 if s(x) must be locally concave at t(j+3),
      !          to (-1) if s(x) must be locally convex at t(j+3) and to 0
      !          if no convexity constraint is imposed at t(j+3),j=1,2,..,n-6.
      !          e(n-5),...,e(n) are not used. unchanged on exit.
      !  maxtr : integer. on entry maxtr must contain an over-estimate of the
      !          total number of records in the used tree structure, to indic-
      !          ate the storage space available to the routine. maxtr >=1
      !          in most practical situation maxtr=100 will be sufficient.
      !          always large enough is
      !                         n-5       n-6
      !              maxtr =  (     ) + (     )  with l the greatest
      !                          l        l+1
      !          integer <= (n-6)/2 . unchanged on exit.
      !  maxbin: integer. on entry maxbin must contain an over-estimate of the
      !          number of knots where s(x) will have a zero second derivative
      !          maxbin >=1. in most practical situation maxbin = 10 will be
      !          sufficient. always large enough is maxbin=n-6.
      !          unchanged on exit.
      !    c   : real array of dimension at least (n).
      !          on successful exit, this array will contain the coefficients
      !          c(1),c(2),..,c(n-4) in the b-spline representation of s(x)
      !    sq  : real. on successful exit, sq contains the weighted sum of
      !          squared residuals of the spline approximation returned.
      !    sx  : real array of dimension at least m. on successful exit
      !          this array will contain the spline values s(x(i)),i=1,...,m
      !   bind : logical array of dimension at least (n). on successful exit
      !          this array will indicate the knots where s''(x)=0, i.e.
      !                s''(t(j+3)) == 0 if  bind(j) = .true.
      !                s''(t(j+3)) /= 0 if  bind(j) = .false., j=1,2,...,n-6
      !   wrk  : real array of dimension at least  m*4+n*7+maxbin*(maxbin+n+1)
      !          used as working space.
      !   lwrk : integer. on entry,lwrk must specify the actual dimension of
      !          the array wrk as declared in the calling (sub)program.lwrk
      !          must not be too small (see wrk). unchanged on exit.
      !   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1))
      !          used as working space.
      !   kwrk : integer. on entry,kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program. kwrk
      !          must not be too small (see iwrk). unchanged on exit.
      !   ier   : integer. error flag
      !      ier=0 : successful exit.
      !      ier>0 : abnormal termination: no approximation is returned
      !        ier=1  : the number of knots where s''(x)=0 exceeds maxbin.
      !                 probably causes : maxbin too small.
      !        ier=2  : the number of records in the tree structure exceeds
      !                 maxtr.
      !                 probably causes : maxtr too small.
      !        ier=3  : the algorithm finds no solution to the posed quadratic
      !                 programming problem.
      !                 probably causes : rounding errors.
      !        ier=10 : on entry, the input data are controlled on validity.
      !                 the following restrictions must be satisfied
      !                   m>3, maxtr>=1, maxbin>=1, 8<=n<=m+4,w(i) > 0,
      !                   x(1)<x(2)<...<x(m), t(1)<=t(2)<=t(3)<=t(4)<=x(1),
      !                   x(1)<t(5)<t(6)<...<t(n-4)<x(m)<=t(n-3)<=...<=t(n),
      !                   kwrk>=maxtr*4+2*(maxbin+1),
      !                   lwrk>=m*4+n*7+maxbin*(maxbin+n+1),
      !                   the schoenberg-whitney conditions, i.e. there must
      !                   be a subset of data points xx(j) such that
      !                     t(j) < xx(j) < t(j+4), j=1,2,...,n-4
      !                 if one of these restrictions is found to be violated
      !                 control is immediately repassed to the calling program
      !
      !
      !  other subroutines required:
      !    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno,fpchec
      !
      !  references:
      !   dierckx p. : an algorithm for cubic spline fitting with convexity
      !                constraints, computing 24 (1980) 349-371.
      !   dierckx p. : an algorithm for least-squares cubic spline fitting
      !                with convexity and concavity constraints, report tw39,
      !                dept. computer science, k.u.leuven, 1978.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !   p. dierckx
      !   dept. computer science, k.u.leuven
      !   celestijnenlaan 200a, b-3001 heverlee, belgium.
      !   e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : march 1978
      !  latest update : march 1987.
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) sq
      integer m,n,maxtr,maxbin,lwrk,kwrk,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),w(m),t(n),e(n),c(n),sx(m),wrk(lwrk)
      integer iwrk(kwrk)
      logical bind(n)
      !  ..local scalars..
      integer i,ia,ib,ic,iq,iu,iz,izz,ji,jib,jjb,jl,jr,ju,kwest, &
       lwest,mb,nm,n6

      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if(m<4 .or. n<8) go to 40
      if(maxtr<1 .or. maxbin<1) go to 40
      lwest = 7*n+m*4+maxbin*(1+n+maxbin)
      kwest = 4*maxtr+2*(maxbin+1)
      if(lwrk<lwest .or. kwrk<kwest) go to 40
      if(w(1)<=0.) go to 40
      do 10 i=2,m
         if(x(i-1)>=x(i) .or. w(i)<=0.) go to 40
  10  continue
      call fpchec(x,m,t,n,3,ier)
      if (ier==FITPACK_OK) go to 20
      go to 40
      !  set numbers e(i)
  20  n6 = n-6
      do 30 i=1,n6
        if(e(i)>0.) e(i) = one
        if(e(i)<0.) e(i) = -one
  30  continue
      !  we partition the working space and determine the spline approximation
      nm = n+maxbin
      mb = maxbin+1
      ia = 1
      ib = ia+4*n
      ic = ib+nm*maxbin
      iz = ic+n
      izz = iz+n
      iu = izz+n
      iq = iu+maxbin
      ji = 1
      ju = ji+maxtr
      jl = ju+maxtr
      jr = jl+maxtr
      jjb = jr+maxtr
      jib = jjb+mb
      call fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,wrk(ia), &
       wrk(ib),wrk(ic),wrk(iz),wrk(izz),wrk(iu),wrk(iq),iwrk(ji), &
       iwrk(ju),iwrk(jl),iwrk(jr),iwrk(jjb),iwrk(jib),ier)
  40  return
      end subroutine cocosp


      recursive subroutine concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin, &
       n,t,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier)

      !  given the set of data points (x(i),y(i)) and the set of positive
      !  numbers w(i), i=1,2,...,m,subroutine concon determines a cubic spline
      !  approximation s(x) which satisfies the following local convexity
      !  constraints  s''(x(i))*v(i) <= 0, i=1,2,...,m.
      !  the number of knots n and the position t(j),j=1,2,...n is chosen
      !  automatically by the routine in a way that
      !       sq = sum((w(i)*(y(i)-s(x(i))))**2) be <= s.
      !  the fit is given in the b-spline representation (b-spline coef-
      !  ficients c(j),j=1,2,...n-4) and can be evaluated by means of
      !  subroutine splev.
      !
      !  calling sequence:
      !
      !     call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,
      !    * sx,bind,wrk,lwrk,iwrk,kwrk,ier)
      !
      !  parameters:
      !    iopt: integer flag.
      !          if iopt=0, the routine will start with the minimal number of
      !          knots to guarantee that the convexity conditions will be
      !          satisfied. if iopt=1, the routine will continue with the set
      !          of knots found at the last call of the routine.
      !          attention: a call with iopt=1 must always be immediately
      !          preceded by another call with iopt=1 or iopt=zero
      !          unchanged on exit.
      !    m   : integer. on entry m must specify the number of data points.
      !          m > 3. unchanged on exit.
      !    x   : real array of dimension at least (m). before entry, x(i)
      !          must be set to the i-th value of the independent variable x,
      !          for i=1,2,...,m. these values must be supplied in strictly
      !          ascending order. unchanged on exit.
      !    y   : real array of dimension at least (m). before entry, y(i)
      !          must be set to the i-th value of the dependent variable y,
      !          for i=1,2,...,m. unchanged on exit.
      !    w   : real array of dimension at least (m). before entry, w(i)
      !          must be set to the i-th value in the set of weights. the
      !          w(i) must be strictly positive. unchanged on exit.
      !    v   : real array of dimension at least (m). before entry, v(i)
      !          must be set to 1 if s(x) must be locally concave at x(i),
      !          to (-1) if s(x) must be locally convex at x(i) and to 0
      !          if no convexity constraint is imposed at x(i).
      !    s   : real. on entry s must specify an over-estimate for the
      !          the weighted sum of squared residuals sq of the requested
      !          spline. s >=zero unchanged on exit.
      !   nest : integer. on entry nest must contain an over-estimate of the
      !          total number of knots of the spline returned, to indicate
      !          the storage space available to the routine. nest >=8.
      !          in most practical situation nest=m/2 will be sufficient.
      !          always large enough is  nest=m+4. unchanged on exit.
      !  maxtr : integer. on entry maxtr must contain an over-estimate of the
      !          total number of records in the used tree structure, to indic-
      !          ate the storage space available to the routine. maxtr >=1
      !          in most practical situation maxtr=100 will be sufficient.
      !          always large enough is
      !                         nest-5      nest-6
      !              maxtr =  (       ) + (        )  with l the greatest
      !                           l          l+1
      !          integer <= (nest-6)/2 . unchanged on exit.
      !  maxbin: integer. on entry maxbin must contain an over-estimate of the
      !          number of knots where s(x) will have a zero second derivative
      !          maxbin >=1. in most practical situation maxbin = 10 will be
      !          sufficient. always large enough is maxbin=nest-6.
      !          unchanged on exit.
      !    n   : integer.
      !          on exit with ier <=0, n will contain the total number of
      !          knots of the spline approximation returned. if the comput-
      !          ation mode iopt=1 is used this value of n should be left
      !          unchanged between subsequent calls.
      !    t   : real array of dimension at least (nest).
      !          on exit with ier<=0, this array will contain the knots of the
      !          spline,i.e. the position of the interior knots t(5),t(6),...,
      !          t(n-4) as well as the position of the additional knots
      !          t(1)=t(2)=t(3)=t(4)=x(1) and t(n-3)=t(n-2)=t(n-1)=t(n)=x(m)
      !          needed for the the b-spline representation.
      !          if the computation mode iopt=1 is used, the values of t(1),
      !          t(2),...,t(n) should be left unchanged between subsequent
      !          calls.
      !    c   : real array of dimension at least (nest).
      !          on successful exit, this array will contain the coefficients
      !          c(1),c(2),..,c(n-4) in the b-spline representation of s(x)
      !    sq  : real. unless ier>0 , sq contains the weighted sum of
      !          squared residuals of the spline approximation returned.
      !    sx  : real array of dimension at least m. on exit with ier<=0
      !          this array will contain the spline values s(x(i)),i=1,...,m
      !          if the computation mode iopt=1 is used, the values of sx(1),
      !          sx(2),...,sx(m) should be left unchanged between subsequent
      !          calls.
      !    bind: logical array of dimension at least nest. on exit with ier<=0
      !          this array will indicate the knots where s''(x)=0, i.e.
      !                s''(t(j+3)) == 0 if  bind(j) = .true.
      !                s''(t(j+3)) /= 0 if  bind(j) = .false., j=1,2,...,n-6
      !          if the computation mode iopt=1 is used, the values of bind(1)
      !          ,...,bind(n-6) should be left unchanged between subsequent
      !          calls.
      !   wrk  : real array of dimension at least (m*4+nest*8+maxbin*(maxbin+
      !          nest+1)). used as working space.
      !   lwrk : integer. on entry,lwrk must specify the actual dimension of
      !          the array wrk as declared in the calling (sub)program.lwrk
      !          must not be too small (see wrk). unchanged on exit.
      !   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1))
      !          used as working space.
      !   kwrk : integer. on entry,kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program. kwrk
      !          must not be too small (see iwrk). unchanged on exit.
      !   ier   : integer. error flag
      !      ier=0 : normal return, s(x) satisfies the concavity/convexity
      !              constraints and sq <= s.
      !      ier<0 : abnormal termination: s(x) satisfies the concavity/
      !              convexity constraints but sq > s.
      !        ier=-3 : the requested storage space exceeds the available
      !                 storage space as specified by the parameter nest.
      !                 probably causes: nest too small. if nest is already
      !                 large (say nest > m/2), it may also indicate that s
      !                 is too small.
      !                 the approximation returned is the least-squares cubic
      !                 spline according to the knots t(1),...,t(n) (n=nest)
      !                 which satisfies the convexity constraints.
      !        ier=-2 : the maximal number of knots n=m+4 has been reached.
      !                 probably causes: s too small.
      !        ier=-1 : the number of knots n is less than the maximal number
      !                 m+4 but concon finds that adding one or more knots
      !                 will not further reduce the value of sq.
      !                 probably causes : s too small.
      !      ier>0 : abnormal termination: no approximation is returned
      !        ier=1  : the number of knots where s''(x)=0 exceeds maxbin.
      !                 probably causes : maxbin too small.
      !        ier=2  : the number of records in the tree structure exceeds
      !                 maxtr.
      !                 probably causes : maxtr too small.
      !        ier=3  : the algorithm finds no solution to the posed quadratic
      !                 programming problem.
      !                 probably causes : rounding errors.
      !        ier=4  : the minimum number of knots (given by n) to guarantee
      !                 that the concavity/convexity conditions will be
      !                 satisfied is greater than nest.
      !                 probably causes: nest too small.
      !        ier=5  : the minimum number of knots (given by n) to guarantee
      !                 that the concavity/convexity conditions will be
      !                 satisfied is greater than m+4.
      !                 probably causes: strongly alternating convexity and
      !                 concavity conditions. normally the situation can be
      !                 coped with by adding n-m-4 extra data points (found
      !                 by linear interpolation e.g.) with a small weight w(i)
      !                 and a v(i) number equal to zero.
      !        ier=10 : on entry, the input data are controlled on validity.
      !                 the following restrictions must be satisfied
      !                   0<=iopt<=1, m>3, nest>=8, s>=0, maxtr>=1, maxbin>=1,
      !                   kwrk>=maxtr*4+2*(maxbin+1), w(i)>0, x(i) < x(i+1),
      !                   lwrk>=m*4+nest*8+maxbin*(maxbin+nest+1)
      !                 if one of these restrictions is found to be violated
      !                 control is immediately repassed to the calling program
      !
      !  further comments:
      !    as an example of the use of the computation mode iopt=1, the
      !    following program segment will cause concon to return control
      !    each time a spline with a new set of knots has been computed.
      !     .............
      !     iopt = 0
      !     s = 0.1e+60  (s very large)
      !     do 10 i=1,m
      !       call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,
      !    *  bind,wrk,lwrk,iwrk,kwrk,ier)
      !       ......
      !       s = sq
      !       iopt=1
      ! 10  continue
      !     .............
      !
      !  other subroutines required:
      !    fpcoco,fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno
      !
      !  references:
      !   dierckx p. : an algorithm for cubic spline fitting with convexity
      !                constraints, computing 24 (1980) 349-371.
      !   dierckx p. : an algorithm for least-squares cubic spline fitting
      !                with convexity and concavity constraints, report tw39,
      !                dept. computer science, k.u.leuven, 1978.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !   p. dierckx
      !   dept. computer science, k.u.leuven
      !   celestijnenlaan 200a, b-3001 heverlee, belgium.
      !   e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : march 1978
      !  latest update : march 1987.
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) s,sq
      integer iopt,m,nest,maxtr,maxbin,n,lwrk,kwrk,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),w(m),v(m),t(nest),c(nest),sx(m),wrk(lwrk)
      integer iwrk(kwrk)
      logical bind(nest)
      !  ..local scalars..
      integer :: lwest,kwest,ie,iw,lww

      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if(iopt<0 .or. iopt>1) go to 30
      if(m<4 .or. nest<8) go to 30
      if(s<0.) go to 30
      if(maxtr<1 .or. maxbin<1) go to 30
      lwest = 8*nest+m*4+maxbin*(1+nest+maxbin)
      kwest = 4*maxtr+2*(maxbin+1)
      if(lwrk<lwest .or. kwrk<kwest) go to 30
      if(iopt>0) go to 20

      ! Zero weights
      if (any(w<=zero)) goto 30

      ! Non-monotonic x
      if (any(x(2:m)<=x(1:m-1))) goto 30

      v = sign(one,v)

  20  ier = 0
      !  we partition the working space and determine the spline approximation
      ie = 1
      iw = ie+nest
      lww = lwrk-nest
      call fpcoco(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx, &
       bind,wrk(ie),wrk(iw),lww,iwrk,kwrk,ier)
  30  return
      end subroutine concon


      recursive subroutine concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb, &
       ie,de,ne,k,s,nest,n,t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)

      !  given the ordered set of m points x(i) in the idim-dimensional space
      !  and given also a corresponding set of strictly increasing values u(i)
      !  and the set of positive numbers w(i),i=1,2,...,m, subroutine concur
      !  determines a smooth approximating spline curve s(u), i.e.
      !      x1 = s1(u)
      !      x2 = s2(u)      ub = u(1) <= u <= u(m) = ue
      !      .........
      !      xidim = sidim(u)
      !  with sj(u),j=1,2,...,idim spline functions of odd degree k with
      !  common knots t(j),j=1,2,...,n.
      !  in addition these splines will satisfy the following boundary
      !  constraints        (l)
      !      if ib > 0 :  sj   (u(1)) = db(idim*l+j) ,l=0,1,...,ib-1
      !  and                (l)
      !      if ie > 0 :  sj   (u(m)) = de(idim*l+j) ,l=0,1,...,ie-1.
      !  if iopt=-1 concur calculates the weighted least-squares spline curve
      !  according to a given set of knots.
      !  if iopt>=0 the number of knots of the splines sj(u) and the position
      !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
      !  ness of s(u) is then achieved by minimalizing the discontinuity
      !  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
      !  n-k-1. the amount of smoothness is determined by the condition that
      !  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
      !  negative constant, called the smoothing factor.
      !  the fit s(u) is given in the b-spline representation and can be
      !  evaluated by means of subroutine curev.
      !
      !  calling sequence:
      !     call concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,nest,n,
      !    * t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
      !
      !  parameters:
      !   iopt  : integer flag. on entry iopt must specify whether a weighted
      !           least-squares spline curve (iopt=-1) or a smoothing spline
      !           curve (iopt=0 or 1) must be determined.if iopt=0 the routine
      !           will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
      !           i=1,2,...,k+1. if iopt=1 the routine will continue with the
      !           knots found at the last call of the routine.
      !           attention: a call with iopt=1 must always be immediately
      !           preceded by another call with iopt=1 or iopt=zero
      !           unchanged on exit.
      !   idim  : integer. on entry idim must specify the dimension of the
      !           curve. 0 < idim < 11.
      !           unchanged on exit.
      !   m     : integer. on entry m must specify the number of data points.
      !           m > k-max(ib-1,0)-max(ie-1,0). unchanged on exit.
      !   u     : real array of dimension at least (m). before entry,
      !           u(i) must be set to the i-th value of the parameter variable
      !           u for i=1,2,...,m. these values must be supplied in
      !           strictly ascending order and will be unchanged on exit.
      !   mx    : integer. on entry mx must specify the actual dimension of
      !           the arrays x and xx as declared in the calling (sub)program
      !           mx must not be too small (see x). unchanged on exit.
      !   x     : real array of dimension at least idim*m.
      !           before entry, x(idim*(i-1)+j) must contain the j-th coord-
      !           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
      !           idim. unchanged on exit.
      !   xx    : real array of dimension at least idim*m.
      !           used as working space. on exit xx contains the coordinates
      !           of the data points to which a spline curve with zero deriv-
      !           ative constraints has been determined.
      !           if the computation mode iopt =1 is used xx should be left
      !           unchanged between calls.
      !   w     : real array of dimension at least (m). before entry, w(i)
      !           must be set to the i-th value in the set of weights. the
      !           w(i) must be strictly positive. unchanged on exit.
      !           see also further comments.
      !   ib    : integer. on entry ib must specify the number of derivative
      !           constraints for the curve at the begin point. 0<=ib<=(k+1)/2
      !           unchanged on exit.
      !   db    : real array of dimension nb. before entry db(idim*l+j) must
      !           contain the l-th order derivative of sj(u) at u=u(1) for
      !           j=1,2,...,idim and l=0,1,...,ib-1 (if ib>0).
      !           unchanged on exit.
      !   nb    : integer, specifying the dimension of db. nb>=max(1,idim*ib)
      !           unchanged on exit.
      !   ie    : integer. on entry ie must specify the number of derivative
      !           constraints for the curve at the end point. 0<=ie<=(k+1)/2
      !           unchanged on exit.
      !   de    : real array of dimension ne. before entry de(idim*l+j) must
      !           contain the l-th order derivative of sj(u) at u=u(m) for
      !           j=1,2,...,idim and l=0,1,...,ie-1 (if ie>0).
      !           unchanged on exit.
      !   ne    : integer, specifying the dimension of de. ne>=max(1,idim*ie)
      !           unchanged on exit.
      !   k     : integer. on entry k must specify the degree of the splines.
      !           k=1,3 or 5.
      !           unchanged on exit.
      !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
      !           factor. s >=zero unchanged on exit.
      !           for advice on the choice of s see further comments.
      !   nest  : integer. on entry nest must contain an over-estimate of the
      !           total number of knots of the splines returned, to indicate
      !           the storage space available to the routine. nest >=2*k+2.
      !           in most practical situation nest=m/2 will be sufficient.
      !           always large enough is nest=m+k+1+max(0,ib-1)+max(0,ie-1),
      !           the number of knots needed for interpolation (s=0).
      !           unchanged on exit.
      !   n     : integer.
      !           unless ier = 10 (in case iopt >=0), n will contain the
      !           total number of knots of the smoothing spline curve returned
      !           if the computation mode iopt=1 is used this value of n
      !           should be left unchanged between subsequent calls.
      !           in case iopt=-1, the value of n must be specified on entry.
      !   t     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the knots of the
      !           spline curve,i.e. the position of the interior knots t(k+2),
      !           t(k+3),..,t(n-k-1) as well as the position of the additional
      !           t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
      !           the b-spline representation.
      !           if the computation mode iopt=1 is used, the values of t(1),
      !           t(2),...,t(n) should be left unchanged between subsequent
      !           calls. if the computation mode iopt=-1 is used, the values
      !           t(k+2),...,t(n-k-1) must be supplied by the user, before
      !           entry. see also the restrictions (ier=10).
      !   nc    : integer. on entry nc must specify the actual dimension of
      !           the array c as declared in the calling (sub)program. nc
      !           must not be too small (see c). unchanged on exit.
      !   c     : real array of dimension at least (nest*idim).
      !           on successful exit, this array will contain the coefficients
      !           in the b-spline representation of the spline curve s(u),i.e.
      !           the b-spline coefficients of the spline sj(u) will be given
      !           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
      !   cp    : real array of dimension at least 2*(k+1)*idim.
      !           on exit cp will contain the b-spline coefficients of a
      !           polynomial curve which satisfies the boundary constraints.
      !           if the computation mode iopt =1 is used cp should be left
      !           unchanged between calls.
      !   np    : integer. on entry np must specify the actual dimension of
      !           the array cp as declared in the calling (sub)program. np
      !           must not be too small (see cp). unchanged on exit.
      !   fp    : real. unless ier = 10, fp contains the weighted sum of
      !           squared residuals of the spline curve returned.
      !   wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
      !           used as working space. if the computation mode iopt=1 is
      !           used, the values wrk(1),...,wrk(n) should be left unchanged
      !           between subsequent calls.
      !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
      !           the array wrk as declared in the calling (sub)program. lwrk
      !           must not be too small (see wrk). unchanged on exit.
      !   iwrk  : integer array of dimension at least (nest).
      !           used as working space. if the computation mode iopt=1 is
      !           used,the values iwrk(1),...,iwrk(n) should be left unchanged
      !           between subsequent calls.
      !   ier   : integer. unless the routine detects an error, ier contains a
      !           non-positive value on exit, i.e.
      !    ier=0  : normal return. the curve returned has a residual sum of
      !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !             ive tolerance set to 0.001 by the program.
      !    ier=-1 : normal return. the curve returned is an interpolating
      !             spline curve, satisfying the constraints (fp=0).
      !    ier=-2 : normal return. the curve returned is the weighted least-
      !             squares polynomial curve of degree k, satisfying the
      !             constraints. in this extreme case fp gives the upper
      !             bound fp0 for the smoothing factor s.
      !    ier=1  : error. the required storage space exceeds the available
      !             storage space, as specified by the parameter nest.
      !             probably causes : nest too small. if nest is already
      !             large (say nest > m/2), it may also indicate that s is
      !             too small
      !             the approximation returned is the least-squares spline
      !             curve according to the knots t(1),t(2),...,t(n). (n=nest)
      !             the parameter fp gives the corresponding weighted sum of
      !             squared residuals (fp>s).
      !    ier=2  : error. a theoretically impossible result was found during
      !             the iteration process for finding a smoothing spline curve
      !             with fp = s. probably causes : s too small.
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=3  : error. the maximal number of iterations maxit (set to 20
      !             by the program) allowed for finding a smoothing curve
      !             with fp=s has been reached. probably causes : s too small
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=10 : error. on entry, the input data are controlled on validity
      !             the following restrictions must be satisfied.
      !             -1<=iopt<=1, k = 1,3 or 5, m>k-max(0,ib-1)-max(0,ie-1),
      !             nest>=2k+2, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
      !             nc >=nest*idim ,u(1)<u(2)<...<u(m),w(i)>0 i=1,2,...,m,
      !             mx>=idim*m,0<=ib<=(k+1)/2,0<=ie<=(k+1)/2,nb>=1,ne>=1,
      !             nb>=ib*idim,ne>=ib*idim,np>=2*(k+1)*idim,
      !             if iopt=-1:2*k+2<=n<=min(nest,mmax) with mmax = m+k+1+
      !                        max(0,ib-1)+max(0,ie-1)
      !                        u(1)<t(k+2)<t(k+3)<...<t(n-k-1)<u(m)
      !                       the schoenberg-whitney conditions, i.e. there
      !                       must be a subset of data points uu(j) such that
      !                         t(j) < uu(j) < t(j+k+1), j=1+max(0,ib-1),...
      !                                                   ,n+k-1-max(0,ie-1)
      !             if iopt>=0: s>=0
      !                         if s=0 : nest >=mmax (see above)
      !             if one of these conditions is found to be violated,control
      !             is immediately repassed to the calling program. in that
      !             case there is no approximation returned.
      !
      !  further comments:
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the curve will be too smooth and signal will be
      !   lost ; if s is too small the curve will pick up too much noise. in
      !   the extreme cases the program will return an interpolating curve if
      !   s=0 and the least-squares polynomial curve of degree k if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the weights w(i). if these are
      !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
      !   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
      !   sqrt(2*m)). if nothing is known about the statistical error in x(i)
      !   each w(i) can be set equal to one and s determined by trial and
      !   error, taking account of the comments above. the best is then to
      !   start with a very large value of s ( to determine the least-squares
      !   polynomial curve and the upper bound fp0 for s) and then to
      !   progressively decrease the value of s ( say by a factor 10 in the
      !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
      !   approximating curve shows more detail) to obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=zero
      !   if iopt=1 the program will continue with the set of knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if concur is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   curve underlying the data. but, if the computation mode iopt=1 is
      !   used, the knots returned may also depend on the s-values at previous
      !   calls (if these were smaller). therefore, if after a number of
      !   trials with different s-values and iopt=1, the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   concur once more with the selected value for s but now with iopt=zero
      !   indeed, concur may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !
      !   the form of the approximating curve can strongly be affected by
      !   the choice of the parameter values u(i). if there is no physical
      !   reason for choosing a particular parameter u, often good results
      !   will be obtained with the choice
      !        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
      !   where
      !        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
      !   other possibilities for q(i) are
      !        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
      !        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
      !        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
      !        q(i)= 1
      !
      !  other subroutines required:
      !    fpback,fpbspl,fpched,fpcons,fpdisc,fpgivs,fpknot,fprati,fprota
      !    curev,fppocu,fpadpo,fpinst
      !
      !  references:
      !   dierckx p. : algorithms for smoothing data with periodic and
      !                parametric splines, computer graphics and image
      !                processing 20 (1982) 171-184.
      !   dierckx p. : algorithms for smoothing data with periodic and param-
      !                etric splines, report tw55, dept. computer science,
      !                k.u.leuven, 1981.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !  latest update : march 1987
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) s,fp
      integer iopt,idim,m,mx,ib,nb,ie,ne,k,nest,n,nc,np,lwrk,ier
      !  ..array arguments..
      real(RKIND) u(m),x(mx),xx(mx),db(nb),de(ne),w(m),t(nest),c(nc),wrk(lwrk &
      )
      real(RKIND) cp(np)
      integer iwrk(nest)
      !  ..local scalars..
      real(RKIND) tol
      integer i,ib1,ie1,ja,jb,jfp,jg,jq,jz,j,k1,k2,lwest,maxit,nmin, &
       ncc,kk,mmin,nmax,mxx
      !  ..
      !  we set up the parameters tol and maxit
      maxit = 20
      tol = smallnum03
      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if(iopt<(-1) .or. iopt>1) go to 90
      if(idim<=0 .or. idim>10) go to 90
      if(k<=0 .or. k>5) go to 90
      k1 = k+1
      kk = k1/2
      if(kk*2/=k1) go to 90
      k2 = k1+1
      if(ib<0 .or. ib>kk) go to 90
      if(ie<0 .or. ie>kk) go to 90
      nmin = 2*k1
      ib1 = max0(0,ib-1)
      ie1 = max0(0,ie-1)
      mmin = k1-ib1-ie1
      if(m<mmin .or. nest<nmin) go to 90
      if(nb<(idim*ib) .or. ne<(idim*ie)) go to 90
      if(np<(2*k1*idim)) go to 90
      mxx = m*idim
      ncc = nest*idim
      if(mx<mxx .or. nc<ncc) go to 90
      lwest = m*k1+nest*(6+idim+3*k)
      if (lwrk<lwest) go to 90
      if (any(w<=zero)) go to 90
      if (any(u(1:m-1)>=u(2:m))) goto 90
      if(iopt>=0) go to 30
      if(n<nmin .or. n>nest) go to 90
      j = n
      do 20 i=1,k1
         t(i) = u(1)
         t(j) = u(m)
         j = j-1
  20  continue
      call fpched(u,m,t,n,k,ib,ie,ier)
      if (ier==FITPACK_OK) go to 40
      go to 90
  30  if(s<zero) go to 90
      nmax = m+k1+ib1+ie1
      if(s==zero .and. nest<nmax) go to 90
      ier = 0
      if(iopt>0) go to 70
      !  we determine a polynomial curve satisfying the boundary constraints.
  40  call fppocu(idim,k,u(1),u(m),ib,db,nb,ie,de,ne,cp,np)
      !  we generate new data points which will be approximated by a spline
      !  with zero derivative constraints.
      wrk(1:k1) = u(1)
      wrk(nmin-k1+1:nmin) = u(m)
      !  evaluate the polynomial curve
      call curev(idim,wrk,nmin,cp,np,k,u,m,xx,mxx,ier)
      !  subtract from the old data, the values of the polynomial curve
      xx(1:mxx) = x(1:mxx)-xx(1:mxx)
      ! we partition the working space and determine the spline curve.
  70  jfp = 1
      jz = jfp+nest
      ja = jz+ncc
      jb = ja+nest*k1
      jg = jb+nest*k2
      jq = jg+nest*k2
      call fpcons(iopt,idim,m,u,mxx,xx,w,ib,ie,k,s,nest,tol,maxit,k1, &
       k2,n,t,ncc,c,fp,wrk(jfp),wrk(jz),wrk(ja),wrk(jb),wrk(jg),wrk(jq),iwrk,ier)
      !  add the polynomial curve to the calculated spline.
      call fpadpo(idim,t,n,c,ncc,k,cp,np,wrk(jz),wrk(ja),wrk(jb))
  90  return
      end subroutine concur


      !  subroutine cualde evaluates at the point u all the derivatives
      !                     (l)
      !     d(idim*l+j) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim
      !  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim
      !  given in its b-spline representation.

      pure subroutine cualde(idim,t,n,c,nc,k1,u,d,nd,ier)

      !
      !  input parameters:
      !    idim : integer, giving the dimension of the spline curve.
      !    t    : array,length n, which contains the position of the knots.
      !    n    : integer, giving the total number of knots of s(u).
      !    c    : array,length nc, which contains the b-spline coefficients.
      !    nc   : integer, giving the total number of coefficients of s(u).
      !    k1   : integer, giving the order of s(u) (order=degree+1).
      !    u    : real, which contains the point where the derivatives must be evaluated.
      !    nd   : integer, giving the dimension of the array d. nd >= k1*idim
      !
      !  output parameters:
      !    d    : array,length nd,giving the different curve derivatives. d(idim*l+j) will contain the
      !           j-th coordinate of the l-th derivative of the curve at the point u.
      !    ier  : error flag
      !      ier = 0 : normal return
      !      ier =10 : invalid input data (see restrictions)
      !
      !  restrictions:
      !    nd >= k1*idim
      !    t(k1) <= u <= t(n-k1+1)
      !
      !  further comments:
      !    if u coincides with a knot, right derivatives are computed
      !    ( left derivatives if u = t(n-k1+1) ).
      !
      !  other subroutines required: fpader.
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory
      !                6 (1972) 50-62.
      !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
      !                applics 10 (1972) 134-149.
      !    dierckx p. : curve and surface fitting with splines, monographs on
      !                 numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1987
      !
      !  ..scalar arguments..
      integer, intent(in)  :: idim,n,nc,k1,nd
      integer, intent(out) :: ier
      real(RKIND), intent(in) :: u
      !  ..array arguments..
      real(RKIND), intent(in) :: t(n),c(nc)
      real(RKIND), intent(out) :: d(nd)
      !  ..local scalars..
      integer :: i,j,kk,l,m,nk1
      !  ..local array..
      real(RKIND) :: h(SIZ_K+1)
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if(nd<(k1*idim)) return
      nk1 = n-k1
      if(u<t(k1) .or. u>t(nk1+1)) return

      !  search for knot interval t(l) <= u < t(l+1)
      l = k1
      do while (.not.(u<t(l+1) .or. l==nk1))
         l = l+1
      end do
      if(t(l)>=t(l+1)) return

      ier = FITPACK_OK

      !  calculate the derivatives.
      j = 1
      do i=1,idim
         call fpader(t,n,c(j),k1,u,l,h)
         m = i
         do kk=1,k1
            d(m) = h(kk)
            m = m+idim
         end do
         j = j+n
      end do
      return
      end subroutine cualde



      !  subroutine curev evaluates in a number of points u(i),i=1,2,...,m a spline curve s(u) of degree k
      !  and dimension idim, given in its b-spline representation.
      pure subroutine curev(idim,t,n,c,nc,k,u,m,x,mx,ier)

      !  calling sequence:
      !     call curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
      !
      !  input parameters:
      !    idim : integer, giving the dimension of the spline curve.
      !    t    : array,length n, which contains the position of the knots.
      !    n    : integer, giving the total number of knots of s(u).
      !    c    : array,length nc, which contains the b-spline coefficients.
      !    nc   : integer, giving the total number of coefficients of s(u).
      !    k    : integer, giving the degree of s(u).
      !    u    : array,length m, which contains the points where s(u) must be evaluated.
      !    m    : integer, giving the number of points where s(u) must be evaluated.
      !    mx   : integer, giving the dimension of the array x. mx >= m*idim
      !
      !  output parameters:
      !    x    : array,length mx,giving the value of s(u) at the different points. x(idim*(i-1)+j) will
      !           contain the j-th coordinate of the i-th point on the curve.
      !    ier  : error flag
      !      ier = 0 : normal return
      !      ier =10 : invalid input data (see restrictions)
      !
      !  restrictions:
      !    m >= 1
      !    mx >= m*idim
      !    t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1.
      !
      !  other subroutines required: fpbspl.
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory 6 (1972) 50-62.
      !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths applics 10 (1972) 134-149.
      !    dierckx p. : curve and surface fitting with splines, monographs on numerical analysis, oxford
      !                 university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer, intent(in) :: idim,n,nc,k,m,mx
      integer, intent(out) :: ier
      !  ..array arguments..
      real(RKIND), intent(in) :: t(n),c(nc),u(m)
      real(RKIND), intent(out) :: x(idim,m) ! x has size (mx), assume 2d (idim,m)
      !  ..local scalars..
      integer :: i,j1,k1,l,ll,l1,nk1
      real(RKIND) :: arg,tb,te
      !  ..local array..
      real(RKIND) h(SIZ_K+1)
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (m<1) return

      ! Check monotonic
      if (m>1 .and. any(u(2:m)<u(1:m-1))) return

      ! Check enough output space
      if (mx<(m*idim)) return

      ier = FITPACK_OK

      !  fetch tb and te, the boundaries of the approximation interval.
      k1  = k+1
      nk1 = n-k1
      tb  = t(k1)
      te  = t(nk1+1)
      l   = k1
      l1  = l+1
      !  main loop for the different points.
      eval_points: do i=1,m

        ! fetch a new u-value arg.
        arg = min(max(u(i),tb),te)

        ! search for knot interval t(l) <= arg < t(l+1)
        do while (.not.(arg<t(l1) .or. l==nk1))
          l = l1
          l1 = l+1
        end do

        ! evaluate the non-zero b-splines at arg.
        call fpbspl(t,n,k,arg,l,h)

        ! find the value of s(u) at u=arg.
        ll = l-k1
        do j1=1,idim
          x(idim,m) = dot_product(h(1:k1),c(ll+1:ll+k1))
          ll = ll+n
        end do
      end do eval_points

      return
      end subroutine curev


      recursive subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)

      !  given the set of data points (x(i),y(i)) and the set of positive
      !  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
      !  approximation of degree k on the interval xb <= x <= xe.
      !  if iopt=-1 curfit calculates the weighted least-squares spline
      !  according to a given set of knots.
      !  if iopt>=0 the number of knots of the spline s(x) and the position
      !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
      !  ness of s(x) is then achieved by minimalizing the discontinuity
      !  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
      !  n-k-1. the amount of smoothness is determined by the condition that
      !  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
      !  negative constant, called the smoothing factor.
      !  the fit s(x) is given in the b-spline representation (b-spline coef-
      !  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
      !  subroutine splev.
      !
      !  calling sequence:
      !     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
      !    * lwrk,iwrk,ier)
      !
      !  parameters:
      !   iopt  : integer flag. on entry iopt must specify whether a weighted
      !           least-squares spline (iopt=-1) or a smoothing spline (iopt=
      !           0 or 1) must be determined. if iopt=0 the routine will start
      !           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
      !           k+1. if iopt=1 the routine will continue with the knots
      !           found at the last call of the routine.
      !           attention: a call with iopt=1 must always be immediately
      !           preceded by another call with iopt=1 or iopt=zero
      !           unchanged on exit.
      !   m     : integer. on entry m must specify the number of data points.
      !           m > k. unchanged on exit.
      !   x     : real array of dimension at least (m). before entry, x(i)
      !           must be set to the i-th value of the independent variable x,
      !           for i=1,2,...,m. these values must be supplied in strictly
      !           ascending order. unchanged on exit.
      !   y     : real array of dimension at least (m). before entry, y(i)
      !           must be set to the i-th value of the dependent variable y,
      !           for i=1,2,...,m. unchanged on exit.
      !   w     : real array of dimension at least (m). before entry, w(i)
      !           must be set to the i-th value in the set of weights. the
      !           w(i) must be strictly positive. unchanged on exit.
      !           see also further comments.
      !   xb,xe : real values. on entry xb and xe must specify the boundaries
      !           of the approximation interval. xb<=x(1), xe>=x(m).
      !           unchanged on exit.
      !   k     : integer. on entry k must specify the degree of the spline.
      !           1<=k<=5. it is recommended to use cubic splines (k=3).
      !           the user is strongly dissuaded from choosing k even,together
      !           with a small s-value. unchanged on exit.
      !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
      !           factor. s >=zero unchanged on exit.
      !           for advice on the choice of s see further comments.
      !   nest  : integer. on entry nest must contain an over-estimate of the
      !           total number of knots of the spline returned, to indicate
      !           the storage space available to the routine. nest >=2*k+2.
      !           in most practical situation nest=m/2 will be sufficient.
      !           always large enough is  nest=m+k+1, the number of knots
      !           needed for interpolation (s=0). unchanged on exit.
      !   n     : integer.
      !           unless ier =10 (in case iopt >=0), n will contain the
      !           total number of knots of the spline approximation returned.
      !           if the computation mode iopt=1 is used this value of n
      !           should be left unchanged between subsequent calls.
      !           in case iopt=-1, the value of n must be specified on entry.
      !   t     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the knots of the
      !           spline,i.e. the position of the interior knots t(k+2),t(k+3)
      !           ...,t(n-k-1) as well as the position of the additional knots
      !           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
      !           the b-spline representation.
      !           if the computation mode iopt=1 is used, the values of t(1),
      !           t(2),...,t(n) should be left unchanged between subsequent
      !           calls. if the computation mode iopt=-1 is used, the values
      !           t(k+2),...,t(n-k-1) must be supplied by the user, before
      !           entry. see also the restrictions (ier=10).
      !   c     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the coefficients
      !           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
      !   fp    : real. unless ier=10, fp contains the weighted sum of
      !           squared residuals of the spline approximation returned.
      !   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
      !           used as working space. if the computation mode iopt=1 is
      !           used, the values wrk(1),...,wrk(n) should be left unchanged
      !           between subsequent calls.
      !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
      !           the array wrk as declared in the calling (sub)program.lwrk
      !           must not be too small (see wrk). unchanged on exit.
      !   iwrk  : integer array of dimension at least (nest).
      !           used as working space. if the computation mode iopt=1 is
      !           used,the values iwrk(1),...,iwrk(n) should be left unchanged
      !           between subsequent calls.
      !   ier   : integer. unless the routine detects an error, ier contains a
      !           non-positive value on exit, i.e.
      !    ier=0  : normal return. the spline returned has a residual sum of
      !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !             ive tolerance set to 0.001 by the program.
      !    ier=-1 : normal return. the spline returned is an interpolating
      !             spline (fp=0).
      !    ier=-2 : normal return. the spline returned is the weighted least-
      !             squares polynomial of degree k. in this extreme case fp
      !             gives the upper bound fp0 for the smoothing factor s.
      !    ier=1  : error. the required storage space exceeds the available
      !             storage space, as specified by the parameter nest.
      !             probably causes : nest too small. if nest is already
      !             large (say nest > m/2), it may also indicate that s is
      !             too small
      !             the approximation returned is the weighted least-squares
      !             spline according to the knots t(1),t(2),...,t(n). (n=nest)
      !             the parameter fp gives the corresponding weighted sum of
      !             squared residuals (fp>s).
      !    ier=2  : error. a theoretically impossible result was found during
      !             the iteration process for finding a smoothing spline with
      !             fp = s. probably causes : s too small.
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=3  : error. the maximal number of iterations maxit (set to 20
      !             by the program) allowed for finding a smoothing spline
      !             with fp=s has been reached. probably causes : s too small
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=10 : error. on entry, the input data are controlled on validity
      !             the following restrictions must be satisfied.
      !             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
      !             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
      !             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
      !                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
      !                       the schoenberg-whitney conditions, i.e. there
      !                       must be a subset of data points xx(j) such that
      !                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
      !             if iopt>=0: s>=0
      !                         if s=0 : nest >= m+k+1
      !             if one of these conditions is found to be violated,control
      !             is immediately repassed to the calling program. in that
      !             case there is no approximation returned.
      !
      !  further comments:
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating spline if
      !   s=0 and the weighted least-squares polynomial of degree k if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the weights w(i). if these are
      !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
      !   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
      !   sqrt(2*m)). if nothing is known about the statistical error in y(i)
      !   each w(i) can be set equal to one and s determined by trial and
      !   error, taking account of the comments above. the best is then to
      !   start with a very large value of s ( to determine the least-squares
      !   polynomial and the corresponding upper bound fp0 for s) and then to
      !   progressively decrease the value of s ( say by a factor 10 in the
      !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
      !   approximation shows more detail) to obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=zero
      !   if iopt=1 the program will continue with the set of knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if curfit is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. but, if the computation mode iopt=1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt=1, the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   curfit once more with the selected value for s but now with iopt=0.
      !   indeed, curfit may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !
      !  other subroutines required:
      !    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
      !
      !  references:
      !   dierckx p. : an algorithm for smoothing, differentiation and integ-
      !                ration of experimental data using spline functions,
      !                j.comp.appl.maths 1 (1975) 165-184.
      !   dierckx p. : a fast algorithm for smoothing data on a rectangular
      !                grid while using spline functions, siam j.numer.anal.
      !                19 (1982) 1286-1304.
      !   dierckx p. : an improved algorithm for curve fitting with spline
      !                functions, report tw54, dept. computer science,k.u.
      !                leuven, 1981.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !  latest update : march 1987
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) xb,xe,s,fp
      integer iopt,m,k,nest,n,lwrk,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
      !  ..local scalars..
      integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,nmin
      !  ..
      !  we set up the parameters tol and maxit
      real(RKIND), parameter :: tol = smallnum03
      integer    , parameter :: maxit = 20

      k1   = k+1
      k2   = k1+1
      nmin = 2*k1

      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (k<=0 .or. k>5)         return
      if (iopt<(-1) .or. iopt>1) return
      if (m<k1 .or. nest<nmin)   return
      lwest = m*k1+nest*(7+3*k)
      if (lwrk<lwest)             return
      if (xb>x(1) .or. xe<x(m))  return
      if (any(x(1:m-1)>x(2:m)))  return

      if (iopt>=0) then
          if (s<zero .or. (s==zero .and. nest<(m+k1))) return
      else
          if (n<nmin .or. n>nest) return
          j = n
          do i=1,k1
             t(i) = xb
             t(j) = xe
             j = j-1
          end do
          call fpchec(x,m,t,n,k,ier); if (ier/=0) return
      endif

      ier = FITPACK_OK

      ! we partition the working space and determine the spline approximation.
      ifp = 1
      iz = ifp+nest
      ia = iz+nest
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp, &
                  wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)

      end subroutine curfit


      !  function dblint calculates the double integral
      !         / xe  / ye
      !        |     |      s(x,y) dx dy
      !    xb /  yb /
      !  with s(x,y) a bivariate spline of degrees kx and ky, given in the b-spline representation.
      real(RKIND) function dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk) result(dblint_res)

      !
      !  calling sequence:
      !     aint = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
      !
      !  input parameters:
      !   tx    : real array, length nx, which contains the position of the knots in the x-direction.
      !   nx    : integer, giving the total number of knots in the x-direction
      !   ty    : real array, length ny, which contains the position of the knots in the y-direction.
      !   ny    : integer, giving the total number of knots in the y-direction
      !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the b-spline coefficients.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   xb,xe : real values, containing the boundaries of the integration
      !   yb,ye   domain. s(x,y) is considered to be identically zero outside the rectangle
      !           (tx(kx+1),tx(nx-kx))*(ty(ky+1),ty(ny-ky))
      !
      !  output parameters:
      !   aint  : real , containing the double integral of s(x,y).
      !   wrk   : real array of dimension at least (nx+ny-kx-ky-2). used as working space.
      !           on exit, wrk(i) will contain the integral
      !                / xe
      !               | ni,kx+1(x) dx , i=1,2,...,nx-kx-1
      !           xb /
      !           with ni,kx+1(x) the normalized b-spline defined on the knots tx(i),...,tx(i+kx+1)
      !           wrk(j+nx-kx-1) will contain the integral
      !                / ye
      !               | nj,ky+1(y) dy , j=1,2,...,ny-ky-1
      !           yb /
      !           with nj,ky+1(y) the normalized b-spline defined on the knots ty(j),...,ty(j+ky+1)
      !
      !  other subroutines required: fpintb
      !
      !  references :
      !    gaffney p.w. : the calculation of indefinite integrals of b-splines
      !                   j. inst. maths applics 17 (1976) 37-41.
      !    dierckx p. : curve and surface fitting with splines, monographs on
      !                 numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer, intent(in) :: nx,ny,kx,ky
      real(RKIND), intent(in) :: xb,xe,yb,ye
      !  ..array arguments..
      real(RKIND), intent(in) :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1))
      real(RKIND), intent(out) :: wrk(nx+ny-kx-ky-2)
      !  ..local scalars..
      integer :: i,j,l,m,nkx1,nky1
      real(RKIND) :: res

      !  ..
      nkx1 = nx-kx-1
      nky1 = ny-ky-1

      !  we calculate the integrals of the normalized b-splines ni,kx+1(x)
      call fpintb(tx,nx,wrk,nkx1,xb,xe)

      !  we calculate the integrals of the normalized b-splines nj,ky+1(y)
      call fpintb(ty,ny,wrk(nkx1+1),nky1,yb,ye)

      !  calculate the integral of s(x,y)
      dblint_res = zero
      x_dim: do i=1,nkx1
        res = wrk(i)
        if (res==zero) cycle x_dim
        m = (i-1)*nky1
        l = nkx1
        y_dim: do j=1,nky1
          m = m+1
          l = l+1
          dblint_res = dblint_res + res*wrk(l)*c(m)
        end do y_dim
      end do x_dim
      return
      end function dblint

      !  function program evapol evaluates the function f(x,y) = s(u,v), defined through the transformation
      !      x = u*rad(v)*cos(v)    y = u*rad(v)*sin(v)
      !  and where s(u,v) is a bicubic spline ( 0<=u<=1 , -pi<=v<=pi ), given in its standard b-spline
      !  representation.

      pure real(RKIND) function evapol(tu,nu,tv,nv,c,rad,x,y) result(e_res)

      !  calling sequence:
      !     f = evapol(tu,nu,tv,nv,c,rad,x,y)
      !
      !  input parameters:
      !   tu    : real array, length nu, which contains the position of the knots in the u-direction.
      !   nu    : integer, giving the total number of knots in the u-direction
      !   tv    : real array, length nv, which contains the position of the knots in the v-direction.
      !   nv    : integer, giving the total number of knots in the v-direction
      !   c     : real array, length (nu-4)*(nv-4), which contains the b-spline coefficients.
      !   rad   : real function subprogram, defining the boundary of the approximation domain. must be
      !           declared external in the calling (sub)-program
      !   x,y   : the co-ordinates of the point where f(x,y) must be evaluated.
      !
      !  output parameter:
      !   e_res : the value of f(x,y)
      !
      !  other subroutines required:
      !    bispev,fpbisp,fpbspl
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory 6 (1972) 50-62.
      !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths applics 10 (1972) 134-149.
      !    dierckx p. : curve and surface fitting with splines, monographs on numerical analysis, oxford
      !                 university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer, intent(in) :: nu,nv
      real(RKIND), intent(in) :: x,y
      !  ..array arguments..
      real(RKIND), intent(in) :: tu(nu),tv(nv),c((nu-4)*(nv-4))
      !  ..user specified function
      procedure(boundary) :: rad
      !  ..local scalars..
      integer :: ier
      integer, parameter :: liwrk = 2, lwrk = 8
      real(RKIND) :: u(1),v(1),r,f(1),dist
      !  ..local arrays
      real(RKIND) :: wrk(lwrk)
      integer :: iwrk(liwrk)
      !  ..
      !  calculate the (u,v)-coordinates of the given point.
      u    = zero
      v    = zero
      dist = x**2+y**2
      if (dist>zero) then
         v(1) = atan2(y,x)
         r    = rad(v(1))
         if (r>zero) &
         u(1) = min(sqrt(dist)/r,one)
      endif
      ! evaluate s(u,v)
      call bispev(tu,nu,tv,nv,c,3,3,u,1,v,1,f,wrk,lwrk,iwrk,liwrk,ier)

      ! Return scalar result
      e_res = f(1)
      return
      end function evapol

      !  subroutine fourco calculates the integrals
      !                    /t(n-3)
      !    ress(i) =      !        s(x)*sin(alfa(i)*x) dx    and
      !              t(4)/
      !                    /t(n-3)
      !    resc(i) =      !        s(x)*cos(alfa(i)*x) dx, i=1,...,m,
      !              t(4)/
      !  where s(x) denotes a cubic spline which is given in its b-spline representation.
      pure subroutine fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)

      !  calling sequence:
      !     call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
      !
      !  input parameters:
      !    t    : real array,length n, containing the knots of s(x).
      !    n    : integer, containing the total number of knots. n>=10.
      !    c    : real array,length n, containing the b-spline coefficients.
      !    alfa : real array,length m, containing the parameters alfa(i).
      !    m    : integer, specifying the number of integrals to be computed.
      !    wrk1 : real array,length n. used as working space
      !    wrk2 : real array,length n. used as working space
      !
      !  output parameters:
      !    ress : real array,length m, containing the integrals ress(i).
      !    resc : real array,length m, containing the integrals resc(i).
      !    ier  : error flag:
      !      ier=0 : normal return.
      !      ier=10: invalid input data (see restrictions).
      !
      !  restrictions:
      !    n >= 10
      !    t(4) < t(5) < ... < t(n-4) < t(n-3).
      !    t(1) <= t(2) <= t(3) <= t(4).
      !    t(n-3) <= t(n-2) <= t(n-1) <= t(n).
      !
      !  other subroutines required: fpbfou,fpcsin
      !
      !  references :
      !    dierckx p. : calculation of fourier coefficients of discrete functions using cubic splines.
      !                 j. computational and applied mathematics 3 (1977) 207-209.
      !    dierckx p. : curve and surface fitting with splines, monographs on
      !                 numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer, intent(in) :: n,m
      integer, intent(out) :: ier
      !  ..array arguments..
      real(RKIND), intent(in) :: t(n),c(n),alfa(m)
      real(RKIND), intent(inout) :: wrk1(n),wrk2(n)
      real(RKIND), intent(out) :: ress(m),resc(m)
      !  ..local scalars..
      integer :: i,n4
      !  ..
      n4 = n-4
      !  before starting computations a data check is made. in the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR

      ! Not enough points
      if (n<10) return

      ! Ends of the support: knots must be monotonic
      if (any(t(1:3)>t(2:4)))       return
      if (any(t(n-2:n)<t(n-3:n-1))) return

      ! Interior: knots must be strictly monotonic
      if (any(t(4:n4)>=t(5:n-3)))   return

      ier = FITPACK_OK
      !  main loop for the different alfa(i).
      alphas: do i=1,m

         !  calculate the integrals
         !    wrk1(j) = integral(nj,4(x)*sin(alfa*x))    and
         !    wrk2(j) = integral(nj,4(x)*cos(alfa*x)),  j=1,2,...,n-4,
         !  where nj,4(x) denotes the normalised cubic b-spline defined on the knots t(j),t(j+1),...,t(j+4).
         call fpbfou(t,n,alfa(i),wrk1,wrk2)

         ! calculate the integrals ress(i) and resc(i).
         ress(i) = dot_product(c(1:n4),wrk1(1:n4))
         resc(i) = dot_product(c(1:n4),wrk2(1:n4))
      end do alphas

      end subroutine fourco


      !  subroutine fpader calculates the derivatives
      !             (j-1)
      !     d(j) = s     (x) , j=1,2,...,k1
      !  of a spline of order k1 at the point t(l)<=x<t(l+1), using the
      !  stable recurrence scheme of de boor
      !  ..
      pure subroutine fpader(t,n,c,k1,x,l,d)
      !  ..scalar arguments..
      real(RKIND), intent(in)  :: x
      integer, intent(in)      :: n,k1,l
      !  ..array arguments..
      real(RKIND), intent(in)  :: t(n),c(n)
      real(RKIND), intent(out) :: d(k1)

      !  ..local scalars..
      integer :: i,ik,j,jj,j1,j2,ki,kj,li,lj,lk
      real(RKIND) :: ak,fac
      !  ..local array..
      real(RKIND) :: h(20)
      !  ..
      lk = l-k1
      do i=1,k1
        ik = i+lk
        h(i) = c(ik)
      end do

      kj  = k1
      fac = one

      order_loop: do j=1,k1
        ki = kj
        j1 = j+1
        if(j>1) then
           i = k1
           do jj=j,k1
              li = i+lk
              lj = li+kj
              h(i) = (h(i)-h(i-1))/(t(lj)-t(li))
              i = i-1
           end do
        endif

        d(j:k1) = h(j:k1)

        if(j<k1) then
            do jj=j1,k1
              ki = ki-1
              i = k1
              do j2=jj,k1
                li = i+lk
                lj = li+ki
                d(i) = ((x-t(li))*d(i)+(t(lj)-x)*d(i-1))/(t(lj)-t(li))
                i = i-1
              end do
            end do
        endif
        d(j) = d(k1)*fac
        ak = k1-j
        fac = fac*ak
        kj = kj-1
      end do order_loop

      return
      end subroutine fpader


      pure subroutine fpadno(maxtr,up,left,right,info,count,merk,jbind,n1,ier)

      !  subroutine fpadno adds a branch of length n1 to the triply linked tree,the information of
      !  which is kept in the arrays up,left,right and info. the information field of the nodes of
      !  this new branch is given in the array jbind. in linking the new branch fpadno takes account
      !  of the property of the tree that info(k) < info(right(k)) ; info(k) < info(left(k))
      !  if necessary the subroutine calls subroutine fpfrno to collect the free nodes of the tree.
      !  if no computer words are available at that moment, the error parameter ier is set to 1.
      !  ..
      !  ..scalar arguments..
      integer, intent(in)    :: maxtr,n1
      integer, intent(inout) :: count,merk
      integer, intent(out)   :: ier
      !  ..array arguments..
      integer, intent(inout) :: up(maxtr),left(maxtr),right(maxtr),info(maxtr)
      integer, intent(in)    :: jbind(n1)
      !  ..local scalars..
      integer :: k,level,point
      logical :: is_left

      !  ..
      point   = 1
      level  = 1
      k       = left(point)
      is_left = .true.
      loop: do while (k/=0 .and. info(k)-jbind(level)<=0)
          point = k
          if (info(k)-jbind(level)<0) then
              k       = right(point)
              is_left = .false.
          else ! info(k)-jbind(level)==0
              level  = level+1
              k       = left(point)
              is_left = .true.
          endif
      end do loop

      loop2: do while (level<=n1)
          count = count+1
          if (count>maxtr) then
             call fpfrno(maxtr,up,left,right,info,point,merk,n1,count,ier)
             if (ier/=FITPACK_OK) return
          endif

          info (count) = jbind(level)
          left (count) = 0
          right(count) = k

          if(is_left) then
              up   (count) = point
              left (point) = count
          else
              is_left      = .true.
              right(point) = count
              up   (count) = up(point)
          endif

          point        = count
          level       = level+1
          k = 0
      end do loop2

      ! Success!
      ier = FITPACK_OK

      end subroutine fpadno


      !  given a idim-dimensional spline curve of degree k, in its b-spline representation ( knots t(j),
      !  j=1,...,n , b-spline coefficients c(j), j=1,...,nc) and given also a polynomial curve in its
      !  b-spline representation ( coefficients cp(j), j=1,...,np), subroutine fpadpo calculates the b-spline
      !  representation (coefficients c(j),j=1,...,nc) of the sum of the two curves.
      pure subroutine fpadpo(idim,t,n,c,nc,k,cp,np,cc,t1,t2)

      !  ..
      !  ..scalar arguments..
      integer, intent(in) :: idim,k,n,nc,np
      !  ..array arguments..
      real(RKIND), intent(in) :: t(n),cp(np)
      real(RKIND), intent(out) :: c(nc),cc(nc),t1(n),t2(n)
      !  ..local scalars..
      integer :: i,ii,j,jj,k1,l,l1,n1,n2,nk1,nk2
      !  ..
      k1  = k+1
      nk1 = n-k1
      !  initialization
      j = 1
      l = 1
      do jj=1,idim
        l1 = j
        do ii=1,k1
          cc(l1) = cp(l)
          l1 = l1+1
          l = l+1
        end do
        j = j+n
        l = l+k1
      end do

      if (nk1/=k1) then
          n1 = k1*2
          j  = n
          l  = n1
          do i=1,k1
            t1(i) = t(i)
            t1(l) = t(j)
            l = l-1
            j = j-1
          end do
          !  find the b-spline representation of the given polynomial curve
          !  according to the given set of knots.
          nk2 = nk1-1
          do l=k1,nk2
            l1 = l+1
            j = 1
            do i=1,idim
              call fpinst(0,t1,n1,cc(j),k,t(l1),l,t2,n2,cc(j),n)
              j = j+n
            end do
            t1(1:n2) = t2(1:n2)
            n1 = n2
          end do
      endif

      !  find the b-spline representation of the resulting curve.
      j = 1
      do jj=1,idim
         l = j
         do i=1,nk1
            c(l) = cc(l)+c(l)
            l = l+1
         end do
         j = j+n
      end do
      return
      end subroutine fpadpo

      !  function fpback calculates the solution of the system of equations a*c = z with
      !  a a n x n upper triangular matrix of bandwidth k.
      pure function fpback(a,z,n,k,nest) result(c)

      !  ..scalar arguments..
      integer, intent(in) :: n,k,nest
      !  ..array arguments..
      real(RKIND), intent(in)  :: a(nest,k),z(n)
      real(RKIND)              :: c(n)
      !  ..local scalars..
      real(RKIND) :: store
      integer :: i,i1,j,k1,l,m
      !  ..
      k1   = k-1
      c(n) = z(n)/a(n,1)
      i    = n-1
      if (i==0) return

      rows: do j=2,n
        store = z(i)
        i1 = merge(j-1,k1,j<=k1)
        m = i
        do l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
        end do
        c(i) = store/a(i,1)
        i = i-1
      end do rows

      end function fpback

      !  function fpbacp calculates the solution of the system of equations g * c = z
      !  with g  a n x n upper triangular matrix of the form
      !            ! a '   !
      !        g = !   ' b !
      !            ! 0 '   !
      !  with b a n x k matrix and a a (n-k) x (n-k) upper triangular matrix of bandwidth k1.
      pure function fpbacp(a,b,z,n,k,k1,nest) result(c)

      !  ..scalar arguments..
      integer, intent(in) :: n,k,k1,nest
      !  ..array arguments..
      real(RKIND), intent(in) :: a(nest,k1),b(nest,k),z(n)
      real(RKIND) :: c(n)
      !  ..local scalars..
      integer :: i,i1,j,l,l0,l1,n2
      real(RKIND) :: store
      !  ..
      n2 = n-k
      l  = n
      do i=1,k
         store = z(l)
         j = k+2-i
         if (i/=1) then
             l0 = l
             do l1=j,k
               l0 = l0+1
               store = store-c(l0)*b(l,l1)
             end do
         endif
         c(l) = store/b(l,j-1)
         l = l-1
         if (l==0) return
      end do
      do i=1,n2
         store = z(i)
         l = n2
         do j=1,k
           l = l+1
           store = store-c(l)*b(i,j)
         end do
         c(i) = store
      end do
      i = n2
      c(i) = c(i)/a(i,1)
      if (i==1) return
      do j=2,n2
         i = i-1
         store = c(i)
         i1 = k
         if(j<=k) i1=j-1
         l = i
         do l0=1,i1
           l = l+1
           store = store-c(l)*a(i,l0+1)
         end do
         c(i) = store/a(i,1)
      end do
      return
      end function fpbacp


      !  subroutine fpbfou calculates the integrals
      !                    /t(n-3)
      !    ress(j) =      !        nj,4(x)*sin(par*x) dx    and
      !              t(4)/
      !                    /t(n-3)
      !    resc(j) =      !        nj,4(x)*cos(par*x) dx ,  j=1,2,...n-4
      !              t(4)/
      !  where nj,4(x) denotes the cubic b-spline defined on the knots t(j),t(j+1),...,t(j+4).
      pure subroutine fpbfou(t,n,par,ress,resc)

      !  calling sequence:
      !     call fpbfou(t,n,par,ress,resc)
      !
      !  input parameters:
      !    t    : real array,length n, containing the knots.
      !    n    : integer, containing the number of knots.
      !    par  : real, containing the value of the parameter par.
      !
      !  output parameters:
      !    ress  : real array,length n, containing the integrals ress(j).
      !    resc  : real array,length n, containing the integrals resc(j).
      !
      !  restrictions:
      !    n >= 10, t(4) < t(5) < ... < t(n-4) < t(n-3).
      !  ..
      !  ..scalar arguments..
      integer, intent(in) :: n
      real(RKIND), intent(in) :: par
      !  ..array arguments..
      real(RKIND), intent(in) :: t(n)
      real(RKIND), intent(out) :: ress(n),resc(n)
      !  ..local scalars..
      integer :: i,ic,ipj,is,j,jj,jp1,jp4,k,li,lj,ll,nmj,nm3,nm7
      real(RKIND) :: ak,beta,c1,c2,delta,fac,f1,f2,f3,sign,s1,s2,term
      !  ..local arrays..
      real(RKIND) :: co(5),si(5),hs(5),hc(5),rs(3),rc(3)
      !  ..
      !  initialization.
      real(RKIND), parameter ::  eps = 0.1e-07_RKIND
      real(RKIND), parameter :: con1 = 0.5e-01_RKIND
      real(RKIND), parameter :: con2 = 0.12e+03_RKIND
      nm3 = n-3
      nm7 = n-7

      term = merge(six/par,zero,par/=zero)
      beta = par*t(4)
      co(1) = cos(beta)
      si(1) = sin(beta)

      !  calculate the integrals ress(j) and resc(j), j=1,2,3 by setting up a divided difference table.
      left: do j=1,3
        jp1 = j+1
        jp4 = j+4
        beta = par*t(jp4)
        co(jp1) = cos(beta)
        si(jp1) = sin(beta)
        call fpcsin(t(4),t(jp4),par,si(1),co(1),si(jp1),co(jp1),rs(j),rc(j))
        i = 5-j
        hs(i) = zero
        hc(i) = zero
        do jj=1,j
          ipj = i+jj
          hs(ipj) = rs(jj)
          hc(ipj) = rc(jj)
        end do
        do jj=1,3
          if (i<jj) i = jj
          k = 5
          li = jp4
          do ll=i,4
            lj = li-jj
            fac = t(li)-t(lj)
            hs(k) = (hs(k)-hs(k-1))/fac
            hc(k) = (hc(k)-hc(k-1))/fac
            k = k-1
            li = li-1
          end do
        end do
        ress(j) = hs(5)-hs(4)
        resc(j) = hc(5)-hc(4)
      end do left

      !  calculate the integrals ress(j) and resc(j),j=4,5,...,n-7.
      center: do j=4,nm7
        jp4 = j+4
        beta = par*t(jp4)
        co(5) = cos(beta)
        si(5) = sin(beta)
        delta = t(jp4)-t(j)

        ! the way of computing ress(j) and resc(j) depends on the value of beta = par*(t(j+4)-t(j)).
        beta = delta*par

        ! if |beta|>1 the integrals are calculated by setting up a divided difference table.
        if(abs(beta)>one) then

            hs(1:5) = si(1:5)
            hc(1:5) = co(1:5)
            do jj=1,3
              k = 5
              li = jp4
              do ll=jj,4
                lj = li-jj
                fac = par*(t(li)-t(lj))
                hs(k) = (hs(k)-hs(k-1))/fac
                hc(k) = (hc(k)-hc(k-1))/fac
                k = k-1
                li = li-1
              end do
            end do
            s2 = (hs(5)-hs(4))*term
            c2 = (hc(5)-hc(4))*term
        else
            ! if |beta|<=1 the integrals are calculated by evaluating a series expansion.
            hs(:4) = par*(t(j+1:j+4)-t(j))
            hc(:4) = hs(:4)
            f3 = con1*sum(hs(1:4))
            c1 = fourth
            s1 = f3
            if(abs(f3)>eps) then
                sign = one
                fac  = con2
                k    = 5
                is   = 0

                series_coef: do ic=1,20
                   k = k+1
                   ak = k
                   fac = fac*ak
                   f1 = zero
                   f3 = zero
                   do i=1,4
                      f1 = f1+hc(i)
                      f2 = f1*hs(i)
                      hc(i) = f2
                      f3 = f3+f2
                   end do
                   f3 = f3*six/fac
                   if (is==0) then
                       sign = -sign
                       is = 1
                       c1 = c1+f3*sign
                   else
                       is = 0
                       s1 = s1+f3*sign
                   end if
                   if(abs(f3)<=eps) exit series_coef
                end do series_coef
            endif
            s2 = delta*(co(1)*s1+si(1)*c1)
            c2 = delta*(co(1)*c1-si(1)*s1)
        endif

        ress(j) = s2
        resc(j) = c2
        co(1:4) = co(2:5)
        si(1:4) = si(2:5)
      end do center
      !  calculate the integrals ress(j) and resc(j),j=n-6,n-5,n-4 by setting
      !  up a divided difference table.
      right: do j=1,3
        nmj = nm3-j
        i = 5-j
        call fpcsin(t(nm3),t(nmj),par,si(4),co(4),si(i-1),co(i-1),rs(j),rc(j))
        hc(i:i+j) = [zero,rc(1:j)]
        hs(i:i+j) = [zero,rs(1:j)]
        do jj=1,3
          if(i<jj) i = jj
          k = 5
          li = nmj
          do ll=i,4
            lj = li+jj
            fac = t(lj)-t(li)
            hs(k) = (hs(k-1)-hs(k))/fac
            hc(k) = (hc(k-1)-hc(k))/fac
            k = k-1
            li = li+1
          end do
        end do
        ress(nmj) = hs(4)-hs(5)
        resc(nmj) = hc(4)-hc(5)
      end do right

      return
      end subroutine fpbfou


      pure subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)

      !  ..scalar arguments..
      integer    , intent(in)  :: nx,ny,kx,ky,mx,my
      !  ..array arguments..
      real(RKIND), intent(in)  :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my)
      integer    , intent(out) :: lx(mx),ly(my)
      real(RKIND), intent(out) :: wx(mx,kx+1),wy(my,ky+1),z(mx*my)

      !  ..local variables..
      integer :: kx1,ky1,l,l1,m,nkx1,nky1,i,i1,j
      real(RKIND) :: arg,sp,tb,te,h(SIZ_K+1)

      ! X
      kx1  = kx+1
      nkx1 = nx-kx1
      tb   = tx(kx1)
      te   = tx(nkx1+1)
      l = kx1
      l1 = l+1
      x_array: do i=1,mx
        arg = x(i)
        if(arg<tb) arg = tb
        if(arg>te) arg = te
        do while (.not.(arg<tx(l1) .or. l==nkx1))
          l = l1
          l1 = l+1
        end do
        call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        wx(i,1:kx1) = h(1:kx1)
      end do x_array

      ! Y
      ky1  = ky+1
      nky1 = ny-ky1
      tb   = ty(ky1)
      te   = ty(nky1+1)
      l    = ky1
      l1   = l+1
      y_array: do i=1,my
        arg = y(i)
        if(arg<tb) arg = tb
        if(arg>te) arg = te
        do while (.not.(arg<ty(l1) .or. l==nky1))
          l  = l1
          l1 = l+1
        end do
        call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        wy(i,1:ky1) = h(1:ky1)
      end do y_array

      m = 0
      do i=1,mx
        l = lx(i)*nky1
        h(1:kx1) = wx(i,1:kx1)
        do j=1,my
          l1 = l+ly(j)
          sp = zero
          do i1=1,kx1
            sp = sp+h(i1)*dot_product(c(l1+1:l1+ky1),wy(j,1:ky1))
            l1 = l1+nky1
          end do
          m = m+1
          z(m) = sp
         end do
      end do
      return
      end subroutine fpbisp


      !  subroutine fpbspl evaluates the (k+1) non-zero b-splines of degree k at t(l) <= x < t(l+1) using
      !  the stable recurrence relation of de boor and cox.
      !  Travis Oliphant 2007 changed so that weighting of 0 is used when knots with multiplicity are present.
      !   Also, notice that l+k <= n and 1 <= l+1-k or else the routine will be accessing memory outside t
      !   Thus it is imperative that that k <= l <= n-k but this is not checked.
      pure subroutine fpbspl(t,n,k,x,l,h)
         integer    , intent(in)  :: n,k,l
         real(RKIND), intent(in)  :: x,t(n)
         real(RKIND), intent(out) :: h(SIZ_K+1)

         ! Local variables
         real(RKIND) :: f,hh(SIZ_K+1)
         integer     :: i,j,li,lj

         h(1) = one
         do j=1,k
           hh(1:j) = h(1:j)
           h(1) = zero
           do i=1,j
             li = l+i
             lj = li-j
             if (t(li)/=t(lj)) then
                f = hh(i)/(t(li)-t(lj))
                h(i) = h(i)+f*(t(li)-x)
                h(i+1) = f*(x-t(lj))
             else
                h(i+1) = zero
             endif
           end do
         end do

      end subroutine fpbspl

      !  subroutine fpchec verifies the number and the position of the knots t(j),j=1,2,...,n of a spline
      !  of degree k, in relation to the number and the position of the data points x(i),i=1,2,...,m.
      !  If all of the following conditions are fulfilled, the error parameter ier is set to zero. if one
      !  of the conditions is violated, an error flag is returned.
      pure subroutine fpchec(x,m,t,n,k,ier)
         integer,     intent(in)  :: m,n,k
         real(RKIND), intent(in) :: x(m),t(n)
         integer,     intent(out) :: ier

         ! Local variables
         integer :: i,j,k1,k2,l,nk1,nk2,nk3
         real(RKIND) :: tj,tl

         ! Init sizes
         k1 = k+1
         k2 = k1+1
         nk1 = n-k1
         nk2 = nk1+1

         ier = FITPACK_INPUT_ERROR

         ! 1) k+1 <= n-k-1 <= m
         if(nk1<k1 .or. nk1>m) return

         ! 2) monotonicity
         !    t(1) <= t(2) <= ... <= t(k+1)
         !    t(n-k) <= t(n-k+1) <= ... <= t(n)

         j = n
         monotonic: do i=1,k
            if(t(i)>t(i+1)) return
            if(t(j)<t(j-1)) return
            j = j-1
         end do monotonic

         ! 3) t(k+1) < t(k+2) < ... < t(n-k)
         do i=k2,nk2
            if(t(i)<=t(i-1)) return
         end do

         ! 4) t(k+1) <= x(i) <= t(n-k)
         ! 5) schoenberg and whitney conditions: they must hold for at least one subset of data points, i.e.
         !    there must be a subset of data points y(j) such that
         !         t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1

         if(x(1)<t(k1) .or. x(m)>t(nk2)) return
         if(x(1)>=t(k2) .or. x(m)<=t(nk1)) return

         i   = 1
         l   = k2
         nk3 = nk1-1

         if (nk3>=2) then
             do j=2,nk3
                tj = t(j)
                l  = l+1
                tl = t(l)
                do while (i<m .and. x(i)<=tj)
                    i = i+1
                    if(i>=m) return
                end do
                if (x(i)>=tl) return
             end do
         endif

         ! All checks passed
         ier = FITPACK_OK

      end subroutine fpchec


      pure subroutine fpched(x,m,t,n,k,ib,ie,ier)

      !  subroutine fpched verifies the number and the position of the knots t(j),j=1,2,...,n of a spline
      !  of degree k,with ib derative constraints at x(1) and ie constraints at x(m), in relation to the
      !  number and the position of the data points x(i),i=1,2,...,m. if all of the following conditions
      !  are fulfilled, the error parameter ier is set to zero. if one of the conditions is violated ier
      !  is set to ten.
      !  ..
      !  ..scalar arguments..
      integer, intent(in) :: m,n,k,ib,ie
      integer, intent(out) :: ier
      !  ..array arguments..
      real(RKIND), intent(in) :: x(m),t(n)
      !  ..local scalars..
      integer :: i,ib1,ie1,j,jj,k1,k2,l,nk1,nk2,nk3
      real(RKIND) :: tj,tl
      !  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ib1 = max(ib-1,0)
      ie1 = max(ie-1,0)
      ier = FITPACK_INPUT_ERROR

      ! 1) k+1 <= n-k-1 <= m + max(0,ib-1) + max(0,ie-1)
      if(nk1<k1 .or. nk1>(m+ib1+ie1)) return

      ! 2) t(1)   <= t(2)     <= ... <= t(k+1)
      !    t(n-k) <= t(n-k+1) <= ... <= t(n)
      j = n
      do i=1,k
        if(t(i)>t(i+1)) return
        if(t(j)<t(j-1)) return
        j = j-1
      end do

      ! 3) t(k+1) < t(k+2) < ... < t(n-k)
      do i=k2,nk2
        if(t(i)<=t(i-1)) return
      end do

      ! 4) t(k+1) <= x(i) <= t(n-k)
      if(x(1)<t(k1) .or. x(m)>t(nk2)) return

      ! 5) the conditions specified by schoenberg and whitney must hold for at least one subset
      !    of data points, i.e. there must be a subset of data points y(j) such that
      !    t(j) < y(j) < t(j+k+1), j=1+ib1,2+ib1,...,n-k-1-ie1
      !       with ib1 = max(0,ib-1), ie1 = max(0,ie-1)
      if (x(1)>=t(k2) .or. x(m)<=t(nk1)) return
      i  = 1
      jj = 2+ib1
      l = jj+k
      nk3 = nk1-1-ie1
      if(nk3>=jj) then
          do j=jj,nk3
            tj = t(j)
            l = l+1
            tl = t(l)
            do
              i = i+1
              if (i>=m) return
              if (x(i)>tj) exit
            end do
            if(x(i)>=tl) return
          end do
      endif

      ! Success! all checks passed
      ier = FITPACK_OK
      return
      end subroutine fpched

      ! subroutine fpchep verifies the number and the position of the knots t(j),j=1,2,...,n of a
      ! periodic spline of degree k, in relation to the number and the position of the data points
      ! x(i),i=1,2,...,m.
      ! if all of the following conditions are fulfilled, ier is set to zero.
      ! if one of the conditions is violated ier is set to ten.
      pure integer function fpchep(x,m,t,n,k) result(ier)
          !  ..scalar arguments..
          integer, intent(in) :: m,n,k
          !  ..array arguments..
          real(RKIND), intent(in) :: x(m),t(n)
          !  ..local scalars..
          integer :: i,i1,i2,j,j1,k1,k2,l,l1,l2,mm,m1,nk1,nk2
          real(RKIND) :: per,tj,tl,xi
          !  ..
          k1  = k+1
          k2  = k1+1
          nk1 = n-k1
          nk2 = nk1+1
          m1  = m-1
          ier = FITPACK_INPUT_ERROR

          ! 1) k+1 <= n-k-1 <= m+k-1
          if(nk1<k1 .or. n>m+2*k) return

          ! 2) t(1) <= t(2) <= ... <= t(k+1)
          !    t(n-k) <= t(n-k+1) <= ... <= t(n)
          j = n
          do i=1,k
             if (t(i)>t(i+1)) return
             if (t(j)<t(j-1)) return
             j = j-1
          end do

          ! 3) t(k+1) < t(k+2) < ... < t(n-k)
          do i=k2,nk2
             if (t(i)<=t(i-1)) return
          end do

          ! 4) t(k+1) <= x(i) <= t(n-k)
          if(x(1)<t(k1) .or. x(m)>t(nk2)) return

          ! 5) the conditions specified by schoenberg and whitney must hold for at least one subset
          !    of data points, i.e. there must be a subset of data points y(j) such that
          !        t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1
          l1 = k1
          l2 = 1
          outer: do l=1,m
             xi = x(l)
             do while (.not.(xi<t(l1+1) .or. l==nk1))
                l1 = l1+1
                l2 = l2+1
                if (l2>k1) exit outer
             end do
          end do outer
          if (l2<=k1) l = m

          per = t(nk2)-t(k1)
          subset_start: do i1=2,l
             i = i1-1
             mm = i+m1
             subset_inner: do j=k1,nk1
                tj = t(j)
                j1 = j+k1
                tl = t(j1)
                xi = -huge(xi)
                do while (xi<=tj)
                    i = i+1
                    if (i>mm) cycle subset_start
                    i2 = i-m1

                    if (i2<=0) then
                      xi = x(i)
                    else
                      xi = x(i2)+per
                    endif
                end do
                if(xi>=tl) cycle subset_start
             end do subset_inner

             ! A full Shoenberg-Whitney subset is found
             ier = FITPACK_OK
             return
          end do subset_start

          ! No subsets found: return with error

      end function fpchep


      recursive subroutine fpclos(iopt,idim,m,u,mx,x,w,k,s,nest,tol, &
        maxit,k1,k2,n,t,nc,c,fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)

      !  ..
      !  ..scalar arguments..
      real(RKIND) s,tol,fp
      integer iopt,idim,m,mx,k,nest,maxit,k1,k2,n,nc,ier
      !  ..array arguments..
      real(RKIND) u(m),x(mx),w(m),t(nest),c(nc),fpint(nest),z(nc),a1(nest,k1), &
       a2(nest,k),b(nest,k2),g1(nest,k2),g2(nest,k1),q(m,k1)
      integer nrdata(nest)
      !  ..local scalars..
      real(RKIND) acc,cos,d1,fac,fpart,fpms,fpold,fp0,f1,f2,f3,p,per,pinv,piv, &
       p1,p2,p3,sin,store,term,ui,wi,rn
      integer i,ich1,ich3,ij,ik,it,iter,i1,i2,i3,j,jj,jk,jper,j1,j2,kk, &
       kk1,k3,l,l0,l1,l5,mm,m1,new,nk1,nk2,nmax,nmin,nplus,npl1, &
       nrint,n10,n11,n7,n8
      !  ..local arrays..
      real(RKIND) h(SIZ_K+1),h1(7),h2(6),xi(SIZ_IDIM)
      !  set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 1: determination of the number of knots and their position     c
      !  **************************************************************      c
      !  given a set of knots we compute the least-squares closed curve      c
      !  sinf(u). if the sum f(p=inf) <= s we accept the choice of knots.    c
      !  if iopt=-1 sinf(u) is the requested curve                           c
      !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
      !    if fp <=s we will continue with the current set of knots.         c
      !    if fp > s we will increase the number of knots and compute the    c
      !       corresponding least-squares curve until finally fp<=s.         c
      !  the initial choice of knots depends on the value of s and iopt.     c
      !    if s=0 we have spline interpolation; in that case the number of   c
      !    knots equals nmax = m+2*k.                                        c
      !    if s > 0 and                                                      c
      !      iopt=0 we first compute the least-squares polynomial curve of   c
      !      degree k; n = nmin = 2*k+2. since s(u) must be periodic we      c
      !      find that s(u) reduces to a fixed point.                        c
      !      iopt=1 we start with the set of knots found at the last         c
      !      call of the routine, except for the case that s > fp0; then     c
      !      we compute directly the least-squares polynomial curve.         c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m1 = m-1
      kk = k
      kk1 = k1
      k3 = 3*k+1
      nmin = 2*k1
      !  determine the length of the period of the splines.
      per = u(m)-u(1)
      if(iopt<0) go to 50
      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      !  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k
      if(s>zero .or. nmax==nmin) go to 30
      !  if s=0, s(u) is an interpolating curve.
      n = nmax
      !  test whether the required storage space exceeds the available one.
      if(n>nest) go to 620
      !  find the position of the interior knots in case of interpolation.
   5  if((k/2)*2 ==k) go to 20
      do 10 i=2,m1
        j = i+k
        t(j) = u(i)
  10  continue
      if(s>0.) go to 50
      kk = k-1
      kk1 = k
      if(kk>0) go to 50
      t(1) = t(m)-per
      t(2) = u(1)
      t(m+1) = u(m)
      t(m+2) = t(3)+per
      jj = 0
      do 15 i=1,m1
        j = i
        do 12 j1=1,idim
          jj = jj+1
          c(j) = x(jj)
          j = j+n
  12    continue
  15  continue
      jj = 1
      j = m
      do 17 j1=1,idim
        c(j) = c(jj)
        j = j+n
        jj = jj+n
  17  continue
      fp = zero
      fpint(n) = fp0
      fpint(n-1) = zero
      nrdata(n) = 0
      go to 630
  20  do 25 i=2,m1
         j = i+k
         t(j) = (u(i)+u(i-1))*half
  25  continue
      go to 50
      !  if s > 0 our initial choice depends on the value of iopt.
      !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
      !  polynomial curve. (i.e. a constant point).
      !  if iopt=1 and fp0>s we start computing the least-squares closed
      !  curve according the set of knots found at the last call of the
      !  routine.
  30  if(iopt==0) go to 35
      if(n==nmin) go to 35
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0>s) go to 50
      !  the case that s(u) is a fixed point is treated separetely.
      !  fp0 denotes the corresponding sum of squared residuals.
  35  fp0 = zero
      d1 = zero
      z(1:idim) = zero
      jj = 0
      do 45 it=1,m1
        wi = w(it)
        call fpgivs(wi,d1,cos,sin)
        do 40 j=1,idim
          jj = jj+1
          fac = wi*x(jj)
          call fprota(cos,sin,fac,z(j))
          fp0 = fp0+fac**2
  40    continue
  45  continue
      z(1:idim) = z(1:idim)/d1
      !  test whether that fixed point is a solution of our problem.
      fpms = fp0-s
      if(fpms<acc .or. nmax==nmin) go to 640
      fpold = fp0
      !  test whether the required storage space exceeds the available one.
      if(n>=nest) go to 620
      !  start computing the least-squares closed curve with one
      !  interior knot.
      nplus = 1
      n = nmin+1
      mm = (m+1)/2
      t(k2) = u(mm)
      nrdata(1) = mm-2
      nrdata(2) = m1-mm
      !  main loop for the different sets of knots. m is a save upper
      !  bound for the number of trials.
  50  do 340 iter=1,m
      !  find nrint, the number of knot intervals.
        nrint = n-nmin+1
      !  find the position of the additional knots which are needed for
      !  the b-spline representation of s(u). if we take
      !      t(k+1) = u(1), t(n-k) = u(m)
      !      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
      !      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
      !  then s(u) will be a smooth closed curve if the b-spline
      !  coefficients satisfy the following conditions
      !      c((i-1)*n+n7+j) = c((i-1)*n+j), j=1,...k,i=1,2,...,idim (**)
      !  with n7=n-2*k-1.
        t(k1) = u(1)
        nk1 = n-k1
        nk2 = nk1+1
        t(nk2) = u(m)
        do 60 j=1,k
          i1 = nk2+j
          i2 = nk2-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
  60    continue
      !  compute the b-spline coefficients of the least-squares closed curve
      !  sinf(u). the observation matrix a is built up row by row while
      !  taking into account condition (**) and is reduced to triangular
      !  form by givens transformations .
      !  at the same time fp=f(p=inf) is computed.
      !  the n7 x n7 triangularised upper matrix a has the form
      !            ! a1 '    !
      !        a = !    ' a2 !
      !            ! 0  '    !
      !  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
      !  matrix of bandwidth k+1 ( n10 = n7-k).
      !  initialization.
        z(1:nc) = zero
        a1(1:nk1,1:kk1) = zero
        n7 = nk1-k
        n10 = n7-kk
        jper = 0
        fp = zero
        l = k1
        jj = 0
        do 290 it=1,m1
      !  fetch the current data point u(it),x(it)
          ui = u(it)
          wi = w(it)
          do 75 j=1,idim
            jj = jj+1
            xi(j) = x(jj)*wi
  75      continue
      !  search for knot interval t(l) <= ui < t(l+1).
  80      if(ui<t(l+1)) go to 85
          l = l+1
          go to 80
      !  evaluate the (k+1) non-zero b-splines at ui and store them in q.
  85      call fpbspl(t,n,k,ui,l,h)
          do 90 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  90      continue
          l5 = l-k1
      !  test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui
          if(l5<n10) go to 285
          if(jper/=0) go to 160
      !  initialize the matrix a2.
          a2(1:n7,1:kk) = zero
          jk = n10+1
          do 110 i=1,kk
            ik = jk
            do 100 j=1,kk1
              if(ik<=0) go to 105
              a2(ik,i) = a1(ik,j)
              ik = ik-1
 100        continue
 105        jk = jk+1
 110      continue
          jper = 1
      !  if one of the b-splines nj,k+1(u),j=n7+1,...nk1 is not zero at ui
      !  we take account of condition (**) for setting up the new row
      !  of the observation matrix a. this row is stored in the arrays h1
      !  (the part with respect to a1) and h2 (the part with
      !  respect to a2).
 160      do 170 i=1,kk
            h1(i) = zero
            h2(i) = zero
 170      continue
          h1(kk1) = zero
          j = l5-n10
          do 210 i=1,kk1
            j = j+1
            l0 = j
 180        l1 = l0-kk
            if(l1<=0) go to 200
            if(l1<=n10) go to 190
            l0 = l1-n10
            go to 180
 190        h1(l1) = h(i)
            go to 210
 200        h2(l0) = h2(l0)+h(i)
 210      continue
      !  rotate the new row of the observation matrix into triangle
      !  by givens transformations.
          if(n10<=0) go to 250
      !  rotation with the rows 1,2,...n10 of matrix a.
          do 240 j=1,n10
            piv = h1(1)
            if(piv/=zero) go to 214
            do 212 i=1,kk
              h1(i) = h1(i+1)
 212        continue
            h1(kk1) = zero
            go to 240
      !  calculate the parameters of the givens transformation.
 214        call fpgivs(piv,a1(j,1),cos,sin)
      !  transformation to the right hand side.
            j1 = j
            do 217 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 217        continue
      !  transformations to the left hand side with respect to a2.
            call fprota(cos,sin,h2(1:kk),a2(j,1:kk))
            if(j==n10) go to 250
            i2 = min0(n10-j,kk)
      !  transformations to the left hand side with respect to a1.
            do 230 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),a1(j,i1))
              h1(i) = h1(i1)
 230        continue
            h1(i1) = zero
 240      continue
      !  rotation with the rows n10+1,...n7 of matrix a.
 250      do 270 j=1,kk
            ij = n10+j
            if(ij<=0) go to 270
            piv = h2(j)
            if (piv==zero) go to 270
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a2(ij,j),cos,sin)
      !  transformations to right hand side.
            j1 = ij
            do 255 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 255        continue
            if(j==kk) go to 280
            j1 = j+1
      !  transformations to left hand side.
            call fprota(cos,sin,h2(j1:kk),a2(ij,j1:kk))
 270      continue
      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
 280      fp = fp+sum(xi(1:idim)**2)
          go to 290
      !  rotation of the new row of the observation matrix into
      !  triangle in case the b-splines nj,k+1(u),j=n7+1,...n-k-1 are all zero
      !  at ui.
 285      j = l5
          do 140 i=1,kk1
            j = j+1
            piv = h(i)
            if (piv==zero) go to 140
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a1(j,1),cos,sin)
      !  transformations to right hand side.
            j1 = j
            do 125 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 125        continue
            if(i==kk1) go to 150
            i2 = 1
            i3 = i+1
      !  transformations to left hand side.
            do 130 i1=i3,kk1
              i2 = i2+1
              call fprota(cos,sin,h(i1),a1(j,i2))
 130        continue
 140      continue
      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
 150      fp = fp + sum(xi(1:idim)**2)
 290    continue
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
      !  backward substitution to obtain the b-spline coefficients .
        j1 = 1
        do j2=1,idim
           c(j1:j1+n-1) = fpbacp(a1,a2,z(j1),n7,kk,kk1,nest)
           j1 = j1+n
        end do
      !  calculate from condition (**) the remaining coefficients.
        do 297 i=1,k
          j1 = i
          do 295 j=1,idim
            j2 = j1+n7
            c(j2) = c(j1)
            j1 = j1+n
 295      continue
 297    continue
        if(iopt<0) go to 660
      !  test whether the approximation sinf(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 660
      !  if f(p=inf) < s accept the choice of knots.
        if(fpms<zero) go to 350
      !  if n=nmax, sinf(u) is an interpolating curve.
        if(n==nmax) go to 630
      !  increase the number of knots.
      !  if n=nest we cannot increase the number of knots because of the
      !  storage capacity limitation.
        if(n==nest) go to 620
      !  determine the number of knots nplus we are going to add.
        npl1 = nplus*2
        rn = nplus
        if(fpold-fp>acc) npl1 = int(rn*fpms/(fpold-fp))
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
        fpold = fp
      !  compute the sum of squared residuals for each knot interval
      !  t(j+k) <= ui <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = zero
        i = 1
        l = k1
        jj = 0
        do 320 it=1,m1
          if(u(it)<t(l)) go to 300
          new = 1
          l = l+1
 300      term = zero
          l0 = l-k2
          do 310 j2=1,idim
            fac = zero
            j1 = l0
            do 305 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 305        continue
            jj = jj+1
            term = term+(w(it)*(fac-x(jj)))**2
            l0 = l0+n
 310      continue
          fpart = fpart+term
          if(new==0) go to 320
          if(l>k2) go to 315
          fpint(nrint) = term
          new = 0
          go to 320
 315      store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 320    continue
        fpint(nrint) = fpint(nrint)+fpart
        do 330 l=1,nplus
      !  add a new knot
          call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1)
      !  if n=nmax we locate the knots as for interpolation
          if(n==nmax) go to 5
      !  test whether we cannot further increase the number of knots.
          if(n==nest) go to 340
 330    continue
      !  restart the computations with the new set of knots.
 340  continue
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 2: determination of the smoothing closed curve sp(u).          c
      !  **********************************************************          c
      !  we have determined the number of knots and their position.          c
      !  we now compute the b-spline coefficients of the smoothing curve     c
      !  sp(u). the observation matrix a is extended by the rows of matrix   c
      !  b expressing that the kth derivative discontinuities of sp(u) at    c
      !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
      !  ponding weights of these additional rows are set to 1/p.            c
      !  iteratively we then have to determine the value of p such that f(p),c
      !  the sum of squared residuals be = s. we already know that the least-c
      !  squares polynomial curve corresponds to p=0, and that the least-    c
      !  squares periodic spline curve corresponds to p=infinity. the        c
      !  iteration process which is proposed here, makes use of rational     c
      !  interpolation. since f(p) is a convex and strictly decreasing       c
      !  function of p, it can be approximated by a rational function        c
      !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
      !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
      !  to calculate the new value of p such that r(p)=s. convergence is    c
      !  guaranteed by taking f1>0 and f3<zero                                 c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  evaluate the discontinuity jump of the kth derivative of the
      !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
 350  call fpdisc(t,n,k2,b,nest)
      !  initial value for p.
      p1 = zero
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      n11 = n10-1
      n8 = n7-1
      p = zero
      l = n7
      do 352 i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if(l==0) go to 356
 352  continue
      p = p + sum(a1(1:n10,1))
 356  rn = n7
      p = rn/p
      ich1 = 0
      ich3 = 0
      !  iteration process to find the root of f(p) = s.
      do 595 iter=1,maxit
      !  form the matrix g  as the matrix a extended by the rows of matrix b.
      !  the rows of matrix b with weight 1/p are rotated into
      !  the triangularised observation matrix a.
      !  after triangularisation our n7 x n7 matrix g takes the form
      !            ! g1 '    !
      !        g = !    ' g2 !
      !            ! 0  '    !
      !  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
      !  matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one/p
      !  store matrix a into g
        do 358 i=1,nc
          c(i) = z(i)
 358    continue
        do i=1,n7
          g1(i,k1) = a1(i,k1)
          g1(i,k2) = zero
          g2(i,1) = zero
          do j=1,k
            g1(i,j) = a1(i,j)
            g2(i,j+1) = a2(i,j)
          end do
        end do
        l = n10
        do j=1,k1
          if(l<=0) exit
          g2(l,1) = a1(l,j)
          l = l-1
        end do
        do 540 it=1,n8
      !  fetch a new row of matrix b and store it in the arrays h1 (the part
      !  with respect to g1) and h2 (the part with respect to g2).
          xi(:idim) = zero
          h1(:k1) = zero
          h2(:k1) = zero
          h1(k2) = zero
          if(it>n11) go to 420
          l = it
          l0 = it
          do 390 j=1,k2
            if(l0==n10) go to 400
            h1(j) = b(it,j)*pinv
            l0 = l0+1
 390      continue
          go to 470
 400      l0 = 1
          do 410 l1=j,k2
            h2(l0) = b(it,l1)*pinv
            l0 = l0+1
 410      continue
          go to 470
 420      l = 1
          i = it-n10
          do 460 j=1,k2
            i = i+1
            l0 = i
 430        l1 = l0-k1
            if(l1<=0) go to 450
            if(l1<=n11) go to 440
            l0 = l1-n11
            go to 430
 440        h1(l1) = b(it,j)*pinv
            go to 460
 450        h2(l0) = h2(l0)+b(it,j)*pinv
 460      continue
          if(n11<=0) go to 510
      !  rotate this row into triangle by givens transformations
      !  rotation with the rows l,l+1,...n11.
 470      do 500 j=l,n11
            piv = h1(1)
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,g1(j,1),cos,sin)
      !  transformation to right hand side.
            j1 = j
            do 475 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 475        continue
      !  transformation to the left hand side with respect to g2.
            do 480 i=1,k1
              call fprota(cos,sin,h2(i),g2(j,i))
 480        continue
            if(j==n11) go to 510
            i2 = min0(n11-j,k1)
      !  transformation to the left hand side with respect to g1.
            do 490 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),g1(j,i1))
              h1(i) = h1(i1)
 490        continue
            h1(i1) = zero
 500      continue
      !  rotation with the rows n11+1,...n7
 510      do 530 j=1,k1
            ij = n11+j
            if(ij<=0) go to 530
            piv = h2(j)
      !  calculate the parameters of the givens transformation
            call fpgivs(piv,g2(ij,j),cos,sin)
      !  transformation to the right hand side.
            j1 = ij
            do 515 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 515        continue
            if(j==k1) go to 540
            j1 = j+1
      !  transformation to the left hand side.
            do 520 i=j1,k1
              call fprota(cos,sin,h2(i),g2(ij,i))
 520        continue
 530      continue
 540    continue
      !  backward substitution to obtain the b-spline coefficients
        j1 = 1
        do 542 j2=1,idim
          c(j1:j1+n-1) = fpbacp(g1,g2,c(j1),n7,k1,k2,nest)
          j1 = j1+n
 542    continue
      !  calculate from condition (**) the remaining b-spline coefficients.
        do 547 i=1,k
          j1 = i
          do 545 j=1,idim
            j2 = j1+n7
            c(j2) = c(j1)
            j1 = j1+n
 545      continue
 547    continue
      !  computation of f(p).
        fp = zero
        l = k1
        jj = 0
        do 570 it=1,m1
          if(u(it)<t(l)) go to 550
          l = l+1
 550      l0 = l-k2
          term = zero
          do 565 j2=1,idim
            fac = zero
            j1 = l0
            do 560 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 560        continue
            jj = jj+1
            term = term+(fac-x(jj))**2
            l0 = l0+n
 565      continue
          fp = fp+term*w(it)**2
 570    continue
      !  test whether the approximation sp(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 660
      !  test whether the maximal number of iterations is reached.
        if (iter==maxit) go to 600
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 580
        if((f2-f3) > acc) go to 575
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 +p2*con1
        go to 595
 575    if(f2<0.) ich3 = 1
 580    if(ich1/=0) go to 590
        if((f1-f2) > acc) go to 585
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 595
        if(p>=p3) p = p2*con1 +p3*con9
        go to 595
 585    if(f2>0.) ich1 = 1
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 590    if(f2>=f1 .or. f2<=f3) go to 610
      !  find the new value for p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 595  continue
      !  error codes and messages.
 600  ier = FITPACK_MAXIT
      go to 660
 610  ier = FITPACK_S_TOO_SMALL
      go to 660
 620  ier = FITPACK_INSUFFICIENT_STORAGE
      go to 660
 630  ier = FITPACK_INTERPOLATING_OK
      go to 660
 640  ier = FITPACK_LEASTSQUARES_OK
      !  the point (z(1),z(2),...,z(idim)) is a solution of our problem.
      !  a constant function is a spline of degree k with all b-spline
      !  coefficients equal to that constant.
      do 650 i=1,k1
        rn = k1-i
        t(i) = u(1)-rn*per
        j = i+k1
        rn = i-1
        t(j) = u(m)+rn*per
 650  continue

      n = nmin
      forall(j=1:idim,i=1:k1) c(n*(j-1)+i) = z(j)

      fp = fp0
      fpint(n) = fp0
      fpint(n-1) = zero
      nrdata(n) = 0
 660  return
      end subroutine fpclos


      recursive subroutine fpcoco(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,bind,e,wrk,lwrk,iwrk,kwrk,ier)

      !  ..scalar arguments..
      real(RKIND) s,sq
      integer iopt,m,nest,maxtr,maxbin,n,lwrk,kwrk,ier
      !  ..array arguments..
      integer iwrk(kwrk)
      real(RKIND) x(m),y(m),w(m),v(m),t(nest),c(nest),sx(m),e(nest),wrk(lwrk)

      logical bind(nest)
      !  ..local scalars..
      integer i,ia,ib,ic,iq,iu,iz,izz,i1,j,k,l,l1,m1,nmax,nr,n4,n6,n8,ji,jib,jjb,jl,jr,ju,mb,nm
      real(RKIND) sql,sqmax,term,tj,xi,half
      !  ..subroutine references..
      !    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno
      !  ..
      !  set constant
      half = 0.5e0
      !  determine the maximal admissible number of knots.
      nmax = m+4
      !  the initial choice of knots depends on the value of iopt.
      !    if iopt=0 the program starts with the minimal number of knots
      !    so that can be guarantied that the concavity/convexity constraints
      !    will be satisfied.
      !    if iopt = 1 the program will continue from the point on where she
      !    left at the foregoing call.
      if(iopt>0) go to 80
      !  find the minimal number of knots.
      !  a knot is located at the data point x(i), i=2,3,...m-1 if
      !    1) v(i) ^= 0    and
      !    2) v(i)*v(i-1) <= 0  or  v(i)*v(i+1) <= zero
      m1 = m-1
      n = 4
      do 20 i=2,m1
        if(v(i)==zero .or. (v(i)*v(i-1)>zero .and. &
        v(i)*v(i+1)>0.)) go to 20
        n = n+1
      !  test whether the required storage space exceeds the available one.
        if(n+4>nest) go to 200
        t(n) = x(i)
  20  continue
      !  find the position of the knots t(1),...t(4) and t(n-3),...t(n) which
      !  are needed for the b-spline representation of s(x).
      do 30 i=1,4
        t(i) = x(1)
        n = n+1
        t(n) = x(m)
  30  continue
      !  test whether the minimum number of knots exceeds the maximum number.
      if(n>nmax) go to 210
      !  main loop for the different sets of knots.
      !  find corresponding values e(j) to the knots t(j+3),j=1,2,...n-6
      !    e(j) will take the value -1,1, or 0 according to the requirement
      !    that s(x) must be locally convex or concave at t(j+3) or that the
      !    sign of s''(x) is unrestricted at that point.
  40  i= 1
      xi = x(1)
      j = 4
      tj = t(4)
      n6 = n-6
      do 70 l=1,n6
  50    if(xi==tj) go to 60
        i = i+1
        xi = x(i)
        go to 50
  60    e(l) = v(i)
        j = j+1
        tj = t(j)
  70  continue
      !  we partition the working space
      nm = n+maxbin
      mb = maxbin+1
      ia = 1
      ib = ia+4*n
      ic = ib+nm*maxbin
      iz = ic+n
      izz = iz+n
      iu = izz+n
      iq = iu+maxbin
      ji = 1
      ju = ji+maxtr
      jl = ju+maxtr
      jr = jl+maxtr
      jjb = jr+maxtr
      jib = jjb+mb
      !  given the set of knots t(j),j=1,2,...n, find the least-squares cubic
      !  spline which satisfies the imposed concavity/convexity constraints.
      call fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,wrk(ia), &
                  wrk(ib),wrk(ic),wrk(iz),wrk(izz),wrk(iu),wrk(iq),iwrk(ji), &
                  iwrk(ju),iwrk(jl),iwrk(jr),iwrk(jjb),iwrk(jib),ier)
      !  if sq <= s or in case of abnormal exit from fpcosp, control is
      !  repassed to the driver program.
      if(sq<=s .or. ier>0) go to 300
      !  calculate for each knot interval t(l-1) <= xi <= t(l) the
      !  sum((wi*(yi-s(xi)))**2).
      !  find the interval t(k-1) <= x <= t(k) for which this sum is maximal
      !  on the condition that this interval contains at least one interior
      !  data point x(nr) and that s(x) is not given there by a straight line.
  80  sqmax = zero
      sql = zero
      l = 5
      nr = 0
      i1 = 1
      n4 = n-4
      do 110 i=1,m
        term = (w(i)*(sx(i)-y(i)))**2
        if(x(i)<t(l) .or. l>n4) go to 100
        term = term*half
        sql = sql+term
        if(i-i1<=1 .or. (bind(l-4).and.bind(l-3))) go to 90
        if(sql<=sqmax) go to 90
        k = l
        sqmax = sql
        nr = i1+(i-i1)/2
  90    l = l+1
        i1 = i
        sql = zero
 100    sql = sql+term
 110  continue
      if(m-i1<=1 .or. (bind(l-4).and.bind(l-3))) go to 120
      if(sql<=sqmax) go to 120
      k = l
      nr = i1+(m-i1)/2
      !  if no such interval is found, control is repassed to the driver
      !  program (ier = -1).
 120  if(nr==0) go to 190
      !  if s(x) is given by the same straight line in two succeeding knot
      !  intervals t(l-1) <= x <= t(l) and t(l) <= x <= t(l+1),delete t(l)
      n8 = n-8
      l1 = 0
      if(n8<=0) go to 150
      do 140 i=1,n8
        if(.not. (bind(i).and.bind(i+1).and.bind(i+2))) go to 140
        l = i+4-l1
        if(k>l) k = k-1
        n = n-1
        l1 = l1+1
        do 130 j=l,n
          t(j) = t(j+1)
 130    continue
 140  continue
      !  test whether we cannot further increase the number of knots.
 150  if(n==nmax) go to 180
      if(n==nest) go to 170
      !  locate an additional knot at the point x(nr).
      j = n
      do 160 i=k,n
        t(j+1) = t(j)
        j = j-1
 160  continue
      t(k) = x(nr)
      n = n+1
      !  restart the computations with the new set of knots.
      go to 40
      !  error codes and messages.
 170  ier = -3
      go to 300
 180  ier = -2
      go to 300
 190  ier = -1
      go to 300
 200  ier = 4
      go to 300
 210  ier = 5
 300  return
      end subroutine fpcoco


      recursive subroutine fpcons(iopt,idim,m,u,mx,x,w,ib,ie,k,s,nest, &
        tol,maxit,k1,k2,n,t,nc,c,fp,fpint,z,a,b,g,q,nrdata,ier)
      !cc         c XXX: mmnin/nmin variables on line 61
      !  ..
      !  ..scalar arguments..
      real(RKIND) s,tol,fp
      integer mx,ib,ie,k,nest,maxit,k1,k2,n,nc
      integer, intent(in)    :: iopt,idim,m
      integer, intent(inout) :: ier
      !  ..array arguments..
      real(RKIND) :: u(m),x(mx),w(m),t(nest),c(nc),fpint(nest),z(nc),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
      !  ..local scalars..
      real(RKIND) :: acc,cos,fac,fpart,fpms,fpold,fp0,f1,f2,f3,p,pinv,piv,p1,p2,p3,rn,sin,store,term,ui,wi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,jb,je,jj,j1,j2,j3,kbe, &
       l,li,lj,l0,mb,me,mm,new,nk1,nmax,nmin,nn,nplus,npl1,nrint,n8,mmin
      !  ..local arrays..
      real(RKIND) h(SIZ_K+1),xi(SIZ_IDIM)

      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 1: determination of the number of knots and their position     c
      !  **************************************************************      c
      !  given a set of knots we compute the least-squares curve sinf(u),    c
      !  and the corresponding sum of squared residuals fp=f(p=inf).         c
      !  if iopt=-1 sinf(u) is the requested curve.                          c
      !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
      !    if fp <=s we will continue with the current set of knots.         c
      !    if fp > s we will increase the number of knots and compute the    c
      !       corresponding least-squares curve until finally fp<=s.         c
      !    the initial choice of knots depends on the value of s and iopt.   c
      !    if s=0 we have spline interpolation; in that case the number of   c
      !    knots equals nmax = m+k+1-max(0,ib-1)-max(0,ie-1)                 c
      !    if s > 0 and                                                      c
      !      iopt=0 we first compute the least-squares polynomial curve of   c
      !      degree k; n = nmin = 2*k+2                                      c
      !      iopt=1 we start with the set of knots found at the last         c
      !      call of the routine, except for the case that s > fp0; then     c
      !      we compute directly the polynomial curve of degree k.           c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      !  find which data points are to be considered.
      mb = merge(2,1,ib>0)
      jb = merge(ib,1,ib>0)
      me = merge(m-1,m,ie>0)
      je = merge(ie,1,ie>0)

      if(iopt<0) go to 60
      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      !  determine nmax, the number of knots for spline interpolation.
      kbe = k1-jb-je
      mmin = kbe+2
      mm = m-mmin
      nmax = nmin+mm
      if(s>zero) go to 40
      !  if s=0, s(u) is an interpolating curve.
      !  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax>nest) go to 420
      !  find the position of the interior knots in case of interpolation.
      if(mm==0) go to 60
  25  i = k2
      j = 3-jb+k/2
      do 30 l=1,mm
        t(i) = u(j)
        i = i+1
        j = j+1
  30  continue
      go to 60
      !  if s>0 our initial choice of knots depends on the value of iopt.
      !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
      !  polynomial curve which is a spline curve without interior knots.
      !  if iopt=1 and fp0>s we start computing the least squares spline curve
      !  according to the set of knots found at the last call of the routine.
  40  if(iopt==0) go to 50
      if(n==nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0>s) go to 60
  50  n = nmin
      fpold = zero
      nplus = 0
      nrdata(1) = m-2
      !  main loop for the different sets of knots. m is a save upper bound
      !  for the number of trials.
  60  do 200 iter = 1,m
        if(n==nmin) ier = -2
      !  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
      !  find the position of the additional knots which are needed for
      !  the b-spline representation of s(u).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = u(1)
          t(i) = u(m)
          i = i-1
  70    continue
      !  compute the b-spline coefficients of the least-squares spline curve
      !  sinf(u). the observation matrix a is built up row by row and
      !  reduced to upper triangular form by givens transformations.
      !  at the same time fp=f(p=inf) is computed.
        fp = zero
      !  nn denotes the dimension of the splines
        nn = nk1-ib-ie
      !  initialize the b-spline coefficients and the observation matrix a.
        do 75 i=1,nc
          z(i) = zero
          c(i) = zero
  75    continue
        if(me<mb) go to 134
        if(nn==0) go to 82
        a(1:nn,1:k1) = zero
  82    l = k1
        jj = (mb-1)*idim
        do 130 it=mb,me
      !  fetch the current data point u(it),x(it).
          ui = u(it)
          wi = w(it)
          do 84 j=1,idim
             jj = jj+1
             xi(j) = x(jj)*wi
  84      continue
      !  search for knot interval t(l) <= ui < t(l+1).
  86      if(ui<t(l+1) .or. l==nk1) go to 90
          l = l+1
          go to 86
      !  evaluate the (k+1) non-zero b-splines at ui and store them in q.
  90      call fpbspl(t,n,k,ui,l,h)
          do 92 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  92      continue
      !  take into account that certain b-spline coefficients must be zero.
          lj = k1
          j = nk1-l-ie
          if(j>=0) go to 94
          lj = lj+j
  94      li = 1
          j = l-k1-ib
          if(j>=0) go to 96
          li = li-j
          j = 0
  96      if(li>lj) go to 120
      !  rotate the new row of the observation matrix into triangle.
          do 110 i=li,lj
            j = j+1
            piv = h(i)
            if (piv==zero) go to 110
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
      !  transformations to right hand side.
            j1 = j
            do 98 j2 =1,idim
               call fprota(cos,sin,xi(j2),z(j1))
               j1 = j1+n
  98        continue
            if(i==lj) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,lj
              i2 = i2+1
      !  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
 120      do 125 j2=1,idim
             fp  = fp+xi(j2)**2
 125      continue
 130    continue
        if(ier==(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
      !  backward substitution to obtain the b-spline coefficients.
        if(nn==0) go to 134
        j1 = 1
        do 132 j2=1,idim
           j3 = j1+ib
           c(j3:j3+nn-1) = fpback(a,z(j1),nn,k1,nest)
           j1 = j1+n
 132    continue
      !  test whether the approximation sinf(u) is an acceptable solution.
 134    if(iopt<0) go to 440
        fpms = fp-s
        if(abs(fpms)<acc) go to 440
      !  if f(p=inf) < s accept the choice of knots.
        if(fpms<0.) go to 250
      !  if n = nmax, sinf(u) is an interpolating spline curve.
        if(n==nmax) go to 430
      !  increase the number of knots.
      !  if n=nest we cannot increase the number of knots because of
      !  the storage capacity limitation.
        if(n==nest) go to 420
      !  determine the number of knots nplus we are going to add.
        if(ier==0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp>acc) npl1 = int(rn*fpms/(fpold-fp))
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
      !  compute the sum of squared residuals for each knot interval
      !  t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = zero
        i = 1
        l = k2
        new = 0
        jj = (mb-1)*idim
        do 180 it=mb,me
          if(u(it)<t(l) .or. l>nk1) go to 160
          new = 1
          l = l+1
 160      term = zero
          l0 = l-k2
          do 175 j2=1,idim
            fac = zero
            j1 = l0
            do 170 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 170        continue
            jj = jj+1
            term = term+(w(it)*(fac-x(jj)))**2
            l0 = l0+n
 175      continue
          fpart = fpart+term
          if(new==0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
      !  add a new knot.
          call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1)
      !  if n=nmax we locate the knots as for interpolation
          if(n==nmax) go to 25
      !  test whether we cannot further increase the number of knots.
          if(n==nest) go to 200
 190    continue
      !  restart the computations with the new set of knots.
 200  continue
      !  test whether the least-squares kth degree polynomial curve is a
      !  solution of our approximation problem.
 250  if(ier==(-2)) go to 440
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 2: determination of the smoothing spline curve sp(u).          c
      !  **********************************************************          c
      !  we have determined the number of knots and their position.          c
      !  we now compute the b-spline coefficients of the smoothing curve     c
      !  sp(u). the observation matrix a is extended by the rows of matrix   c
      !  b expressing that the kth derivative discontinuities of sp(u) at    c
      !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
      !  ponding weights of these additional rows are set to 1/p.            c
      !  iteratively we then have to determine the value of p such that f(p),c
      !  the sum of squared residuals be = s. we already know that the least c
      !  squares kth degree polynomial curve corresponds to p=0, and that    c
      !  the least-squares spline curve corresponds to p=infinity. the       c
      !  iteration process which is proposed here, makes use of rational     c
      !  interpolation. since f(p) is a convex and strictly decreasing       c
      !  function of p, it can be approximated by a rational function        c
      !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
      !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
      !  to calculate the new value of p such that r(p)=s. convergence is    c
      !  guaranteed by taking f1>0 and f3<zero                                 c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  evaluate the discontinuity jump of the kth derivative of the
      !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
      !  initial value for p.
      p1 = zero
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = zero
      do 252 i=1,nn
         p = p+a(i,1)
 252  continue
      rn = nn
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
      !  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
      !  the rows of matrix b with weight 1/p are rotated into the
      !  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 255 i=1,nc
          c(i) = z(i)
 255    continue
        g(1:nn,1:k1) = a(1:nn,1:k1)
        g(1:nn,  k2) = zero
        do 300 it=1,n8
      !  the row of matrix b is rotated into triangle by givens transformation
          do 264 i=1,k2
            h(i) = b(it,i)*pinv
 264      continue
          do 268 j=1,idim
            xi(j) = zero
 268      continue
      !  take into account that certain b-spline coefficients must be zero.
          if(it>ib) go to 274
          j1 = ib-it+2
          j2 = 1
          do 270 i=j1,k2
            h(j2) = h(i)
            j2 = j2+1
 270      continue
          do 272 i=j2,k2
            h(i) = zero
 272      continue
 274      jj = max0(1,it-ib)
          do 290 j=jj,nn
            piv = h(1)
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
      !  transformations to right hand side.
            j1 = j
            do 277 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 277        continue
            if(j==nn) go to 300
            i2 = min0(nn-j,k1)
            do 280 i=1,i2
      !  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = zero
 290      continue
 300    continue
      !  backward substitution to obtain the b-spline coefficients.
        j1 = 1
        do 308 j2=1,idim
          j3 = j1+ib
          c(j3:j3+nn-1) = fpback(g,c(j1),nn,k2,nest)
          if(ib==0) go to 306
          j3 = j1
          do 304 i=1,ib
            c(j3) = zero
            j3 = j3+1
 304      continue
 306      j1 =j1+n
 308    continue
      !  computation of f(p).
        fp = zero
        l = k2
        jj = (mb-1)*idim
        do 330 it=mb,me
          if(u(it)<t(l) .or. l>nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = zero
          do 325 j2=1,idim
            fac = zero
            j1 = l0
            do 320 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 320        continue
            jj = jj+1
            term = term+(fac-x(jj))**2
            l0 = l0+n
 325      continue
          fp = fp+term*w(it)**2
 330    continue
      !  test whether the approximation sp(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 440
      !  test whether the maximal number of iterations is reached.
        if(iter==maxit) go to 400
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 340
        if((f2-f3)>acc) go to 335
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2<0.) ich3=1
 340    if(ich1/=0) go to 350
        if((f1-f2)>acc) go to 345
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 360
        if(p>=p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2>0.) ich1=1
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 350    if(f2>=f1 .or. f2<=f3) go to 410
      !  find the new value for p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 360  continue
      !  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end subroutine fpcons


      pure subroutine fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,a, &
                             b,const,z,zz,u,q,info,up,left,right,jbind,ibind,ier)

      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(out) :: sq
      integer,     intent(in)  :: m,n,maxtr,maxbin,nm,mb
      integer,     intent(out) :: ier
      !  ..array arguments..
      real(RKIND), intent(in)  :: x(m),y(m),w(m)
      real(RKIND), intent(inout) :: t(n),e(n),c(n),sx(m),a(n,4),b(nm,maxbin),const(n),z(n),zz(n),u(maxbin),q(m,4)
      integer, intent(inout)   :: info(maxtr),up(maxtr),left(maxtr),right(maxtr)
      integer, intent(out)     :: ibind(mb),jbind(mb)
      logical, intent(out)     :: bind(n)
      !  ..local scalars..
      integer :: count,i,i1,j,j1,j2,j3,k,kdim,k1,k2,k3,k4,k5,k6,l,l1,l2,l3,merk,nbind,violated,n1,n4,n6
      real(RKIND) :: f,wi,xi
      !  ..local array..
      real(RKIND) :: h(SIZ_K+1)
      !  ..subroutine references..
      !    fpbspl,fpadno,fpdeno,fpfrno,fpseno
      !  ..

      !  *****
      !  if we use the b-spline representation of s(x) our approximation problem results in a quadratic
      !  programming problem: find the b-spline coefficients c(j),j=1,2,...n-4 such that
      !        (1) sumi((wi*(yi-sumj(cj*nj(xi))))**2),i=1,2,...m is minimal
      !        (2) sumj(cj*n''j(t(l+3)))*e(l) <= 0, l=1,2,...n-6.
      !  to solve this problem we use the theil-van de panne procedure.
      !  if the inequality constraints (2) are numbered from 1 to n-6, this algorithm finds a subset of
      !  constraints ibind(1)..ibind(nbind) such that the solution of the minimization problem (1) with
      !  these constraints in equality form, satisfies all constraints. such a feasible solution is optimal
      !  if the lagrange parameters associated with that problem with equality constraints are all positive.
      !  *****

      !  determine n6, the number of inequality constraints.
      n6 = n-6
      n4 = n-4

      !  fix the parameters which determine these constraints.
      forall (i=1:n6) const(i) = e(i)*(t(i+4)-t(i+1))/(t(i+5)-t(i+2))

      !  initialize the triply linked tree which is used to find the subset of constraints
      !  ibind(1),...ibind(nbind).
      count    = 1
      info(1)  = 0
      left(1)  = 0
      right(1) = 0
      up(1)    = 1
      merk     = 1

      !  set up the normal equations  n'nc=n'y  where n denotes the m x (n-4) observation matrix with
      !  elements ni,j = wi*nj(xi)  and y is the column vector with elements yi*wi.
      !  from the properties of the b-splines nj(x),j=1,2,...n-4, it follows that  n'n  is a (n-4) x (n-4)
      !  positive definite bandmatrix of bandwidth 7. the matrices n'n and n'y are built up in a and z.

      !  initialization
      z(1:n4) = zero
      a(1:n4,1:4) = zero
      l   = 4

      rows: do i=1,m
          ! fetch the current row of the observation matrix.
          xi = x(i)
          wi = w(i)**2

          ! search for knot interval  t(l) <= xi < t(l+1)
          do while (xi>=t(l+1) .and. l/=n4)
            l  = l+1
          end do

          ! evaluate the four non-zero cubic b-splines nj(xi),j=l-3,...l.
          call fpbspl(t,n,3,xi,l,h) ! h is a temporary

          ! store in q these values h(1),h(2),...h(4).
          q(i,1:4) = h(1:4)

          ! add the contribution of the current row of the observation matrix n to the normal equations.
          l3 = l-3
          k1 = 0
          add_row: do j1 = l3,l
            k1 = k1+1
            f = h(k1)
            z(j1) = z(j1)+f*wi*y(i)
            k2 = k1
            j2 = 4
            do j3 = j1,l
              a(j3,j2) = a(j3,j2)+f*wi*h(k2)
              k2 = k2+1
              j2 = j2-1
            end do
          end do add_row
      end do rows

      !  since n'n is a symmetric matrix it can be factorized as
      !        (3)  n'n = (r1)'(d1)(r1)
      !  with d1 a diagonal matrix and r1 an (n-4) x (n-4)  unit upper triangular matrix of bandwidth 4.
      !  matrices r1 and d1 are built up in a. at the same time we solve the systems of equations
      !        (4)  (r1)'(z2) = n'y
      !        (5)  (d1) (z1) = (z2)
      !  the vectors z2 and z1 are kept in zz and z.
      solve: do i=1,n4
         k1 = max(5-i,1)
         k2 = i-4+k1
         k3 = k2

         hundred: do j=k1,4
            k4 = j-1
            k5 = 4-j+k1
            f  = a(i,j)
            if (k1<=k4) then
                k6 = k2
                do k=k1,k4
                  f = f-a(i,k)*a(k3,k5)*a(k6,4)
                  k5 = k5+1
                  k6 = k6+1
                end do
            end if
            if (j<4) then
               a(i,j) = f/a(k3,4)
               k3 = k3+1
            endif
         end do hundred
         a(i,4) = f
         f = z(i)
         if(i>1) then
            k4 = i
            do j=k1,3
               k = k1+3-j
               k4 = k4-1
               f = f-a(i,k)*z(k4)*a(k4,4)
            end do
         endif
         z (i) = f/a(i,4)
         zz(i) = f
      end do solve

      !  start computing the least-squares cubic spline without taking account
      !  of any constraint.
      nbind    = 0
      n1       = 1
      ibind(1) = 0

      !  main loop for the least-squares problems with different subsets of the constraints (2) in equality
      !  form. the resulting b-spline coeff. c and lagrange parameters u are the solution of the system
      !            ! n'n  b' ! ! c !   ! n'y !
      !        (6) !         ! !   ! = !     !
      !            !  b   0  ! ! u !   !  0  !
      !  z1 is stored into array c.
      least_squares: do

          c(1:n4) = z(1:n4)

          !  if there are no equality constraints, compute the coeff. c directly.
          has_constraints: if (nbind>0) then
              !  initialization
              kdim = n4+nbind
              b(1:kdim,1:nbind) = zero

              !  matrix b is built up,expressing that the constraints nrs ibind(1),...
              !  ibind(nbind) must be satisfied in equality form.
              do i=1,nbind
                 l        = ibind(i)
                 b(l,i)   = e(l)
                 b(l+1,i) = -(e(l)+const(l))
                 b(l+2,i) = const(l)
              end do

              !  find the matrix (b1) as the solution of the system of equations
              !        (7)  (r1)'(d1)(b1) = b'
              !  (b1) is built up in the upper part of the array b(rows 1,...n-4).
              make_b1: do k1=1,nbind
                 l = ibind(k1)
                 do i=l,n4
                    f = b(i,k1)
                    if (i/=1) then
                        k2 = min(i-1,3)
                        do k3=1,k2
                          l1 = i-k3
                          l2 = 4-k3
                          f = f-b(l1,k1)*a(i,l2)*a(l1,4)
                        end do
                    endif
                    b(i,k1) = f/a(i,4)
                 end do
              end do make_b1

              !  factorization of the symmetric matrix  -(b1)'(d1)(b1)
              !        (8)  -(b1)'(d1)(b1) = (r2)'(d2)(r2)
              !  with (d2) a diagonal matrix and (r2) an nbind x nbind unit upper triangular matrix.
              !  the matrices r2 and d2 are built up in the lower part of the array b (rows n-3,n-2,
              !  ...n-4+nbind).
              factor: do i=1,nbind
                 i1 = i-1
                 do j=i,nbind
                    f = zero
                    do k=1,n4
                      f = f+b(k,i)*b(k,j)*a(k,4)
                    end do
                    k1 = n4+1
                    if (i1/=0) then
                       do k=1,i1
                          f = f+b(k1,i)*b(k1,j)*b(k1,k)
                          k1 = k1+1
                       end do
                    endif
                    b(k1,j) = -f*merge(one/b(k1,i),one,i/=j)
                 end do
              end do factor

              !  according to (3),(7) and (8) the system of equations (6) becomes
              !         ! (r1)'    0  ! ! (d1)    0  ! ! (r1)  (b1) ! ! c !   ! n'y !
              !    (9)  !             ! !            ! !            ! !   ! = !     !
              !         ! (b1)'  (r2)'! !   0   (d2) ! !   0   (r2) ! ! u !   !  0  !
              !  backward substitution to obtain the b-spline coefficients c(j),j=1,..n-4 and the
              !  lagrange parameters u(j),j=1,2,...nbind. first step of the backward substitution:
              !  solve the system
              !             ! (r1)'(d1)      0     ! ! (c1) !   ! n'y !
              !        (10) !                      ! !      ! = !     !
              !             ! (b1)'(d1)  (r2)'(d2) ! ! (u1) !   !  0  !
              !  from (4) and (5) we know that this is equivalent to
              !        (11)  (c1) = (z1)
              !        (12)  (r2)'(d2)(u1) = -(b1)'(z2)

              solve2: do i=1,nbind
                 f  = dot_product(b(1:n4,i),zz(1:n4))
                 i1 = i-1
                 k1 = n4+1
                 if (i1/=0) then ! [FP] this is probably unnecessary
                    do j=1,i1
                       f = f+u(j)*b(k1,i)*b(k1,j)
                       k1 = k1+1
                    end do
                 end if
                 u(i) = -f/b(k1,i)
              end do solve2

              !  second step of the backward substitution: solve the system
              !             ! (r1)  (b1) ! ! c !   ! c1 !
              !        (13) !            ! !   ! = !    !
              !             !   0   (r2) ! ! u !   ! u1 !
              k1 = nbind
              k2 = kdim
              !  find the lagrange parameters u.
              finish2: do i=1,nbind
                f = u(k1)
                if (i/=1) then
                   k3 = k1+1
                   do j=k3,nbind
                     f = f-u(j)*b(k2,j)
                   end do
                endif
                u(k1) = f
                k1 = k1-1
                k2 = k2-1
              end do finish2

              !  find the b-spline coefficients c.
              c(1:n4) = c(1:n4) - matmul(b(1:n4,1:nbind),u(1:nbind))

          end if has_constraints

          k1 = n4
          do i=2,n4
             k1 = k1-1
             f  = c(k1)
             k2 = max(5-i,1)
             k3 = k1
             l = 3
             do j=k2,3
               k3 = k3+1
               f = f-a(k3,l)*c(k3)
               l = l-1
             end do
             c(k1) = f
          end do

          !  test whether the solution of the least-squares problem with the constraints ibind(1),...
          !  ibind(nbind) in equality form, satisfies all of the constraints (2).
          k = 1

          !  number counts the number of violated inequality constraints.
          violated = 0

          test_constraints: do j=1,n6
            l = ibind(k)
            k = k+1
            if (j==l) cycle test_constraints
            k = k-1
            ! test whether constraint j is satisfied

            f = e(j)*(c(j)-c(j+1))+const(j)*(c(j+2)-c(j+1))
            if(f<=zero) cycle test_constraints

            ! if constraint j is not satisfied, add a branch of length nbind+1 to the tree. the nodes
            ! of this branch contain in their information field the number of the constraints ibind(1),
            ! ...ibind(nbind) and j, arranged in increasing order.

            violated = violated+1
            k1 = k-1
            if (k1>0) jbind(1:k1) = ibind(1:k1)
            jbind(k) = j
            if (l>0) jbind(k+1:nbind+1) = ibind(k:nbind)

            call fpadno(maxtr,up,left,right,info,count,merk,jbind,n1,ier)

            ! test whether the storage space which is required for the tree,
            ! exceeds the available storage space.
            if (ier/=0) then
                ier = FITPACK_S_TOO_SMALL
                return
            end if
          end do test_constraints

          !  test whether the solution of the least-squares problem with equality
          !  constraints is a feasible solution.
          if (violated==0) then

              !  test whether the feasible solution is optimal.
              bind(1:n6) = .false.
              if (nbind>0) bind(ibind(1:nbind)) = u(1:nbind)>zero

              if (all(.not.bind(1:n6))) exit least_squares

          end if

          !  test whether there are still cases with nbind constraints in
          !  equality form to be considered.
          if (merk<=1) then
              nbind = n1
              !  test whether the number of knots where s''(x)=0 exceeds maxbin.
              if(nbind>maxbin) then
                 ier = 1
                 return
              end if
              n1 = n1+1
              ibind(n1) = 0
              !  search which cases with nbind constraints in equality form
              !  are going to be considered.
              call fpdeno(maxtr,up,left,right,nbind,merk)
              !  test whether the quadratic programming problem has a solution.
              if (merk==1) then
                 ier = 3
                 return
              end if
          endif

          ! find a new case with nbind constraints in equality form.
          call fpseno(maxtr,up,left,right,info,merk,ibind,nbind)

      end do least_squares

      !  test whether the feasible solution is optimal.
      ier = FITPACK_OK
      bind(1:n6) = .false.
      if (nbind>0) bind(ibind(1:nbind)) = u(1:nbind)>zero

      !  evaluate s(x) at the data points x(i) and calculate the weighted
      !  sum of squared residual right hand sides sq.
      sq = zero
      l   = 4
      evaluate_error: do i=1,m
         do while (x(i)>=t(l+1) .and. l/=n4)
            l = l+1
         end do
         sx(i) = c(l-3)*q(i,1)+c(l-2)*q(i,2)+c(l-1)*q(i,3)+c(l)*q(i,4)
         sq = sq+(w(i)*(y(i)-sx(i)))**2
      end do evaluate_error

      return

      end subroutine fpcosp


      pure elemental subroutine fpcsin(a,b,par,sia,coa,sib,cob,ress,resc)

      !  fpcsin calculates the integrals ress=integral((b-x)**3*sin(par*x))
      !  and resc=integral((b-x)**3*cos(par*x)) over the interval (a,b),
      !  given sia=sin(par*a),coa=cos(par*a),sib=sin(par*b) and cob=cos(par*b)
      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(in)  :: a,b,par,sia,coa,sib,cob
      real(RKIND), intent(out) :: ress,resc
      !  ..local scalars..
      integer :: i,j
      real(RKIND) :: ab,ab4,ai,alfa,beta,b2,b4,fac,f1,f2
      real(RKIND), parameter :: eps = smallnum10

      ab   = b-a
      ab4  = ab**4
      alfa = ab*par
      ! the way of calculating the integrals ress and resc depends on
      ! the value of alfa = (b-a)*par.

      if(abs(alfa)<=one) then
         ! ress and resc are found by evaluating a series expansion.
         fac = fourth
         f1  = fac
         f2  = zero
         i   = 4
         series: do j=1,5
           i   = i+1
           ai  = i
           fac = fac*alfa/ai
           f2  = f2+fac
           if (abs(fac)<=eps) exit series
           i   = i+1
           ai  = i
           fac = -fac*alfa/ai
           f1  = f1+fac
           if (abs(fac)<=eps) exit series
         end do series
         ress = ab4*(coa*f2+sia*f1)
         resc = ab4*(coa*f1-sia*f2)

      else
         ! integration by parts.
         beta = one/alfa
         b2   = beta**2
         b4   = six*b2**2
         f1   = three*b2*(one-two*b2)
         f2   = beta*(one-six*b2)
         ress = ab4*(coa*f2+sia*f1+sib*b4)
         resc = ab4*(coa*f1-sia*f2+cob*b4)
         return
      endif

      end subroutine fpcsin


      recursive subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol, &
         maxit,k1,k2,n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)

      !  ..
      !  ..scalar arguments..
      real(RKIND) xb,xe,s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),w(m),t(nest),c(nest),fpint(nest), &
       z(nest),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
      !  ..local scalars..
      real(RKIND) acc,cos,fpart,fpms,fpold,fp0,f1,f2,f3, &
       p,pinv,piv,p1,p2,p3,rn,sin,store,term,wi,xi,yi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,k3,l,l0, &
       mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
      !  ..local arrays..
      real(RKIND) h(SIZ_K+1)

      !  set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 1: determination of the number of knots and their position     c
      !  **************************************************************      c
      !  given a set of knots we compute the least-squares spline sinf(x),   c
      !  and the corresponding sum of squared residuals fp=f(p=inf).         c
      !  if iopt=-1 sinf(x) is the requested approximation.                  c
      !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
      !    if fp <=s we will continue with the current set of knots.         c
      !    if fp > s we will increase the number of knots and compute the    c
      !       corresponding least-squares spline until finally fp<=s.        c
      !    the initial choice of knots depends on the value of s and iopt.   c
      !    if s=0 we have spline interpolation; in that case the number of   c
      !    knots equals nmax = m+k+1.                                        c
      !    if s > 0 and                                                      c
      !      iopt=0 we first compute the least-squares polynomial of         c
      !      degree k; n = nmin = 2*k+2                                      c
      !      iopt=1 we start with the set of knots found at the last         c
      !      call of the routine, except for the case that s > fp0; then     c
      !      we compute directly the least-squares polynomial of degree k.   c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      if(iopt<0) go to 60
      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      !  determine nmax, the number of knots for spline interpolation.
      nmax = m+k1
      if(s>zero) go to 45
      !  if s=0, s(x) is an interpolating spline.
      !  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax>nest) go to 420
      !  find the position of the interior knots in case of interpolation.
  10  mk1 = m-k1
      if(mk1==0) go to 60
      k3 = k/2
      i = k2
      j = k3+2
      if(k3*2==k) go to 30
      do 20 l=1,mk1
        t(i) = x(j)
        i = i+1
        j = j+1
  20  continue
      go to 60
  30  do 40 l=1,mk1
        t(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  40  continue
      go to 60
      !  if s>0 our initial choice of knots depends on the value of iopt.
      !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
      !  polynomial of degree k which is a spline without interior knots.
      !  if iopt=1 and fp0>s we start computing the least squares spline
      !  according to the set of knots found at the last call of the routine.
  45  if(iopt==0) go to 50
      if(n==nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0>s) go to 60
  50  n = nmin
      fpold = zero
      nplus = 0
      nrdata(1) = m-2
      !  main loop for the different sets of knots. m is a save upper bound
      !  for the number of trials.
  60  do 200 iter = 1,m
        if(n==nmin) ier = -2
      !  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
      !  find the position of the additional knots which are needed for
      !  the b-spline representation of s(x).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = xb
          t(i) = xe
          i = i-1
  70    continue
      !  compute the b-spline coefficients of the least-squares spline
      !  sinf(x). the observation matrix a is built up row by row and
      !  reduced to upper triangular form by givens transformations.
      !  at the same time fp=f(p=inf) is computed.
        fp = zero
      !  initialize the observation matrix a.
        z(1:nk1) = zero
        a(1:nk1,1:k1) = zero
        l = k1
        do 130 it=1,m
      !  fetch the current data point x(it),y(it).
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
      !  search for knot interval t(l) <= xi < t(l+1).
  85      if(xi<t(l+1) .or. l==nk1) go to 90
          l = l+1
          go to 85
      !  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  90      call fpbspl(t,n,k,xi,l,h)
          do 95 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  95      continue
      !  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
            j = j+1
            piv = h(i)
            if(piv==zero) go to 110
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
      !  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i==k1) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,k1
              i2 = i2+1
      !  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
 120      fp = fp+yi*yi
 130    continue
        if(ier==(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
      !  backward substitution to obtain the b-spline coefficients.
        c(:nk1) = fpback(a,z,nk1,k1,nest)
      !  test whether the approximation sinf(x) is an acceptable solution.
        if(iopt<0) go to 440
        fpms = fp-s
        if(abs(fpms)<acc) go to 440
      !  if f(p=inf) < s accept the choice of knots.
        if(fpms<zero) go to 250
      !  if n = nmax, sinf(x) is an interpolating spline.
        if(n==nmax) go to 430
      !  increase the number of knots.
      !  if n=nest we cannot increase the number of knots because of
      !  the storage capacity limitation.
        if(n==nest) go to 420
      !  determine the number of knots nplus we are going to add.
        if(ier==0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp>acc) npl1 = int(rn*fpms/(fpold-fp))
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
      !  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
      !  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = zero
        i = 1
        l = k2
        new = 0
        do 180 it=1,m
          if(x(it)<t(l) .or. l>nk1) go to 160
          new = 1
          l = l+1
 160      term = zero
          l0 = l-k2
          do 170 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 170      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new==0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
      !  add a new knot.
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
      !  if n=nmax we locate the knots as for interpolation.
          if(n==nmax) go to 10
      !  test whether we cannot further increase the number of knots.
          if(n==nest) go to 200
 190    continue
      !  restart the computations with the new set of knots.
 200  continue
      !  test whether the least-squares kth degree polynomial is a solution
      !  of our approximation problem.
 250  if(ier==(-2)) go to 440
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 2: determination of the smoothing spline sp(x).                c
      !  ***************************************************                 c
      !  we have determined the number of knots and their position.          c
      !  we now compute the b-spline coefficients of the smoothing spline    c
      !  sp(x). the observation matrix a is extended by the rows of matrix   c
      !  b expressing that the kth derivative discontinuities of sp(x) at    c
      !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
      !  ponding weights of these additional rows are set to 1/p.            c
      !  iteratively we then have to determine the value of p such that      c
      !  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
      !  the least-squares kth degree polynomial corresponds to p=0, and     c
      !  that the least-squares spline corresponds to p=infinity. the        c
      !  iteration process which is proposed here, makes use of rational     c
      !  interpolation. since f(p) is a convex and strictly decreasing       c
      !  function of p, it can be approximated by a rational function        c
      !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
      !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
      !  to calculate the new value of p such that r(p)=s. convergence is    c
      !  guaranteed by taking f1>0 and f3<zero                                 c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  evaluate the discontinuity jump of the kth derivative of the
      !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
      !  initial value for p.
      p1 = zero
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 255 i=1,nk1
         p = p+a(i,1)
 255  continue
      rn = nk1
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
      !  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
      !  the rows of matrix b with weight 1/p are rotated into the
      !  triangularised observation matrix a which is stored in g.
        pinv = one/p
        c(1:nk1) = z(1:nk1)
        g(1:nk1,1:k1) = a(1:nk1,1:k1)
        g(1:nk1,  k2) = zero
        do 300 it=1,n8
      !  the row of matrix b is rotated into triangle by givens transformation
          h(1:k2) = b(it,1:k2)*pinv
          yi = zero
          do 290 j=it,nk1
            piv = h(1)
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
      !  transformations to right hand side.
            call fprota(cos,sin,yi,c(j))
            if(j==nk1) go to 300
            i2 = k1
            if(j>n8) i2 = nk1-j
            do 280 i=1,i2
      !  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = zero
 290      continue
 300    continue
      !  backward substitution to obtain the b-spline coefficients.
        c(:nk1) = fpback(g,c,nk1,k2,nest)
      !  computation of f(p).
        fp = zero
        l = k2
        do 330 it=1,m
          if(x(it)<t(l) .or. l>nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = zero
          do 320 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 320      continue
          fp = fp+(w(it)*(term-y(it)))**2
 330    continue
      !  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 440
      !  test whether the maximal number of iterations is reached.
        if(iter==maxit) go to 400
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 340
        if((f2-f3)>acc) go to 335
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2<zero) ich3=1
 340    if(ich1/=0) go to 350
        if((f1-f2)>acc) go to 345
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 360
        if(p>=p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2>zero) ich1=1
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 350    if(f2>=f1 .or. f2<=f3) go to 410
      !  find the new value for p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 360  continue
      !  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end subroutine fpcurf


      !  subroutine fpcuro finds the real zeros of a cubic polynomial p(x) = a*x**3+b*x**2+c*x+d.
      pure subroutine fpcuro(a,b,c,d,x,n)

      !
      !  calling sequence:
      !     call fpcuro(a,b,c,d,x,n)
      !
      !  input parameters:
      !    a,b,c,d: real values, containing the coefficients of p(x).
      !
      !  output parameters:
      !    x      : real array,length 3, which contains the real zeros of p(x)
      !    n      : integer, giving the number of real zeros of p(x).
      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(in) :: a,b,c,d
      integer,     intent(out) :: n
      !  ..array argument..
      real(RKIND), intent(out) :: x(3)
      !  ..local scalars..
      integer :: i
      real(RKIND) :: a1,b1,c1,df,disc,d1,f,p3,q,r,step,u,u1,u2,y

      !  set constants
      real(RKIND), parameter :: ovfl = 1.0e4_RKIND
      real(RKIND), parameter :: tent = 0.1_RKIND
      real(RKIND), parameter :: pi3  = datan(one)/0.75_RKIND
      real(RKIND), parameter :: e3   = tent/0.3_RKIND

      a1 = abs(a)
      b1 = abs(b)
      c1 = abs(c)
      d1 = abs(d)

      if (max(b1,c1,d1)<a1*ovfl) then

         !  p(x) is a third degree polynomial.
         b1 = b/a*e3
         c1 = c/a
         d1 = d/a
         q = c1*e3-b1*b1
         r = b1*b1*b1+(d1-b1*c1)*half
         disc = q*q*q+r*r
         if (disc>zero) then

            u = sqrt(disc)
            u1 = -r+u
            u2 = -r-u
            n = 1
            x(1) = sign(abs(u1)**e3,u1)+sign(abs(u2)**e3,u2)-b1

         else

            u = sign(sqrt(abs(q)),r)
            p3 = atan2(sqrt(-disc),abs(r))*e3
            u2 = u+u
            n = 3
            x(1) = -u2*cos(p3)-b1
            x(2) = u2*cos(pi3-p3)-b1
            x(3) = u2*cos(pi3+p3)-b1

         end if

      elseif(max(c1,d1)<b1*ovfl) then

         !  p(x) is a second degree polynomial.
         disc = c*c-four*b*d

         if (disc<zero) then
            n = 0
            return
         else
            n = 2
            u = sqrt(disc)
            b1 = b+b
            x(1) = (-c+u)/b1
            x(2) = (-c-u)/b1
         endif

      elseif (d1<c1*ovfl) then

         !  p(x) is a first degree polynomial.
         n = 1
         x(1) = -d/c

      else

         !  p(x) is a constant function.
         n = 0
         return

      end if

      ! apply a newton iteration to improve the accuracy of the roots.
      do i=1,n
        y    = x(i)
        f    = ((a*y+b)*y+c)*y+d
        df   = (three*a*y+two*b)*y+c
        step = merge(f/df,zero,abs(f)<abs(df)*tent)
        x(i) = y-step
      end do

      end subroutine fpcuro

      ! (l u)-decomposition of a cyclic tridiagonal matrix with the non-zero
      ! elements stored as follows
      !
      !    | a(1,2) a(1,3)                                    a(1,1)  |
      !    | a(2,1) a(2,2) a(2,3)                                     |
      !    |        a(3,1) a(3,2) a(3,3)                              |
      !    |               ...............                            |
      !    |                               a(n-1,1) a(n-1,2) a(n-1,3) |
      !    | a(n,3)                                  a(n,1)   a(n,2)  |
      pure subroutine fpcyt1(a,n,nn)

          !  ..scalar arguments..
          integer, intent(in) :: n,nn
          !  ..array arguments..
          real(RKIND), intent(inout) :: a(nn,6)
          !  ..local scalars..
          real(RKIND) aa,beta,gamma,sum,teta,v
          integer i,n1,n2
          !  ..
          n2   = n-2
          beta = one/a(1,2)
          gamma = a(n,3)
          teta = a(1,1)*beta
          a(1,4) = beta
          a(1,5) = gamma
          a(1,6) = teta
          sum = gamma*teta

          internal_rows: do i=2,n2
             v = a(i-1,3)*beta
             aa = a(i,1)
             beta = one/(a(i,2)-aa*v)
             gamma = -gamma*v
             teta = -teta*aa*beta
             a(i,4) = beta
             a(i,5) = gamma
             a(i,6) = teta
             sum = sum+gamma*teta
          end do internal_rows

          n1 = n-1
          v = a(n2,3)*beta
          aa = a(n1,1)
          beta = one/(a(n1,2)-aa*v)
          gamma = a(n,1)-gamma*v
          teta = (a(n1,3)-teta*aa)*beta
          a(n1,4) = beta
          a(n1,5) = gamma
          a(n1,6) = teta
          a(n,4) = one/(a(n,2)-(sum+gamma*teta))

      return
      end subroutine fpcyt1


      ! subroutine fpcyt2 solves a linear n x n system
      !     a * c = b
      ! where matrix a is a cyclic tridiagonal matrix, decomposed using subroutine fpsyt1.
      pure subroutine fpcyt2(a,n,b,c,nn)
          !  ..scalar arguments..
          integer, intent(in) :: n,nn
          !  ..array arguments..
          real(RKIND), intent(in)  :: a(nn,6),b(n)
          real(RKIND), intent(out) :: c(n)
          !  ..local scalars..
          real(RKIND) :: cc,sum
          integer :: i,j,j1,n1
          !  ..
          c(1) = b(1)*a(1,4)
          sum  = c(1)*a(1,5)
          n1   = n-1
          do i=2,n1
             c(i) = (b(i)-a(i,1)*c(i-1))*a(i,4)
             sum = sum+c(i)*a(i,5)
          end do
          cc    = (b(n)-sum)*a(n,4)
          c(n)  = cc
          c(n1) = c(n1)-cc*a(n1,6)
          j = n1
          do i=3,n
             j1 = j-1
             c(j1) = c(j1)-c(j)*a(j1,3)*a(j1,4)-cc*a(j1,6)
             j = j1
          end do
          return

      end subroutine fpcyt2


      ! subroutine fpdeno frees the nodes of all branches of a triply linked tree with length < nbind
      ! by putting to zero their up field. on exit the parameter merk points to the terminal node of
      ! the most left branch of length nbind or takes the value 1 if there is no such branch.
      pure subroutine fpdeno(maxtr,up,left,right,nbind,merk)
      !  ..
      !  ..scalar arguments..
      integer, intent(in)  :: maxtr,nbind
      integer, intent(inout) :: merk
      !  ..array arguments..
      integer, intent(inout) :: up(maxtr),left(maxtr),right(maxtr)
      !  ..local scalars ..
      integer i,j,k,l,level,point
      !  ..
      i = 1
      level = 0
  10  point = i
      i = left(point)
      if (i==0) go to 20
      level = level+1
      go to 10
  20  if(level==nbind) go to 70
  30  i = right(point)
      j = up(point)
      up(point) = 0
      k = left(j)
      if (point/=k) go to 50
      if (i==0) then
          level = level-1
          if (level==0) go to 80
          point = j
          go to 30
      endif
  40  left(j) = i
      go to 10
  50  l = right(k)
      if (point/=l) then
          k = l
          go to 50
      endif
  60  right(k) = i
      point = k
  70  i = right(point)
      if(i/=0) go to 10
      i = up(point)
      level = level-1
      if (level/=0) then
         point = i
         go to 70
      endif
  80  k = 1
      l = left(k)
      if (up(l)/=0) then
          do while (k/=0)
             merk = k
             k = left(k)
          end do
      endif
      return
      end subroutine fpdeno


      !  subroutine fpdisc calculates the discontinuity jumps of the kth
      !  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
      pure subroutine fpdisc(t,n,k2,b,nest)

      !  ..scalar arguments..
      integer, intent(in) :: n,k2,nest
      !  ..array arguments..
      real(RKIND), intent(in) :: t(n)
      real(RKIND), intent(out) :: b(nest,k2)
      !  ..local scalars..
      real(RKIND) :: an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
      !  ..local array..
      real(RKIND) :: h(12)
      !  ..
      k1    = k2-1
      k     = k1-1
      nk1   = n-k1
      nrint = nk1-k
      an    = nrint
      fac   = an/(t(nk1+1)-t(k1))
      do l=k2,nk1
         lmk = l-k1
         do j=1,k1
            ik = j+k1
            lj = l+j
            lk = lj-k2
            h(j) = t(l)-t(lk)
            h(ik) = t(l)-t(lj)
         end do
         lp = lmk
         do j=1,k2
           jk = j
           prod = h(j)
           do i=1,k
             jk = jk+1
             prod = prod*h(jk)*fac
           end do
           lk = lp+k1
           b(lmk,j) = (t(lk)-t(lp))/prod
           lp = lp+1
         end do
      end do
      return
      end subroutine fpdisc


      !  subroutine fpfrno collects the free nodes (up field zero) of the triply linked tree the
      !  information of which is kept in the arrays up,left,right and info. the maximal length of the
      !  branches of the tree is given by n1. if no free nodes are found, the error flag ier is set
      !  to 1.
      pure subroutine fpfrno(maxtr,up,left,right,info,point,merk,n1,count,ier)

      !  ..scalar arguments..
      integer, intent(in)    :: maxtr,n1
      integer, intent(inout) :: point,merk
      integer, intent(out)   :: count,ier
      !  ..array arguments..
      integer, intent(inout) :: up(maxtr),left(maxtr),right(maxtr),info(maxtr)
      !  ..local scalars
      integer :: i,j,k,l,n,level
      !  ..
      ier    = 1
      level  = 1
      count  = 2
      if (n1==2) return ! No free nodes

      free_nodes_left: do while (level<=n1)
         j = 0
         i = 1
         k = 0
         l = 0
      20 march_left: do while (j/=level)
            l = left(i)
            if (l==0) exit march_left
            i = l
            j = j+1
         end do march_left
      30 if (i<count .or. l==0) then
             ! continue
         elseif (i==count) then
             count = count+1
             ! continue
         elseif (up(count)==0) then
             up   (count) = up   (i)
             left (count) = left (i)
             right(count) = right(i)
             info (count) = info (i)
             if( merk==i) merk  = count
             if(point==i) point = count
             if (k==0) then
                n = up(i)
                left(n) = count
             else
                right(k) = count
             endif
             l = left(i)
             do while (l/=0)
               up(l) = count
               l = right(l)
             end do
             up(i) = 0
             i = count
             count = count+1
         else
             count = count+1
             go to 30
         endif

         hundred10: do
             l = right(i)
             k = i

             if (l==0) then
                 l = up(i)
                 j = j-1
                 if (j==0) then
                     exit hundred10
                 else
                     i = l
                     cycle hundred10
                 endif
             else
                 i = l
                 go to 20
             endif
         end do hundred10

         level = level+1

      end do free_nodes_left

      if(count<=maxtr) ier = FITPACK_OK
      return
      end subroutine fpfrno


      !  subroutine fpgivs calculates the parameters of a givens transformation .
      elemental subroutine fpgivs(piv,ww,cos,sin)
          real(RKIND), intent(in)    :: piv
          real(RKIND), intent(inout) :: ww
          real(RKIND), intent(out)   :: cos,sin
          !  ..local scalars..
          real(RKIND) :: dd,store

          store = abs(piv)
          dd  = merge(store*sqrt(one+(ww/piv)**2), &
                      ww   *sqrt(one+(piv/ww)**2), store>=ww)
          cos = ww/dd
          sin = piv/dd
          ww  = dd
          return
      end subroutine fpgivs


      ! Compute spline coefficients on a rectangular grid
      recursive subroutine fpgrdi(ifsu,ifsv,ifbu,ifbv,iback,u,mu,v, &
       mv,z,mz,dz,iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm, &
       mvnu,spu,spv,right,q,au,av1,av2,bu,bv,aa,bb,cc,cosi,nru,nrv)

      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(out)   :: sq
      real(RKIND), intent(inout) :: fp  ! computed if iback==0
      real(RKIND), intent(out)   :: p
      integer :: ifsu,ifsv,ifbu,ifbv,iback,mu,mv,mz,iop0,iop1,nu,nv,nc,mm,mvnu
      !  ..array arguments..
      real(RKIND), intent(inout) :: fpu(nu),fpv(nv) ! if iback==0
      real(RKIND) :: u(mu),v(mv),z(mz),dz(3),tu(nu),tv(nv),c(nc), &
                     spu(mu,4),spv(mv,4),right(mm),q(mvnu),au(nu,5),av1(nv,6), &
                     av2(nv,4),aa(2,mv),bb(2,nv),cc(nv),cosi(2,nv),bu(nu,5),bv(nv,5)
      integer, intent(inout) :: nru(mu),nrv(mv)
      !  ..local scalars..
      real(RKIND) :: arg,co,dz1,dz2,dz3,fac,fac0,pinv,piv,si,term

      integer :: i,ic,ii,ij,ik,iq,irot,it,iz,i0,i1,i2,i3,j,jj,jk,jper,j0,k,k1,&
                 l,l0,l1,mvv,ncof,nrold,nroldu,nroldv,number,numu,numu1,numv,&
                 numv1,nuu,nu4,nu7,nu8,nu9,nv11,nv4,nv7,nv8,n1

      !  ..local arrays..
      real(RKIND) :: h(SIZ_K+1),h1(5),h2(4)

      !  let
      !               |   (spu)    |            |   (spv)    |
      !        (au) = | ---------- |     (av) = | ---------- |
      !               | (1/p) (bu) |            | (1/p) (bv) |
      !
      !                                | z  ' 0 |
      !                            q = | ------ |
      !                                | 0  ' 0 |
      !
      !  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
      !                coefficients.
      !       z      : the mu x mv matrix which contains the function values.
      !       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
      !                according to the least-squares problems in the u-,resp.
      !                v-direction.
      !       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
      !                containing the discontinuity jumps of the derivatives
      !                of the b-splines in the u-,resp.v-variable at the knots
      !  the b-spline coefficients of the smoothing spline are then calculated
      !  as the least-squares solution of the following over-determined linear
      !  system of equations
      !
      !    (1)  (av) c (au)' = q
      !
      !  subject to the constraints
      !
      !    (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
      !
      !    (3)  if iop0 = 0  c(1,j) = dz(1)
      !            iop0 = 1  c(1,j) = dz(1)
      !                      c(2,j) = dz(1)+(dz(2)*cosi(1,j)+dz(3)*cosi(2,j))*
      !                               tu(5)/3. = cc(j) , j=1,2,...nv-4
      !
      !    (4)  if iop1 = 1  c(nu-4,j) = 0, j=1,2,...,nv-4.
      !
      !  initialization
      nu4  = nu-4
      nu7  = nu-7
      nu8  = nu-8
      nu9  = nu-9
      nv4  = nv-4
      nv7  = nv-7
      nv8  = nv-8
      nv11 = nv-11
      nuu  = nu4-iop0-iop1-1
      pinv = merge(one/p,one,p>zero)

      !  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv and iop0 and
      !  on the value of p whether the matrices (spu), (spv), (bu), (bv) and
      !  (cosi) still must be determined.
      if (ifsu==0) then

          !  calculate the non-zero elements of the matrix (spu) which is the ob-
          !  servation matrix according to the least-squares spline approximation
          !  problem in the u-direction.
          l  = 4
          l1 = 5
          number = 0

          do it=1,mu
            arg = u(it)
            do while (.not.(arg<tu(l1) .or. l==nu4))
               l  = l1
               l1 = l+1
               number = number+1
            end do
            call fpbspl(tu,nu,3,arg,l,h)
            spu(it,1:4) = h(1:4)
            nru(it) = number

          end do
          ifsu = 1

      endif

      !  calculate the non-zero elements of the matrix (spv) which is the ob-
      !  servation matrix according to the least-squares spline approximation
      !  problem in the v-direction.
      if (ifsv==0) then
          l  = 4
          l1 = 5
          number = 0
          do it=1,mv
             arg = v(it)
             do while (.not.(arg<tv(l1) .or. l==nv4))
                l = l1
                l1 = l+1
                number = number+1
             end do
             call fpbspl(tv,nv,3,arg,l,h)
             spv(it,1:4) = h(1:4)
             nrv(it) = number
          end do
          ifsv = 1

          if (iop0==0) then
              !  calculate the coefficients of the interpolating splines for cos(v) and sin(v).
              cosi(:,:nv4) = zero
              if (nv7>=4) then
                  do i=1,nv7
                     l = i+3
                     arg = tv(l)
                     call fpbspl(tv,nv,3,arg,l,h)
                     av1(i,1:3) = h(1:3)
                     cosi(1,i) = cos(arg)
                     cosi(2,i) = sin(arg)
                  end do
                  call fpcyt1(av1,nv7,nv)
                  do j=1,2
                     call fpcyt2(av1,nv7,cosi(j,1:nv7),right,nv)
                     cosi(j,2:nv7+1) = right(1:nv7)
                     cosi(j,1)       = cosi(j,nv7+1)
                     cosi(j,nv7+2)   = cosi(j,2)
                     cosi(j,nv4)     = cosi(j,3)
                  end do
              endif
          endif
      end if

      if (p>zero) then
          !  calculate the non-zero elements of the matrix (bu).
          if (ifbu==0 .and. nu8/=0) then
             call fpdisc(tu,nu,5,bu,nu)
             ifbu = 1
          endif
          !  calculate the non-zero elements of the matrix (bv).
          if (ifbv==0 .and. nv8/=0) then
             call fpdisc(tv,nv,5,bv,nv)
             ifbv = 1
          endif
      endif
      !  substituting (2),(3) and (4) into (1), we obtain the overdetermined system
      !         (5)  (avv) (cr) (auu)' = (qq)
      !  from which the nuu*nv7 remaining coefficients
      !         c(i,j) , i=2+iop0,3+iop0,...,nu-4-iop1 ; j=1,2,...,nv-7 ,
      !  the elements of (cr), are then determined in the least-squares sense.
      !  simultaneously, we compute the resulting sum of squared residuals sq.

      dz1 = dz(1)
      aa (1,1:mv) = dz1
      if (nv8/=0 .and. p>zero) bb(1,1:nv8) = zero
      mvv = mv
      if (iop0/=0) then
          fac = tu(5)/three
          dz2 = dz(2)*fac
          dz3 = dz(3)*fac
          cc(1:nv4) = dz1+dz2*cosi(1,1:nv4)+dz3*cosi(2,1:nv4)
          forall (i=1:mv) aa(2,i) = dot_product(spv(i,1:4),cc(nrv(i)+1:nrv(i)+4))

          if (nv8/=0 .and. p>zero) then
             forall(i=1:nv8) bb(2,i) = pinv*dot_product(cc(i:i+4),bv(i,1:5))
             mvv = mvv+nv8
          endif
      endif

      !  we first determine the matrices (auu) and (qq). then we reduce the
      !  matrix (auu) to upper triangular form (ru) using givens rotations.
      !  we apply the same transformations to the rows of matrix qq to obtain
      !  the (mv+nv8) x nuu matrix g.
      !  we store matrix (ru) into au and g into q.
      l = mvv*nuu
      !  initialization.
      sq = zero
      q(1:l) = zero
      au(1:nuu,1:5) = zero
      l     = 0
      nrold = 0
      n1    = nrold+1
      auu_iterations: do it=1,mu
         number = nru(it)

         ! find the appropriate column of q.
         inner: do
            right(1:mvv) = zero
            if (nrold/=number) then
               if (p<=zero) then
                  nrold = n1
                  n1 = n1+1
                  cycle inner
               endif

               ! fetch a new row of matrix (bu).
               h(1:5) = bu(n1,1:5)*pinv
               i0 = 1
               i1 = 5
            else
               ! fetch a new row of matrix (spu).
               h(1:4) = spu(it,1:4)

               ! find the appropriate column of q.
               right(1:mv) = z(l+1:l+mv)
               l = l+mv

               i0 = 1
               i1 = 4
            endif
            if (nu7-number == iop1) i1 = i1-1
            j0 = n1

            ! take into account that we eliminate the constraints (3)
            do while (j0-1<=iop0)
               fac0 = h(i0)
               right(1:mv) = right(1:mv)-fac0*aa(j0,1:mv)
               if (mv/=mvv) forall(jj=1:nv8) right(mv+jj) = right(mv+jj)-fac0*bb(j0,jj)
               j0 = j0+1
               i0 = i0+1
            end do

            irot = max(0,nrold-iop0-1)
            if(irot<0) irot = 0

            ! rotate the new row of matrix (auu) into triangle.
            rotate_auu: do i=i0,i1
               irot = irot+1
               piv = h(i)
               if (piv/=zero) then
                  ! calculate the parameters of the givens transformation.
                  call fpgivs(piv,au(irot,1),co,si)
                  ! apply that transformation to the rows of matrix (qq).
                  iq = (irot-1)*mvv
                  call fprota(co,si,right(1:mvv),q(iq+1:iq+mvv))

                  ! apply that transformation to the columns of (auu).
                  if (i==i1) cycle rotate_auu
                  i2 = 1
                  i3 = i+1
                  do j=i3,i1
                     i2 = i2+1
                     call fprota(co,si,h(j),au(irot,i2))
                  end do
               endif
            end do rotate_auu

            ! we update the sum of squared residuals
            sq = sq+sum(right(:mvv)**2)

            if (nrold==number) cycle auu_iterations
            nrold = n1
            n1 = n1+1
        end do inner
      end do auu_iterations

      !  we determine the matrix (avv) and then we reduce it to upper triangular form (rv) using
      !  givens rotations. we apply the same transformations to the columns of matrix g to obtain
      !  the (nv-7) x (nu-5-iop0-iop1) matrix h. we store matrix (rv) into av1 and av2, h into c.
      !  the nv7 x nv7 upper triangular matrix (rv) has the form
      !              | av1 '     |
      !       (rv) = |     ' av2 |
      !              |  0  '     |
      !  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 upper triangular matrix of bandwidth 5.
      ncof = nuu*nv7

      !  initialization.
      c  (1:ncof) = zero
      av1(1:nv4,1:5) = zero
      av2(1:nv4,1:4) = zero
      jper = 0
      nrold = 0
      avv_iterations: do it=1,mv
         number = nrv(it)

         avv_inner: do

         if (nrold/=number) then

            if (p<=zero) then
               nrold = nrold+1
               cycle avv_inner
            end if

            ! fetch a new row of matrix (bv).
            n1 = nrold+1
            h(1:5) = bv(n1,1:5)*pinv

            ! find the appropriate row of g.
            if (mv/=mvv) then
               l = mv+n1
               right(1:nuu) = q(l:l+(nuu-1)*mvv:mvv)
            else
               right(1:nuu) = zero
            endif

         else

            ! fetch a new row of matrix (spv)
            h(1:5) = [spv(it,1:4),zero]

            ! find the appropriate row of g.
            right(1:nuu) = q(it:it+(nuu-1)*mvv:mvv)

         endif

         ! test whether there are non-zero values in the new row of (avv)
         ! corresponding to the b-splines n(j,v),j=nv7+1,...,nv4.
         if (nrold>=nv11) then

             if (jper==0) then
                ! initialize the matrix (av2).
                jk = nv11+1
                do i=1,4
                   ik = jk
                   do j=1,5
                      if (ik<=0) exit
                      av2(ik,i) = av1(ik,j)
                      ik = ik-1
                   end do
                   jk = jk+1
                end do
                jper = 1
             endif

             ! if one of the non-zero elements of the new row corresponds to one of the b-splines n(j;v),
             ! j=nv7+1,...,nv4, we take account of condition (2) for setting up this row of (avv). the row
             ! is stored in h1  the part with respect to av1) and h2 (the part with respect to av2).
             h1 = zero
             h2 = zero
             do i=1,5
                j  = nrold-nv11+i
                l0 = j
                l1 = l0-4
                do while (l1>nv11)
                   l0 = l1-nv11
                   l1 = l0-4
                end do
                if (l1<=0) then
                   h2(l0) = h2(l0) + h(i)
                else ! (l1<=nv11)
                   h1(l1) = h(i)
                endif
             end do

             ! rotate the new row of (avv) into triangle.
             if (nv11>0) then
                ! rotations with the rows 1,2,...,nv11 of (avv).
                avv_rot: do j=1,nv11
                   piv = h1(1)
                   i2 = min0(nv11-j,4)

                   if (piv/=zero) then

                      ! calculate the parameters of the givens transformation.
                      call fpgivs(piv,av1(j,1),co,si)

                      ! apply that transformation to the columns of matrix g.
                      ic = j
                      do i=1,nuu
                         call fprota(co,si,right(i),c(ic))
                         ic = ic+nv7
                      end do

                      ! apply that transformation to the rows of (avv) with respect to av2.
                      call fprota(co,si,h2(1:4),av2(j,1:4))

                      ! apply that transformation to the rows of (avv) with respect to av1.
                      if(i2==0) exit avv_rot

                      call fprota(co,si,h1(2:i2+1),av1(j,2:i2+1))
                   endif

                   h1(1:i2+1) = [h1(2:i2+1),zero]

                end do avv_rot
             endif

             ! rotations with the rows nv11+1,...,nv7 of avv.
             avv_rot2: do j=1,4
                ij  = nv11+j
                piv = h2(j)
                if (ij>0 .and. piv/=zero) then

                   ! calculate the parameters of the givens transformation.
                   call fpgivs(piv,av2(ij,j),co,si)

                   ! apply that transformation to the columns of matrix g.
                   ic = ij
                   do i=1,nuu
                      call fprota(co,si,right(i),c(ic))
                      ic = ic+nv7
                   end do

                   ! apply that transformation to the rows of (avv) with respect to av2.
                   if (j<4) call fprota(co,si,h2(j+1:4),av2(ij,j+1:4))

                endif
             end do avv_rot2

         else

             ! rotation into triangle of the new row of (avv), in case the elements
             ! corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
             irot = nrold
             rot_avv_3: do i=1,5
                irot = irot+1
                piv = h(i)
                if (piv/=zero) then
                   ! calculate the parameters of the givens transformation.
                   call fpgivs(piv,av1(irot,1),co,si)

                   ! apply that transformation to the columns of matrix g.
                   ic = irot
                   do j=1,nuu
                      call fprota(co,si,right(j),c(ic))
                      ic = ic+nv7
                   end do

                   ! apply that transformation to the rows of (avv).
                   if (i<5) then
                      i2 = 1
                      i3 = i+1
                      do j=i3,5
                         i2 = i2+1
                         call fprota(co,si,h(j),av1(irot,i2))
                      end do
                   endif
                endif
             end do rot_avv_3

         endif

         ! we update the sum of squared residuals
         sq = sq+sum(right(:nuu)**2)

         if (nrold==number) exit avv_inner
         nrold = nrold+1
         end do avv_inner
      end do avv_iterations

      !  test whether the b-spline coefficients must be determined.
      spline_coefs: if (iback==0) then

          !  backward substitution to obtain the b-spline coefficients as the solution of the linear
          !  system    (rv) (cr) (ru)' = h.

          !  first step: solve the system  (rv) (c1) = h.
          k = 1
          do i=1,nuu
             c(k:k+nv7-1) = fpbacp(av1,av2,c(k),nv7,4,5,nv)
             k = k+nv7
          end do

          !  second step: solve the system  (cr) (ru)' = (c1).
          k = 0
          do j=1,nv7
            k = k+1
            right(:nuu) = c(k:k+(nuu-1)*nv7:nv7)
            right(:nuu) = fpback(au,right,nuu,5,nu)
            c(k:k+(nuu-1)*nv7:nv7) = right(:nuu)
          end do

          !  calculate from the conditions (2)-(3)-(4), the remaining b-spline coefficients.
          ncof = nu4*nv4

          j = 0

          i = nv4 ! i = "last" element pointer
          q (1:nv4) = dz1
          if (iop0/=0) then
             i=2*nv4
             q(nv4+1:2*nv4) = cc(1:nv4)
          end if


          if (nuu/=0) then
             do l=1,nuu
                ii = i
                do k=1,nv7
                   i = i+1
                   j = j+1
                   q(i) = c(j)
                end do
                do k=1,3
                   ii = ii+1
                   i  = i+1
                   q(i) = q(ii)
                end do
             end do
          endif

          if (iop1/=0) q(i+1:i+nv4) = zero

          c(1:ncof) = q(1:ncof)

          !  calculate the quantities
          !    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu; j=1,2,..,mv
          !    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
          !    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7, tu(r+3) <= u(i) <= tu(r+4)
          !    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7, tv(r+3) <= v(j) <= tv(r+4)
          fp  = zero
          fpu = zero
          fpv = zero
          iz = 0
          nroldu = 0
          !  main loop for the different grid points.
          do i1=1,mu
             numu = nru(i1)
             numu1 = numu+1
             nroldv = 0
             do i2=1,mv
                numv = nrv(i2)
                numv1 = numv+1
                iz = iz+1
                !  evaluate s(u,v) at the current grid point by making the sum of the
                !  cross products of the non-zero b-splines at (u,v), multiplied with
                !  the appropriate b-spline coefficients.
                term = zero
                k1 = numu*nv4+numv
                do l1=1,4
                   fac  = spu(i1,l1)
                   term = term + fac*dot_product(spv(i2,1:4),c(k1+1:k1+4))
                   k1 = k1+nv4
                end do
                !  calculate the squared residual at the current grid point.
                term = (z(iz)-term)**2
                !  adjust the different parameters.
                fp = fp+term
                fpu(numu1) = fpu(numu1)+term
                fpv(numv1) = fpv(numv1)+term
                fac = term*half

                if(numv/=nroldv) then
                   fpv(numv1) = fpv(numv1)-fac
                   fpv(numv)  = fpv(numv)+fac
                endif
                nroldv = numv
                if (numu/=nroldu) then
                   fpu(numu1) = fpu(numu1)-fac
                   fpu(numu)  = fpu(numu)+fac
                endif
              end do
              nroldu = numu
          end do
      endif spline_coefs

      end subroutine fpgrdi


      !  find the least-squares spline sinf(u,v) and calculate for each knot interval tu(j+3)<=u<=tu(j+4)
      !  (tv(j+3)<=v<=tv(j+4)) the sum of squared residuals fpintu(j),j=1,2,...,nu-7 (fpintv(j),j=1,2,...
      !  ,nv-7) for the data points having their absciss (ordinate)-value belonging to that interval.
      !  fp gives the total sum of squared residuals.
      subroutine fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,tu,nu,tv,nv,p,c,nc,fp,fpu,fpv, &
                        mm,mvnu,spu,spv,right,q,au,au1,av,av1,bu,bv,nru,nrv)

      !  ..
      !  ..scalar arguments..
      integer, intent(inout) :: ifsu ! (spu) needs to be determined if /=0 [todo: make logical]
      integer, intent(inout) :: ifsv ! (spv) needs to be determined
      integer, intent(inout) :: ifbu ! (bu)  needs to be determined
      integer, intent(inout) :: ifbv ! (bv)  needs to be determined
      integer, intent(in) :: idim,mu,mv,mz,nu,nv,nc,mm,mvnu ! sizes

      !  ..array arguments..
      real(RKIND), intent(inout) :: c(nc*idim),right(mm*idim),q(mvnu),au(nu,5),av(nv,5)
      real(RKIND), intent(in) :: u(mu),v(mv),z(mz*idim),tu(nu),tv(nv),au1(nu,4),av1(nv,4)
      integer,     intent(in) :: ipar(2)

      real(RKIND), intent(in)  :: p
      real(RKIND), intent(out) :: fp,fpu(nu),fpv(nv)

      ! Only modified if ifsu, ifsv
      real(RKIND), intent(inout) :: spu(mu,4),spv(mv,4)
      integer,     intent(inout) :: nru(mu),nrv(mv)

      ! Only modified if ifbu,ifbv
      real(RKIND), intent(inout) :: bu(nu,5),bv(nv,5)

      !  ..local scalars..
      real(RKIND) :: arg,fac,term,value
      integer :: i,id,ii,it,iz,i1,i2,j,jz,k,k1,k2,l,l1,l2,mvv,k0,muu,ncof,nroldu,nroldv,number,nmd,numu, &
                 numu1,numv,numv1,nuu,nvv,nu4,nu7,nu8,nv4,nv7,nv8,n33
      !  ..local arrays..
      real(RKIND) :: h(SIZ_K+1)
      !  ..subroutine references..
      !    fpback,fpbspl,fpdisc,fpbacp,fptrnp,fptrpe
      !  ..
      !  let
      !               |   (spu)    |            |   (spv)    |
      !        (au) = | ---------- |     (av) = | ---------- |
      !               | (1/p) (bu) |            | (1/p) (bv) |
      !
      !                                | z  ' 0 |
      !                            q = | ------ |
      !                                | 0  ' 0 |
      !
      !  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline coefficients.
      !       z      : the mu x mv matrix which contains the function values.
      !       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices according to the
      !                least-squares problems in the u-,resp. v-direction.
      !       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices containing the discontinuity
      !                jumps of the derivatives of the b-splines in the u-,resp.v-variable at the knots
      !
      !  the b-spline coefficients of the smoothing spline are then calculated as the least-squares
      !  solution of the following over-determined linear system of equations
      !
      !    (1)  (av) c (au)' = q
      !
      !  subject to the constraints
      !
      !    (2)  c(nu-3+i,j) = c(i,j), i=1,2,3 ; j=1,2,...,nv-4
      !            if(ipar(1)/=0)
      !
      !    (3)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
      !            if(ipar(2)/=0)
      !

      !  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      muu = mu - merge(1,0,ipar(1)/=0)
      mvv = mv - merge(1,0,ipar(2)/=0)

      !  it depends on the value of the flags ifsu,ifsv,ifbu  and ibvand
      !  on the value of p whether the matrices (spu), (spv), (bu) and (bv)
      !  still must be determined.
      compute_spu: if (ifsu==0) then
          !  calculate the non-zero elements of the matrix (spu) [mu x nu-4] which is the observation
          !  matrix according to the least-squares spline approximation problem in the u-direction.
          l  = 4
          l1 = 5
          number = 0
          spu_rows: do it=1,muu
            arg = u(it)
            do while (arg>=tu(l1) .and. l/=nu4)
               l  = l1
               l1 = l+1
               number = number+1
            end do
            call fpbspl(tu,nu,3,arg,l,h)
            spu(it,1:4) = h(1:4)
            nru(it) = number
          end do spu_rows

          ! Set (spu) known
          ifsu = 1
      endif compute_spu

      !  calculate the non-zero elements of the matrix (spv) [mv x nv-4] which is the observation
      !  matrix according to the least-squares spline approximation problem in the v-direction.
      compute_spv: if (ifsv==0) then
          l  = 4
          l1 = 5
          number = 0
          spv_rows: do it=1,mvv
            arg = v(it)
            do while (arg>=tv(l1) .and. l/=nv4)
               l = l1
               l1 = l+1
               number = number+1
            end do
            call fpbspl(tv,nv,3,arg,l,h)
            spv(it,1:4) = h(1:4)
            nrv(it) = number
          end do spv_rows

          ! Set (spv) known
          ifsv = 1
      endif compute_spv

      if (p>zero) then

         ! calculate the non-zero elements of the matrix (bu).
         if (ifbu==0 .and. nu8>0) then
            call fpdisc(tu,nu,5,bu,nu)
            ifbu = 1
         endif

         ! calculate the non-zero elements of the matrix (bv).
         if (ifbv==0 .and. nv8/=0) then
            call fpdisc(tv,nv,5,bv,nv)
            ifbv = 1
         endif

      endif

      !  substituting (2)  and (3) into (1), we obtain the overdetermined system
      !         (4)  (avv) (cr) (auu)' = (qq)
      !  from which the nuu*nvv remaining coefficients
      !         c(i,j) , i=1,...,nu-4-3*ipar(1) ; j=1,...,nv-4-3*ipar(2) ,
      !  the elements of (cr), are then determined in the least-squares sense.
      !  we first determine the matrices (auu) and (qq). then we reduce the matrix (auu) to upper
      !  triangular form (ru) using givens rotations.
      !  we apply the same transformations to the rows of matrix qq to obtain the (mv) x nuu matrix g.
      !  we store matrix (ru) into au (and au1 if ipar(1)=1) and g into q.
      if (ipar(1)==0) then
         nuu = nu4
         call fptrnp(mu,mv,idim,nu,nru,spu,p,bu,z,au,q,right)
      else
         nuu = nu7
         call fptrpe(mu,mv,idim,nu,nru,spu,p,bu,z,au,au1,q,right)
      endif

      !  we determine the matrix (avv) and then we reduce this matrix to upper triangular form (rv)
      !  using givens rotations. we apply the same transformations to the columns of matrix g to obtain
      !  the (nvv) x (nuu) matrix h. we store matrix (rv) into av (and av1 if ipar(2)=1) and h into c.
      if (ipar(2)==0) then
         nvv = nv4
         call fptrnp(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,c,right)
      else
         nvv = nv7
         call fptrpe(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,av1,c,right)
      endif

      !  backward substitution to obtain the b-spline coefficients as the solution of the linear system
      !   (rv) (cr) (ru)' = h.
      !  first step: solve the system  (rv) (c1) = h.
      ncof = nuu*nvv
      k = 1
      do ii=1,idim
          do i=1,nuu
             if (ipar(2)/=0) then
                c(k:k+nvv-1) = fpbacp(av,av1,c(k),nvv,4,5,nv)
             else
                c(k:k+nvv-1) = fpback(av,c(k),nvv,5,nv)
             end if
             k = k+nvv
          end do
      end do

      !  second step: solve the system  (cr) (ru)' = (c1).
      do ii=1,idim
          k = (ii-1)*ncof
          do j=1,nvv
            k = k+1
            l = k
            do i=1,nuu
                right(i) = c(l)
                l = l+nvv
            end do
            if (ipar(1)/=0) then
                right(:nuu) = fpbacp(au,au1,right,nuu,4,5,nu)
            else
                right(:nuu) = fpback(au,right,nuu,5,nu)
            end if
            l = k
            do i=1,nuu
                c(l) = right(i)
                l = l+nvv
            end do
          end do
      end do

      !  calculate from the conditions (2)-(3), the remaining b-spline coefficients.
      if (ipar(2)/=0) then
          i = 0
          j = 0
          do id=1,idim
              do l=1,nuu
                 ii = i
                 do k=1,nvv
                    i = i+1
                    j = j+1
                    q(i) = c(j)
                 end do
                 do k=1,3
                    ii = ii+1
                    i = i+1
                    q(i) = q(ii)
                 end do
              end do
          end do
          ncof = nv4*nuu
          nmd  = ncof*idim
          c(1:nmd) = q(1:nmd)
      endif

      if (ipar(1)/=0) then
          i = 0
          j = 0
          n33 = 3*nv4
          do id=1,idim
             ii = i
             do k=1,ncof
                i = i+1
                j = j+1
                q(i) = c(j)
             end do
             do k=1,n33
                ii = ii+1
                i = i+1
                q(i) = q(ii)
             end do
          end do
          ncof = nv4*nu4
          nmd  = ncof*idim
          c(1:nmd) = q(1:nmd)
      endif

      !  calculate the quantities
      !    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
      !    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
      !    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
      !                  tu(r+3) <= u(i) <= tu(r+4)
      !    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
      !                  tv(r+3) <= v(j) <= tv(r+4)
      fp = zero
      fpu(1:nu) = zero
      fpv(1:nv) = zero
      nroldu = 0

      !  main loop for the different grid points.
      u_points: do i1=1,muu
         numu   = nru(i1)
         numu1  = numu+1
         nroldv = 0
         iz     = (i1-1)*mv

         v_points: do i2=1,mvv
            numv  = nrv(i2)
            numv1 = numv+1
            iz    = iz+1

            ! evaluate s(u,v) at the current grid point by making the sum of the cross products of the
            ! non-zero b-splines at (u,v), multiplied with the appropriate b-spline coefficients.
            term = zero
            k0 = numu*nv4+numv
            jz = iz
            do id=1,idim
               k1 = k0
               value = zero
               do l1=1,4
                  k2 = k1
                  fac = spu(i1,l1)
                  do l2=1,4
                     k2 = k2+1
                     value = value+fac*spv(i2,l2)*c(k2)
                  end do
                  k1 = k1+nv4
               end do
               ! calculate the squared residual at the current grid point.
               term = term+(z(jz)-value)**2
               jz = jz+mz
               k0 = k0+ncof
            end do

            ! adjust the different parameters.
            fp = fp+term
            fpu(numu1) = fpu(numu1)+term
            fpv(numv1) = fpv(numv1)+term
            fac = term*half

            if (numv/=nroldv) then
               fpv(numv1) = fpv(numv1)-fac
               fpv(numv)  = fpv(numv) +fac
            endif
            nroldv = numv
            if (numu/=nroldu) then
               fpu(numu1) = fpu(numu1)-fac
               fpu(numu)  = fpu(numu) +fac
            endif
         end do v_points
         nroldu = numu
      end do u_points
      return
      end subroutine fpgrpa


      !  ..
      !  the b-spline coefficients of the smoothing spline are calculated as the least-squares
      !  solution of the over-determined linear system of equations  (ay) c (ax)' = q  where
      !
      !               |   (spx)    |            |   (spy)    |
      !        (ax) = | ---------- |     (ay) = | ---------- |
      !               | (1/p) (bx) |            | (1/p) (by) |
      !
      !                                | z  ' 0 |
      !                            q = | ------ |
      !                                | 0  ' 0 |
      !
      !  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the b-spline coefficients.
      !       z      : the my x mx matrix which contains the function values.
      !       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation matrices according to the
      !                least-squares problems in the x- and y-direction.
      !       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1) matrices which contain
      !                the discontinuity jumps of the derivatives of the b-splines in the x- and
      !                y-direction.
      recursive subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz, &
                                  kx,ky,tx,nx,ty,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2, &
                                  spx,spy,right,q,ax,ay,bx,by,nrx,nry)

      !  ..
      !  ..scalar arguments..
      real(RKIND) ::p,fp
      integer :: ifsx,ifsy,ifbx,ifby,mx,my,mz,kx,ky,nx,ny,nc,mm,mynx,kx1,kx2,ky1,ky2
      !  ..array arguments..
      real(RKIND) :: x(mx),y(my),z(mz),tx(nx),ty(ny),c(nc),spx(mx,kx1),spy(my,ky1),&
                     right(mm),q(mynx),ax(nx,kx2),bx(nx,kx2),ay(ny,ky2),by(ny,ky2),fpx(nx),fpy(ny)
      integer :: nrx(mx),nry(my)
      !  ..local scalars..
      real(RKIND) :: arg,cos,fac,pinv,piv,sin,term
      integer :: i,ibandx,ibandy,ic,iq,irot,it,iz,i1,i2,i3,j,k,k1,k2,l,l1,l2,ncof,nk1x,nk1y,&
                 nrold,nroldx,nroldy,number,numx,numx1,numy,numy1,n1
      !  ..local arrays..
      real(RKIND) :: h(SIZ_K+1)

      nk1x = nx-kx1
      nk1y = ny-ky1
      pinv = merge(one/p,one,p>zero)

      !  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on the value of p whether
      !  the matrices (spx),(spy),(bx) and (by) still must be determined.
      if (ifsx==0) then

          !  calculate the non-zero elements of the matrix (spx) which is the observation matrix
          !  according to the least-squares spline approximation problem in the x-direction.
          l  = kx1
          l1 = kx2
          number = 0
          get_nrx: do it=1,mx
            arg = x(it)
            do while (arg>=tx(l1) .and. l/=nk1x)
               l  = l1
               l1 = l+1
               number = number+1
            end do
            call fpbspl(tx,nx,kx,arg,l,h)
            spx(it,1:kx1) = h(1:kx1)
            nrx(it) = number
          end do get_nrx

          ifsx = 1
      endif

      if (ifsy==0) then

          ! calculate the non-zero elements of the matrix (spy) which is the observation matrix
          ! according to the least-squares spline approximation problem in the y-direction.
          l  = ky1
          l1 = ky2
          number = 0
          get_nry: do it=1,my
             arg = y(it)
             do while (arg>=ty(l1) .and. l/=nk1y)
                l = l1
                l1 = l+1
                number = number+1
             end do
            call fpbspl(ty,ny,ky,arg,l,h)
            spy(it,1:ky1) = h(1:ky1)
            nry(it) = number
          end do get_nry

          ifsy = 1
      endif

      if (p>zero) then
          !  calculate the non-zero elements of the matrix (bx).
          if (ifbx==0 .and. nx/=2*kx1) then
             call fpdisc(tx,nx,kx2,bx,nx)
             ifbx = 1
          endif
          !  calculate the non-zero el ements of the matrix (by).
          if (ifby==0 .and. ny/=2*ky1) then
             call fpdisc(ty,ny,ky2,by,ny)
             ifby = 1
          endif
      endif

      !  reduce the matrix (ax) to upper triangular form (rx) using givens rotations. apply the
      !  same transformations to the rows of matrix q to obtain the my x (nx-kx-1) matrix g.
      !  store matrix (rx) into (ax) and g into q.
      l = my*nk1x
      !  initialization.
      q(1:l) = zero
      ax(1:nk1x,kx2) = zero
      l = 0
      nrold = 0
      !  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = kx1
      givens_ax: do it=1,mx
         number = nrx(it)
         inner_ax: do
           if(nrold==number) then
              ! fetch a new row of matrix (spx).
              h(ibandx) = zero
              h(1:kx1) = spx(it,1:kx1)
              ! find the appropriate column of q.
              do j=1,my
                 l = l+1
                 right(j) = z(l)
              end do
              irot = number
           elseif (p<=zero) then
              nrold = nrold+1
              cycle inner_ax
           else
              ibandx = kx2
              ! fetch a new row of matrix (bx).
              n1 = nrold+1
              h(1:kx2) = bx(n1,1:kx2)*pinv
              ! find the appropriate column of q.
              right(1:my) = zero
              irot = nrold
           endif

           ! rotate the new row of matrix (ax) into triangle.
           rot_new_row: do i=1,ibandx
              irot = irot+1
              piv = h(i)
              if (piv==zero) cycle rot_new_row

              ! calculate the parameters of the givens transformation.
              call fpgivs(piv,ax(irot,1),cos,sin)
              ! apply that transformation to the rows of matrix q.
              iq = (irot-1)*my
              call fprota(cos,sin,right(1:my),q(iq+1:iq+my))

              ! apply that transformation to the columns of (ax).
              if (i<ibandx) then
                 i2 = 1
                 i3 = i+1
                 do j=i3,ibandx
                    i2 = i2+1
                    call fprota(cos,sin,h(j),ax(irot,i2))
                 end do
              endif
           end do rot_new_row

           if (nrold==number) exit inner_ax

           nrold = nrold+1
         end do inner_ax
      end do givens_ax

      !  reduce the matrix (ay) to upper triangular form (ry) using givens rotations. apply the same
      !  transformations to the columns of matrix g to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
      !  store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y

      !  initialization.
      c(1:ncof) = zero
      ay(1:nk1y,1:ky2) = zero
      nrold = 0

      !  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = ky1
      givens_ay: do it=1,my
         number = nry(it)
         inner_ay: do
            if (nrold==number) then
                ! fetch a new row of matrix (spy)
                h(ibandy) = zero
                h(1:ky1) = spy(it,1:ky1)
                ! find the appropriate row of g.
                l = it
                do j=1,nk1x
                   right(j) = q(l)
                   l = l+my
                end do
                irot = number
            elseif (p<=zero) then
                nrold = nrold+1
                cycle inner_ay
            else
                ibandy = ky2
                ! fetch a new row of matrix (by).
                n1 = nrold+1
                h(1:ky2) = by(n1,1:ky2)*pinv
                ! find the appropriate row of g.
                right(1:nk1x) = zero
                irot = nrold
            endif

            ! rotate the new row of matrix (ay) into triangle.
            rot_new_rowy: do i=1,ibandy
               irot = irot+1
               piv = h(i)
               if (piv==zero) cycle rot_new_rowy

               ! calculate the parameters of the givens transformation.
               call fpgivs(piv,ay(irot,1),cos,sin)
               ! apply that transformation to the columns of matrix g.
               ic = irot
               do j=1,nk1x
                  call fprota(cos,sin,right(j),c(ic))
                  ic = ic+nk1y
               end do
               ! apply that transformation to the columns of matrix (ay).
               if (i<ibandy) then
                   i2 = 1
                   i3 = i+1
                   do j=i3,ibandy
                      i2 = i2+1
                      call fprota(cos,sin,h(j),ay(irot,i2))
                   end do
               endif
            end do rot_new_rowy
            if (nrold==number) exit inner_ay
            nrold = nrold+1
         end do inner_ay
      end do givens_ay

      !  backward substitution to obtain the b-spline coefficients as the
      !  solution of the linear system    (ry) c (rx)' = h.
      !  first step: solve the system  (ry) (c1) = h.
      k = 1
      do i=1,nk1x
        c(k:k+nk1y-1) = fpback(ay,c(k),nk1y,ibandy,ny)
        k = k+nk1y
      end do

      !  second step: solve the system  c (rx)' = (c1).
      k = 0
      do j=1,nk1y
         k = k+1
         l = k
         do i=1,nk1x
            right(i) = c(l)
            l = l+nk1y
         end do
         right(:nk1x) = fpback(ax,right,nk1x,ibandx,nx)
         l = k
         do i=1,nk1x
            c(l) = right(i)
            l = l+nk1y
         end do
      end do

      !  calculate the quantities
      !    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
      !    fp = sumi=1,mx(sumj=1,my(res(i,j)))
      !    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
      !                  tx(r+kx) <= x(i) <= tx(r+kx+1)
      !    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
      !                  ty(r+ky) <= y(j) <= ty(r+ky+1)
      fp     = zero
      fpx    = zero
      fpy    = zero
      nk1y   = ny-ky1
      iz     = 0
      nroldx = 0

      !  main loop for the different grid points.
      grid_x: do i1=1,mx
         numx = nrx(i1)
         numx1 = numx+1
         nroldy = 0
         grid_y: do i2=1,my
            numy = nry(i2)
            numy1 = numy+1
            iz = iz+1
            ! evaluate s(x,y) at the current grid point by making the sum of the
            ! cross products of the non-zero b-splines at (x,y), multiplied with
            ! the appropriate b-spline coefficients.
            term = zero
            k1 = numx*nk1y+numy
            do l1=1,kx1
               term = term+spx(i1,l1)*dot_product(spy(i2,1:ky1),c(k1+1:k1+ky1))
               k1 = k1+nk1y
            end do

            ! calculate the squared residual at the current grid point.
            term = (z(iz)-term)**2
            ! adjust the different parameters.
            fp = fp+term
            fpx(numx1) = fpx(numx1)+term
            fpy(numy1) = fpy(numy1)+term
            fac = term*half
            if (numy/=nroldy) then
               fpy(numy1) = fpy(numy1)-fac
               fpy(numy)  = fpy(numy) +fac
            endif
            nroldy = numy
            if (numx/=nroldx) then
               fpx(numx1) = fpx(numx1)-fac
               fpx(numx)  = fpx(numx) +fac
            endif
         end do grid_y
         nroldx = numx
      end do grid_x
      return
      end subroutine fpgrre

      recursive subroutine fpgrsp(ifsu,ifsv,ifbu,ifbv,iback,u,mu,v, &
       mv,r,mr,dr,iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm, &
       mvnu,spu,spv,right,q,au,av1,av2,bu,bv,a0,a1,b0,b1,c0,c1, &
       cosi,nru,nrv)

      !  ..
      !  ..scalar arguments..
      real(RKIND) :: p,sq,fp
      integer :: ifsu,ifsv,ifbu,ifbv,iback,mu,mv,mr,iop0,iop1,nu,nv,nc,mm,mvnu
      !  ..array arguments..
      real(RKIND) :: u(mu),v(mv),r(mr),dr(6),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv), &
                     spu(mu,4),spv(mv,4),right(mm),q(mvnu),au(nu,5),av1(nv,6),c0(nv), &
                     av2(nv,4),a0(2,mv),b0(2,nv),cosi(2,nv),bu(nu,5),bv(nv,5),c1(nv), &
                     a1(2,mv),b1(2,nv)
      integer :: nru(mu),nrv(mv)
      !  ..local scalars..
      real(RKIND) :: arg,co,dr01,dr11,fac,fac0,fac1,pinv,piv,si,term
      integer :: i,ic,ii,ij,ik,iq,irot,it,ir,i0,i1,i2,i3,j,jj,jk,jper, &
                 j0,j1,k,k1,l,l0,l1,mvv,ncof,nrold,nroldu,nroldv,number, &
                 numu,numu1,numv,numv1,nuu,nu4,nu7,nu8,nu9,nv11,nv4,nv7,nv8,n1
      !  ..local arrays..
      real(RKIND) :: h(SIZ_K+1),h1(5),h2(4)

      !  let
      !               |     (spu)      |            |     (spv)      |
      !        (au) = | -------------- |     (av) = | -------------- |
      !               | sqrt(1/p) (bu) |            | sqrt(1/p) (bv) |
      !
      !                                | r  ' 0 |
      !                            q = | ------ |
      !                                | 0  ' 0 |
      !
      !  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
      !                coefficients.
      !       r      : the mu x mv matrix which contains the function values.
      !       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
      !                according to the least-squares problems in the u-,resp.
      !                v-direction.
      !       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
      !                containing the discontinuity jumps of the derivatives
      !                of the b-splines in the u-,resp.v-variable at the knots
      !  the b-spline coefficients of the smoothing spline are then calculated
      !  as the least-squares solution of the following over-determined linear
      !  system of equations
      !
      !  (1)  (av) c (au)' = q
      !
      !  subject to the constraints
      !
      !  (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
      !
      !  (3)  if iop0 = 0  c(1,j) = dr(1)
      !          iop0 = 1  c(1,j) = dr(1)
      !                    c(2,j) = dr(1)+(dr(2)*cosi(1,j)+dr(3)*cosi(2,j))*
      !                            tu(5)/3. = c0(j) , j=1,2,...nv-4
      !
      !  (4)  if iop1 = 0  c(nu-4,j) = dr(4)
      !          iop1 = 1  c(nu-4,j) = dr(4)
      !                    c(nu-5,j) = dr(4)+(dr(5)*cosi(1,j)+dr(6)*cosi(2,j))
      !                                *(tu(nu-4)-tu(nu-3))/3. = c1(j)
      !

      !  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nu9 = nu-9
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      nv11 = nv-11
      nuu = nu4-iop0-iop1-2
      pinv = merge(one/p,one,p>zero)

      !  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv,iop0,iop1
      !  and on the value of p whether the matrices (spu), (spv), (bu), (bv),
      !  (cosi) still must be determined.
      if (ifsu==0) then
          !  calculate the non-zero elements of the matrix (spu) which is the observation matrix
          !  according to the least-squares spline approximation problem in the u-direction.
          l  = 4
          l1 = 5
          number = 0
          do it=1,mu
            arg = u(it)
            do while (.not.(arg<tu(l1) .or. l==nu4))
                l = l1
                l1 = l+1
                number = number+1
            end do
            call fpbspl(tu,nu,3,arg,l,h)
            spu(it,1:4) = h(1:4)
            nru(it) = number
          end do
          ifsu = 1
      endif

      !  calculate the non-zero elements of the matrix (spv) which is the ob-
      !  servation matrix according to the least-squares spline approximation
      !  problem in the v-direction.
      if(ifsv==0) then

          l  = 4
          l1 = 5
          number = 0
          do it=1,mv
            arg = v(it)
            do while (.not.(arg<tv(l1) .or. l==nv4))
                l = l1
                l1 = l+1
                number = number+1
            end do
            call fpbspl(tv,nv,3,arg,l,h)
            spv(it,1:4) = h(1:4)
            nrv(it) = number
          end do
          ifsv = 1


          !  calculate the coefficients of the interpolating splines for cos(v) and sin(v).
          if (iop0/=0 .or. iop1/=0) then

              cosi(:,1:nv4) = zero
              if (nv7>=4) then
                  do i=1,nv7
                     l = i+3
                     arg = tv(l)
                     call fpbspl(tv,nv,3,arg,l,h)
                     av1(i,1:3) = h(1:3)
                     cosi(:,i)  = [cos(arg),sin(arg)]
                  end do
                  call fpcyt1(av1,nv7,nv)
                  do j=1,2
                     call fpcyt2(av1,nv7,cosi(j,1:nv7),right,nv)
                     cosi(j,1:nv4) = [right(nv7),right(1:nv7),right(1:2)]
                  end do
              endif
          endif

      endif

      if (p>zero) then
          !  calculate the non-zero elements of the matrix (bu).
          if(ifbu==0 .and. nu8/=0) then
             call fpdisc(tu,nu,5,bu,nu)
             ifbu = 1
          endif

          !  calculate the non-zero elements of the matrix (bv).
          if (ifbv==0 .and. nv8/=0) then
             call fpdisc(tv,nv,5,bv,nv)
             ifbv = 1
          endif
      endif
      !  substituting (2),(3) and (4) into (1), we obtain the overdetermined system
      !         (5)  (avv) (cc) (auu)' = (qq)
      !  from which the nuu*nv7 remaining coefficients
      !         c(i,j) , i=2+iop0,3+iop0,...,nu-5-iop1,j=1,2,...,nv-7.
      !  the elements of (cc), are then determined in the least-squares sense.
      !  simultaneously, we compute the resulting sum of squared residuals sq.
      dr01 = dr(1)
      dr11 = dr(4)
      a0(1,1:mv) = dr01
      a1(1,1:mv) = dr11
      if (nv8/=0 .and. p>zero) then
         b0(1,1:nv8) = zero
         b1(1,1:nv8) = zero
      endif
      mvv = mv
      if (iop0/=0) then
          fac       = (tu(5)-tu(4))/three
          c0(1:nv4) = dr01+fac*matmul(dr(2:3),cosi(:,:nv4))

          forall (i=1:mv) a0(2,i) = dot_product(spv(i,1:4),c0(nrv(i)+1:nrv(i)+4))

          if (nv8/=0 .and. p>zero) then
              forall (i=1:nv8) b0(2,i) = pinv*dot_product(bv(i,1:5),c0(i:i+4))*pinv
              mvv = mv+nv8
          endif
      endif
      if (iop1/=0) then
          fac       = (tu(nu4)-tu(nu4+1))/three
          c1(1:nv4) = dr11 + fac*matmul(dr(5:6),cosi(:,:nv4))

          forall (i=1:mv) a1(2,i) = dot_product(c1(nrv(i)+1:nrv(i)+4),spv(i,1:4))

          if (nv8/=0 .and. p>zero) then
              forall (i=1:nv8) b1(2,i) = pinv*dot_product(bv(i,1:5),c1(i:i+4))
              mvv = mv+nv8
          endif
      endif

      !  we first determine the matrices (auu) and (qq). then we reduce the matrix (auu) to an unit
      !  upper triangular form (ru) using givens rotations without square roots. we apply the same
      !  transformations to the rows of matrix qq to obtain the mv x nuu matrix g.
      !  we store matrix (ru) into au and g into q.
      l = mvv*nuu
      !  initialization.
      sq = zero

      if (l/=0) then
         q(1:l) = zero
         au(1:nuu,1:5) = zero
         l = 0
      endif
      nrold = 0
      n1 = nrold+1
      auu_iterations: do it=1,mu
         number = nru(it)

         auu_inner: do

         ! find the appropriate column of q.
         right(1:mvv) = zero

         if (nrold==number) then
            !  fetch a new row of matrix (spu).
            h(1:4) = spu(it,1:4)

            !  find the appropriate column of q.
            right(1:mv) = r(l+1:l+mv)
            l  = l+mv
            i0 = 1
            i1 = 4
         else

            if (p<=zero) then
               nrold = n1
               n1 = n1+1
               cycle auu_inner
            end if
            !  fetch a new row of matrix (bu).
            h(1:5) = bu(n1,1:5)*pinv
            i0 = 1
            i1 = 5
         end if

         j0 = n1
         j1 = nu7-number

         ! take into account that we eliminate the constraints (3)
         do while (j0-1<=iop0)
             fac0 = h(i0)
             right(1:mv) = right(1:mv)-fac0*a0(j0,1:mv)
             if (mv/=mvv) then
                 j = mv
                 do jj=1,nv8
                    j = j+1
                    right(j) = right(j)-fac0*b0(j0,jj)
                 end do
             endif
             j0 = j0+1
             i0 = i0+1
         end do

         ! take into account that we eliminate the constraints (4)
         do while (j1-1<=iop1)
            fac1 = h(i1)
            right(1:mv) = right(1:mv)-fac1*a1(j1,1:mv)
            if (mv/=mvv) then
                j = mv
                do jj=1,nv8
                   j = j+1
                   right(j) = right(j)-fac1*b1(j1,jj)
                end do
            endif
            j1 = j1+1
            i1 = i1-1
         end do

         irot = max(0,nrold-iop0-1)

         ! rotate the new row of matrix (auu) into triangle.
         if (i0<=i1) then
            auu_rot: do i=i0,i1
               irot = irot+1
               piv = h(i)
               if (piv==zero) cycle auu_rot
               ! calculate the parameters of the givens transformation.
               call fpgivs(piv,au(irot,1),co,si)
               ! apply that transformation to the rows of matrix (qq).
               iq = (irot-1)*mvv

               do j=1,mvv
                  iq = iq+1
                  call fprota(co,si,right(j),q(iq))
               end do
               ! apply that transformation to the columns of (auu).
               if (i<i1) then
                  i2 = 1
                  i3 = i+1
                  do j=i3,i1
                     i2 = i2+1
                     call fprota(co,si,h(j),au(irot,i2))
                  end do
               endif
           end do auu_rot
         endif

         ! we update the sum of squared residuals.
         sq = sq+sum(right(:mvv)**2)
         if (nrold==number) exit auu_inner
         nrold = n1
         n1 = n1+1
         end do auu_inner
      end do auu_iterations

      if (nuu/=0) then
          !  we determine the matrix (avv) and then we reduce her to an unit upper triangular form (rv)
          !  using givens rotations without square roots. we apply the same transformations to the
          !  columns of matrix g to obtain the (nv-7) x (nu-6-iop0-iop1) matrix h. we store matrix (rv)
          !  into av1 and av2, h into c. the nv7 x nv7 triangular unit upper matrix (rv) has the form
          !              | av1 '     |
          !       (rv) = |     ' av2 |
          !              |  0  '     |
          !  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 unit upper triangular matrix of
          !  bandwidth 5.
          ncof = nuu*nv7
          !  initialization.
          c(:ncof) = zero
          av1(1:nv4,1:5) = zero
          av2(1:nv4,1:4) = zero
          jper = 0
          nrold = 0
          avv_iterations: do it=1,mv
             number = nrv(it)

             avv_inner: do

             if (nrold==number) then

               ! fetch a new row of matrix (spv)
               h(1:5) = [spv(it,1:4),zero]

               ! find the appropriate row of g.
               right(1:nuu) = q(it:it+(nuu-1)*mvv:mvv)

            else

               if (p<=zero) then
                  nrold = nrold+1
                  cycle avv_inner
               endif
               ! fetch a new row of matrix (bv).
               n1 = nrold+1
               h(1:5) = bv(n1,1:5)*pinv

               ! find the appropriate row of g.
               if(mv/=mvv) then
                  l = mv+n1
                  right(1:nuu) = q(l:l+(nuu-1)*mvv:mvv)
               else
                  right(1:nuu) = zero
               endif

            endif

            ! test whether there are non-zero values in the new row of (avv)
            ! corresponding to the b-splines n(j;v),j=nv7+1,...,nv4.
            if (nrold<nv11) then

               ! rotation into triangle of the new row of (avv), in case the elements
               ! corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
               do i=1,5
                  irot = nrold+i
                  piv  = h(i)

                  if (piv==zero) cycle

                  ! calculate the parameters of the givens transformation.
                  call fpgivs(piv,av1(irot,1),co,si)

                  ! apply that transformation to the columns of matrix g.
                  ic = irot
                  call fprota(co,si,right(1:nuu),c(ic:ic+(nuu-1)*nv7:nv7))

                  ! apply that transformation to the rows of (avv).
                  if (i<5) then
                     i2 = 1
                     i3 = i+1
                     do j=i3,5
                        i2 = i2+1
                        call fprota(co,si,h(j),av1(irot,i2))
                     end do
                  endif
               end do

            else

               ! initialize the matrix (av2).
               if (jper==0) then

                  jk = nv11+1
                  do i=1,4
                     ik = jk
                     do j=1,5
                         if (ik<=0) exit
                         av2(ik,i) = av1(ik,j)
                         ik = ik-1
                     end do
                     jk = jk+1
                  end do
                  jper = 1

               endif

               ! if one of the non-zero elements of the new row corresponds to one of the b-splines
               ! n(j;v),j=nv7+1,...,nv4, we take account of condition (2) for setting up this row of
               ! (avv). the row is stored in h1(the part with respect to av1) and h2 (the part with
               ! respect to av2).
               h1 = zero
               h2 = zero
               j = nrold-nv11
               do i=1,5
                   j = j+1
                   l0 = j
                   l1 = l0-4
                   do while (l1>nv11)
                      l0 = l1-nv11
                      l1 = l0-4
                   end do
                   if (l1<=0) then
                      h2(l0) = h2(l0) + h(i)
                   else ! l1<=nv11
                      h1(l1) = h(i)
                   end if
               end do

               ! rotate the new row of (avv) into triangle.
               if (nv11>0) then

                  ! rotations with the rows 1,2,...,nv11 of (avv).
                  do j=1,nv11
                     piv = h1(1)
                     i2 = min(nv11-j,4)
                     if (piv/=zero) then

                        ! calculate the parameters of the givens transformation.
                        call fpgivs(piv,av1(j,1),co,si)

                        ! apply that transformation to the columns of matrix g.
                        call fprota(co,si,right(1:nuu),c(j:j+nv7*(nuu-1):nv7))

                        ! apply that transformation to the rows of (avv) with respect to av2.
                        call fprota(co,si,h2(1:4),av2(j,1:4))

                        ! apply that transformation to the rows of (avv) with respect to av1.
                        if (i2/=0) call fprota(co,si,h1(2:i2+1),av1(j,2:i2+1))

                     endif
                     h1(1:i2+1) = [h1(2:i2+1),zero]
                  end do

               endif

               ! rotations with the rows nv11+1,...,nv7 of avv.
               avv_rot: do j=1,4
                  ij  = nv11+j
                  j1  = j+1
                  piv = h2(j)

                  if (ij<=0 .or. piv==zero) cycle avv_rot

                  ! calculate the parameters of the givens transformation.
                  call fpgivs(piv,av2(ij,j),co,si)

                  ! apply that transformation to the columns of matrix g.
                  call fprota(co,si,right(1:nuu),c(ij:ij+nv7*(nuu-1):nv7))

                  ! apply that transformation to the rows of (avv) with respect to av2.
                  if (j<4) call fprota(co,si,h2(j1:4),av2(ij,j1:4))

              end do avv_rot

            end if

            ! we update the sum of squared residuals.
            sq = sq+sum(right(:nuu)**2)

            if (nrold==number) exit avv_inner

            nrold = nrold+1
            end do avv_inner

          end do avv_iterations

          !  test whether the b-spline coefficients must be determined.
          if (iback/=0) return
          !  backward substitution to obtain the b-spline coefficients as the
          !  solution of the linear system    (rv) (cr) (ru)' = h.
          !  first step: solve the system  (rv) (c1) = h.
          k = 1
          do i=1,nuu
             c(k:k+nv7-1) = fpbacp(av1,av2,c(k),nv7,4,5,nv)
             k = k+nv7
          end do

          !  second step: solve the system  (cr) (ru)' = (c1).
          do k=1,nv7
            right(:nuu) = c(k:k+(nuu-1)*nv7)
            right(:nuu) = fpback(au,right,nuu,5,nu)
            c(k:k+(nuu-1)*nv7:nv7) = right(:nuu)
          end do

      endif

      !  calculate from the conditions (2)-(3)-(4), the remaining b-spline
      !  coefficients.
      ncof = nu4*nv4

      q(1:nv4) = dr01
      q(ncof-nv4+1:ncof) = dr11

      i = nv4
      j = 0
      if (iop0/=0) then
         q(i+1:i+nv4) = c0(1:nv4)
         i = i+nv4
      endif
      if(nuu/=0) then
         do l=1,nuu
            ii = i
            do k=1,nv7
               i = i+1
               j = j+1
               q(i) = c(j)
            end do
            do k=1,3
               ii = ii+1
               i  = i+1
               q(i) = q(ii)
            end do
         end do
      endif
      if(iop1/=0) then
          do l=1,nv4
             i = i+1
             q(i) = c1(l)
          end do
      endif
      c(1:ncof) = q(1:ncof)
      !  calculate the quantities
      !    res(i,j) = (r(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
      !    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
      !    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
      !                  tu(r+3) <= u(i) <= tu(r+4)
      !    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
      !                  tv(r+3) <= v(j) <= tv(r+4)
      fp = zero
      fpu = zero
      fpv = zero
      ir = 0
      nroldu = 0

      !  main loop for the different grid points.
      grid_u: do i1=1,mu
         numu = nru(i1)
         numu1 = numu+1
         nroldv = 0
         grid_v: do i2=1,mv
            numv = nrv(i2)
            numv1 = numv+1
            ir = ir+1
            ! evaluate s(u,v) at the current grid point by making the sum of the
            ! cross products of the non-zero b-splines at (u,v), multiplied with
            ! the appropriate b-spline coefficients.
            term = zero
            k1 = numu*nv4+numv
            do l1 = 1,4
               fac  = spu(i1,l1)
               term = term+fac*dot_product(spv(i2,1:4),c(k1+1:k1+4))
               k1   = k1+nv4
            end do
            ! calculate the squared residual at the current grid point.
            term = (r(ir)-term)**2
            ! adjust the different parameters.
            fp = fp+term
            fpu(numu1) = fpu(numu1)+term
            fpv(numv1) = fpv(numv1)+term
            fac = term*half
            if (numv/=nroldv) then
               fpv(numv1) = fpv(numv1)-fac
               fpv(numv) = fpv(numv)+fac
            endif
            nroldv = numv
            if (numu/=nroldu) then
               fpu(numu1) = fpu(numu1)-fac
               fpu(numu) = fpu(numu)+fac
            endif
         end do grid_v
         nroldu = numu
      end do grid_u
      return
      end subroutine fpgrsp


      !  given the b-spline representation (knots t(j),j=1,2,...,n, b-spline coefficients c(j),j=1,2,...,
      !  n-k-1) of a spline of degree k, fpinst calculates the b-spline representation (knots
      !  tt(j),j=1,2,...,nn, b-spline coefficients cc(j),j=1,2,...,nn-k-1) of the same spline if an
      !  additional knot is inserted at the point x situated in the inter val t(l)<=x<t(l+1).
      !  iopt/=0: periodic spline; at leas one of the following conditions must be fulfilled: l>2*k or l<n-2*k.
      !  iopt==0: non-periodic spline
      pure subroutine fpinst(iopt,t,n,c,k,x,l,tt,nn,cc,nest)

      !
      !  ..scalar arguments..
      integer, intent(in) :: k,n,l,iopt,nest
      integer, intent(out) :: nn
      real(RKIND), intent(in) :: x
      !  ..array arguments..
      real(RKIND), intent(in)  :: t(nest),c(nest)
      real(RKIND), intent(out) :: tt(nest),cc(nest)
      !  ..local scalars..
      real(RKIND) :: fac,per
      integer :: i,i1,j,k1,m,mk,nk,nk1,nl,ll
      !  ..
      k1  = k+1
      nk1 = n-k1
      !  the new knots
      ll = l+1
      tt(1:n+1) = [t(1:l),x,t(ll:n)]
      !  the new b-spline coefficients
      i = nk1
      do j=l,nk1
         cc(i+1) = c(i)
         i = i-1
      end do
      i = l
      do j=1,k
         m = i+k1
         fac = (x-tt(i))/(tt(m)-tt(i))
         i1 = i-1
         cc(i) = fac*c(i)+(one-fac)*c(i1)
         i = i1
      end do
      cc(1:i) = c(1:i)
      nn = n+1

      ! incorporate the boundary conditions for a periodic spline.
      if (iopt/=0) then

          nk = nn-k
          nl = nk-k1
          per = tt(nk)-tt(k1)
          i = k1
          j = nk
          if (ll>nl) then
              do m=1,k
                 mk = m+nl
                 cc(m) = cc(mk)
                 i = i-1
                 j = j-1
                 tt(i) = tt(j)-per
              end do
          elseif (ll<=(k1+k)) then
             do m=1,k
                mk = m+nl
                cc(mk) = cc(m)
                i = i+1
                j = j+1
                tt(j) = tt(i)+per
             end do
          endif
      endif

      end subroutine fpinst


      pure subroutine fpintb(t,n,bint,nk1,x,y)

      !  subroutine fpintb calculates integrals of the normalized b-splines
      !  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
      !  it makes use of the formulae of gaffney for the calculation of
      !  indefinite integrals of b-splines.
      !
      !  calling sequence:
      !     call fpintb(t,n,bint,nk1,x,y)
      !
      !  input parameters:
      !    t    : real array,length n, containing the position of the knots.
      !    n    : integer value, giving the number of knots.
      !    nk1  : integer value, giving the number of b-splines of degree k,
      !           defined on the set of knots ,i.e. nk1 = n-k-1.
      !    x,y  : real values, containing the end points of the integration
      !           interval.
      !  output parameter:
      !    bint : array,length nk1, containing the integrals of the b-splines.
      !  ..
      !  ..scalars arguments..
      integer,     intent(in)  :: n,nk1
      real(RKIND), intent(in)  :: x,y
      !  ..array arguments..
      real(RKIND), intent(in)  :: t(n)
      real(RKIND), intent(out) :: bint(nk1)

      !  ..local scalars..
      integer :: i,ia,ib,it,j,j1,k,k1,l,li,lj,lk,l0
      logical :: lmin
      real(RKIND) :: a,ak,arg,b,f
      !  ..local arrays..
      real(RKIND) aint(6),h(SIZ_K+1),h1(6)

      integer, parameter :: nit = 2 ! number of iterations

      !  initialization.
      k1   = n-nk1
      ak   = k1
      k    = k1-1
      bint = zero

      !  the integration limits are arranged in increasing order.
      if (x==y) return

      a = max(x,y)
      b = min(x,y)
      lmin = x>y

      a = max(t(k1),a)
      b = min(t(nk1+1),b)
      if (a>b) return

      !  using the expression of gaffney for the indefinite integral of a b-spline we find that
      !  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
      !    where for t(l) <= x < t(l+1)
      !    res(j,x) = 0, j=1,2,...,l-k-1
      !             = 1, j=l+1,l+2,...,nk1
      !             = aint(j+k-l+1), j=l-k,l-k+1,...,l
      !               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
      !                 i=0,1,...,k
      l  = k1
      l0 = l+1
      !  set arg = a.
      arg = a
      ia  = 0
      iterations: do it=1,nit
      !  search for the knot interval t(l) <= arg < t(l+1).
        do while (.not.(arg<t(l0) .or. l==nk1))
           l  = l0
           l0 = l+1
        end do
      !  calculation of aint(j), j=1,2,...,k+1.
      !  initialization.
        aint(1)  = (arg-t(l))/(t(l+1)-t(l))
        aint(2:) = zero
        h1(1) = one
        do j=1,k
      !  evaluation of the non-zero b-splines of degree j at arg,i.e.
      !    h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
          h(1) = zero
          do i=1,j
            li = l+i
            lj = li-j
            f = h1(i)/(t(li)-t(lj))
            h(i) = h(i)+f*(t(li)-arg)
            h(i+1) = f*(arg-t(lj))
          end do
      !  updating of the integrals aint.
          j1 = j+1
          do i=1,j1
            li = l+i
            lj = li-j1
            aint(i) = aint(i)+h(i)*(arg-t(lj))/(t(li)-t(lj))
            h1(i) = h(i)
          end do
        end do

        if (it<nit) then
            ! updating of the integrals bint
            lk = l-k
            ia = lk
            do i=1,k1
              bint(lk) = -aint(i)
              lk = lk+1
            end do

            arg = b
        endif
      end do iterations
      !  updating of the integrals bint.
      lk = l-k
      ib = lk-1
      do i=1,k1
         bint(lk) = bint(lk)+aint(i)
         lk = lk+1
      end do
      if (ib>=ia) bint(ia:ib) = bint(ia:ib)+one

      !  the scaling factors are taken into account.
      f = one/ak
      do i=1,nk1
        j = i+k1
        bint(i) = bint(i)*(t(j)-t(i))*f
      end do
      !  the order of the integration limits is taken into account.
      if (lmin) bint = -bint

      end subroutine fpintb


      !  subroutine fpknot locates an additional knot for a spline of degree k and adjusts the
      !  corresponding parameters,i.e.
      !    t     : the position of the knots.
      !    n     : the number of knots.
      !    nrint : the number of knotintervals.
      !    fpint : the sum of squares of residual right hand sides
      !            for each knot interval.
      !    nrdata: the number of data points inside each knot interval.
      !  istart indicates that the smallest data point at which the new knot may be added is x(istart+1)
      pure subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)

      !  ..
      !  ..scalar arguments..
      integer, intent(in)    :: m,nest,istart
      integer, intent(inout) :: n,nrint
      !  ..array arguments..
      real(RKIND), intent(in)    :: x(m)
      real(RKIND), intent(inout) :: t(nest)
      real(RKIND), intent(inout) :: fpint(nest)
      integer    , intent(inout) :: nrdata(nest)

      !  ..local scalars..
      real(RKIND) :: an,am,fpmax
      integer :: ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,next,nrx,number
      !  ..
      number = 0
      maxpt  = 0
      maxbeg = 0
      k      = (n-nrint-1)/2
      !  search for knot interval t(number+k) <= x <= t(number+k+1) where fpint(number) is maximal on the
      !  condition that nrdata(number)/=0 .
      fpmax = zero
      jbegin = istart
      do j=1,nrint
        jpoint = nrdata(j)

        if (fpmax<fpint(j) .and. jpoint/=0) then
           fpmax = fpint(j)
           number = j
           maxpt = jpoint
           maxbeg = jbegin
        endif

        jbegin = jbegin+jpoint+1
      end do
      !  let coincide the new knot t(number+k+1) with a data point x(nrx)
      !  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx   = maxbeg+ihalf
      next  = number+1

      !  adjust the different parameters.
      if(next<=nrint) then
         do j=next,nrint
            jj = next+nrint-j
            fpint(jj+1) = fpint(jj)
            nrdata(jj+1) = nrdata(jj)
            jk = jj+k
            t(jk+1) = t(jk)
         end do
      endif

      nrdata(number) = ihalf-1
      nrdata(next)   = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n     = n+1
      nrint = nrint+1
      return
      end subroutine fpknot

      ! given the set of function values z(i,j) defined on the rectangular grid (u(i),v(j)),
      ! i=1,2,...,mu;j=1,2,...,mv, fpopdi determines a smooth bicubic spline approximation with
      ! given knots tu(i),i=1,..,nu in the u-direction and tv(j),j=1,2,...,nv in the v-direction.
      ! this spline sp(u,v) will be periodic in the variable v and will satisfy the following
      ! constraints
      !
      !     s(tu(1),v) = dz(1) , tv(4) <=v<= tv(nv-3)
      !
      !  and (if iopt(2) = 1)
      !
      !     d s(tu(1),v)
      !     ------------ =  dz(2)*cos(v)+dz(3)*sin(v) , tv(4) <=v<= tv(nv-3)
      !     d u
      !
      !  and (if iopt(3) = 1)
      !
      !     s(tu(nu),v)  =  0   tv(4) <=v<= tv(nv-3)
      !
      ! where the parameters dz(i) correspond to the derivative values g(i,j) as defined in
      ! subroutine pogrid.
      !
      ! the b-spline coefficients of sp(u,v) are determined as the least-squares solution  of an
      ! overdetermined linear system which depends on the value of p and on the values dz(i),i=1,2,3.
      ! the corresponding sum of squared residuals sq is a simple quadratic function in the variables
      ! dz(i). these may or may not be provided. the values dz(i) which are not given will be
      ! determined so as to minimize the resulting sum of squared residuals sq. in that case the user
      ! must provide some initial guess dz(i) and some estimate (dz(i)-step, dz(i)+step) of the range
      ! of possible values for these latter.
      !
      !  sp(u,v) also depends on the parameter p (p>0) in such a way that
      !    - if p tends to infinity, sp(u,v) becomes the least-squares spline with given knots,
      !      satisfying the constraints.
      !    - if p tends to zero, sp(u,v) becomes the least-squares polynomial, satisfying the
      !      constraints.
      !    - the function  f(p)=sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2) is continuous and
      !      strictly decreasing for p>0.
      !
      recursive subroutine fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,ider,tu,nu,tv,nv,&
                                  nuest,nvest,p,step,c,nc,fp,fpu,fpv,nru,nrv,wrk,lwrk)

      !
      !  ..scalar arguments..
      integer, intent(inout) :: ifsu,ifsv,ifbu,ifbv
      integer, intent(in)    :: mu,mv,mz,nu,nv,nc,lwrk,nuest,nvest
      real(RKIND) :: z0,p,step,fp
      !  ..array arguments..
      integer ider(2),nru(mu),nrv(mv),iopt(3)
      real(RKIND) u(mu),v(mv),z(mz),dz(3),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv)
      real(RKIND), intent(inout) :: wrk(lwrk)
      !  ..local scalars..
      real(RKIND) res,sq,sqq,step1,step2
      integer i,id0,iop0,iop1,i1,j,l,laa,lau,lav1,lav2,lbb,lbu,lbv, &
       lcc,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number
      !  ..local arrays..
      integer nr(3)
      real(RKIND) delta(3),dzz(3),sum(3),a(6,6),g(6)

      !  we partition the working space
      lsu  = 1
      lsv  = lsu+4*mu
      lri  = lsv+4*mv
      mm   = max0(nuest,mv+nvest)
      lq   = lri+mm
      mvnu = nuest*(mv+nvest-8)
      lau  = lq+mvnu
      lav1 = lau+5*nuest
      lav2 = lav1+6*nvest
      lbu  = lav2+4*nvest
      lbv  = lbu+5*nuest
      laa  = lbv+5*nvest
      lbb  = laa+2*mv
      lcc  = lbb+2*nvest
      lcs  = lcc+nvest

      !  we calculate the smoothing spline sp(u,v) according to the input values dz(i),i=1,2,3.
      iop0 = iopt(2)
      iop1 = iopt(3)
      call fpgrdi(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,z,mz,dz,iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,&
                  fp,fpu,fpv,mm,mvnu,wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),&
                  wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),wrk(lcc),wrk(lcs),nru,nrv)

      id0 = ider(1)
      if (id0==0) then
        res = (z0-dz(1))**2
        fp = fp+res
        sq = sq+res
      endif

      ! in case all derivative values dz(i) are given (step<=0) or in case we have spline
      ! interpolation, we accept this spline as a solution.
      if (step<=zero .or. sq<=zero) return

      dzz(1:3) = dz(1:3)

      ! number denotes the number of derivative values dz(i) that still must be optimized.
      ! let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if (id0<=0) then
          number   = 1
          nr   (1) = 1
          delta(1) = step
      endif

      if (iop0/=0 .and. ider(2)==0) then
          step2 = step*three/tu(5)
          nr(number+1) = 2
          nr(number+2) = 3
          delta(number+1) = step2
          delta(number+2) = step2
          number = number+2
      endif

      if(number==0) return

      ! the sum of squared residuals sq is a quadratic polynomial in the parameters g(j). we
      ! determine the unknown coefficients of this polymomial by calculating (number+1)*(number+2)/2
      ! different splines according to specific values for g(j).
      parameter_splines: do i=1,number
         l      = nr(i)
         step1  = delta(i)

         dzz(l) = dz(l)+step1
         call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),&
                     fp,fpu,fpv,mm,mvnu,wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                     wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),wrk(lcc),wrk(lcs),nru,nrv)

         if (id0==0) sum(i) = sum(i)+(z0-dzz(1))**2

         dzz(l) = dz(l)-step1
         call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,&
                     fp,fpu,fpv,mm,mvnu,wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                     wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),wrk(lcc),wrk(lcs),nru,nrv)

         if (id0==0) sqq = sqq+(z0-dzz(1))**2

         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if (a(i,i)<=zero) then
            number = 0
            exit parameter_splines
         endif
         g(i) = (sqq-sum(i))/(step1+step1)
         dzz(l) = dz(l)
      end do parameter_splines

      if (number>1) then
         do i=2,number
            l1      = nr(i)
            step1   = delta(i)
            dzz(l1) = dz(l1)+step1
            i1      = i-1
            do j=1,i1
               l2      = nr(j)
               step2   = delta(j)
               dzz(l2) = dz(l2)+step2

               call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,&
                           fp,fpu,fpv,mm,mvnu,wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                           wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),wrk(lcc),wrk(lcs),nru,nrv)

               if (id0==0) sqq = sqq+(z0-dzz(1))**2
               a (i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
               dzz(l2) = dz(l2)
            end do
            dzz(l1) = dz(l1)
         end do
      endif

      ! the optimal values g(j) are found as the solution of the system
      ! d (sq) / d (g(j)) = 0 , j=1,...,number.
      if (number>0) then
          call fpsysy(a,number,g)
          do i=1,number
             l = nr(i)
             dz(l) = dz(l)+g(i)
          end do
      endif

      ! we determine the spline sp(u,v) according to the optimal values g(j).
      call fpgrdi(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,z,mz,dz,iop0,iop1,tu,nu,tv,nv,p,c,nc,sq, &
                  fp,fpu,fpv,mm,mvnu,wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                  wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),wrk(lcc),wrk(lcs),nru,nrv)
      if (id0==0) fp = fp+(z0-dz(1))**2
      return
      end subroutine fpopdi

      recursive subroutine fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r, &
       mr,r0,r1,dr,iopt,ider,tu,nu,tv,nv,nuest,nvest,p,step,c,nc, &
       fp,fpu,fpv,nru,nrv,wrk,lwrk)

      !  given the set of function values r(i,j) defined on the rectangular
      !  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopsp determines a
      !  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
      !  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
      !  spline sp(u,v) will be periodic in the variable v and will satisfy
      !  the following constraints
      !
      !     s(tu(1),v) = dr(1) , tv(4) <=v<= tv(nv-3)
      !
      !     s(tu(nu),v) = dr(4) , tv(4) <=v<= tv(nv-3)
      !
      !  and (if iopt(2) = 1)
      !
      !     d s(tu(1),v)
      !     ------------ =  dr(2)*cos(v)+dr(3)*sin(v) , tv(4) <=v<= tv(nv-3)
      !     d u
      !
      !  and (if iopt(3) = 1)
      !
      !     d s(tu(nu),v)
      !     ------------- =  dr(5)*cos(v)+dr(6)*sin(v) , tv(4) <=v<= tv(nv-3)
      !     d u
      !
      !  where the parameters dr(i) correspond to the derivative values at the
      !  poles as defined in subroutine spgrid.
      !
      !  the b-spline coefficients of sp(u,v) are determined as the least-
      !  squares solution  of an overdetermined linear system which depends
      !  on the value of p and on the values dr(i),i=1,...,6. the correspond-
      !  ing sum of squared residuals sq is a simple quadratic function in
      !  the variables dr(i). these may or may not be provided. the values
      !  dr(i) which are not given will be determined so as to minimize the
      !  resulting sum of squared residuals sq. in that case the user must
      !  provide some initial guess dr(i) and some estimate (dr(i)-step,
      !  dr(i)+step) of the range of possible values for these latter.
      !
      !  sp(u,v) also depends on the parameter p (p>0) in such a way that
      !    - if p tends to infinity, sp(u,v) becomes the least-squares spline
      !      with given knots, satisfying the constraints.
      !    - if p tends to zero, sp(u,v) becomes the least-squares polynomial,
      !      satisfying the constraints.
      !    - the function  f(p)=sumi=1,mu(sumj=1,mv((r(i,j)-sp(u(i),v(j)))**2)
      !      is continuous and strictly decreasing for p>0.
      !
      !  ..scalar arguments..
      integer :: ifsu,ifsv,ifbu,ifbv,mu,mv,mr,nu,nv,nuest,nvest,nc,lwrk
      real(RKIND) :: r0,r1,fp
      real(RKIND), intent(out) :: p
      !  ..array arguments..
      integer :: ider(4),nru(mu),nrv(mv),iopt(3)
      real(RKIND) :: u(mu),v(mv),r(mr),dr(6),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv),wrk(lwrk),step(2)
      !  ..local scalars..
      real(RKIND) :: sq,sqq,sq0,sq1,step1,step2
      integer :: i,id0,iop0,iop1,i1,j,l,lau,lav1,lav2,la0,la1,lbu,lbv,lb0, &
       lb1,lc0,lc1,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number, id1
      !  ..local arrays..
      integer :: nr(6)
      logical :: zeroed
      real(RKIND) :: delta(6),drr(6),sum(6),a(6,6),g(6)

      !  we partition the working space
      lsu  = 1
      lsv  = lsu+4*mu
      lri  = lsv+4*mv
      mm   = max(nuest,mv+nvest)
      lq   = lri+mm
      mvnu = nuest*(mv+nvest-8)
      lau  = lq+mvnu
      lav1 = lau+5*nuest
      lav2 = lav1+6*nvest
      lbu  = lav2+4*nvest
      lbv  = lbu+5*nuest
      la0  = lbv+5*nvest
      la1  = la0+2*mv
      lb0  = la1+2*mv
      lb1  = lb0+2*nvest
      lc0  = lb1+2*nvest
      lc1  = lc0+nvest
      lcs  = lc1+nvest
      !  we calculate the smoothing spline sp(u,v) according to the input
      !  values dr(i),i=1,...,6.
      iop0 = iopt(2)
      iop1 = iopt(3)
      id0  = ider(1)
      id1  = ider(3)
      call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr,                 &
                  iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,      &
                  wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),   &
                  wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),  &
                  wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)

      sq0 = merge((r0-dr(1))**2,zero,id0==0)
      sq1 = merge((r1-dr(4))**2,zero,id1==0)
      sq  = sq+sq0+sq1
      ! in case all derivative values dr(i) are given (step<=0) or in case
      ! we have spline interpolation, we accept this spline as a solution.
      if (sq<=zero .or. all(step(1:2)<=zero)) return
      drr = dr

      ! number denotes the number of derivative values dr(i) that still must
      ! be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if (id0<=0) then
         number   = 1
         nr   (1) = 1
         delta(1) = step(1)
      endif
      if (iop0/=0 .and. ider(2)==0) then
         step2 = step(1)*three/(tu(5)-tu(4))
         nr   (number+1:number+2) = [2,3]
         delta(number+1:number+2) = step2
         number = number+2
      endif
      if (id1<=0) then
         number = number+1
         nr   (number) = 4
         delta(number) = step(2)
      end if
      if (iop1/=0 .and. ider(4)==0) then
         step2 = step(2)*three/(tu(nu)-tu(nu-4))
         nr   (number+1:number+2) = [5,6]
         delta(number+1:number+2) = step2
         number = number+2
      endif

      if(number==0) return

      ! the sum of squared residulas sq is a quadratic polynomial in the parameters g(j).
      ! we determine the unknown coefficients of this polymomial by calculating (number+1)*(number+2)/2
      ! different splines according to specific values for g(j).
      zeroed = .false.
      do i=1,number
         l      = nr(i)
         step1  = delta(i)
         drr(l) = dr(l)+step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr, &
                     iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),fp,fpu,fpv,mm,mvnu, &
                     wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                     wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                     wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)

         if (id0==0) sq0 = (r0-drr(1))**2
         if (id1==0) sq1 = (r1-drr(4))**2
         sum(i) = sum(i)+sq0+sq1
         drr(l) = dr(l)-step1

         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr, &
                     iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu, &
                     wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                     wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                     wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)

         if (id0==0) sq0 = (r0-drr(1))**2
         if (id1==0) sq1 = (r1-drr(4))**2

         sqq    = sqq+sq0+sq1
         drr(l) = dr(l)
         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if (a(i,i)<=zero) then
             zeroed = .true.
             exit
         end if
         g(i)   = (sqq-sum(i))/(2*step1)
      end do

      if (.not.zeroed) then

          if (number>1) then
             do i=2,number
                l1      = nr(i)
                step1   = delta(i)
                drr(l1) = dr(l1)+step1
                i1 = i-1
                do j=1,i1
                    l2      = nr(j)
                    step2   = delta(j)
                    drr(l2) = dr(l2)+step2
                    call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr, &
                                iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu, &
                                wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                                wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                                wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
                    if(id0==0) sq0 = (r0-drr(1))**2
                    if(id1==0) sq1 = (r1-drr(4))**2
                    sqq = sqq+sq0+sq1
                    a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
                    drr(l2) = dr(l2)
                end do
                drr(l1) = dr(l1)
             end do
          end if

          ! the optimal values g(j) are found as the solution of the system
          ! d (sq) / d (g(j)) = 0 , j=1,...,number.
          call fpsysy(a,number,g)
          forall(i=1:number) dr(nr(i))=dr(nr(i))+g(i)

      endif

      ! we determine the spline sp(u,v) according to the optimal values g(j).
      call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr, &
                  iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu, &
                  wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1), &
                  wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0), &
                  wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      if(id0==0) sq0 = (r0-dr(1))**2
      if(id1==0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
      return
      end subroutine fpopsp


      !  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m according to the panel
      !  tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong to. for each panel a stack is constructed
      !  containing the numbers of data points lying inside; index(j),j=1,2,...,nreg points to the first
      !  data point in the jth panel while nummer(i),i=1,2,...,m gives the number of the next data point
      !  in the panel.
      pure subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
      !  ..scalar arguments..
      integer,     intent(in)  :: m,kx,ky,nx,ny,nreg
      !  ..array arguments..
      real(RKIND), intent(in)  :: x(m),y(m),tx(nx),ty(ny)
      integer,     intent(out) :: nummer(m),index(nreg)
      !  ..local scalars..
      real(RKIND) :: xi,yi
      integer :: im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy

      !  ..
      kx1 = kx+1
      ky1 = ky+1
      nk1x = nx-kx1
      nk1y = ny-ky1
      nyy = nk1y-ky
      index = 0

      points: do im=1,m
         xi = x(im)
         yi = y(im)
         l = kx1
         l1 = l+1
         do while (.not.(xi<tx(l1) .or. l==nk1x))
            l = l1
            l1 = l+1
         end do

         k = ky1
         k1 = k+1
         do while (.not.(yi<ty(k1) .or. k==nk1y))
            k = k1
            k1 = k+1
         end do
         num = (l-kx1)*nyy+k-ky
         nummer(im) = index(num)
         index(num) = im
      end do points
      return
      end subroutine fporde



      subroutine fppara(iopt,idim,m,u,mx,x,w,ub,ue,k,s,nest,tol,maxit, &
                        k1,k2,n,t,nc,c,fp,fpint,z,a,b,g,q,nrdata,ier)
      !  ..
      !  ..scalar arguments..
      real(RKIND) :: ub,ue,s,tol,fp
      integer, intent(in) :: idim,maxit
      integer :: iopt,m,mx,k,nest,k1,k2,n,nc,ier
      !  ..array arguments..
      real(RKIND) :: u(m),x(mx),w(m),t(nest),c(nc),fpint(nest),z(nc),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer :: nrdata(nest)
      !  ..local scalars..
      real(RKIND) :: acc,cos,fac,fpart,fpms,fpold,fp0,f1,f2,f3,p,pinv,piv,p1,p2,p3,rn,sin,store,term,ui,wi
      integer :: i,it,iter,i1,i2,i3,j,jj,j1,j2,k3,l,l0,mk1,nk1,nmax,nmin,nplus,npl1,nrint,n8
      logical :: new,check1,check3
      !  ..local arrays..
      real(RKIND) :: h(SIZ_K+1),xi(idim)

      !  set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND

      ! *****
      !  part 1: determination of the number of knots and their position
      ! *****
      !  given a set of knots we compute the least-squares curve sinf(u), and the corresponding sum
      !  of squared residuals fp=f(p=inf).
      !  if iopt=-1 sinf(u) is the requested curve.
      !  if iopt=0 or iopt=1 we check whether we can accept the knots:
      !    if fp <=s we will continue with the current set of knots.
      !    if fp > s we will increase the number of knots and compute the corresponding least-
      !    squares curve until finally fp<=s.
      !  the initial choice of knots depends on the value of s and iopt.
      !    if s=0 we have spline interpolation; in that case the number of knots equals nmax = m+k+1.
      !    if (s>0 and iopt=0) we first compute the least-squares polynomial curve of degree k;
      !      n = nmin = 2*k+2
      !    iopt=1 we start with the set of knots found at the last call of the routine, except for
      !    the case that s > fp0; then we compute directly the polynomial curve of degree k.
      ! *****

      !  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      bootstrap: if (iopt>=0) then

          !  calculation of acc, the absolute tolerance for the root of f(p)=s.
          acc = tol*s

          !  determine nmax, the number of knots for spline interpolation.
          nmax = m+k1

          interpolating: if (s<=zero) then

              !  if s=0, s(u) is an interpolating curve.
              !  check that the required storage space exceeds the available one.
              n = nmax
              if (nmax>nest) then
                 ier = FITPACK_INSUFFICIENT_STORAGE
                 return
              end if

              !  find the position of the interior knots in case of interpolation.
              mk1 = m-k1
              if (mk1/=0) then
                  k3 = k/2
                  i  = k2
                  j  = k3+2
                  do l=1,mk1
                    t(i) = merge( u(j) , (u(j)+u(j-1))*half , k3*2/=k)
                    i = i+1
                    j = j+1
                  end do
              endif

          else interpolating

              !  if s>0 our initial choice of knots depends on the value of iopt.
              !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
              !  polynomial curve which is a spline curve without interior knots.
              !  if iopt=1 and fp0>s we start computing the least squares spline curve
              !  according to the set of knots found at the last call of the routine.
              use_last_call: if (iopt/=0 .and. n/=nmin)  then
                 fp0   = fpint(n)
                 fpold = fpint(n-1)
                 nplus = nrdata(n)
                 if (fp0<=s) then
                    n         = nmin
                    fpold     = zero
                    nplus     = 0
                    nrdata(1) = m-2
                 end if
              else use_last_call
                 n = nmin
                 fpold = zero
                 nplus = 0
                 nrdata(1) = m-2
              endif use_last_call
          endif interpolating

      endif bootstrap

      !  main loop for the different sets of knots. m is a safe upper bound
      !  for the number of trials.
      iter = 0
      main_loop: do while (iter<=m)

        iter = iter+1

        if (n==nmin) ier = FITPACK_LEASTSQUARES_OK

        ! find nrint, tne number of knot intervals.
        nrint = n-nmin+1
        ! find the position of the additional knots which are needed for
        ! the b-spline representation of s(u).
        nk1 = n-k1
        t(1:k1)    = ub
        t(nk1+1:n) = ue

        ! compute the b-spline coefficients of the least-squares spline curve
        ! sinf(u). the observation matrix a is built up row by row and
        ! reduced to upper triangular form by givens transformations.
        ! at the same time fp=f(p=inf) is computed.
        fp = zero

        ! initialize the b-spline coefficients and the observation matrix a.
        z(1:nc) = zero
        a(1:nk1,1:k1) = zero
        l = k1
        coefs: do it=1,m

           ! fetch the current data point u(it),x(it).
           ui = u(it)
           wi = w(it)
           xi = x((it-1)*idim+1:it*idim)*wi

           ! search for knot interval t(l) <= ui < t(l+1).
           do while (ui>=t(l+1) .and. l/=nk1)
              l = l+1
           end do

           ! evaluate the (k+1) non-zero b-splines at ui and store them in q.
           call fpbspl(t,n,k,ui,l,h)
           q(it,1:k1) = h(1:k1)
           h(:k1) = wi*h(:k1)

           ! rotate the new row of the observation matrix into triangle.
           j = l-k1
           rotate_row: do i=1,k1

              j = j+1

              piv = h(i); if (piv==zero) cycle rotate_row

              ! calculate the parameters of the givens transformation.
              call fpgivs(piv,a(j,1),cos,sin)

              ! transformations to right hand side.
              call fprota(cos,sin,xi,z(j:j+idim*n:n))

              ! transformations to left hand side.
              not_last: if (i<k1) then
                 i2 = 1
                 i3 = i+1
                 do i1 = i3,k1
                   i2 = i2+1
                   call fprota(cos,sin,h(i1),a(j,i2))
                 end do
              endif not_last
           end do rotate_row

           !  add contribution of this row to the sum of squares of residual right hand sides.
           fp = fp + sum(xi**2)

        end do coefs

        if (ier==FITPACK_LEASTSQUARES_OK) fp0 = fp
        fpint(n-1:n) = [fpold,fp0]
        nrdata(n)    = nplus

        ! backward substitution to obtain the b-spline coefficients.
        j1 = 1
        do j2=1,idim
           c(j1:j1+nk1-1) = fpback(a,z(j1),nk1,k1,nest)
           j1 = j1+n
        end do

        ! test whether the approximation sinf(u) is an acceptable solution.
        if (iopt<0) return ! was done already

        fpms = fp-s; if(abs(fpms)<acc) return

        ! if f(p=inf) < s accept the choice of knots.
        if(fpms<zero) exit main_loop

        ! if n = nmax, sinf(u) is an interpolating spline curve.
        if (n==nmax) then
           ier = FITPACK_INTERPOLATING_OK
           return
        endif

        ! increase the number of knots.
        ! if n=nest we cannot increase the number of knots because of the storage capacity limitation.
        if (n==nest) then
            ier = FITPACK_INSUFFICIENT_STORAGE
            return
        end if

        ! determine the number of knots nplus we are going to add.
        if (ier==FITPACK_OK) then
           npl1 = nplus*2
           rn = nplus
           if (fpold-fp>acc) npl1 = int(rn*fpms/(fpold-fp))
           nplus = min(nplus*2,max(npl1,nplus/2,1))
        else
           nplus = 1
           ier   = FITPACK_OK
        endif

        ! Initialize iterate
        fpold = fp

        ! compute the sum of squared residuals for each knot interval
        ! t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = zero
        i = 1
        l = k2
        new = .false.
        jj = 0
        square_residuals: do it=1,m
            if (u(it)>=t(l) .and. l<=nk1) then
               new = .true.
               l = l+1
            endif
            term = zero
            l0 = l-k2
            do j2=1,idim
               fac  = dot_product(c(l0+1:l0+k1),q(it,1:k1))
               jj   = jj+1
               term = term+(w(it)*(fac-x(jj)))**2
               l0   = l0+n
            end do
            fpart = fpart+term
            if (new) then
               store    = term*half
               fpint(i) = fpart-store
               i        = i+1
               fpart    = store
               new      = .false.
            endif
        end do square_residuals

        fpint(nrint) = fpart
        add_new_knots: do l=1,nplus

            ! add a new knot.
            call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1)

            ! if n=nmax we locate the knots as for interpolation
            if (n==nmax) then
               !  find the position of the interior knots in case of interpolation.
               mk1 = m-k1
               if (mk1/=0) then
                  k3 = k/2
                  i  = k2
                  j  = k3+2
                  do jj=1,mk1
                    t(i) = merge( u(j) , (u(j)+u(j-1))*half , k3*2/=k)
                    i = i+1
                    j = j+1
                  end do
               endif

               ! Restart main loop
               iter = 0
               cycle main_loop

            end if

            ! test whether we cannot further increase the number of knots.
            if (n==nest) exit add_new_knots
        end do add_new_knots
      !  restart the computations with the new set of knots.
      end do main_loop

      ! test whether the least-squares kth degree polynomial curve is a solution
      ! of our approximation problem.
      if (ier==FITPACK_LEASTSQUARES_OK) return

      ! *****
      !  part 2: determination of the smoothing spline curve sp(u).
      ! *****
      !  we have determined the number of knots and their position.
      !  we now compute the b-spline coefficients of the smoothing curve sp(u). the observation matrix a
      !  is extended by the rows of matrix b expressing that the kth derivative discontinuities of sp(u)
      !  at the interior knots t(k+2),...t(n-k-1) must be zero. the corresponding weights of these
      !  additional rows are set to 1/p.
      !  iteratively we then have to determine the value of p such that f(p), the sum of squared
      !  residuals be = s. we already know that the least squares kth degree polynomial curve corresponds
      !  to p=0, and that the least-squares spline curve corresponds to p=infinity. the iteration process
      !  which is proposed here, makes use of rational interpolation. since f(p) is a convex and strictly
      !  decreasing function of p, it can be approximated by a rational function r(p) = (u*p+v)/(p+w).
      !  three values of p(p1,p2,p3) with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)
      !  are used to calculate the new value of p such that r(p)=s. convergence is guaranteed by taking
      !  f1>0 and f3<zero
      ! *****

      ! evaluate the discontinuity jump of the kth derivative of the b-splines at the knots
      ! t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)

      !  initial value for p.
      p1 = zero
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p  = sum(a(1:nk1,1))
      rn = nk1
      p  = rn/p
      check1 = .false.
      check3 = .false.
      n8 = n-nmin

      !  iteration process to find the root of f(p) = s.
      iter = 0
      find_root: do while (iter<maxit)

         iter = iter+1

         ! the rows of matrix b with weight 1/p are rotated into the
         ! triangularised observation matrix a which is stored in g.
         pinv = one/p
         c(1:nc)       = z(1:nc)
         g(1:nk1,k2)   = zero
         g(1:nk1,1:k1) = a(1:nk1,1:k1)
         b_rows: do it=1,n8
            ! the row of matrix b is rotated into triangle by givens transformation
            h(1:k2) = b(it,1:k2)*pinv
            xi = zero
            b_cols: do j=it,nk1
               piv = h(1)

               ! calculate the parameters of the givens transformation.
               call fpgivs(piv,g(j,1),cos,sin)

               ! transformations to right hand side.
               call fprota(cos,sin,xi,c(j:j+idim*n:n))
               if (j==nk1) cycle b_rows

               !  transformations to left hand side.
               i2 = merge(nk1-j,k1,j>n8)
               do i=1,i2
                  call fprota(cos,sin,h(i+1),g(j,i+1))
                  h(i) = h(i+1)
               end do
               h(i2+1) = zero
            end do b_cols
         end do b_rows

         ! backward substitution to obtain the b-spline coefficients.
         j1 = 1
         do j2=1,idim
            c(j1:j1+nk1-1) = fpback(g,c(j1),nk1,k2,nest)
            j1 =j1+n
         end do

         ! computation of f(p).
         fp = zero
         l = k2
         jj = 0
         get_fp: do it=1,m
            if (u(it)>=t(l) .and. l<=nk1) l = l+1
            l0 = l-k2
            term = zero
            do j2=1,idim
              fac  = dot_product(c(l0+1:l0+k1),q(it,1:k1))
              jj   = jj+1
              term = term+(fac-x(jj))**2
              l0   = l0+n
            end do
            fp = fp+term*w(it)**2
         end do get_fp

         ! SUCCESS! the approximation sp(u) is an acceptable solution.
         fpms = fp-s; if(abs(fpms)<acc) return

         ! Reinitialize p to carry out one more iteration
         p2 = p
         f2 = fpms
         if (.not.check3) then
            if((f2-f3)>acc) then
               if (f2<zero) check3=.true.
            else
               ! our initial choice of p is too large.
               p3 = p2
               f3 = f2
               p  = p*con4
               if (p<=p1) p=p1*con9 + p2*con1
               cycle find_root
            endif
         endif

         if (.not.check1) then
            if ((f1-f2)>acc) then
                if (f2>zero) check1 = .true.
            else
                ! our initial choice of p is too small
                p1 = p2
                f1 = f2
                p = p/con4
                if (p3<zero) cycle find_root
                if (p>=p3) p = p2*con1 + p3*con9
                cycle find_root
            endif
         endif

         ! test whether the iteration process proceeds as theoretically expected.
         if (f2>=f1 .or. f2<=f3) then
            ier = FITPACK_S_TOO_SMALL
            return
         endif

         ! find the new value for p.
         call fprati(p1,f1,p2,f2,p3,f3,p)

      end do find_root

      ! Maximum number of iterations reached
      ier = FITPACK_MAXIT
      return

      end subroutine fppara


      subroutine fppasu(iopt,ipar,idim,u,mu,v,mv,z,mz,s,nuest,nvest, &
       tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu,reducv,fpintu, &
       fpintv,lastdi,nplusu,nplusv,nru,nrv,nrdatu,nrdatv,wrk,lwrk,ier)

      !  ..
      !  ..scalar arguments..
      real(RKIND) s,tol,fp,fp0,fpold,reducu,reducv
      integer iopt,idim,mu,mv,mz,nuest,nvest,maxit,nc,nu,nv,lastdi, &
       nplusu,nplusv,lwrk,ier
      !  ..array arguments..
      real(RKIND) u(mu),v(mv),z(mz*idim),tu(nuest),tv(nvest),c(nc*idim), &
       fpintu(nuest),fpintv(nvest),wrk(lwrk)
      integer ipar(2),nrdatu(nuest),nrdatv(nvest),nru(mu),nrv(mv)
      !  ..local scalars
      real(RKIND) acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,peru,perv,ub,ue,vb,ve
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,iter,j,lau1,lav1,laa, &
       l,lau,lav,lbu,lbv,lq,lri,lsu,lsv,l1,l2,l3,l4,mm,mpm,mvnu,ncof, &
       nk1u,nk1v,nmaxu,nmaxv,nminu,nminv,nplu,nplv,npl1,nrintu, &
       nrintv,nue,nuk,nve,nuu,nvv

      !   set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND

      !  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s

      !  find nmaxu and nmaxv which denote the number of knots in u- and v-
      !  direction in case of spline interpolation.
      nmaxu = mu+4+2*ipar(1)
      nmaxv = mv+4+2*ipar(2)

      !  find nue and nve which denote the maximum number of knots
      !  allowed in each direction
      nue = min0(nmaxu,nuest)
      nve = min0(nmaxv,nvest)

      !  set boundaries of the approximation domain
      ub = u(1)
      ue = u(mu)
      vb = v(1)
      ve = v(mv)
      !  we partition the working space.
      lsu = 1
      lsv = lsu+mu*4
      lri = lsv+mv*4
      mm = max0(nuest,mv)
      lq = lri+mm*idim
      mvnu = nuest*mv*idim
      lau = lq+mvnu
      nuk = nuest*5
      lbu = lau+nuk
      lav = lbu+nuk
      nuk = nvest*5
      lbv = lav+nuk
      laa = lbv+nuk
      lau1 = lau
      if(ipar(1)==0) go to 10
      peru = ue-ub
      lau1 = laa
      laa = laa+4*nuest
  10  lav1 = lav
      if(ipar(2)==0) go to 20
      perv = ve-vb
      lav1 = laa
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 1: determination of the number of knots and their position.     c
      ! ****************************************************************     c
      !  given a set of knots we compute the least-squares spline sinf(u,v), c
      !  and the corresponding sum of squared residuals fp=f(p=inf).         c
      !  if iopt=-1  sinf(u,v) is the requested approximation.               c
      !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
      !    if fp <=s we will continue with the current set of knots.         c
      !    if fp > s we will increase the number of knots and compute the    c
      !       corresponding least-squares spline until finally fp<=s.        c
      !    the initial choice of knots depends on the value of s and iopt.   c
      !    if s=0 we have spline interpolation; in that case the number of   c
      !    knots equals nmaxu = mu+4+2*ipar(1) and  nmaxv = mv+4+2*ipar(2)   c
      !    if s>0 and                                                        c
      !     *iopt=0 we first compute the least-squares polynomial            c
      !          nu=nminu=8 and nv=nminv=8                                   c
      !     *iopt=1 we start with the knots found at the last call of the    c
      !      routine, except for the case that s > fp0; then we can compute  c
      !      the least-squares polynomial directly.                          c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  determine the number of knots for polynomial approximation.
  20  nminu = 8
      nminv = 8
      if(iopt<0) go to 100

      if(s>zero) go to 60
      !  if s = 0, s(u,v) is an interpolating spline.
      nu = nmaxu
      nv = nmaxv
      !  test whether the required storage space exceeds the available one.
      if(nv>nvest .or. nu>nuest) go to 420
      !  find the position of the interior knots in case of interpolation.
      !  the knots in the u-direction.
      nuu = nu-8
      if(nuu==0) go to 40
      i = 5
      j = 3-ipar(1)
      do 30 l=1,nuu
        tu(i) = u(j)
        i = i+1
        j = j+1
  30  continue
      !  the knots in the v-direction.
  40  nvv = nv-8
      if(nvv==0) go to 60
      i = 5
      j = 3-ipar(2)
      do 50 l=1,nvv
        tv(i) = v(j)
        i = i+1
        j = j+1
  50  continue
      go to 100
      !  if s > 0 our initial choice of knots depends on the value of iopt.
  60  if(iopt==0) go to 90
      if(fp0<=s) go to 90
      !  if iopt=1 and fp0 > s we start computing the least- squares spline
      !  according to the set of knots found at the last call of the routine.
      !  we determine the number of grid coordinates u(i) inside each knot
      !  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mpm = mu-1
      do 70 i=2,mpm
        nrdatu(j) = nrdatu(j)+1
        if(u(i)<tu(l)) go to 70
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
  70  continue
      !  we determine the number of grid coordinates v(i) inside each knot
      !  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      mpm = mv-1
      do 80 i=2,mpm
        nrdatv(j) = nrdatv(j)+1
        if(v(i)<tv(l)) go to 80
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
  80  continue
      go to 100
      !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
      !  polynomial (which is a spline without interior knots).
  90  nu = nminu
      nv = nminv
      nrdatu(1) = mu-2
      nrdatv(1) = mv-2
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
 100  mpm = mu+mv
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
      !  main loop for the different sets of knots.mpm=mu+mv is a save upper
      !  bound for the number of trials.
      do 250 iter=1,mpm
        if(nu==nminu .and. nv==nminv) ier = -2
      !  find nrintu (nrintv) which is the number of knot intervals in the
      !  u-direction (v-direction).
        nrintu = nu-nminu+1
        nrintv = nv-nminv+1
      !  find ncof, the number of b-spline coefficients for the current set
      !  of knots.
        nk1u = nu-4
        nk1v = nv-4
        ncof = nk1u*nk1v
      !  find the position of the additional knots which are needed for the
      !  b-spline representation of s(u,v).
        if(ipar(1)/=0) go to 110
        i = nu
        do 105 j=1,4
          tu(j) = ub
          tu(i) = ue
          i = i-1
 105    continue
        go to 120
 110    l1 = 4
        l2 = l1
        l3 = nu-3
        l4 = l3
        tu(l2) = ub
        tu(l3) = ue
        do 115 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tu(l2) = tu(l4)-peru
          tu(l3) = tu(l1)+peru
 115    continue
 120    if(ipar(2)/=0) go to 130
        i = nv
        do 125 j=1,4
          tv(j) = vb
          tv(i) = ve
          i = i-1
 125    continue
        go to 140
 130    l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 135 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-perv
          tv(l3) = tv(l1)+perv
 135    continue
      !  find the least-squares spline sinf(u,v) and calculate for each knot
      !  interval tu(j+3)<=u<=tu(j+4) (tv(j+3)<=v<=tv(j+4)) the sum
      !  of squared residuals fpintu(j),j=1,2,...,nu-7 (fpintv(j),j=1,2,...
      !  ,nv-7) for the data points having their absciss (ordinate)-value
      !  belonging to that interval.
      !  fp gives the total sum of squared residuals.
 140    call fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,tu, &
        nu,tv,nv,p,c,nc,fp,fpintu,fpintv,mm,mvnu,wrk(lsu),wrk(lsv), &
        wrk(lri),wrk(lq),wrk(lau),wrk(lau1),wrk(lav),wrk(lav1), &
        wrk(lbu),wrk(lbv),nru,nrv)
        if(ier==(-2)) fp0 = fp
      !  test whether the least-squares spline is an acceptable solution.
        if(iopt<0) go to 440
        fpms = fp-s
        if(abs(fpms) < acc) go to 440
      !  if f(p=inf) < s, we accept the choice of knots.
        if(fpms<0.) go to 300
      !  if nu=nmaxu and nv=nmaxv, sinf(u,v) is an interpolating spline.
        if(nu==nmaxu .and. nv==nmaxv) go to 430
      !  increase the number of knots.
      !  if nu=nue and nv=nve we cannot further increase the number of knots
      !  because of the storage capacity limitation.
        if(nu==nue .and. nv==nve) go to 420
        ier = 0
      !  adjust the parameter reducu or reducv according to the direction
      !  in which the last added knots were located.
        if (lastdi<0) go to 150
        if (lastdi==0) go to 170
        go to 160
 150    reducu = fpold-fp
        go to 170
 160    reducv = fpold-fp
      !  store the sum of squared residuals for the current set of knots.
 170    fpold = fp
      !  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu==nminu) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu>acc) npl1 = int(rn*fpms/reducu)
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
      !  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 1
        if(nv==nminv) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv>acc) npl1 = int(rn*fpms/reducv)
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
 190    if (nplu<nplv) go to 210
        if (nplu==nplv) go to 200
        go to 230
 200    if(lastdi<0) go to 230
 210    if(nu==nue) go to 230
      !  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        do 220 l=1,nplusu
      !  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,1)
      !  test whether we cannot further increase the number of knots in the
      !  u-direction.
          if(nu==nue) go to 250
 220    continue
        go to 250
 230    if(nv==nve) go to 210
      !  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
      !  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
      !  test whether we cannot further increase the number of knots in the
      !  v-direction.
          if(nv==nve) go to 250
 240    continue
      !  restart the computations with the new set of knots.
 250  continue
      !  test whether the least-squares polynomial is a solution of our
      !  approximation problem.
 300  if(ier==(-2)) go to 440
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 2: determination of the smoothing spline sp(u,v)                c
      ! *****************************************************                c
      !  we have determined the number of knots and their position. we now   c
      !  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
      !  this smoothing spline varies with the parameter p in such a way thatc
      !  f(p)=suml=1,idim(sumi=1,mu(sumj=1,mv((z(i,j,l)-sp(u(i),v(j),l))**2) c
      !  is a continuous, strictly decreasing function of p. moreover the    c
      !  least-squares polynomial corresponds to p=0 and the least-squares   c
      !  spline to p=infinity. iteratively we then have to determine the     c
      !  positive value of p such that f(p)=s. the process which is proposed c
      !  here makes use of rational interpolation. f(p) is approximated by a c
      !  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
      !  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
      !  are used to calculate the new value of p such that r(p)=s.          c
      !  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      ich1 = 0
      ich3 = 0
      !  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
      !  find the smoothing spline sp(u,v) and the corresponding sum of
      !  squared residuals fp.
        call fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,tu, &
        nu,tv,nv,p,c,nc,fp,fpintu,fpintv,mm,mvnu,wrk(lsu),wrk(lsv), &
        wrk(lri),wrk(lq),wrk(lau),wrk(lau1),wrk(lav),wrk(lav1), &
        wrk(lbu),wrk(lbv),nru,nrv)
      !  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 440
      !  test whether the maximum allowable number of iterations has been
      !  reached.
        if(iter==maxit) go to 400
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 320
        if((f2-f3)>acc) go to 310
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2<0.) ich3 = 1
 320    if(ich1/=0) go to 340
        if((f1-f2)>acc) go to 330
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 350
        if(p>=p3) p = p2*con1 + p3*con9
        go to 350
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 330    if(f2>zero) ich1 = 1
 340    if(f2>=f1 .or. f2<=f3) go to 410
      !  find the new value of p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 350  continue
      !  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end subroutine fppasu


      recursive subroutine fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit, &
         k1,k2,n,t,c,fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)

      !  ..
      !  ..scalar arguments..
      real(RKIND) s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),w(m),t(nest),c(nest),fpint(nest),z(nest), &
       a1(nest,k1),a2(nest,k),b(nest,k2),g1(nest,k2),g2(nest,k1), &
       q(m,k1)
      integer nrdata(nest)
      !  ..local scalars..
      real(RKIND) acc,cos,c1,d1,fpart,fpms,fpold,fp0,f1,f2,f3,p,per,pinv,piv, &
       p1,p2,p3,sin,store,term,wi,xi,yi,rn
      integer i,ich1,ich3,ij,ik,it,iter,i1,i2,i3,j,jk,jper,j1,j2,kk, &
       kk1,k3,l,l0,l1,l5,mm,m1,new,nk1,nk2,nmax,nmin,nplus,npl1, &
       nrint,n10,n11,n7,n8
      !  ..local arrays..
      real(RKIND) h(SIZ_K+1),h1(7),h2(6)

      !  set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 1: determination of the number of knots and their position     c
      !  **************************************************************      c
      !  given a set of knots we compute the least-squares periodic spline   c
      !  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c
      !  the initial choice of knots depends on the value of s and iopt.     c
      !    if s=0 we have spline interpolation; in that case the number of   c
      !    knots equals nmax = m+2*k.                                        c
      !    if s > 0 and                                                      c
      !      iopt=0 we first compute the least-squares polynomial of         c
      !      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c
      !      find that s(x) is a constant function.                          c
      !      iopt=1 we start with the set of knots found at the last         c
      !      call of the routine, except for the case that s > fp0; then     c
      !      we compute directly the least-squares periodic polynomial.      c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m1 = m-1
      kk = k
      kk1 = k1
      k3 = 3*k+1
      nmin = 2*k1
      !  determine the length of the period of s(x).
      per = x(m)-x(1)
      if(iopt<0) go to 50
      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      !  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k
      if(s>0. .or. nmax==nmin) go to 30
      !  if s=0, s(x) is an interpolating spline.
      n = nmax
      !  test whether the required storage space exceeds the available one.
      if(n>nest) go to 620
      !  find the position of the interior knots in case of interpolation.
   5  if((k/2)*2 == k) go to 20
      do 10 i=2,m1
        j = i+k
        t(j) = x(i)
  10  continue
      if(s>0.) go to 50
      kk = k-1
      kk1 = k
      if(kk>0) go to 50
      t(1) = t(m)-per
      t(2) = x(1)
      t(m+1) = x(m)
      t(m+2) = t(3)+per
      do 15 i=1,m1
        c(i) = y(i)
  15  continue
      c(m) = c(1)
      fp = 0.
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
      go to 630
  20  do 25 i=2,m1
        j = i+k
        t(j) = (x(i)+x(i-1))*half
  25  continue
      go to 50
      !  if s > 0 our initial choice depends on the value of iopt.
      !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
      !  periodic polynomial. (i.e. a constant function).
      !  if iopt=1 and fp0>s we start computing the least-squares periodic
      !  spline according the set of knots found at the last call of the
      !  routine.
  30  if(iopt==0) go to 35
      if(n==nmin) go to 35
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0>s) go to 50
      !  the case that s(x) is a constant function is treated separetely.
      !  find the least-squares constant c1 and compute fp0 at the same time.
  35  fp0 = 0.
      d1 = 0.
      c1 = 0.
      do 40 it=1,m1
        wi = w(it)
        yi = y(it)*wi
        call fpgivs(wi,d1,cos,sin)
        call fprota(cos,sin,yi,c1)
        fp0 = fp0+yi**2
  40  continue
      c1 = c1/d1
      !  test whether that constant function is a solution of our problem.
      fpms = fp0-s
      if(fpms<acc .or. nmax==nmin) go to 640
      fpold = fp0
      !  test whether the required storage space exceeds the available one.
      if(nmin>=nest) go to 620
      !  start computing the least-squares periodic spline with one
      !  interior knot.
      nplus = 1
      n = nmin+1
      mm = (m+1)/2
      t(k2) = x(mm)
      nrdata(1) = mm-2
      nrdata(2) = m1-mm
      !  main loop for the different sets of knots. m is a save upper
      !  bound for the number of trials.
  50  do 340 iter=1,m
      !  find nrint, the number of knot intervals.
        nrint = n-nmin+1
      !  find the position of the additional knots which are needed for
      !  the b-spline representation of s(x). if we take
      !      t(k+1) = x(1), t(n-k) = x(m)
      !      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
      !      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
      !  then s(x) is a periodic spline with period per if the b-spline
      !  coefficients satisfy the following conditions
      !      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
        t(k1) = x(1)
        nk1 = n-k1
        nk2 = nk1+1
        t(nk2) = x(m)
        do 60 j=1,k
          i1 = nk2+j
          i2 = nk2-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
  60    continue
      !  compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
      !  periodic spline sinf(x). the observation matrix a is built up row
      !  by row while taking into account condition (**) and is reduced to
      !  triangular form by givens transformations .
      !  at the same time fp=f(p=inf) is computed.
      !  the n7 x n7 triangularised upper matrix a has the form
      !            ! a1 '    !
      !        a = !    ' a2 !
      !            ! 0  '    !
      !  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
      !  matrix of bandwidth k+1 ( n10 = n7-k).
      !  initialization.
        z (1:nk1) = zero
        a1(1:nk1,1:kk1) = zero
        n7 = nk1-k
        n10 = n7-kk
        jper = 0
        fp = zero
        l = k1
        do 290 it=1,m1
      !  fetch the current data point x(it),y(it)
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
      !  search for knot interval t(l) <= xi < t(l+1).
  80      if(xi<t(l+1)) go to 85
          l = l+1
          go to 80
      !  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  85      call fpbspl(t,n,k,xi,l,h)
          do 90 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  90      continue
          l5 = l-k1
      !  test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
          if(l5<n10) go to 285
          if(jper/=0) go to 160
      !  initialize the matrix a2.
          a2(1:n7,1:kk) = zero
          jk = n10+1
          do 110 i=1,kk
            ik = jk
            do 100 j=1,kk1
              if(ik<=0) go to 105
              a2(ik,i) = a1(ik,j)
              ik = ik-1
 100        continue
 105        jk = jk+1
 110      continue
          jper = 1
      !  if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
      !  we take account of condition (**) for setting up the new row
      !  of the observation matrix a. this row is stored in the arrays h1
      !  (the part with respect to a1) and h2 (the part with
      !  respect to a2).
 160      h1 = zero
          h2 = zero
          j = l5-n10
          do 210 i=1,kk1
            j = j+1
            l0 = j
 180        l1 = l0-kk
            if(l1<=0) go to 200
            if(l1<=n10) go to 190
            l0 = l1-n10
            go to 180
 190        h1(l1) = h(i)
            go to 210
 200        h2(l0) = h2(l0)+h(i)
 210      continue
      !  rotate the new row of the observation matrix into triangle
      !  by givens transformations.
          if(n10<=0) go to 250
      !  rotation with the rows 1,2,...n10 of matrix a.
          do 240 j=1,n10
            piv = h1(1)
            if(piv/=zero) go to 214
            do 212 i=1,kk
              h1(i) = h1(i+1)
 212        continue
            h1(kk1) = 0.
            go to 240
      !  calculate the parameters of the givens transformation.
 214        call fpgivs(piv,a1(j,1),cos,sin)
      !  transformation to the right hand side.
            call fprota(cos,sin,yi,z(j))
      !  transformations to the left hand side with respect to a2.
            do 220 i=1,kk
              call fprota(cos,sin,h2(i),a2(j,i))
 220        continue
            if(j==n10) go to 250
            i2 = min0(n10-j,kk)
      !  transformations to the left hand side with respect to a1.
            do 230 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),a1(j,i1))
              h1(i) = h1(i1)
 230        continue
            h1(i1) = zero
 240      continue
      !  rotation with the rows n10+1,...n7 of matrix a.
 250      do 270 j=1,kk
            ij = n10+j
            if(ij<=0) go to 270
            piv = h2(j)
            if (piv==zero) go to 270
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a2(ij,j),cos,sin)
      !  transformations to right hand side.
            call fprota(cos,sin,yi,z(ij))
            if(j==kk) go to 280
            j1 = j+1
      !  transformations to left hand side.
            do 260 i=j1,kk
              call fprota(cos,sin,h2(i),a2(ij,i))
 260        continue
 270      continue
      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
 280      fp = fp+yi**2
          go to 290
      !  rotation of the new row of the observation matrix into
      !  triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero
      !  at xi.
 285      j = l5
          do 140 i=1,kk1
            j = j+1
            piv = h(i)
            if (piv==zero) go to 140
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a1(j,1),cos,sin)
      !  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i==kk1) go to 150
            i2 = 1
            i3 = i+1
      !  transformations to left hand side.
            do 130 i1=i3,kk1
              i2 = i2+1
              call fprota(cos,sin,h(i1),a1(j,i2))
 130        continue
 140      continue
      !  add contribution of this row to the sum of squares of residual
      !  right hand sides.
 150      fp = fp+yi**2
 290    continue
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus

        ! backward substitution to obtain the b-spline coefficients c(j),j=1,.n
        c(1:n7) = fpbacp(a1,a2,z,n7,kk,kk1,nest)

        ! calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
        do 295 i=1,k
          j    = i+n7
          c(j) = c(i)
 295    continue
        if(iopt<0) go to 660
      !  test whether the approximation sinf(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 660
      !  if f(p=inf) < s accept the choice of knots.
        if(fpms<0.) go to 350
      !  if n=nmax, sinf(x) is an interpolating spline.
        if(n==nmax) go to 630
      !  increase the number of knots.
      !  if n=nest we cannot increase the number of knots because of the
      !  storage capacity limitation.
        if(n==nest) go to 620
      !  determine the number of knots nplus we are going to add.
        npl1 = nplus*2
        rn = nplus
        if(fpold-fp>acc) npl1 = int(rn*fpms/(fpold-fp))
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
        fpold = fp
      !  compute the sum(wi*(yi-s(xi))**2) for each knot interval
      !  t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = zero
        i = 1
        l = k1
        do 320 it=1,m1
          if(x(it)<t(l)) go to 300
          new = 1
          l = l+1
 300      term = 0.
          l0 = l-k2
          do 310 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 310      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new==0) go to 320
          if(l>k2) go to 315
          fpint(nrint) = term
          new = 0
          go to 320
 315      store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 320    continue
        fpint(nrint) = fpint(nrint)+fpart
        do 330 l=1,nplus
      !  add a new knot
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
      !  if n=nmax we locate the knots as for interpolation.
          if(n==nmax) go to 5
      !  test whether we cannot further increase the number of knots.
          if(n==nest) go to 340
 330    continue
      !  restart the computations with the new set of knots.
 340  continue
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 2: determination of the smoothing periodic spline sp(x).       c
      !  *************************************************************       c
      !  we have determined the number of knots and their position.          c
      !  we now compute the b-spline coefficients of the smoothing spline    c
      !  sp(x). the observation matrix a is extended by the rows of matrix   c
      !  b expressing that the kth derivative discontinuities of sp(x) at    c
      !  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
      !  ponding weights of these additional rows are set to 1/sqrt(p).      c
      !  iteratively we then have to determine the value of p such that      c
      !  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c
      !  the least-squares constant function corresponds to p=0, and that    c
      !  the least-squares periodic spline corresponds to p=infinity. the    c
      !  iteration process which is proposed here, makes use of rational     c
      !  interpolation. since f(p) is a convex and strictly decreasing       c
      !  function of p, it can be approximated by a rational function        c
      !  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
      !  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
      !  to calculate the new value of p such that r(p)=s. convergence is    c
      !  guaranteed by taking f1>0 and f3<zero                                 c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  evaluate the discontinuity jump of the kth derivative of the
      !  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
 350  call fpdisc(t,n,k2,b,nest)
      !  initial value for p.
      p1 = zero
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      n11 = n10-1
      n8 = n7-1
      p = 0.
      l = n7
      do 352 i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if(l==0) go to 356
 352  continue
      do 354 i=1,n10
         p = p+a1(i,1)
 354  continue
 356  rn = n7
      p = rn/p
      ich1 = 0
      ich3 = 0
      !  iteration process to find the root of f(p) = s.
      do 595 iter=1,maxit
      !  form the matrix g  as the matrix a extended by the rows of matrix b.
      !  the rows of matrix b with weight 1/p are rotated into
      !  the triangularised observation matrix a.
      !  after triangularisation our n7 x n7 matrix g takes the form
      !            ! g1 '    !
      !        g = !    ' g2 !
      !            ! 0  '    !
      !  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
      !  matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one/p
      !  store matrix a into g
        c(1:n7)        = z(1:n7)
        g1(1:n7,1:k)   = a1(1:n7,1:k)
        g1(1:n7,k1)    = a1(1:n7,k1)
        g1(1:n7,k2)    = zero
        g2(1:n7,1)     = zero
        g2(1:n7,2:k+1) = a2(1:n7,1:k)

        l = n10
        do 370 j=1,k1
          if(l<=0) go to 375
          g2(l,1) = a1(l,j)
          l = l-1
 370    continue
 375    do 540 it=1,n8
      !  fetch a new row of matrix b and store it in the arrays h1 (the part
      !  with respect to g1) and h2 (the part with respect to g2).
          yi = 0.
          do 380 i=1,k1
            h1(i) = 0.
            h2(i) = 0.
 380      continue
          h1(k2) = 0.
          if(it>n11) go to 420
          l = it
          l0 = it
          do 390 j=1,k2
            if(l0==n10) go to 400
            h1(j) = b(it,j)*pinv
            l0 = l0+1
 390      continue
          go to 470
 400      l0 = 1
          do 410 l1=j,k2
            h2(l0) = b(it,l1)*pinv
            l0 = l0+1
 410      continue
          go to 470
 420      l = 1
          i = it-n10
          do 460 j=1,k2
            i = i+1
            l0 = i
 430        l1 = l0-k1
            if(l1<=0) go to 450
            if(l1<=n11) go to 440
            l0 = l1-n11
            go to 430
 440        h1(l1) = b(it,j)*pinv
            go to 460
 450        h2(l0) = h2(l0)+b(it,j)*pinv
 460      continue
          if(n11<=0) go to 510
      !  rotate this row into triangle by givens transformations without
      !  square roots.
      !  rotation with the rows l,l+1,...n11.
 470      do 500 j=l,n11
            piv = h1(1)
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,g1(j,1),cos,sin)
      !  transformation to right hand side.
            call fprota(cos,sin,yi,c(j))
      !  transformation to the left hand side with respect to g2.
            do 480 i=1,k1
              call fprota(cos,sin,h2(i),g2(j,i))
 480        continue
            if(j==n11) go to 510
            i2 = min0(n11-j,k1)
      !  transformation to the left hand side with respect to g1.
            do 490 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),g1(j,i1))
              h1(i) = h1(i1)
 490        continue
            h1(i1) = 0.
 500      continue
      !  rotation with the rows n11+1,...n7
 510      do 530 j=1,k1
            ij = n11+j
            if(ij<=0) go to 530
            piv = h2(j)
      !  calculate the parameters of the givens transformation
            call fpgivs(piv,g2(ij,j),cos,sin)
      !  transformation to the right hand side.
            call fprota(cos,sin,yi,c(ij))
            if(j==k1) go to 540
            j1 = j+1
      !  transformation to the left hand side.
            do 520 i=j1,k1
              call fprota(cos,sin,h2(i),g2(ij,i))
 520        continue
 530      continue
 540    continue
      !  backward substitution to obtain the b-spline coefficients
      !  c(j),j=1,2,...n7 of sp(x).
        c(:n7) = fpbacp(g1,g2,c,n7,k1,k2,nest)
      !  calculate from condition (**) the b-spline coefficients c(n7+j),j=1,.
        do 545 i=1,k
          j = i+n7
          c(j) = c(i)
 545    continue
      !  computation of f(p).
        fp = 0.
        l = k1
        do 570 it=1,m1
          if(x(it)<t(l)) go to 550
          l = l+1
 550      l0 = l-k2
          term = 0.
          do 560 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 560      continue
          fp = fp+(w(it)*(term-y(it)))**2
 570    continue
      !  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 660
      !  test whether the maximal number of iterations is reached.
        if(iter==maxit) go to 600
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 580
        if((f2-f3) > acc) go to 575
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 +p2*con1
        go to 595
 575    if(f2<0.) ich3 = 1
 580    if(ich1/=0) go to 590
        if((f1-f2) > acc) go to 585
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 595
        if(p>=p3) p = p2*con1 +p3*con9
        go to 595
 585    if(f2>0.) ich1 = 1
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 590    if(f2>=f1 .or. f2<=f3) go to 610
      !  find the new value for p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 595  continue
      !  error codes and messages.
 600  ier = 3
      go to 660
 610  ier = 2
      go to 660
 620  ier = 1
      go to 660
 630  ier = -1
      go to 660
 640  ier = -2
      !  the least-squares constant function c1 is a solution of our problem.
      !  a constant function is a spline of degree k with all b-spline
      !  coefficients equal to that constant c1.
      do 650 i=1,k1
        rn = k1-i
        t(i) = x(1)-rn*per
        c(i) = c1
        j = i+k1
        rn = i-1
        t(j) = x(m)+rn*per
 650  continue
      n = nmin
      fp = fp0
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
 660  return
      end subroutine fpperi


      recursive subroutine fppocu(idim,k,a,b,ib,db,nb,ie,de,ne,cp,np)

      !  subroutine fppocu finds a idim-dimensional polynomial curve p(u) =
      !  (p1(u),p2(u),...,pidim(u)) of degree k, satisfying certain derivative
      !  constraints at the end points a and b, i.e.
      !                  (l)
      !    if ib > 0 : pj   (a) = db(idim*l+j), l=0,1,...,ib-1
      !                  (l)
      !    if ie > 0 : pj   (b) = de(idim*l+j), l=0,1,...,ie-1
      !
      !  the polynomial curve is returned in its b-spline representation
      !  ( coefficients cp(j), j=1,2,...,np )
      !  ..
      !  ..scalar arguments..
      integer idim,k,ib,nb,ie,ne,np
      real(RKIND) a,b
      !  ..array arguments..
      real(RKIND) db(nb),de(ne),cp(np)
      !  ..local scalars..
      real(RKIND) ab,aki
      integer i,id,j,jj,l,ll,k1,k2
      !  ..local array..
      real(RKIND) work(6,6)
      !  ..
      k1 = k+1
      k2 = 2*k1
      ab = b-a
      do 110 id=1,idim
        do 10 j=1,k1
          work(j,1) = 0.
  10    continue
        if(ib==0) go to 50
        l = id
        do 20 i=1,ib
          work(1,i) = db(l)
          l = l+idim
  20    continue
        if(ib==1) go to 50
        ll = ib
        do 40 j=2,ib
          ll =  ll-1
          do 30 i=1,ll
            aki = k1-i
            work(j,i) = ab*work(j-1,i+1)/aki + work(j-1,i)
  30      continue
  40    continue
  50    if(ie==0) go to 90
        l = id
        j = k1
        do 60 i=1,ie
          work(j,i) = de(l)
          l = l+idim
          j = j-1
  60    continue
        if(ie==1) go to 90
        ll = ie
        do 80 jj=2,ie
          ll =  ll-1
          j = k1+1-jj
          do 70 i=1,ll
            aki = k1-i
            work(j,i) = work(j+1,i) - ab*work(j,i+1)/aki
            j = j-1
  70      continue
  80    continue
  90    l = (id-1)*k2
        do 100 j=1,k1
          l = l+1
          cp(l) = work(j,1)
 100    continue
 110  continue
      return
      end subroutine fppocu


      recursive subroutine fppogr(iopt,ider,u,mu,v,mv,z,mz,z0,r,s, &
       nuest,nvest,tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu, &
       reducv,fpintu,fpintv,dz,step,lastdi,nplusu,nplusv,lasttu,nru, &
       nrv,nrdatu,nrdatv,wrk,lwrk,ier)

      !  ..
      !  ..scalar arguments..
      integer mu,mv,mz,nuest,nvest,maxit,nc,nu,nv,lastdi,nplusu,nplusv, &
       lasttu,lwrk,ier
      real(RKIND) z0,r,s,tol,fp,fp0,fpold,reducu,reducv,step
      !  ..array arguments..
      integer iopt(3),ider(2),nrdatu(nuest),nrdatv(nvest),nru(mu), &
       nrv(mv)
      real(RKIND) u(mu),v(mv),z(mz),tu(nuest),tv(nvest),c(nc),fpintu(nuest), &
       fpintv(nvest),dz(3),wrk(lwrk)
      !  ..local scalars..
      real(RKIND) acc,fpms,f1,f2,f3,p,per,pi,p1,p2,p3,vb,ve,zmax,zmin,rn
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,istart,iter,i1,i2,j,ju, &
       ktu,l,l1,l2,l3,l4,mpm,mumin,mu0,mu1,nn,nplu,nplv,npl1,nrintu, &
       nrintv,nue,numax,nve,nvmax
      !  ..local arrays..
      integer idd(2)
      real(RKIND) dzz(3)

      !   set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND
      !   initialization
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
      mumin = 4-iopt(3)
      if(ider(1)>=0) mumin = mumin-1
      if(iopt(2)==1 .and. ider(2)==1) mumin = mumin-1
      pi = datan2(0d0,-one)
      per = pi+pi
      vb = v(1)
      ve = vb+per

      !  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s

      !  numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+5+iopt(2)+iopt(3)
      nvmax = mv+7
      nue = min0(numax,nuest)
      nve = min0(nvmax,nvest)

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 1: determination of the number of knots and their position.     c
      ! ****************************************************************     c
      !  given a set of knots we compute the least-squares spline sinf(u,v)  c
      !  and the corresponding sum of squared residuals fp = f(p=inf).       c
      !  if iopt(1)=-1  sinf(u,v) is the requested approximation.            c
      !  if iopt(1)>=0  we check whether we can accept the knots:            c
      !    if fp <= s we will continue with the current set of knots.        c
      !    if fp >  s we will increase the number of knots and compute the   c
      !       corresponding least-squares spline until finally fp <= s.      c
      !    the initial choice of knots depends on the value of s and iopt.   c
      !    if s=0 we have spline interpolation; in that case the number of   c
      !     knots in the u-direction equals nu=numax=mu+5+iopt(2)+iopt(3)    c
      !     and in the v-direction nv=nvmax=mv+7.                            c
      !    if s>0 and                                                        c
      !      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  c
      !       spline without interior knots : nu=8 ; nv=8.                   c
      !      iopt(1)=1 we start with the set of knots found at the last call c
      !       of the routine, except for the case that s > fp0; then we      c
      !       compute the least-squares polynomial directly.                 c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(iopt(1)<0) go to 120

      if(s>0.) go to 100
      !  if s = 0, s(u,v) is an interpolating spline.
      nu = numax
      nv = nvmax
      !  test whether the required storage space exceeds the available one.
      if(nu>nuest .or. nv>nvest) go to 420
      !  find the position of the knots in the v-direction.
      do 10 l=1,mv
        tv(l+3) = v(l)
  10  continue
      tv(mv+4) = ve
      l1 = mv-2
      l2 = mv+5
      do 20 i=1,3
         tv(i) = v(l1)-per
         tv(l2) = v(i+1)+per
         l1 = l1+1
         l2 = l2+1
  20  continue
      !  if not all the derivative values g(i,j) are given, we will first
      !  estimate these values by computing a least-squares spline
      idd(1) = ider(1)
      if(idd(1)==0) idd(1) = 1
      if(idd(1)>0) dz(1) = z0
      idd(2) = ider(2)
      if(ider(1)<0) go to 30
      if(iopt(2)==0 .or. ider(2)/=0) go to 70
      ! we set up the knots in the u-direction for computing the least-squares
      ! spline.
  30  i1 = 3
      i2 = mu-2
      nu = 4
      do 40 i=1,mu
         if(i1>i2) go to 50
         nu = nu+1
         tu(nu) = u(i1)
         i1 = i1+2
  40  continue
  50  do 60 i=1,4
         tu(i) = 0.
         nu = nu+1
         tu(nu) = r
  60  continue
      ! we compute the least-squares spline for estimating the derivatives.
      call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd, &
        tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
        wrk,lwrk)
      ifsu = 0
      ! if all the derivatives at the origin are known, we compute the
      ! interpolating spline.
      ! we set up the knots in the u-direction, needed for interpolation.
  70  nn = numax-8
      if(nn==0) go to 95
      ju = 2-iopt(2)
      do 80 l=1,nn
        tu(l+4) = u(ju)
        ju = ju+1
  80  continue
      nu = numax
      l = nu
      do 90 i=1,4
         tu(i) = 0.
         tu(l) = r
         l = l-1
  90  continue
      ! we compute the interpolating spline.
  95  call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd, &
        tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
        wrk,lwrk)
      go to 430
      !  if s>0 our initial choice of knots depends on the value of iopt(1).
 100  ier = 0
      if(iopt(1)==0) go to 115
      step = -step
      if(fp0<=s) go to 115
      !  if iopt(1)=1 and fp0 > s we start computing the least-squares spline
      !  according to the set of knots found at the last call of the routine.
      !  we determine the number of grid coordinates u(i) inside each knot
      !  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mu0 = 2-iopt(2)
      mu1 = mu-2+iopt(3)
      do 105 i=mu0,mu1
        nrdatu(j) = nrdatu(j)+1
        if(u(i)<tu(l)) go to 105
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
 105  continue
      !  we determine the number of grid coordinates v(i) inside each knot
      !  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      do 110 i=2,mv
        nrdatv(j) = nrdatv(j)+1
        if(v(i)<tv(l)) go to 110
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
 110  continue
      idd(1) = ider(1)
      idd(2) = ider(2)
      go to 120
      !  if iopt(1)=0 or iopt(1)=1 and s >= fp0,we start computing the least-
      !  squares polynomial (which is a spline without interior knots).
 115  ier = -2
      idd(1) = ider(1)
      idd(2) = 1
      nu = 8
      nv = 8
      nrdatu(1) = mu-3+iopt(2)+iopt(3)
      nrdatv(1) = mv-1
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
      !  main loop for the different sets of knots.mpm=mu+mv is a save upper
      !  bound for the number of trials.
 120  mpm = mu+mv
      do 270 iter=1,mpm
      !  find nrintu (nrintv) which is the number of knot intervals in the
      !  u-direction (v-direction).
        nrintu = nu-7
        nrintv = nv-7
      !  find the position of the additional knots which are needed for the
      !  b-spline representation of s(u,v).
        i = nu
        do 130 j=1,4
          tu(j) = 0.
          tu(i) = r
          i = i-1
 130    continue
        l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 140 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-per
          tv(l3) = tv(l1)+per
 140    continue
      !  find an estimate of the range of possible values for the optimal
      !  derivatives at the origin.
        ktu = nrdatu(1)+2-iopt(2)
        if(nrintu==1) ktu = mu
        if(ktu<mumin) ktu = mumin
        if(ktu==lasttu) go to 150
         zmin = z0
         zmax = z0
         l = mv*ktu
         do 145 i=1,l
            if(z(i)<zmin) zmin = z(i)
            if(z(i)>zmax) zmax = z(i)
 145     continue
         step = zmax-zmin
         lasttu = ktu
      !  find the least-squares spline sinf(u,v).
 150    call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd, &
         tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
         wrk,lwrk)
        if(step<0.) step = -step
        if(ier==(-2)) fp0 = fp
      !  test whether the least-squares spline is an acceptable solution.
        if(iopt(1)<0) go to 440
        fpms = fp-s
        if(abs(fpms) < acc) go to 440
      !  if f(p=inf) < s, we accept the choice of knots.
        if(fpms<0.) go to 300
      !  if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
        if(nu==numax .and. nv==nvmax) go to 430
      !  increase the number of knots.
      !  if nu=nue and nv=nve we cannot further increase the number of knots
      !  because of the storage capacity limitation.
        if(nu==nue .and. nv==nve) go to 420
        if(ider(1)==0) fpintu(1) = fpintu(1)+(z0-c(1))**2
        ier = 0
      !  adjust the parameter reducu or reducv according to the direction
      !  in which the last added knots were located.
        if (lastdi<0) go to 160
        if (lastdi==0) go to 155
        go to 170
 155     nplv = 3
         idd(2) = ider(2)
         fpold = fp
         go to 230
 160    reducu = fpold-fp
        go to 175
 170    reducv = fpold-fp
      !  store the sum of squared residuals for the current set of knots.
 175    fpold = fp
      !  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu==8) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu>acc) npl1 = int(rn*fpms/reducu)
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
      !  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 3
        if(nv==8) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv>acc) npl1 = int(rn*fpms/reducv)
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
      !  test whether we are going to add knots in the u- or v-direction.
 190    if (nplu<nplv) go to 210
        if (nplu==nplv) go to 200
        go to 230
 200    if(lastdi<0) go to 230
 210    if(nu==nue) go to 230
      !  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        istart = 0
        if(iopt(2)==0) istart = 1
        do 220 l=1,nplusu
      !  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,istart)
      !  test whether we cannot further increase the number of knots in the
      !  u-direction.
          if(nu==nue) go to 270
 220    continue
        go to 270
 230    if(nv==nve) go to 210
      !  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
      !  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
      !  test whether we cannot further increase the number of knots in the
      !  v-direction.
          if(nv==nve) go to 270
 240    continue
      !  restart the computations with the new set of knots.
 270  continue
      !  test whether the least-squares polynomial is a solution of our
      !  approximation problem.
 300  if(ier==(-2)) go to 440
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 2: determination of the smoothing spline sp(u,v)                c
      ! *****************************************************                c
      !  we have determined the number of knots and their position. we now   c
      !  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
      !  this smoothing spline depends on the parameter p in such a way that c
      !    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             c
      !  is a continuous, strictly decreasing function of p. moreover the    c
      !  least-squares polynomial corresponds to p=0 and the least-squares   c
      !  spline to p=infinity. then iteratively we have to determine the     c
      !  positive value of p such that f(p)=s. the process which is proposed c
      !  here makes use of rational interpolation. f(p) is approximated by a c
      !  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
      !  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
      !  are used to calculate the new value of p such that r(p)=s.          c
      !  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      dzz(1) = dz(1)
      dzz(2) = dz(2)
      dzz(3) = dz(3)
      ich1 = 0
      ich3 = 0
      !  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
      !  find the smoothing spline sp(u,v) and the corresponding sum f(p).
        call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dzz,iopt,idd, &
         tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
         wrk,lwrk)
      !  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 440
      !  test whether the maximum allowable number of iterations has been
      !  reached.
        if(iter==maxit) go to 400
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 320
        if((f2-f3)>acc) go to 310
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2<zero) ich3 = 1
 320    if(ich1/=0) go to 340
        if((f1-f2)>acc) go to 330
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 350
        if(p>=p3) p = p2*con1 + p3*con9
        go to 350
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 330    if(f2>0.) ich1 = 1
 340    if(f2>=f1 .or. f2<=f3) go to 410
      !  find the new value of p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 350  continue
      !  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end subroutine fppogr


      recursive subroutine fppola(iopt1,iopt2,iopt3,m,u,v,z,w,rad,s, &
       nuest,nvest,eta,tol,maxit,ib1,ib3,nc,ncc,intest,nrest,nu,tu,nv, &
       tv,c,fp,sup,fpint,coord,f,ff,row,cs,cosi,a,q,bu,bv,spu,spv,h, &
       index,nummer,wrk,lwrk,ier)

      !  ..scalar arguments..
      integer iopt1,iopt2,iopt3,m,nuest,nvest,maxit,ib1,ib3,nc,ncc, &
       intest,nrest,lwrk,ier
      integer, intent(inout) :: nu,nv
      real(RKIND) s,eta,tol,fp,sup
      !  ..array arguments..
      integer :: index(nrest),nummer(m)
      real(RKIND) :: u(m),v(m),z(m),w(m),tu(nuest),tv(nvest),c(nc),fpint(intest),coord(intest),f(ncc),ff(nc),row(nvest), &
                     cs(nvest),cosi(5,nvest),a(ncc,ib1),q(ncc,ib3),bu(nuest,5),bv(nvest,5),spu(m,4),spv(m,4),h(ib3),wrk(lwrk)
      !  ..user supplied function..
      real(RKIND) :: rad
      !  ..local scalars..
      real(RKIND) :: acc,arg,co,c1,c2,c3,c4,dmax,eps,fac,fac1,fac2,fpmax,fpms,f1,f2,f3,hui,huj,p,pinv,piv,p1,p2,p3, &
                     r,ratio,si,sigma,sq,store,uu,u2,u3,wi,zi,rn
      integer :: i,iband,iband3,iband4,ich1,ich3,ii,il,in,ipar,ipar1,irot,iter,i1,i2,j,jrot,j1,j2,l,la,lf,lh,ll,&
                 lu,lv,lwest,l1,l2,l3,l4,ncof,ncoff,nvv,nv4,nreg,nrint,nrr,nr1,nuu,nu4,num,num1,numin,nvmin,rank,iband1,jlu
      !  ..local arrays..
      real(RKIND), dimension(SIZ_K+1) :: hu,hv

      !  set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND

      fpms = zero
      ipar = iopt2*(iopt2+3)/2
      ipar1 = ipar+1
      eps = sqrt(eta)
      iband1 = 0

      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s

      if(iopt1<0) go to 90
      numin = 9
      nvmin = 9+iopt2*(iopt2+1)

      if(iopt1==0) go to 10
      if(s<sup) then
        if (nv<nvmin) go to 70
        go to 90
      endif
      !  if iopt1 = 0 we begin by computing the weighted least-squares
      !  polymomial of the form
      !     s(u,v) = f(1)*(1-u**3)+f(2)*u**3+f(3)*(u**2-u**3)+f(4)*(u-u**3)
      !  where f(4) = 0 if iopt2> 0 , f(3) = 0 if iopt2 > 1 and
      !        f(2) = 0 if iopt3> 0.
      !  the corresponding weighted sum of squared residuals gives the upper
      !  bound sup for the smoothing factor s.
  10  sup = zero
      f(1:4) = zero
      a(1:4,1:4) = zero
      do 50 i=1,m
         wi = w(i)
         zi = z(i)*wi
         uu = u(i)
         u2 = uu*uu
         u3 = uu*u2
         h(1) = (one-u3)*wi
         h(2) = u3*wi
         h(3) = u2*(one-uu)*wi
         h(4) = uu*(one-u2)*wi
         if(iopt3/=0) h(2) = zero
         if(iopt2>1) h(3) = zero
         if(iopt2>0) h(4) = zero
         do 40 j=1,4
            piv = h(j)
            if (piv==zero) go to 40
            call fpgivs(piv,a(j,1),co,si)
            call fprota(co,si,zi,f(j))
            if(j==4) go to 40
            j1 = j+1
            j2 = 1
            do 30 l=j1,4
               j2 = j2+1
               call fprota(co,si,h(l),a(j,j2))
  30        continue
  40     continue
         sup = sup+zi*zi
  50  continue
      if(a(4,1)/=zero) f(4) = f(4)/a(4,1)
      if(a(3,1)/=zero) f(3) = (f(3)-a(3,2)*f(4))/a(3,1)
      if(a(2,1)/=zero) f(2) = (f(2)-a(2,2)*f(3)-a(2,3)*f(4))/a(2,1)
      if(a(1,1)/=zero) f(1) = (f(1)-a(1,2)*f(2)-a(1,3)*f(3)-a(1,4)*f(4))/a(1,1)

      !  find the b-spline representation of this least-squares polynomial
      c1 = f(1)
      c4 = f(2)
      c2 = f(4)/three+c1
      c3 = (f(3)+two*f(4))/three+c1
      nu = 8
      nv = 8
      do 60 i=1,4
         c(i)    = c1
         c(i+4)  = c2
         c(i+8)  = c3
         c(i+12) = c4
         tu(i)   = zero
         tu(i+4) = one
         rn      = 2*i-9
         tv(i)   = rn*pi
         rn      = 2*i-1
         tv(i+4) = rn*pi
  60  continue
      fp = sup
      !  test whether the least-squares polynomial is an acceptable solution
      fpms = sup-s
      if(fpms<acc) go to 960
      !  test whether we cannot further increase the number of knots.
  70  if(nuest<numin .or. nvest<nvmin) go to 950
      !  find the initial set of interior knots of the spline in case iopt1=0.
      nu = numin
      nv = nvmin
      tu(5) = half
      nvv = nv-8
      rn = nvv+1
      fac = pi2/rn
      forall (i=1:nvv) tv(i+4) = i*fac-pi

      !  ************************************************************************************************************
      !  part 1 : computation of least-squares bicubic splines.
      !  ************************************************************************************************************
      !  if iopt1<0 we compute the least-squares bicubic spline according to the given set of knots.
      !  if iopt1>=0 we compute least-squares bicubic splines with increasing numbers of knots until the
      !  corresponding sum f(p=inf)<=s.
      !  the initial set of knots then depends on the value of iopt1:
      !    if iopt1=0 we start with one interior knot in the u-direction (0.5) and 1+iopt2*(iopt2+1) in the
      !               v-direction.
      !    if iopt1>0 we start with the set of knots found at the last call of the routine.
      !  ************************************************************************************************************
      !  main loop for the different sets of knots. m is a save upper bound for the number of trials.
  90  do 570 iter=1,m
      !  find the position of the additional knots which are needed for the
      !  b-spline representation of s(u,v).
         l1 = 4
         l2 = l1
         l3 = nv-3
         l4 = l3
         tv(l2) = -pi
         tv(l3) = pi
         do 120 i=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tv(l2) = tv(l4)-pi2
            tv(l3) = tv(l1)+pi2
 120     continue
        l = nu
        do 130 i=1,4
          tu(i) = zero
          tu(l) = one
          l = l-1
 130    continue
      !  find nrint, the total number of knot intervals and nreg, the number
      !  of panels in which the approximation domain is subdivided by the
      !  intersection of knots.
        nuu = nu-7
        nvv = nv-7
        nrr = nvv/2
        nr1 = nrr+1
        nrint = nuu+nvv
        nreg = nuu*nvv
      !  arrange the data points according to the panel they belong to.
        call fporde(u,v,m,3,3,tu,nu,tv,nv,nummer,index,nreg)
        if(iopt2==0) go to 195
      !  find the b-spline coefficients cosi of the cubic spline
      !  approximations for cr(v)=rad(v)*cos(v) and sr(v) = rad(v)*sin(v)
      !  if iopt2=1, and additionally also for cr(v)**2,sr(v)**2 and
      !  2*cr(v)*sr(v) if iopt2=2
        a(1:nvv,1:nvv) = zero
        cosi(1:ipar,1:nvv) = zero
      !  the coefficients cosi are obtained from interpolation conditions
      !  at the knots tv(i),i=4,5,...nv-4.
        do 175 i=1,nvv
           l2 = i+3
           arg = tv(l2)
           call fpbspl(tv,nv,3,arg,l2,hv)
           do 145 j=1,nvv
              row(j) = 0.
 145       continue
           ll = i
           do 150 j=1,3
              if(ll>nvv) ll= 1
              row(ll) = row(ll)+hv(j)
              ll = ll+1
 150       continue
           co = cos(arg)
           si = sin(arg)
           r = rad(arg)
           cs(1) = co*r
           cs(2) = si*r
           if(iopt2==1) go to 155
           cs(3) = cs(1)*cs(1)
           cs(4) = cs(2)*cs(2)
           cs(5) = cs(1)*cs(2)
 155       do 170 j=1,nvv
              piv = row(j)
              if (piv==zero) go to 170
              call fpgivs(piv,a(j,1),co,si)
              do 160 l=1,ipar
                 call fprota(co,si,cs(l),cosi(l,j))
 160          continue
              if(j==nvv) go to 175
              j1 = j+1
              j2 = 1
              do 165 l=j1,nvv
                 j2 = j2+1
                 call fprota(co,si,row(l),a(j,j2))
 165          continue
 170       continue
 175    continue
         do 190 l=1,ipar
            cs(1:nvv) = cosi(l,1:nvv)
            cosi(l,1:nvv) = fpback(a,cs,nvv,nvv,ncc)
 190     continue
      !  find ncof, the dimension of the spline and ncoff, the number
      !  of coefficients in the standard b-spline representation.
 195    nu4 = nu-4
        nv4 = nv-4
        ncoff = nu4*nv4
        ncof = ipar1+nvv*(nu4-1-iopt2-iopt3)
      !  find the bandwidth of the observation matrix a.
        iband = 4*nvv
        if(nuu-iopt2-iopt3<=1) iband = ncof
        iband1 = iband-1
      !  initialize the observation matrix a.
        f(1:ncof) = zero
        a(1:ncof,1:iband) = zero
      !  initialize the sum of squared residuals.
        fp = zero
        ratio = one+tu(6)/tu(5)
      !  fetch the data points in the new order. main loop for the different panels.
        panels: do num=1,nreg
           !  fix certain constants for the current panel; jrot records the column number of the first
           ! non-zero element in a row of the observation matrix according to a data point of the panel.
           num1 = num-1
           lu = num1/nvv
           l1 = lu+4
           lv = num1-lu*nvv+1
           l2 = lv+3
           jrot = 0
           if(lu>iopt2) jrot = ipar1+(lu-iopt2-1)*nvv
           lu = lu+1

           !  test whether there are still data points in the current panel.
           in = index(num)
           points_left: do while (in/=0)

              ! fetch a new data point.
              wi = w(in)
              zi = z(in)*wi

              ! evaluate for the u-direction, the 4 non-zero b-splines at u(in)
              call fpbspl(tu,nu,3,u(in),l1,hu)

              ! evaluate for the v-direction, the 4 non-zero b-splines at v(in)
              call fpbspl(tv,nv,3,v(in),l2,hv)

              ! store the value of these b-splines in spu and spv resp.
              spu(in,:) = hu(1:4)
              spv(in,:) = hv(1:4)

              ! initialize the new row of observation matrix.
              h(1:iband) = zero

              ! calculate the non-zero elements of the new row by making the cross
              ! products of the non-zero b-splines in u- and v-direction and
              ! by taking into account the conditions of the splines.
              row(1:nvv) = zero

              ! take into account the periodicity condition of the bicubic splines.
              ll = lv
              do i=1,4
                 if (ll>nvv) ll=1
                 row(ll) = row(ll)+hv(i)
                 ll = ll+1
              end do

              ! take into account the other conditions of the splines.
              if (iopt2/=0 .and. lu<=iopt2+1) &
              cs(1:ipar) = matmul(cosi(1:ipar,1:nvv),row(1:nvv))

              ! fill in the non-zero elements of the new row.
              j1 = 0
              new_row: do j =1,4
                jlu = j+lu
                huj = hu(j)
                if (jlu>iopt2+2) then
                    if (jlu>nu4 .and. iopt3/=0) cycle new_row
                    h(j1+1:j1+nvv) = row(1:nvv)*huj
                    j1 = j1+nvv
                elseif (jlu==1 .or. jlu==2) then
                    h(1) = huj
                    j1 = 1
                elseif (jlu==3) then
                    h(1:3) = [h(1)+huj,huj*cs(1:2)]
                    j1 = 3
                elseif (jlu==4) then
                    h(1)   = h(1)+huj
                    h(2:3) = h(2:3)+huj*ratio*cs(1:2)
                    h(4:6) = huj*cs(3:5)
                    j1 = 6
                endif
              end do new_row
              h(1:iband) = wi*h(1:iband)

              ! rotate the row into triangle by givens transformations.
              irot = jrot
              rotate: do i=1,iband
                irot = irot+1
                piv  = h(i)
                if (piv==zero) cycle rotate

                ! calculate the parameters of the givens transformation.
                call fpgivs(piv,a(irot,1),co,si)

                ! apply that transformation to the right hand side.
                call fprota(co,si,zi,f(irot))

                if (i==iband) exit rotate

                ! apply that transformation to the left hand side.
                i2 = 1
                do j=i+1,iband
                  i2 = i2+1
                  call fprota(co,si,h(j),a(irot,i2))
                end do
              end do rotate

              ! add the contribution of the row to the sum of squares of residual right hand sides.
              fp = fp+zi**2

              ! find the number of the next data point in the panel.
              in = nummer(in)
           end do points_left
        end do panels

        ! find dmax, the maximum value for the diagonal elements in the reduced triangle.
        dmax = max(zero,maxval(a(1:ncof,1)))

        ! check whether the observation matrix is rank deficient.
        sigma = eps*dmax

        if (all(a(i,1:ncof)>sigma)) then

           ! backward substitution in case of full rank.
           c(:ncof) = fpback(a,f,ncof,iband,ncc)
           rank = ncof
           q(1:ncof,1) = a(1:ncof,1)/dmax

        else

           ! in case of rank deficiency, find the minimum norm solution.
           lwest = ncof*iband+ncof+iband
           if(lwest>lwrk) go to 925

           lf = 1
           lh = lf+ncof
           la = lh+iband
           ff(1:ncof) = f(1:ncof)
           q(1:ncof,1:iband) = a(1:ncof,1:iband)
           call fprank(q,ff,ncof,iband,ncc,sigma,c,sq,rank,wrk(la),wrk(lf),wrk(lh))
           q(1:ncof,1) = q(1:ncof,1)/dmax

           ! add to the sum of squared residuals, the contribution of reducing the rank.
           fp = fp+sq

        endif

        ! find the coefficients in the standard b-spline representation of the spline.
        call fprppo(nu,nv,iopt2,iopt3,cosi,ratio,c,ff,ncoff)

        ! test whether the least-squares spline is an acceptable solution.
        if(iopt1<0) then
          if (fp<=0) go to 970
          go to 980
        endif
        fpms = fp-s
        if(abs(fpms)<=acc) then
            if (fp<=0) go to 970
            go to 980
        endif
      !  if f(p=inf) < s, accept the choice of knots.
        if(fpms<zero) go to 580
      !  test whether we cannot further increase the number of knots
        if(m<ncof) go to 935
      !  search where to add a new knot.
      !  find for each interval the sum of squared residuals fpint for the
      !  data points having the coordinate belonging to that knot interval.
      !  calculate also coord which is the same sum, weighted by the position
      !  of the data points considered.
        fpint(1:nrint) = zero
        coord(1:nrint) = zero

        do 490 num=1,nreg
          num1 = num-1
          lu = num1/nvv
          l1 = lu+1
          lv = num1-lu*nvv
          l2 = lv+1+nuu
          jrot = lu*nv4+lv
          in = index(num)
 460      if(in==0) go to 490
          store = zero
          i1 = jrot
          do 480 i=1,4
            hui = spu(in,i)
            j1 = i1
            do 470 j=1,4
              j1 = j1+1
              store = store+hui*spv(in,j)*c(j1)
 470        continue
            i1 = i1+nv4
 480      continue
          store = (w(in)*(z(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*u(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*v(in)
          in = nummer(in)
          go to 460
 490    continue
      ! bring together the information concerning knot panels which are
      ! symmetric with respect to the origin.
        do 495 i=1,nrr
          l1 = nuu+i
          l2 = l1+nrr
          fpint(l1) = fpint(l1)+fpint(l2)
          coord(l1) = coord(l1)+coord(l2)-pi*fpint(l2)
 495    continue
      !  find the interval for which fpint is maximal on the condition that
      !  there still can be added a knot.
        l1 = 1
        l2 = nuu+nrr
        if(nuest<nu+1) l1=nuu+1
        if(nvest<nv+2) l2=nuu
      !  test whether we cannot further increase the number of knots.
        if(l1>l2) go to 950
 500    fpmax = 0.
        l = 0
        do 510 i=l1,l2
          if(fpmax>=fpint(i)) go to 510
          l = i
          fpmax = fpint(i)
 510    continue
        if(l==0) go to 930
      !  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
      !  test in what direction the new knot is going to be added.
        if(l>nuu) go to 530
      !  addition in the u-direction
        l4 = l+4
        fpint(l) = zero
        fac1 = tu(l4)-arg
        fac2 = arg-tu(l4-1)
        if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) go to 500
        j = nu
        do 520 i=l4,nu
          tu(j+1) = tu(j)
          j = j-1
 520    continue
        tu(l4) = arg
        nu = nu+1
        go to 570
      !  addition in the v-direction
 530    l4 = l+4-nuu
        fpint(l) = zero
        fac1 = tv(l4)-arg
        fac2 = arg-tv(l4-1)
        if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) go to 500
        ll = nrr+4
        j = ll
        do 550 i=l4,ll
          tv(j+1) = tv(j)
          j = j-1
 550    continue
        tv(l4) = arg
        nv = nv+2
        nrr = nrr+1
        do 560 i=5,ll
          j = i+nrr
          tv(j) = tv(i)+pi
 560    continue
      !  restart the computations with the new set of knots.
 570  continue
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 2: determination of the smoothing bicubic spline.               c
      ! ******************************************************               c
      ! we have determined the number of knots and their position. we now    c
      ! compute the coefficients of the smoothing spline sp(u,v).            c
      ! the observation matrix a is extended by the rows of a matrix, expres-c
      ! sing that sp(u,v) must be a constant function in the variable        c
      ! v and a cubic polynomial in the variable u. the corresponding        c
      ! weights of these additional rows are set to 1/(p). iteratively       c
      ! we than have to determine the value of p such that f(p) = sum((w(i)* c
      ! (z(i)-sp(u(i),v(i))))**2)  be = s.                                   c
      ! we already know that the least-squares polynomial corresponds to p=0,c
      ! and that the least-squares bicubic spline corresponds to p=infin.    c
      ! the iteration process makes use of rational interpolation. since f(p)c
      ! is a convex and strictly decreasing function of p, it can be approx- c
      ! imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c
      ! three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c
      ! f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c
      ! of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<zeroc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  evaluate the discontinuity jumps of the 3-th order derivative of
      !  the b-splines at the knots tu(l),l=5,...,nu-4.
 580  call fpdisc(tu,nu,5,bu,nuest)
      !  evaluate the discontinuity jumps of the 3-th order derivative of
      !  the b-splines at the knots tv(l),l=5,...,nv-4.
      call fpdisc(tv,nv,5,bv,nvest)
      !  initial value for p.
      p1 = 0.
      f1 = sup-s
      p3 = -one
      f3 = fpms
      p = zero
      do 590 i=1,ncof
        p = p+a(i,1)
 590  continue
      rn = ncof
      p = rn/p
      !  find the bandwidth of the extended observation matrix.
      iband4 = iband+ipar1
      if(iband4>ncof) iband4 = ncof
      iband3 = iband4 -1
      ich1 = 0
      ich3 = 0
      nuu = nu4-iopt3-1
      !  iteration process to find the root of f(p)=s.
      do 920 iter=1,maxit
        pinv = one/p
      !  store the triangularized observation matrix into q.
        ff(1:ncof) = f(1:ncof)
        q(1:ncof,1:iband4) = zero
        q(1:ncof,1:iband)  = a(1:ncof,1:iband)

      !  extend the observation matrix with the rows of a matrix, expressing
      !  that for u=constant sp(u,v) must be a constant function.
        do 720 i=5,nv4
          ii = i-4
          do 635 l=1,nvv
             row(l) = 0.
 635      continue
          ll = ii
          do 640  l=1,5
             if(ll>nvv) ll=1
             row(ll) = row(ll)+bv(ii,l)
             ll = ll+1
 640      continue
          do 721 j=1,nuu
      !  initialize the new row.
            h(1:iband) = zero
      !  fill in the non-zero elements of the row. jrot records the column
      !  number of the first non-zero element in the row.
            if(j>iopt2) go to 665
            if(j==2) go to 655
            cs(1:2) = matmul(cosi(1:2,1:nvv),row(1:nvv))
            h(1) = cs(1)
            h(2) = cs(2)
            jrot = 2
            go to 675
 655        cs(3:5) = matmul(cosi(3:5,1:nvv),row(1:nvv))
            h(1) = cs(1)*ratio
            h(2) = cs(2)*ratio
            h(3) = cs(3)
            h(4) = cs(4)
            h(5) = cs(5)
            jrot = 2
            go to 675
 665        h(1:nvv) = row(1:nvv)
             jrot = ipar1+1+(j-iopt2-1)*nvv
 675        h(1:iband) = h(1:iband)*pinv
            zi = 0.
      !  rotate the new row into triangle by givens transformations.
            do 710 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if (piv==zero) then
                 if (i2<=0) go to 721
                 go to 690
              endif
      !  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
      !  apply that givens transformation to the right hand side.
              call fprota(co,si,zi,ff(irot))
              if(i2==0) go to 721
      !  apply that givens transformation to the left hand side.
              do 680 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 680          continue
 690          do 700 l=1,i2
                h(l) = h(l+1)
 700          continue
              h(i2+1) = zero
 710        continue
 721      continue
 720    continue
      !  extend the observation matrix with the rows of a matrix expressing
      !  that for v=constant. sp(u,v) must be a cubic polynomial.
        do 810 i=5,nu4
          ii = i-4
          do 811 j=1,nvv
      !  initialize the new row
            h(1:iband4) = zero
      !  fill in the non-zero elements of the row. jrot records the column
      !  number of the first non-zero element in the row.
            j1 = 1
            do 760 l=1,5
               il = ii+l-1
               if(il==nu4 .and. iopt3/=0) go to 760
               if(il>iopt2+1) go to 750
               go to (735,740,745),il
 735           h(1) = bu(ii,l)
               j1 = j+1
               go to 760
 740           h(1) = h(1)+bu(ii,l)
               h(2) = bu(ii,l)*cosi(1,j)
               h(3) = bu(ii,l)*cosi(2,j)
               j1 = j+3
               go to 760
 745           h(1) = h(1)+bu(ii,l)
               h(2) = bu(ii,l)*cosi(1,j)*ratio
               h(3) = bu(ii,l)*cosi(2,j)*ratio
               h(4) = bu(ii,l)*cosi(3,j)
               h(5) = bu(ii,l)*cosi(4,j)
               h(6) = bu(ii,l)*cosi(5,j)
               j1 = j+6
               go to 760
 750           h(j1) = bu(ii,l)
               j1 = j1+nvv
 760        continue
            do 765 l=1,iband4
              h(l) = h(l)*pinv
 765        continue
            zi = 0.
            jrot = 1
            if(ii>iopt2+1) jrot = ipar1+(ii-iopt2-2)*nvv+j
      !  rotate the new row into triangle by givens transformations.
            do 800 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if (piv==zero) then
                if (i2<=0) go to 811
                go to 780
              endif
      !  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
      !  apply that givens transformation to the right hand side.
              call fprota(co,si,zi,ff(irot))
              if(i2==0) go to 811
      !  apply that givens transformation to the left hand side.
              do 770 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 770          continue
 780          do 790 l=1,i2
                h(l) = h(l+1)
 790          continue
              h(i2+1) = zero
 800        continue
 811      continue
 810    continue
      !  find dmax, the maximum value for the diagonal elements in the
      !  reduced triangle.
        dmax = 0.
        do 820 i=1,ncof
          if(q(i,1)<=dmax) go to 820
          dmax = q(i,1)
 820    continue
      !  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 830 i=1,ncof
          if(q(i,1)<=sigma) go to 840
 830    continue
      !  backward substitution in case of full rank.
        c(:ncof) = fpback(q,ff,ncof,iband4,ncc)
        rank = ncof
        go to 845
      !  in case of rank deficiency, find the minimum norm solution.
 840    lwest = ncof*iband4+ncof+iband4
        if(lwrk<lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,ncc,sigma,c,sq,rank,wrk(la), &
         wrk(lf),wrk(lh))
 845    do 850 i=1,ncof
           q(i,1) = q(i,1)/dmax
 850    continue
      !  find the coefficients in the standard b-spline representation of
      !  the polar spline.
        call fprppo(nu,nv,iopt2,iopt3,cosi,ratio,c,ff,ncoff)
      !  compute f(p).
        fp = 0.
        do 890 num = 1,nreg
          num1 = num-1
          lu = num1/nvv
          lv = num1-lu*nvv
          jrot = lu*nv4+lv
          in = index(num)
 860      if(in==0) go to 890
          store = zero
          i1 = jrot
          do 880 i=1,4
            hui = spu(in,i)
            j1 = i1
            do 870 j=1,4
              j1 = j1+1
              store = store+hui*spv(in,j)*c(j1)
 870        continue
            i1 = i1+nv4
 880      continue
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
          go to 860
 890    continue
      !  test whether the approximation sp(u,v) is an acceptable solution
        fpms = fp-s
        if(abs(fpms)<=acc) go to 980
      !  test whether the maximum allowable number of iterations has been
      !  reached.
        if(iter==maxit) go to 940
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 900
        if((f2-f3)>acc) go to 895
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 + p2*con1
        go to 920
 895    if(f2<0.) ich3 = 1
 900    if(ich1/=0) go to 910
        if((f1-f2)>acc) go to 905
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 920
        if(p>=p3) p = p2*con1 +p3*con9
        go to 920
 905    if(f2>0.) ich1 = 1
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 910    if(f2>=f1 .or. f2<=f3) go to 945
      !  find the new value of p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 920  continue
      !  error codes and messages.
 925  ier = lwest
      go to 990
 930  ier = 5
      go to 990
 935  ier = 4
      go to 990
 940  ier = FITPACK_MAXIT
      go to 990
 945  ier = FITPACK_S_TOO_SMALL
      go to 990
 950  ier = FITPACK_INSUFFICIENT_STORAGE
      go to 990
 960  ier = FITPACK_LEASTSQUARES_OK
      go to 990
 970  ier = FITPACK_INTERPOLATING_OK
      fp  = zero
 980  if (ncof/=rank) ier = -rank
 990  return
      end subroutine fppola


      ! subroutine fprank finds the minimum norm solution of a leastsquares problem in case of
      ! rank deficiency.
      pure subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
      !  ..scalar arguments..
      integer,     intent(in)  :: n         ! the dimension of a.
      integer,     intent(in)  :: m         ! the bandwidth of a.
      integer,     intent(in)  :: na
      integer,     intent(out) :: rank      ! the rank of matrix a.
      real(RKIND), intent(in)  :: tol       ! a threshold to determine the rank of a.
      real(RKIND), intent(out) :: sq        ! the contribution of reducing the rank to the sum of squared residuals.

      !  ..array arguments..
      real(RKIND), intent(out) :: c(n)      ! the minimum norm solution.
      real(RKIND), intent(out) :: a(na,m)   ! the non-zero elements of the observation matrix after
                                            ! triangularization by givens transformations.
      real(RKIND), intent(out)  :: f(n)     ! the transformed right hand side.

      real(RKIND), intent(inout) :: aa(n,m),ff(n),h(m) ! working arrays
      !  ..local scalars..
      integer     :: i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
      real(RKIND) :: cos,fac,piv,sin,yi
      real(RKIND) :: store,stor1,stor2,stor3

      !  ..
      m1 = m-1
      !  the rank deficiency nl is considered to be the number of sufficient small diagonal elements of a.
      nl = 0
      sq = zero
      rows: do i=1,n
         if (a(i,1)>tol) cycle rows

         ! if a sufficient small diagonal element is found, we put it to zero. the remainder of the row
         ! corresponding to that zero diagonal element is then rotated into triangle by givens rotations .
         ! the rank deficiency is increased by one.
         nl = nl+1
         if (i==n) cycle rows
         yi     = f(i)
         h(1:m) = [a(i,2:m),zero]
         i1 = i+1
         in_row: do ii=i1,n
            i2 = min(n-ii,m1)
            piv = h(1)
            if (piv==zero) then
               if (i2==0) exit in_row
               h(1:i2) = h(2:i2+1)
            else
               call fpgivs(piv,a(ii,1),cos,sin)
               call fprota(cos,sin,yi,f(ii))
               if (i2==0) exit in_row
               do j=1,i2
                 j1 = j+1
                 call fprota(cos,sin,h(j1),a(ii,j1))
                 h(j) = h(j1)
               end do
            endif
            h(i2+1) = zero
         end do in_row
         ! add to the sum of squared residuals the contribution of deleting
         ! the row with small diagonal element.
         sq = sq+yi**2
      end do rows

      !  rank denotes the rank of a.
      rank = n-nl

      !  let b denote the (rank*n) upper trapezoidal matrix which can be obtained from the (n*n) upper
      !  triangular matrix a by deleting the rows and interchanging the columns corresponding to a zero
      !  diagonal element. if this matrix is factorized using givens transformations as  b = (r) (u)  where
      !    r is a (rank*rank) upper triangular matrix,
      !    u is a (rank*n) orthonormal matrix
      !  then the minimal least-squares solution c is given by c = b' v, where v is the solution of the
      !  system  (r) (r)' v = g  and g denotes the vector obtained from the old right hand side f, by
      !  removing the elements corresponding to a zero diagonal element of a.

      ! initialization.
      aa(1:rank,1:m) = zero

      !  form in aa the upper triangular matrix obtained from a by removing rows and columns with zero
      !  diagonal elements. form in ff the new right hand side by removing the elements of the old right
      !  hand side corresponding to a deleted row.
      ii = 0
      make_U: do i=1,n
         if (a(i,1)<=tol) cycle make_U
         ii = ii+1
         ff(ii) = f(i)
         aa(ii,1) = a(i,1)
         jj = ii
         kk = 1
         j  = i
         j1 = min(j-1,m1)
         if (j1/=0) then  ! probably unnecessary
            do k=1,j1
              j = j-1
              if (a(j,1)>tol) then
                 kk = kk+1
                 jj = jj-1
                 aa(jj,kk) = a(j,k+1)
              endif
            end do
         endif
      end do make_U

      !  form successively in h the columns of a with a zero diagonal element.
      ii = 0
      make_h: do i=1,n
        ii = ii+1
        if (a(i,1)>tol) cycle make_h
        ii = ii-1
        if (ii==0) cycle make_h
        jj = 1
        j  = i
        j1 = min(j-1,m1)
        do k=1,j1
           j = j-1
           if (a(j,1)>tol) then
             h(jj) = a(j,k+1)
             jj = jj+1
           endif
        end do
        h(jj:m) = zero

        ! rotate this column into aa by givens transformations.
        jj = ii
        rotate_col: do i1=1,ii
           j1 = min(jj-1,m1)
           piv = h(1)
           if (piv/=zero) then
              call fpgivs(piv,aa(jj,1),cos,sin)
              if (j1==0) cycle make_h
              kk = jj
              do j2=1,j1
                 j3 = j2+1
                 kk = kk-1
                 call fprota(cos,sin,h(j3),aa(kk,j3))
                 h(j2) = h(j3)
              end do
           else
              if (j1==0) cycle make_h
              do j2=1,j1
                j3 = j2+1
                h(j2) = h(j3)
              end do
           endif
           jj = jj-1
           h(j3) = zero
        end do rotate_col
      end do make_h

      !  solve the system (aa) (f1) = ff
      ff(rank) = ff(rank)/aa(rank,1)
      i = rank-1
      if (i/=0) then
         do j=2,rank
            store = ff(i)
            i1 = min0(j-1,m1)
            k = i
            do ii=1,i1
               k = k+1
               stor1 = ff(k)
               stor2 = aa(i,ii+1)
               store = store-stor1*stor2
            end do
            stor1 = aa(i,1)
            ff(i) = store/stor1
            i = i-1
         end do
      endif

      !  solve the system  (aa)' (f2) = f1
      ff(1) = ff(1)/aa(1,1)
      if (rank>1) then
          do j=2,rank
            store = ff(j)
            i1 = min0(j-1,m1)
            k = j
            do ii=1,i1
               k = k-1
               stor1 = ff(k)
               stor2 = aa(k,ii+1)
               store = store-stor1*stor2
            end do
            stor1 = aa(j,1)
            ff(j) = store/stor1
         end do
      endif

      !  premultiply f2 by the transpoze of a.
      k = 0
      do i=1,n
        store = zero
        if (a(i,1)>tol) k = k+1
        j1 = min(i,m)
        kk = k
        ij = i+1
        do j=1,j1
           ij = ij-1
           if (a(ij,1)>tol) then
              stor1 = a(ij,j)
              stor2 = ff(kk)
              store = store+stor1*stor2
              kk = kk-1
           endif
        end do
        c(i) = store
      end do

      !  add to the sum of squared residuals the contribution of putting
      !  to zero the small diagonal elements of matrix (a).
      stor3 = zero
      add_to_res: do i=1,n
         if (a(i,1)<=tol) then
            store = f(i)
            i1 = min(n-i,m1)
            if (i1>0) then
               do j=1,i1
                  ij = i+j
                  stor1 = c(ij)
                  stor2 = a(i,j+1)
                  store = store-stor1*stor2
               end do
            endif
            fac = a(i,1)*c(i)
            stor1 = a(i,1)
            stor2 = c(i)
            stor1 = stor1*stor2
            stor3 = stor3+stor1*(stor1-store-store)
         endif
      end do add_to_res

      fac = stor3
      sq  = sq+fac
      return
      end subroutine fprank


      !  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati  gives the value of p such
      !  that the rational interpolating function of the form r(p) = (u*p+v)/(p+w) equals zero at p.
      elemental subroutine fprati(p1,f1,p2,f2,p3,f3,p)
      !  ..scalar arguments..
      real(RKIND), intent(inout) :: p1,f1,p3,f3
      real(RKIND), intent(in)    :: p2,f2
      real(RKIND), intent(out)   :: p

      !  ..local scalars..
      real(RKIND) :: h1,h2,h3
      !  ..
      if (p3>zero) then
         h1 = f1*(f2-f3)
         h2 = f2*(f3-f1)
         h3 = f3*(f1-f2)
         p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
      else
         !  value of p in case p3 = infinity.
         p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      end if

      !  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
      if (f2>=zero) then
        p1 = p2
        f1 = f2
      else
        p3 = p2
        f3 = f2
      endif

      return
      end subroutine fprati


      recursive subroutine fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye, &
       kx,ky,s,nxest,nyest,tol,maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold, &
       reducx,reducy,fpintx,fpinty,lastdi,nplusx,nplusy,nrx,nry, &
       nrdatx,nrdaty,wrk,lwrk,ier)

      !  ..
      !  ..scalar arguments..
      real(RKIND) xb,xe,yb,ye,s,tol,fp,fp0,fpold,reducx,reducy
      integer iopt,mx,my,mz,kx,ky,nxest,nyest,maxit,nc,nx,ny,lastdi, &
       nplusx,nplusy,lwrk,ier
      !  ..array arguments..
      real(RKIND) x(mx),y(my),z(mz),tx(nxest),ty(nyest),c(nc),fpintx(nxest), &
       fpinty(nyest),wrk(lwrk)
      integer nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my)
      !  ..local scalars
      real(RKIND) acc,fpms,f1,f2,f3,p,p1,p2,p3,rn
      integer i,ich1,ich3,ifbx,ifby,ifsx,ifsy,iter,j,kx1,kx2,ky1,ky2, &
       k3,l,lax,lay,lbx,lby,lq,lri,lsx,lsy,mk1,mm,mpm,mynx,ncof, &
       nk1x,nk1y,nmaxx,nmaxy,nminx,nminy,nplx,nply,npl1,nrintx, &
       nrinty,nxe,nxk,nye

      !   set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND

      !  we partition the working space.
      kx1 = kx+1
      ky1 = ky+1
      kx2 = kx1+1
      ky2 = ky1+1
      lsx = 1
      lsy = lsx+mx*kx1
      lri = lsy+my*ky1
      mm = max0(nxest,my)
      lq = lri+mm
      mynx = nxest*my
      lax = lq+mynx
      nxk = nxest*kx2
      lbx = lax+nxk
      lay = lbx+nxk
      lby = lay+nyest*ky2
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 1: determination of the number of knots and their position.     c
      ! ****************************************************************     c
      !  given a set of knots we compute the least-squares spline sinf(x,y), c
      !  and the corresponding sum of squared residuals fp=f(p=inf).         c
      !  if iopt=-1  sinf(x,y) is the requested approximation.               c
      !  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
      !    if fp <=s we will continue with the current set of knots.         c
      !    if fp > s we will increase the number of knots and compute the    c
      !       corresponding least-squares spline until finally fp<=s.        c
      !    the initial choice of knots depends on the value of s and iopt.   c
      !    if s=0 we have spline interpolation; in that case the number of   c
      !    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               c
      !    if s>0 and                                                        c
      !     *iopt=0 we first compute the least-squares polynomial of degree  c
      !      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       c
      !     *iopt=1 we start with the knots found at the last call of the    c
      !      routine, except for the case that s > fp0; then we can compute  c
      !      the least-squares polynomial directly.                          c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  determine the number of knots for polynomial approximation.
      nminx = 2*kx1
      nminy = 2*ky1
      if(iopt<0) go to 120
      !  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      !  find nmaxx and nmaxy which denote the number of knots in x- and y-
      !  direction in case of spline interpolation.
      nmaxx = mx+kx1
      nmaxy = my+ky1
      !  find nxe and nye which denote the maximum number of knots
      !  allowed in each direction
      nxe = min0(nmaxx,nxest)
      nye = min0(nmaxy,nyest)
      if(s>0.) go to 100
      !  if s = 0, s(x,y) is an interpolating spline.
      nx = nmaxx
      ny = nmaxy
      !  test whether the required storage space exceeds the available one.
      if(ny>nyest .or. nx>nxest) go to 420
      !  find the position of the interior knots in case of interpolation.
      !  the knots in the x-direction.
      mk1 = mx-kx1
      if(mk1==0) go to 60
      k3 = kx/2
      i = kx1+1
      j = k3+2
      if(k3*2==kx) go to 40
      do 30 l=1,mk1
        tx(i) = x(j)
        i = i+1
        j = j+1
  30  continue
      go to 60
  40  do 50 l=1,mk1
        tx(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  50  continue
      !  the knots in the y-direction.
  60  mk1 = my-ky1
      if(mk1==0) go to 120
      k3 = ky/2
      i = ky1+1
      j = k3+2
      if(k3*2==ky) go to 80
      do 70 l=1,mk1
        ty(i) = y(j)
        i = i+1
        j = j+1
  70  continue
      go to 120
  80  do 90 l=1,mk1
        ty(i) = (y(j)+y(j-1))*half
        i = i+1
        j = j+1
  90  continue
      go to 120
      !  if s > 0 our initial choice of knots depends on the value of iopt.
 100  if(iopt==0) go to 115
      if(fp0<=s) go to 115
      !  if iopt=1 and fp0 > s we start computing the least- squares spline
      !  according to the set of knots found at the last call of the routine.
      !  we determine the number of grid coordinates x(i) inside each knot
      !  interval (tx(l),tx(l+1)).
      l = kx2
      j = 1
      nrdatx(1) = 0
      mpm = mx-1
      do 105 i=2,mpm
        nrdatx(j) = nrdatx(j)+1
        if(x(i)<tx(l)) go to 105
        nrdatx(j) = nrdatx(j)-1
        l = l+1
        j = j+1
        nrdatx(j) = 0
 105  continue
      !  we determine the number of grid coordinates y(i) inside each knot
      !  interval (ty(l),ty(l+1)).
      l = ky2
      j = 1
      nrdaty(1) = 0
      mpm = my-1
      do 110 i=2,mpm
        nrdaty(j) = nrdaty(j)+1
        if(y(i)<ty(l)) go to 110
        nrdaty(j) = nrdaty(j)-1
        l = l+1
        j = j+1
        nrdaty(j) = 0
 110  continue
      go to 120
      !  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
      !  polynomial of degree kx in x and ky in y (which is a spline without
      !  interior knots).
 115  nx = nminx
      ny = nminy
      nrdatx(1) = mx-2
      nrdaty(1) = my-2
      lastdi = 0
      nplusx = 0
      nplusy = 0
      fp0 = 0.
      fpold = 0.
      reducx = 0.
      reducy = 0.
 120  mpm = mx+my
      ifsx = 0
      ifsy = 0
      ifbx = 0
      ifby = 0
      p = -one
      !  main loop for the different sets of knots.mpm=mx+my is a save upper
      !  bound for the number of trials.
      do 250 iter=1,mpm
        if(nx==nminx .and. ny==nminy) ier = -2
      !  find nrintx (nrinty) which is the number of knot intervals in the
      !  x-direction (y-direction).
        nrintx = nx-nminx+1
        nrinty = ny-nminy+1
      !  find ncof, the number of b-spline coefficients for the current set
      !  of knots.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
      !  find the position of the additional knots which are needed for the
      !  b-spline representation of s(x,y).
        i = nx
        do 130 j=1,kx1
          tx(j) = xb
          tx(i) = xe
          i = i-1
 130    continue
        i = ny
        do 140 j=1,ky1
          ty(j) = yb
          ty(i) = ye
          i = i-1
 140    continue
      !  find the least-squares spline sinf(x,y) and calculate for each knot
      !  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
      !  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
      !  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
      !  value belonging to that interval.
      !  fp gives the total sum of squared residuals.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty, &
        ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx), &
        wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby), &
        nrx,nry)
        if(ier==(-2)) fp0 = fp
      !  test whether the least-squares spline is an acceptable solution.
        if(iopt<0) go to 440
        fpms = fp-s
        if(abs(fpms) < acc) go to 440
      !  if f(p=inf) < s, we accept the choice of knots.
        if(fpms<0.) go to 300
      !  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
        if(nx==nmaxx .and. ny==nmaxy) go to 430
      !  increase the number of knots.
      !  if nx=nxe and ny=nye we cannot further increase the number of knots
      !  because of the storage capacity limitation.
        if(nx==nxe .and. ny==nye) go to 420
        ier = 0
      !  adjust the parameter reducx or reducy according to the direction
      !  in which the last added knots were located.
        if (lastdi<0) go to 150
        if (lastdi==0) go to 170
        go to 160
 150    reducx = fpold-fp
        go to 170
 160    reducy = fpold-fp
      !  store the sum of squared residuals for the current set of knots.
 170    fpold = fp
      !  find nplx, the number of knots we should add in the x-direction.
        nplx = 1
        if(nx==nminx) go to 180
        npl1 = nplusx*2
        rn = nplusx
        if(reducx>acc) npl1 = int(rn*fpms/reducx)
        nplx = min0(nplusx*2,max0(npl1,nplusx/2,1))
      !  find nply, the number of knots we should add in the y-direction.
 180    nply = 1
        if(ny==nminy) go to 190
        npl1 = nplusy*2
        rn = nplusy
        if(reducy>acc) npl1 = int(rn*fpms/reducy)
        nply = min0(nplusy*2,max0(npl1,nplusy/2,1))
 190    if (nplx<nply) go to 210
        if (nplx==nply) go to 200
        go to 230
 200    if(lastdi<0) go to 230
 210    if(nx==nxe) go to 230
      !  addition in the x-direction.
        lastdi = -1
        nplusx = nplx
        ifsx = 0
        do 220 l=1,nplusx
      !  add a new knot in the x-direction
          call fpknot(x,mx,tx,nx,fpintx,nrdatx,nrintx,nxest,1)
      !  test whether we cannot further increase the number of knots in the
      !  x-direction.
          if(nx==nxe) go to 250
 220    continue
        go to 250
 230    if(ny==nye) go to 210
      !  addition in the y-direction.
        lastdi = 1
        nplusy = nply
        ifsy = 0
        do 240 l=1,nplusy
      !  add a new knot in the y-direction.
          call fpknot(y,my,ty,ny,fpinty,nrdaty,nrinty,nyest,1)
      !  test whether we cannot further increase the number of knots in the
      !  y-direction.
          if(ny==nye) go to 250
 240    continue
      !  restart the computations with the new set of knots.
 250  continue
      !  test whether the least-squares polynomial is a solution of our
      !  approximation problem.
 300  if(ier==(-2)) go to 440
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 2: determination of the smoothing spline sp(x,y)                c
      ! *****************************************************                c
      !  we have determined the number of knots and their position. we now   c
      !  compute the b-spline coefficients of the smoothing spline sp(x,y).  c
      !  this smoothing spline varies with the parameter p in such a way thatc
      !    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             c
      !  is a continuous, strictly decreasing function of p. moreover the    c
      !  least-squares polynomial corresponds to p=0 and the least-squares   c
      !  spline to p=infinity. iteratively we then have to determine the     c
      !  positive value of p such that f(p)=s. the process which is proposed c
      !  here makes use of rational interpolation. f(p) is approximated by a c
      !  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
      !  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
      !  are used to calculate the new value of p such that r(p)=s.          c
      !  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      ich1 = 0
      ich3 = 0
      !  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
      !  find the smoothing spline sp(x,y) and the corresponding sum of
      !  squared residuals fp.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty, &
        ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx), &
        wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby), &
        nrx,nry)
      !  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<acc) go to 440
      !  test whether the maximum allowable number of iterations has been
      !  reached.
        if(iter==maxit) go to 400
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 320
        if((f2-f3)>acc) go to 310
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2<0.) ich3 = 1
 320    if(ich1/=0) go to 340
        if((f1-f2)>acc) go to 330
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 350
        if(p>=p3) p = p2*con1 + p3*con9
        go to 350
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 330    if(f2>0.) ich1 = 1
 340    if(f2>=f1 .or. f2<=f3) go to 410
        ! find the new value of p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 350  continue
      !  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end subroutine fpregr

      ! subroutine fprota applies a givens rotation to a and b.
      elemental subroutine fprota(cos,sin,a,b)

          !  ..scalar arguments..
          real(RKIND), intent(in)    :: cos,sin
          real(RKIND), intent(inout) :: a,b

          ! ..local scalars..
          real(RKIND) stor1,stor2

          !  ..
          stor1 = a
          stor2 = b
          b = cos*stor2+sin*stor1
          a = cos*stor1-sin*stor2
          return

      end subroutine fprota


      !  given the coefficients of a constrained bicubic spline, as determined in subroutine fppola,
      !  fprppo calculates the coefficients in the standard b-spline representation of bicubic splines.
      pure subroutine fprppo(nu,nv,if1,if2,cosi,ratio,c,f,ncoff)

      !  ..scalar arguments..
      real(RKIND), intent(in) :: ratio
      integer, intent(in) :: if1,if2,nu,nv,ncoff
      !  ..array arguments
      real(RKIND), intent(inout) :: c(ncoff),f(ncoff)
      real(RKIND), intent(in) :: cosi(5,nv)
      !  ..local scalars..
      integer :: i,iopt,ii,j,k,l,nu4,nvv
      !  ..
      nu4  = nu-4
      nvv  = nv-7
      iopt = if1+1
      f    = zero
      i    = 0
      j    = 0

      main_loop: do l=1,nu4
         ii = i
         if (l==nu4 .and. if2/=0) then

            exit main_loop

         elseif (l>iopt) then

             do k=1,nvv
                i = i+1
                j = j+1
                f(i) = c(j)
             end do

         elseif (l==1) then

             do k=1,nvv
                i = i+1
                f(i) = c(1)
             end do
             j = 1

         elseif (l==2) then

             do k=1,nvv
                i = i+1
                f(i) = c(1)+c(2)*cosi(1,k)+c(3)*cosi(2,k)
             end do
             j = 3


         elseif (l==3) then

             do k=1,nvv
                i = i+1
                f(i) = c(1)+ratio*(c(2)*cosi(1,k)+c(3)*cosi(2,k))+ &
                       c(4)*cosi(3,k)+c(5)*cosi(4,k)+c(6)*cosi(5,k)
             end do
             j = 6

         end if

         do k=1,3
            ii   = ii+1
            i    = i+1
            f(i) = f(ii)
         end do

      end do main_loop

      c = f

      return
      end subroutine fprppo

      ! given the coefficients of a spherical spline function, subroutine fprpsp calculates the
      ! coefficients in the standard b-spline representation of this bicubic spline.
      pure subroutine fprpsp(nt,np,co,si,c,f,ncoff)
          !  ..
          !  ..scalar arguments
          integer,     intent(in)    :: nt,np,ncoff

          !  ..array arguments
          real(RKIND), intent(in)    :: co(np),si(np)
          real(RKIND), intent(inout) :: c(ncoff),f(ncoff) ! f is a working array

          !  ..local scalars
          real(RKIND) :: cn,c1,c2,c3
          integer :: i,ii,j,k,l,ncof,npp,np4,nt4
          !  ..
          nt4 = nt-4
          np4 = np-4
          npp = np4-3
          ncof = 6+npp*(nt4-4)
          c1 = c(1)
          cn = c(ncof)
          f(1:np4)             = c1
          f(ncoff-np4+1:ncoff) = cn

          i = np4
          j = 1
          do l=3,nt4
             ii = i
             if (l==3 .or. l==nt4) then
                if (l==nt4) c1 = cn
                c2 = c(j+1)
                c3 = c(j+2)
                j = j+2
                do k=1,npp
                   i = i+1
                   f(i) = c1+c2*co(k)+c3*si(k)
                end do
             else
                do k=1,npp
                   i = i+1
                   j = j+1
                   f(i) = c(j)
                end do
             endif

             do k=1,3
                ii = ii+1
                i = i+1
                f(i) = f(ii)
             end do
          end do

          c(1:ncoff) = f(1:ncoff)

      return
      end subroutine fprpsp


      ! subroutine fpseno fetches a branch of a triply linked tree the information of which is kept
      ! in the arrays up,left,right and info.
      ! the branch has a specified length nbind and is determined by the parameter merk which points to
      ! its terminal node. the information field of the nodes of this branch is stored in the array
      ! ibind. on exit merk points to a new branch of length nbind or takes the value 1 if no such
      !  branch was found.
      pure subroutine fpseno(maxtr,up,left,right,info,merk,ibind,nbind)

          !  ..scalar arguments..
          integer, intent(in)    :: maxtr   ! Tree array sizes
          integer, intent(inout) :: merk    ! (in) terminal node of the branch
          integer, intent(in)    :: nbind
          !  ..array arguments..
          integer, intent(in)    :: up(maxtr),left(maxtr),right(maxtr),info(maxtr)
          integer, intent(out)   :: ibind(nbind)
          !  ..scalar arguments..
          integer :: i,j,k
          !  ..
          k = merk
          j = nbind
          do i=1,nbind
             ibind(j) = info(k)
             k = up(k)
             j = j-1
          end do

          do
              k = right(merk)
              if (k/=0) exit
              merk = up(merk)
              if (merk<=1) return
          end do

          do
              merk = k
              k = left(merk)
              if (k==0) exit
          end do

      end subroutine fpseno


      recursive subroutine fpspgr(iopt,ider,u,mu,v,mv,r,mr,r0,r1,s, &
       nuest,nvest,tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu, &
       reducv,fpintu,fpintv,dr,step,lastdi,nplusu,nplusv,lastu0, &
       lastu1,nru,nrv,nrdatu,nrdatv,wrk,lwrk,ier)

      !  ..
      !  ..scalar arguments..
      integer, intent(in) :: mu,mv,mr,nuest,nvest,maxit,nc
      integer :: nu,nv,lastdi,nplusu,nplusv,lastu0,lastu1,lwrk,ier
      real(RKIND) :: r0,r1,s,tol,fp,fp0,fpold,reducu,reducv
      !  ..array arguments..
      integer :: iopt(3),ider(4),nrdatu(nuest),nrdatv(nvest),nru(mu),nrv(mv)
      real(RKIND) :: u(mu),v(mv),r(mr),tu(nuest),tv(nvest),c(nc),fpintu(nuest),fpintv(nvest),dr(6), &
                     wrk(lwrk),step(2)
      !  ..local scalars..
      real(RKIND) acc,fpms,f1,f2,f3,p,p1,p2,p3,vb,ve,rmax,rmin,rn
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,istart,iter,i1,i2,j,ju, &
       ktu,l,l1,l2,l3,l4,mpm,mumin,mu0,mu1,nn,nplu,nplv,npl1,nrintu, &
       nrintv,nue,numax,nve,nvmax
      !  ..local arrays..
      integer idd(4)
      real(RKIND) drr(6)

      !   set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND

      !   initialization
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one

      mumin = 4
      if (ider(1)>=0) mumin = mumin-1
      if (iopt(2)==1 .and. ider(2)==1) mumin = mumin-1
      if (ider(3)>=0) mumin = mumin-1
      if (iopt(3)==1 .and. ider(4)==1) mumin = mumin-1
      if (mumin==0) mumin = 1

      vb = v(1)
      ve = vb+pi2

      ! *****
      ! part 1: determination of the number of knots and their position.
      ! *****
      !  given a set of knots we compute the least-squares spline sinf(u,v) and the corresponding sum of
      !  squared residuals fp = f(p=inf).
      !  if iopt(1)=-1  sinf(u,v) is the requested approximation.
      !  if iopt(1)>=0  we check whether we can accept the knots:
      !    if fp <= s we will continue with the current set of knots.
      !    if fp >  s we will increase the number of knots and compute the corresponding least-squares spline
      !               until finally fp <= s.
      !    the initial choice of knots depends on the value of s and iopt.
      !    if s=0 we have spline interpolation; in that case the number of knots in the u-direction equals
      !     nu=numax=mu+6+iopt(2)+iopt(3) and in the v-direction nv=nvmax=mv+7.
      !    if s>0 and
      !      iopt(1)=0 we first compute the least-squares polynomial,i.e. a spline without interior knots:
      !                nu=8 ; nv=8.
      !      iopt(1)=1 we start with the set of knots found at the last call of the routine, except for the
      !                case that s > fp0; then we compute the least-squares polynomial directly.
      ! *****

      if (iopt(1)<0) go to 120

      !  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      !  numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+6+iopt(2)+iopt(3)
      nvmax = mv+7
      nue = min(numax,nuest)
      nve = min(nvmax,nvest)

      if (s>zero) go to 100

      !  if s = 0, s(u,v) is an interpolating spline.
      nu = numax
      nv = nvmax

      !  test whether the required storage space exceeds the available one.
      if (nu>nuest .or. nv>nvest) then
         ier = FITPACK_INSUFFICIENT_STORAGE
         return
      end if

      !  find the position of the knots in the v-direction.
      tv(4:mv+4) = [v(1:mv),ve]
      l1 = mv-2
      l2 = mv+5
      do i=1,3
         tv(i) = v(l1)-pi2
         tv(l2) = v(i+1)+pi2
         l1 = l1+1
         l2 = l2+1
      end do
      !  if not all the derivative values g(i,j) are given, we will first
      !  estimate these values by computing a least-squares spline
      idd = ider
      if(idd(1)==0) idd(1) = 1
      if(idd(1)>0)  dr(1) = r0
      if(idd(3)==0) idd(3) = 1
      if(idd(3)>0)  dr(4) = r1
      if(ider(1)<0 .or. ider(3)<0 .or. &
         (iopt(2)/=0 .and. ider(2)==0)) go to 30
         if(iopt(3)==0 .or. ider(4)/=0) go to 70
      ! we set up the knots in the u-direction for computing the least-squares
      ! spline.
  30  i1 = 3
      i2 = mu-2
      nu = 4
      do i=1,mu
         if (i1>i2) exit
         nu = nu+1
         tu(nu) = u(i1)
         i1 = i1+2
      end do
      do i=1,4
         tu(i) = zero
         nu = nu+1
         tu(nu) = pi
      end do

      ! we compute the least-squares spline for estimating the derivatives.
      call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,         &
                  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
                  wrk,lwrk)
      ifsu = 0
      ! if all the derivatives at the origin are known, we compute the
      ! interpolating spline.
      ! we set up the knots in the u-direction, needed for interpolation.
  70  nn = numax-8
      if(nn==0) go to 95
      ju = 2-iopt(2)
      do l=1,nn
        tu(l+4) = u(ju)
        ju = ju+1
      end do
      nu = numax
      l  = nu
      do i=1,4
         tu(i) = zero
         tu(l) = pi
         l = l-1
      end do

      ! we compute the interpolating spline.
  95  call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,         &
                  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv, &
                  wrk,lwrk)

      ier = FITPACK_INTERPOLATING_OK
      fp  = zero
      return

      !  if s>0 our initial choice of knots depends on the value of iopt(1).
 100  ier = FITPACK_OK
      if (iopt(1)/=0) step = -step

      if (fp0<=s .or. iopt(1)==0) then

          !  if iopt(1)=0 or iopt(1)=1 and fp0 <= s,we start computing the least-
          !  squares polynomial (which is a spline without interior knots).
          ier = FITPACK_LEASTSQUARES_OK
          idd = [ider(1),1,ider(3),1]
          nu = 8
          nv = 8
          nrdatu(1) = mu-2+iopt(2)+iopt(3)
          nrdatv(1) = mv-1
          lastdi = 0
          nplusu = 0
          nplusv = 0
          fp0    = zero
          fpold  = zero
          reducu = zero
          reducv = zero


      else

          !  if iopt(1)=1 and fp0 > s we start computing the least-squares spline
          !  according to the set of knots found at the last call of the routine.
          !  we determine the number of grid coordinates u(i) inside each knot
          !  interval (tu(l),tu(l+1)).
          l = 5
          j = 1
          nrdatu(1) = 0
          mu0 = 2-iopt(2)
          mu1 = mu-1+iopt(3)
          do i=mu0,mu1
            nrdatu(j) = nrdatu(j)+1
            if (u(i)>=tu(l)) then
                nrdatu(j) = nrdatu(j)-1
                l = l+1
                j = j+1
                nrdatu(j) = 0
            endif
          end do
          !  we determine the number of grid coordinates v(i) inside each knot
          !  interval (tv(l),tv(l+1)).
          l = 5
          j = 1
          nrdatv(1) = 0
          do i=2,mv
            nrdatv(j) = nrdatv(j)+1
            if (v(i)>=tv(l)) then
               nrdatv(j) = nrdatv(j)-1
               l = l+1
               j = j+1
               nrdatv(j) = 0
            endif
          end do
          idd = ider

      endif

      !  main loop for the different sets of knots.mpm=mu+mv is a safe upper
      !  bound for the number of iterations
 120  mpm = mu+mv
      do 270 iter=1,mpm
         ! find nrintu (nrintv) which is the number of knot intervals in the
         ! u-direction (v-direction).
         nrintu = nu-7
         nrintv = nv-7

         ! find the position of the additional knots which are needed for the
         ! b-spline representation of s(u,v).
         i = nu
         do j=1,4
            tu(j) = zero
            tu(i) = pi
            i = i-1
         end do
         l1 = 4
         l2 = l1
         l3 = nv-3
         l4 = l3
         tv(l2) = vb
         tv(l3) = ve
         do j=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tv(l2) = tv(l4)-pi2
            tv(l3) = tv(l1)+pi2
         end do

         ! find an estimate of the range of possible values for the optimal derivatives at the origin.
         ktu = max(mumin,nrdatu(1)+2-iopt(2))
         if (ktu/=lastu0) then
             l       = mv*ktu
             rmin    = min(r0,minval(r(:l)))
             rmax    = max(r0,maxval(r(:l)))
             step(1) = rmax-rmin
             lastu0  = ktu
         endif

         ktu = max(mumin,nrdatu(nrintu)+2-iopt(3))
         if (ktu/=lastu1) then
            l = mv*ktu
            rmin = min(r1,minval(r(mr-l+1:mr)))
            rmax = max(r1,maxval(r(mr-l+1:mr)))
            step(2) = rmax-rmin
            lastu1 = ktu
         endif

         ! find the least-squares spline sinf(u,v).
         call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt, &
                     idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru, &
                     nrv,wrk,lwrk)

         step = abs(step)
         if (ier==FITPACK_LEASTSQUARES_OK) fp0 = fp

         ! test whether the least-squares spline is an acceptable solution.
         if (iopt(1)<0) return
         fpms = fp-s; if (abs(fpms) < acc) return ! success

         ! if f(p=inf) < s, we accept the choice of knots.
         if (fpms<zero) go to 300

         ! if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
         if(nu==numax .and. nv==nvmax) then
            ier = FITPACK_INTERPOLATING_OK
            fp  = zero
            return
         endif

         ! increase the number of knots.
         ! if nu=nue and nv=nve we cannot further increase the number of knots
         ! because of the storage capacity limitation.
         if (nu==nue .and. nv==nve) then
            ier = FITPACK_INSUFFICIENT_STORAGE
            return
         endif

         if(ider(1)==0) fpintu(1)      = fpintu(1)     +(r0-dr(1))**2
         if(ider(3)==0) fpintu(nrintu) = fpintu(nrintu)+(r1-dr(4))**2
         ier = FITPACK_OK

         ! adjust the parameter reducu or reducv according to the direction
         ! in which the last added knots were located.
         first_knot: if (lastdi==0) then

            nplv = 3
            idd(2) = ider(2)
            idd(4) = ider(4)
            fpold = fp

         else first_knot

            if (lastdi<0) reducu = fpold-fp
            if (lastdi>0) reducv = fpold-fp

            ! store the sum of squared residuals for the current set of knots.
            fpold = fp

            ! find nplu, the number of knots we should add in the u-direction.
            nplu = 1
            if (nu/=8) then
               npl1 = nplusu*2
               rn = nplusu
               if (reducu>acc) npl1 = int(rn*fpms/reducu)
               nplu = min(nplusu*2,max(npl1,nplusu/2,1))
            endif

            ! find nplv, the number of knots we should add in the v-direction.
            nplv = 3
            if (nv/=8) then
               npl1 = nplusv*2
               rn = nplusv
               if (reducv>acc) npl1 = int(rn*fpms/reducv)
               nplv = min(nplusv*2,max(npl1,nplusv/2,1))
            endif

         endif first_knot



         ! test whether we are going to add knots in the u- or v-direction.
         ! lastdi = last knot direction: lastdi==0  = not yet set
         !                               lastdi==1  = v direction
         !                               lastdi==-1 = u direction
        choose_dir: if ( nv<nve .and. &
                         ((nu==nue .and. nplu<nplv .or. (nplu==nplv .and. lastdi>=0)) &
                          .or. nplu>nplv &
                          .or. (nplu==nplv .and. lastdi<0) )) then

            ! addition in the v-direction.
            lastdi = 1
            nplusv = nplv
            ifsv   = 0

            add_v_knots: do l=1,nplusv

               ! add a new knot in the v-direction.
               call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)

               ! test whether we cannot further increase the number of knots in the v-direction.
               if (nv==nve) exit add_v_knots

            end do add_v_knots

        else choose_dir

            ! addition in the u-direction.
            lastdi = -1
            nplusu = nplu
            ifsu   = 0
            istart = merge(1,0,iopt(2)==0)
            add_u_knots: do l=1,nplusu

               ! add a new knot in the u-direction
               call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,istart)

               ! test whether we cannot further increase the number of knots in the u-direction.
               if (nu==nue) exit add_u_knots

            end do add_u_knots

        endif choose_dir

        ! restart the computations with the new set of knots.
        270  continue

        ! test whether the least-squares polynomial is a solution of our approximation problem.
 300    if(ier==FITPACK_LEASTSQUARES_OK) return

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 2: determination of the smoothing spline sp(u,v)                c
      ! *****************************************************                c
      !  we have determined the number of knots and their position. we now   c
      !  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
      !  this smoothing spline depends on the parameter p in such a way that c
      !    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             c
      !  is a continuous, strictly decreasing function of p. moreover the    c
      !  least-squares polynomial corresponds to p=0 and the least-squares   c
      !  spline to p=infinity. then iteratively we have to determine the     c
      !  positive value of p such that f(p)=s. the process which is proposed c
      !  here makes use of rational interpolation. f(p) is approximated by a c
      !  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
      !  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
      !  are used to calculate the new value of p such that r(p)=s.          c
      !  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  initial value for p.
      p1   = zero
      f1   = fp0-s
      p3   = -one
      f3   = fpms
      p    = one
      drr  = dr
      ich1 = 0
      ich3 = 0


      !  iteration process to find the root of f(p)=s.
      find_root: do iter = 1,maxit

         ! find the smoothing spline sp(u,v) and the corresponding sum f(p).
         call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,drr,iopt, &
                     idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu, &
                     fpintv,nru,nrv,wrk,lwrk)

         ! test whether the approximation sp(u,v) is an acceptable solution.
         fpms = fp-s; if(abs(fpms)<acc) return ! success

         ! carry out one more step of the iteration process.
         p2 = p
         f2 = fpms
         if (ich3==0) then
            if ((f2-f3)<=acc) then
               ! our initial choice of p is too large.
               p3 = p2
               f3 = f2
               p = p*con4
               if (p<=p1) p = p1*con9 + p2*con1
               cycle find_root
            elseif (f2<zero) then
                ich3 = 1
            endif
        endif

        if (ich1==0) then
           if ((f1-f2)<=acc) then
              ! our initial choice of p is too small
              p1 = p2
              f1 = f2
              p = p/con4
              if (p3<zero) cycle find_root
              if (p>=p3) p = p2*con1 + p3*con9
              cycle find_root
           elseif (f2>zero) then
              ! test whether the iteration process proceeds as theoretically expected.
              ich1 = 1
           endif
        endif

        if (f2>=f1 .or. f2<=f3) then
           ier = FITPACK_S_TOO_SMALL
           return
        end if

        ! find the new value of p.
        call fprati(p1,f1,p2,f2,p3,f3,p)

      end do find_root

      ! max number of iterations exceeded
      ier = FITPACK_MAXIT
      return

      end subroutine fpspgr


      recursive subroutine fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eta,tol,maxit, &
       ib1,ib3,nc,ncc,intest,nrest,nt,tt,np,tp,c,fp,sup,fpint,coord,f, &
       ff,row,coco,cosi,a,q,bt,bp,spt,spp,h,index,nummer,wrk,lwrk,ier)

      !  ..
      !  ..scalar arguments..
      integer iopt,m,ntest,npest,maxit,ib1,ib3,nc,ncc,intest,nrest, &
       nt,np,lwrk,ier
      real(RKIND) s,eta,tol,fp,sup
      !  ..array arguments..
      real(RKIND) teta(m),phi(m),r(m),w(m),tt(ntest),tp(npest),c(nc), &
       fpint(intest),coord(intest),f(ncc),ff(nc),row(npest),coco(npest), &
       cosi(npest),a(ncc,ib1),q(ncc,ib3),bt(ntest,5),bp(npest,5), &
       spt(m,4),spp(m,4),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
      !  ..local scalars..
      real(RKIND) aa,acc,arg,cn,co,c1,dmax,d1,d2,eps,facc,facs,fac1,fac2,fn, &
       fpmax,fpms,f1,f2,f3,hti,htj,p,pinv,piv,p1,p2,p3,ri,si, &
       sigma,sq,store,wi,rn
      integer i,iband,iband1,iband3,iband4,ich1,ich3,ii,ij,il,in,irot, &
       iter,i1,i2,i3,j,jlt,jrot,j1,j2,l,la,lf,lh,ll,lp,lt,lwest,l1,l2, &
       l3,l4,ncof,ncoff,npp,np4,nreg,nrint,nrr,nr1,ntt,nt4,nt6,num, &
       num1,rank
      !  ..local arrays..
      real(RKIND), dimension(SIZ_K+1) :: hp,ht

      !  set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND
      eps = sqrt(eta)

      ! Initializations
      lwest = 0
      ntt   = 0
      iband1 = 0


      if(iopt<0) go to 70
      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt==0) go to 10
      if(s<sup) then
        if (np<11) go to 60
        go to 70
      endif
      !  if iopt=0 we begin by computing the weighted least-squares polynomial
      !  of the form
      !     s(teta,phi) = c1*f1(teta) + cn*fn(teta)
      !  where f1(teta) and fn(teta) are the cubic polynomials satisfying
      !     f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta).
      !  the corresponding weighted sum of squared residuals gives the upper
      !  bound sup for the smoothing factor s.
  10  sup = zero
      d1 = zero
      d2 = zero
      c1 = zero
      cn = zero
      fac1 = pi*(one + half)
      fac2 = (one + one)/pi**3
      aa = zero
      do 40 i=1,m
         wi = w(i)
         ri = r(i)*wi
         arg = teta(i)
         fn = fac2*arg*arg*(fac1-arg)
         f1 = (one-fn)*wi
         fn = fn*wi
         if(fn==zero) go to 20
         call fpgivs(fn,d1,co,si)
         call fprota(co,si,f1,aa)
         call fprota(co,si,ri,cn)
 20      if(f1==zero) go to 30
         call fpgivs(f1,d2,co,si)
         call fprota(co,si,ri,c1)
 30      sup = sup+ri*ri
 40   continue
      if(d2/=zero) c1 = c1/d2
      if(d1/=zero) cn = (cn-aa*c1)/d1
      !  find the b-spline representation of this least-squares polynomial
      nt = 8
      np = 8
      c(1:8)   = c1
      c(9:16)  = cn
      tt(1:4)  = zero
      tt(5:8)  = pi
      tp(1:4)  = zero
      tp(5:8)  = pi2
      fp = sup
      !  test whether the least-squares polynomial is an acceptable solution
      fpms = sup-s
      if(fpms<acc) go to 960
      !  test whether we cannot further increase the number of knots.
  60  if(npest<11 .or. ntest<9) go to 950
      !  find the initial set of interior knots of the spherical spline in
      !  case iopt = 0.
      np = 11
      tp(5:7) = pi*[half,one,onep5]
      nt = 9
      tt(5) = tp(5)
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  part 1 : computation of least-squares spherical splines.            c
      !  ********************************************************            c
      !  if iopt < 0 we compute the least-squares spherical spline according c
      !  to the given set of knots.                                          c
      !  if iopt >=0 we compute least-squares spherical splines with increas-c
      !  ing numbers of knots until the corresponding sum f(p=inf)<=s.       c
      !  the initial set of knots then depends on the value of iopt:         c
      !    if iopt=0 we start with one interior knot in the teta-direction   c
      !              (pi/2) and three in the phi-direction (pi/2,pi,3*pi/2). c
      !    if iopt>0 we start with the set of knots found at the last call   c
      !              of the routine.                                         c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  main loop for the different sets of knots. m is a save upper bound
      !  for the number of trials.
  70  do 570 iter=1,m
      !  find the position of the additional knots which are needed for the
      !  b-spline representation of s(teta,phi).
         l1 = 4
         l2 = l1
         l3 = np-3
         l4 = l3
         tp(l2) = 0.
         tp(l3) = pi2
         do 80 i=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tp(l2) = tp(l4)-pi2
            tp(l3) = tp(l1)+pi2
  80     continue
        tt(1:4)     = zero
        tt(nt-3:nt) = pi
      !  find nrint, the total number of knot intervals and nreg, the number
      !  of panels in which the approximation domain is subdivided by the
      !  intersection of knots.
        ntt = nt-7
        npp = np-7
        nrr = npp/2
        nr1 = nrr+1
        nrint = ntt+npp
        nreg = ntt*npp
      !  arrange the data points according to the panel they belong to.
        call fporde(teta,phi,m,3,3,tt,nt,tp,np,nummer,index,nreg)
      !  find the b-spline coefficients coco and cosi of the cubic spline
      !  approximations sc(phi) and ss(phi) for cos(phi) and sin(phi).
        coco(1:npp) = zero
        cosi(1:npp) = zero
        a(1:npp,1:npp) = zero
      !  the coefficients coco and cosi are obtained from the conditions
      !  sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4.
        do 150 i=1,npp
           l2 = i+3
           arg = tp(l2)
           call fpbspl(tp,np,3,arg,l2,hp)
           row(1:npp) = zero
           ll = i
           do 120 j=1,3
              if(ll>npp) ll= 1
              row(ll) = row(ll)+hp(j)
              ll = ll+1
 120       continue
           facc = cos(arg)
           facs = sin(arg)
           do 140 j=1,npp
              piv = row(j)
              if (piv==zero) go to 140
              call fpgivs(piv,a(j,1),co,si)
              call fprota(co,si,facc,coco(j))
              call fprota(co,si,facs,cosi(j))
              if(j==npp) go to 150
              j1 = j+1
              i2 = 1
              do 130 l=j1,npp
                 i2 = i2+1
                 call fprota(co,si,row(l),a(j,i2))
 130          continue
 140       continue
 150    continue
        coco(:npp) = fpback(a,coco,npp,npp,ncc)
        cosi(:npp) = fpback(a,cosi,npp,npp,ncc)
      !  find ncof, the dimension of the spherical spline and ncoff, the
      !  number of coefficients in the standard b-spline representation.
        nt4 = nt-4
        np4 = np-4
        ncoff = nt4*np4
        ncof = 6+npp*(ntt-1)
      !  find the bandwidth of the observation matrix a.
        iband = 4*npp
        if(ntt==4) iband = 3*(npp+1)
        if(ntt<4) iband = ncof
        iband1 = iband-1
      !  initialize the observation matrix a.
        f(1:ncof) = zero
        a(1:ncof,1:iband) = zero
      !  initialize the sum of squared residuals.
        fp = zero
      !  fetch the data points in the new order. main loop for the
      !  different panels.
        do 340 num=1,nreg
      !  fix certain constants for the current panel; jrot records the column
      !  number of the first non-zero element in a row of the observation
      !  matrix according to a data point of the panel.
          num1 = num-1
          lt = num1/npp
          l1 = lt+4
          lp = num1-lt*npp+1
          l2 = lp+3
          lt = lt+1
          jrot = 0
          if(lt>2) jrot = 3+(lt-3)*npp
      !  test whether there are still data points in the current panel.
          in = index(num)
 170      if(in==0) go to 340
      !  fetch a new data point.
          wi = w(in)
          ri = r(in)*wi
      !  evaluate for the teta-direction, the 4 non-zero b-splines at teta(in)
          call fpbspl(tt,nt,3,teta(in),l1,ht)
      !  evaluate for the phi-direction, the 4 non-zero b-splines at phi(in)
          call fpbspl(tp,np,3,phi(in),l2,hp)
      !  store the value of these b-splines in spt and spp resp.
          spp(in,1:4) = hp(1:4)
          spt(in,1:4) = ht(1:4)
      !  initialize the new row of observation matrix.
         h(1:iband) = zero
      !  calculate the non-zero elements of the new row by making the cross
      !  products of the non-zero b-splines in teta- and phi-direction and
      !  by taking into account the conditions of the spherical splines.
         row(1:npp) = zero
      !  take into account the condition (3) of the spherical splines.
          ll = lp
          do 210 i=1,4
             if(ll>npp) ll=1
             row(ll) = row(ll)+hp(i)
             ll = ll+1
 210      continue
          ! take into account the other conditions of the spherical splines.
          if(lt<=2 .or. lt>=(ntt-1)) then
             facc = dot_product(row(:npp),coco(:npp))
             facs = dot_product(row(:npp),cosi(:npp))
          else
             facc = zero
             facs = zero
          endif
         ! fill in the non-zero elements of the new row.
         j1 = 0
         new_row: do j =1,4
            jlt = j+lt
            htj = ht(j)
            if (jlt==3) then
                h(1:3) = [h(1)+htj,facc*htj,facs*htj]
                j1 = 3
            elseif (jlt==nt4) then
                h(j1+1:j1+3) = htj*[facc,facs,one]
                j1 = j1+2
            elseif (jlt>2 .and. jlt<=nt4) then
                h(j1+1:j1+1:npp) = row(1:npp)*htj
                j1 = j1+npp
            else
                j1 = j1+1
                h(j1) = h(j1)+htj
            endif
          end do new_row
          h(:iband) = h(:iband1)*wi
      !  rotate the row into triangle by givens transformations.
          irot = jrot
          do 310 i=1,iband
            irot = irot+1
            piv = h(i)
            if (piv==zero) go to 310
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
      !  apply that transformation to the right hand side.
            call fprota(co,si,ri,f(irot))
            if(i==iband) go to 320
      !  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 300 j=i3,iband
              i2 = i2+1
              call fprota(co,si,h(j),a(irot,i2))
 300        continue
 310      continue
      !  add the contribution of the row to the sum of squares of residual
      !  right hand sides.
 320      fp = fp+ri**2
      !  find the number of the next data point in the panel.
          in = nummer(in)
          go to 170
 340    continue
      !  find dmax, the maximum value for the diagonal elements in the reduced
      !  triangle.
        dmax = max(zero,maxval(a(:ncof,1),a(:ncof,1)>=zero))
      !  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        if (any(a(1:ncof,1)<=sigma)) then
            lwest = ncof*iband+ncof+iband
            if(lwrk<lwest) go to 925
            lf = 1
            lh = lf+ncof
            la = lh+iband
            ff(1:ncof) = zero
            q(1:ncof,1:iband) = a(1:ncof,1:iband)
            call fprank(q,ff,ncof,iband,ncc,sigma,c,sq,rank,wrk(la), &
             wrk(lf),wrk(lh))
             q(1:ncof,1) = q(1:ncof,1)/dmax
          !  add to the sum of squared residuals, the contribution of reducing
          !  the rank.
            fp = fp+sq
        else
          !  backward substitution in case of full rank.
            c(:ncof) = fpback(a,f,ncof,iband,ncc)
            rank = ncof
            q(1:ncof,1) = a(1:ncof,1)/dmax
        endif
      !  in case of rank deficiency, find the minimum norm solution.

      !  find the coefficients in the standard b-spline representation of
      !  the spherical spline.
        call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
      !  test whether the least-squares spline is an acceptable solution.
        if(iopt<0) then
          if (fp<=0) go to 970
          go to 980
        endif
        fpms = fp-s
        if(abs(fpms)<=acc) then
          if (fp<=0) go to 970
          go to 980
        endif
      !  if f(p=inf) < s, accept the choice of knots.
        if(fpms<zero) go to 580
      !  test whether we cannot further increase the number of knots.
        if(ncof>m) go to 935
      !  search where to add a new knot.
      !  find for each interval the sum of squared residuals fpint for the
      !  data points having the coordinate belonging to that knot interval.
      !  calculate also coord which is the same sum, weighted by the position
      !  of the data points considered.
        fpint(:nrint) = zero
        coord(:nrint) = zero
        do 490 num=1,nreg
          num1 = num-1
          lt = num1/npp
          l1 = lt+1
          lp = num1-lt*npp
          l2 = lp+1+ntt
          jrot = lt*np4+lp
          in = index(num)
 460      if(in==0) go to 490
          store = 0.
          i1 = jrot
          do 480 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 470 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 470        continue
            i1 = i1+np4
 480      continue
          store = (w(in)*(r(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*teta(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*phi(in)
          in = nummer(in)
          go to 460
 490    continue
      !  find the interval for which fpint is maximal on the condition that
      !  there still can be added a knot.
        l1 = 1
        l2 = nrint
        if(ntest<nt+1) l1=ntt+1
        if(npest<np+2) l2=ntt
      !  test whether we cannot further increase the number of knots.
        if(l1>l2) go to 950
 500    fpmax = zero
        l = 0
        do 510 i=l1,l2
          if(fpmax>=fpint(i)) go to 510
          l = i
          fpmax = fpint(i)
 510    continue
        if(l==0) go to 930
      !  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
      !  test in what direction the new knot is going to be added.
        if(l>ntt) go to 530
      !  addition in the teta-direction
        l4 = l+4
        fpint(l) = zero
        fac1 = tt(l4)-arg
        fac2 = arg-tt(l4-1)
        if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) go to 500
        j = nt
        do 520 i=l4,nt
          tt(j+1) = tt(j)
          j = j-1
 520    continue
        tt(l4) = arg
        nt = nt+1
        go to 570
      !  addition in the phi-direction
 530    l4 = l+4-ntt
        if(arg<pi) go to 540
        arg = arg-pi
        l4 = l4-nrr
 540    fpint(l) = zero
        fac1 = tp(l4)-arg
        fac2 = arg-tp(l4-1)
        if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) go to 500
        ll = nrr+4
        j = ll
        do 550 i=l4,ll
          tp(j+1) = tp(j)
          j = j-1
 550    continue
        tp(l4) = arg
        np = np+2
        nrr = nrr+1
        do 560 i=5,ll
          j = i+nrr
          tp(j) = tp(i)+pi
 560    continue
      !  restart the computations with the new set of knots.
 570  continue
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 2: determination of the smoothing spherical spline.             c
      ! ********************************************************             c
      ! we have determined the number of knots and their position. we now    c
      ! compute the coefficients of the smoothing spline sp(teta,phi).       c
      ! the observation matrix a is extended by the rows of a matrix, expres-c
      ! sing that sp(teta,phi) must be a constant function in the variable   c
      ! phi and a cubic polynomial in the variable teta. the corresponding   c
      ! weights of these additional rows are set to 1/(p). iteratively       c
      ! we than have to determine the value of p such that f(p) = sum((w(i)* c
      ! (r(i)-sp(teta(i),phi(i))))**2)  be = s.                              c
      ! we already know that the least-squares polynomial corresponds to p=0,c
      ! and that the least-squares spherical spline corresponds to p=infin.  c
      ! the iteration process makes use of rational interpolation. since f(p)c
      ! is a convex and strictly decreasing function of p, it can be approx- c
      ! imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c
      ! three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c
      ! f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c
      ! of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<zeroc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !  evaluate the discontinuity jumps of the 3-th order derivative of
      !  the b-splines at the knots tt(l),l=5,...,nt-4.
 580  call fpdisc(tt,nt,5,bt,ntest)
      !  evaluate the discontinuity jumps of the 3-th order derivative of
      !  the b-splines at the knots tp(l),l=5,...,np-4.
      call fpdisc(tp,np,5,bp,npest)
      !  initial value for p.
      p1 = zero
      f1 = sup-s
      p3 = -one
      f3 = fpms
      p  = sum(a(1:ncof,1))
      rn = ncof
      p = rn/p
      !  find the bandwidth of the extended observation matrix.
      iband4 = iband+3
      if (ntt<=4) iband4 = ncof
      iband3 = iband4 -1
      ich1 = 0
      ich3 = 0
      !  iteration process to find the root of f(p)=s.
      do 920 iter=1,maxit
        pinv = one/p
      !  store the triangularized observation matrix into q.
        ff(1:ncof) = f(1:ncof)
        q(1:ncof,1:iband4) = zero
        q(1:ncof,1:iband) = q(1:ncof,1:iband)
      !  extend the observation matrix with the rows of a matrix, expressing
      !  that for teta=cst. sp(teta,phi) must be a constant function.
        nt6 = nt-6
        do 720 i=5,np4
          ii = i-4
          row(1:npp) = zero
          ll = ii
          do 620  l=1,5
             if(ll>npp) ll=1
             row(ll) = row(ll)+bp(ii,l)
             ll = ll+1
 620      continue
          facc = dot_product(row(1:npp),coco(1:npp))
          facs = dot_product(row(1:npp),cosi(1:npp))
          do 721 j=1,nt6
      !  initialize the new row.
            h(1:iband) = zero
      !  fill in the non-zero elements of the row. jrot records the column
      !  number of the first non-zero element in the row.
            jrot = 4+(j-2)*npp
            if(j>1 .and. j<nt6) then
                h(1:npp) = row(1:npp)
            else
                h(1) = facc
                h(2) = facs
                if(j==1) jrot = 2
            endif
            h(1:iband) = h(1:iband)*pinv
            ri = zero
      !  rotate the new row into triangle by givens transformations.
            do 710 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if (piv==zero) then
                if (i2<=0) go to 721
                go to 690
              endif
      !  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
      !  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2==0) go to 721
      !  apply that givens transformation to the left hand side.
              do 680 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 680          continue
 690          h(1:i2+1) = [h(2:i2+1),zero]
 710        continue
 721      continue
 720    continue
      !  extend the observation matrix with the rows of a matrix expressing
      !  that for phi=cst. sp(teta,phi) must be a cubic polynomial.
        do 810 i=5,nt4
          ii = i-4
          do 811 j=1,npp
      !  initialize the new row
          h(1:iband4) = zero
      !  fill in the non-zero elements of the row. jrot records the column
      !  number of the first non-zero element in the row.
            j1 = 1
            do 760 l=1,5
               il = ii+l
               ij = npp
               if(il/=3 .and. il/=nt4) go to 750
               j1 = j1+3-j
               j2 = j1-2
               ij = 0
               if(il/=3) go to 740
               j1 = 1
               j2 = 2
               ij = j+2
 740           h(j2) = bt(ii,l)*coco(j)
               h(j2+1) = bt(ii,l)*cosi(j)
 750           h(j1) = h(j1)+bt(ii,l)
               j1 = j1+ij
 760        continue
            do 765 l=1,iband4
               h(l) = h(l)*pinv
 765        continue
            ri = zero
            jrot = 1
            if(ii>2) jrot = 3+j+(ii-3)*npp
      !  rotate the new row into triangle by givens transformations.
            do 800 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if (piv==zero) then
                if (i2<=0) go to 811
                go to 780
              endif
      !  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
      !  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2==0) go to 811
      !  apply that givens transformation to the left hand side.
              do 770 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 770          continue
 780          h(1:i2+1) = [h(2:i2+1),zero]
 800        continue
 811      continue
 810    continue
      !  find dmax, the maximum value for the diagonal elements in the
      !  reduced triangle.
        dmax  = max(zero,maxval(q(1:ncof,1)))
      !  check whether the matrix is rank deficient.
        sigma = max(eps*dmax,maxval(q(1:ncof,1)))
      !  backward substitution in case of full rank.
        c(:ncof) = fpback(q,ff,ncof,iband4,ncc)
        rank = ncof
        go to 845
      !  in case of rank deficiency, find the minimum norm solution.
        lwest = ncof*iband4+ncof+iband4
        if(lwrk<lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,ncc,sigma,c,sq,rank,wrk(la), &
         wrk(lf),wrk(lh))
 845    q(1:ncof,1) = q(1:ncof,1)/dmax
      !  find the coefficients in the standard b-spline representation of
      !  the spherical spline.
        call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
      !  compute f(p).
        fp = zero
        do 890 num = 1,nreg
          num1 = num-1
          lt = num1/npp
          lp = num1-lt*npp
          jrot = lt*np4+lp
          in = index(num)
 860      if(in==0) go to 890
          store = zero
          i1 = jrot
          do 880 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 870 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 870        continue
            i1 = i1+np4
 880      continue
          fp = fp+(w(in)*(r(in)-store))**2
          in = nummer(in)
          go to 860
 890    continue
      !  test whether the approximation sp(teta,phi) is an acceptable solution
        fpms = fp-s
        if(abs(fpms)<=acc) go to 980
      !  test whether the maximum allowable number of iterations has been
      !  reached.
        if(iter==maxit) go to 940
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 900
        if((f2-f3)>acc) go to 895
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 + p2*con1
        go to 920
 895    if(f2<zero) ich3 = 1
 900    if(ich1/=0) go to 910
        if((f1-f2)>acc) go to 905
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<0.) go to 920
        if(p>=p3) p = p2*con1 +p3*con9
        go to 920
 905    if(f2>0.) ich1 = 1
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 910    if(f2>=f1 .or. f2<=f3) go to 945
      !  find the new value of p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 920  continue
      !  error codes and messages.
 925  ier = lwest
      go to 990
 930  ier = 5
      go to 990
 935  ier = 4
      go to 990
 940  ier = 3
      go to 990
 945  ier = 2
      go to 990
 950  ier = 1
      go to 990
 960  ier = -2
      go to 990
 970  ier = -1
      fp = zero
 980  if(ncof/=rank) ier = -rank
 990  return
      end subroutine fpsphe


      ! Once all inputs checked, do the actual b-spline surface evaluation
      pure subroutine fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,wu,wv,lu,lv)

      !  ..scalar arguments..
      integer, intent(in) :: idim,nu,nv,mu,mv
      !  ..array arguments..
      integer, intent(out)    :: lu(mu),lv(mv)
      real(RKIND), intent(in) :: tu(nu),tv(nv),c((nu-4)*(nv-4)*idim),u(mu),v(mv)
      real(RKIND), intent(out) :: wu(mu,4),wv(mv,4),f(mu*mv*idim)
      !  ..local scalars..
      integer :: i,i1,j,j1,k,l,l1,l3,m,nuv,nu4,nv4
      real(RKIND) :: arg,sp,tb,te
      !  ..local arrays..
      real(RKIND) :: h(SIZ_K+1)

      ! Process u
      nu4 = nu-4
      tb = tu(4)
      te = tu(nu4+1)
      l  = 4
      l1 = l+1
      do i=1,mu
        arg = min(max(tb,u(i)),te)
        do while (.not.(arg<tu(l1) .or. l==nu4))
           l = l1
           l1 = l+1
        end do
        call fpbspl(tu,nu,3,arg,l,h)
        lu(i) = l-4
        wu(i,1:4) = h(1:4)
      end do

      ! Process v
      nv4 = nv-4
      tb = tv(4)
      te = tv(nv4+1)
      l = 4
      l1 = l+1
      do i=1,mv
         arg = min(max(v(i),tb),te)
         do while (.not.(arg<tv(l1) .or. l==nv4))
           l = l1
           l1 = l+1
         end do
         call fpbspl(tv,nv,3,arg,l,h)
         lv(i) = l-4
         wv(i,1:4) = h(1:4)
      end do

      m = 0
      nuv = nu4*nv4
      dims: do k=1,idim
        l3 = (k-1)*nuv
        do i=1,mu
          l = lu(i)*nv4+l3
          h(1:4) = wu(i,1:4)
          do j=1,mv
            l1 = l+lv(j)
            sp = zero
            do i1=1,4
              do j1=1,4
                sp = sp+c(l1+j1)*h(i1)*wv(j,j1)
              end do
              l1 = l1+nv4
            end do
            m = m+1
            f(m) = sp
          end do
         end do
      end do dims
      return
      end subroutine fpsuev


      recursive subroutine fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kxx,kyy, &
       s,nxest, nyest,eta,tol,maxit,nmax,km1,km2,ib1,ib3,nc,intest, &
       nrest,nx0,tx,ny0,ty,c,fp,fp0,fpint,coord,f,ff,a,q,bx,by,spx, &
       spy,h,index,nummer,wrk,lwrk,ier)

      !  ..
      !  ..scalar arguments..
      real(RKIND) xb,xe,yb,ye,s,eta,tol,fp,fp0
      integer iopt,m,kxx,kyy,nxest,nyest,maxit,nmax,km1,km2,ib1,ib3, &
       nc,intest,nrest,nx0,ny0,lwrk,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c(nc),fpint(intest), &
       coord(intest),f(nc),ff(nc),a(nc,ib1),q(nc,ib3),bx(nmax,km2), &
       by(nmax,km2),spx(m,km1),spy(m,km1),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
      !  ..local scalars..
      real(RKIND) acc,arg,cos,dmax,fac1,fac2,fpmax,fpms,f1,f2,f3,hxi,p,pinv, &
       piv,p1,p2,p3,sigma,sin,sq,store,wi,x0,x1,y0,y1,zi,eps, &
       rn
      integer i,iband,iband1,iband3,iband4,ibb,ichang,ich1,ich3,ii, &
       in,irot,iter,i1,i2,i3,j,jrot,jxy,j1,kx,kx1,kx2,ky,ky1,ky2,l, &
       la,lf,lh,lwest,lx,ly,l1,l2,n,ncof,nk1x,nk1y,nminx,nminy,nreg, &
       nrint,num,num1,nx,nxe,nxx,ny,nye,nyy,n1,rank
      !  ..local arrays..
      real(RKIND), dimension(SIZ_K+1) :: hx,hy

      !  set constants
      real(RKIND), parameter :: con1 = 0.1e0_RKIND
      real(RKIND), parameter :: con9 = 0.9e0_RKIND
      real(RKIND), parameter :: con4 = 0.4e-01_RKIND
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 1: determination of the number of knots and their position.     c
      ! ****************************************************************     c
      ! given a set of knots we compute the least-squares spline sinf(x,y),  c
      ! and the corresponding weighted sum of squared residuals fp=f(p=inf). c
      ! if iopt=-1  sinf(x,y) is the requested approximation.                c
      ! if iopt=0 or iopt=1 we check whether we can accept the knots:        c
      !   if fp <=s we will continue with the current set of knots.          c
      !   if fp > s we will increase the number of knots and compute the     c
      !      corresponding least-squares spline until finally  fp<=s.        c
      ! the initial choice of knots depends on the value of s and iopt.      c
      !   if iopt=0 we first compute the least-squares polynomial of degree  c
      !     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c
      !     fp0=f(0) denotes the corresponding weighted sum of squared       c
      !     residuals                                                        c
      !   if iopt=1 we start with the knots found at the last call of the    c
      !     routine, except for the case that s>=fp0; then we can compute    c
      !     the least-squares polynomial directly.                           c
      ! eventually the independent variables x and y (and the corresponding  c
      ! parameters) will be switched if this can reduce the bandwidth of the c
      ! system to be solved.                                                 c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! ichang denotes whether(1) or not(-1) the directions have been inter changed.
      ichang = -1
      x0 = xb
      x1 = xe
      y0 = yb
      y1 = ye
      kx = kxx
      ky = kyy
      kx1 = kx+1
      ky1 = ky+1
      nxe = nxest
      nye = nyest
      eps = sqrt(eta)
      if(iopt<0) go to 20
      !  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt==0) go to 10
      if(fp0>s) go to 20
      !  initialization for the least-squares polynomial.
  10  nminx = 2*kx1
      nminy = 2*ky1
      nx = nminx
      ny = nminy
      ier = -2
      go to 30
  20  nx = nx0
      ny = ny0
      !  main loop for the different sets of knots. m is a save upper bound
      !  for the number of trials.
  30  do 420 iter=1,m
      !  find the position of the additional knots which are needed for the
      !  b-spline representation of s(x,y).
        l = nx
        do 40 i=1,kx1
          tx(i) = x0 ! 1, 2, ..., kx1
          tx(l) = x1 ! nx, nx-1, ..., nx+1-kx1
          l = l-1
  40    continue
        l = ny
        do 50 i=1,ky1
          ty(i) = y0
          ty(l) = y1
          l = l-1
  50    continue
      !  find nrint, the total number of knot intervals and nreg, the number
      !  of panels in which the approximation domain is subdivided by the
      !  intersection of knots.
        nxx = nx-2*kx1+1
        nyy = ny-2*ky1+1
        nrint = nxx+nyy
        nreg = nxx*nyy
      !  find the bandwidth of the observation matrix a.
      !  if necessary, interchange the variables x and y, in order to obtain
      !  a minimal bandwidth.
        iband1 = kx*(ny-ky1)+ky
        l = ky*(nx-kx1)+kx
        if(iband1<=l) go to 130
        iband1 = l
        ichang = -ichang
        call swap_RKIND(x,y)
        call swap_RKIND(x0,y0)
        call swap_RKIND(x1,y1)
        n  = min(nx,ny)
        n1 = n+1
        call swap_RKIND(tx ,ty)
        call swap_int  (nx ,ny)
        call swap_int  (nxe,nye)
        call swap_int  (nxx,nyy)
        call swap_int  (kx ,ky)
        kx1 = kx+1
        ky1 = ky+1
 130    iband = iband1+1
      !  arrange the data points according to the panel they belong to.
        call fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
      !  find ncof, the number of b-spline coefficients.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
      !  initialize the observation matrix a.
        f(1:ncof) = zero
        a(1:ncof,1:iband) = zero
      !  initialize the sum of squared residuals.
        fp = zero
      !  fetch the data points in the new order. main loop for the
      !  different panels.
        do 250 num=1,nreg
      !  fix certain constants for the current panel; jrot records the column
      !  number of the first non-zero element in a row of the observation
      !  matrix according to a data point of the panel.
          num1 = num-1
          lx = num1/nyy
          l1 = lx+kx1
          ly = num1-lx*nyy
          l2 = ly+ky1
          jrot = lx*nk1y+ly
      !  test whether there are still data points in the panel.
          in = index(num)
 150      if(in==0) go to 250
      !  fetch a new data point.
          wi = w(in)
          zi = z(in)*wi
      !  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
          call fpbspl(tx,nx,kx,x(in),l1,hx)
      !  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
          call fpbspl(ty,ny,ky,y(in),l2,hy)
      !  store the value of these b-splines in spx and spy respectively.
          spx(in,1:kx1) = hx(1:kx1)
          spy(in,1:ky1) = hy(1:ky1)
      !  initialize the new row of observation matrix.
          h(1:iband) = zero

      !  calculate the non-zero elements of the new row by making the cross
      !  products of the non-zero b-splines in x- and y-direction.
          i1 = 0
          do 200 i=1,kx1
            hxi = hx(i)
            j1 = i1
            do 190 j=1,ky1
              j1 = j1+1
              h(j1) = hxi*hy(j)*wi
 190        continue
            i1 = i1+nk1y
 200      continue
      !  rotate the row into triangle by givens transformations .
          irot = jrot
          do 220 i=1,iband
            irot = irot+1
            piv = h(i)
            if (piv==zero) go to 220
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),cos,sin)
      !  apply that transformation to the right hand side.
            call fprota(cos,sin,zi,f(irot))
            if(i==iband) go to 230
      !  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 210 j=i3,iband
              i2 = i2+1
              call fprota(cos,sin,h(j),a(irot,i2))
 210        continue
 220      continue
      !  add the contribution of the row to the sum of squares of residual
      !  right hand sides.
 230      fp = fp+zi**2
      !  find the number of the next data point in the panel.
          in = nummer(in)
          go to 150
 250    continue
      !  find dmax, the maximum value for the diagonal elements in the reduced triangle.
        dmax = max(zero,maxval(a(1:ncof,1)))
      !  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        if (all(a(1:ncof,1)>sigma)) then
            ! backward substitution in case of full rank.
            c(:ncof) = fpback(a,f,ncof,iband,nc)
            rank = ncof
            q(:ncof,1) = a(:ncof,1)/dmax
        else
            ! in case of rank deficiency, find the minimum norm solution.
            !  check whether there is sufficient working space
            lwest = ncof*iband+ncof+iband
            if(lwrk<lwest) go to 780
            ff(1:ncof) = f(1:ncof)
            q(1:ncof,1:iband)=a(1:ncof,1:iband)
            lf =1
            lh = lf+ncof
            la = lh+iband
            call fprank(q,ff,ncof,iband,nc,sigma,c,sq,rank,wrk(la),wrk(lf),wrk(lh))
            q(1:ncof,1)=q(1:ncof,1)/dmax
            !  add to the sum of squared residuals, the contribution of reducing the rank.
            fp = fp+sq
        endif
        if(ier==(-2)) fp0 = fp
      !  test whether the least-squares spline is an acceptable solution.
        if(iopt<0) go to 820
        fpms = fp-s
        if(abs(fpms)<=acc) then
          if (fp<=0) go to 815
          go to 820
        endif
      !  test whether we can accept the choice of knots.
        if(fpms<0.) go to 430
      !  test whether we cannot further increase the number of knots.
        if(ncof>m) go to 790
        ier = 0
      !  search where to add a new knot.
      !  find for each interval the sum of squared residuals fpint for the
      !  data points having the coordinate belonging to that knot interval.
      !  calculate also coord which is the same sum, weighted by the position
      !  of the data points considered.
        fpint(:nrint) = zero
        coord(:nrint) = zero
        do 360 num=1,nreg
          num1 = num-1
          lx = num1/nyy
          l1 = lx+1
          ly = num1-lx*nyy
          l2 = ly+1+nxx
          jrot = lx*nk1y+ly
          in = index(num)
 330      if(in==0) go to 360
          store = zero
          i1 = jrot
          do 350 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 340 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 340        continue
            i1 = i1+nk1y
 350      continue
          store = (w(in)*(z(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*x(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*y(in)
          in = nummer(in)
          go to 330
 360    continue
      !  find the interval for which fpint is maximal on the condition that
      !  there still can be added a knot.
 370    l = 0
        fpmax = zero
        l1 = 1
        l2 = nrint
        if(nx==nxe) l1 = nxx+1
        if(ny==nye) l2 = nxx
        if(l1>l2) go to 810
        do 380 i=l1,l2
          if(fpmax>=fpint(i)) go to 380
          l = i
          fpmax = fpint(i)
 380    continue
      !  test whether we cannot further increase the number of knots.
        if(l==0) go to 785
      !  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
      !  test in what direction the new knot is going to be added.
        if(l>nxx) go to 400
      !  addition in the x-direction.
        jxy = l+kx1
        fpint(l) = zero
        fac1 = tx(jxy)-arg
        fac2 = arg-tx(jxy-1)
        if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) go to 370
        j = nx
        do 390 i=jxy,nx
          tx(j+1) = tx(j)
          j = j-1
 390    continue
        tx(jxy) = arg
        nx = nx+1
        go to 420
      !  addition in the y-direction.
 400    jxy = l+ky1-nxx
        fpint(l) = zero
        fac1 = ty(jxy)-arg
        fac2 = arg-ty(jxy-1)
        if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) go to 370
        j = ny
        do 410 i=jxy,ny
          ty(j+1) = ty(j)
          j = j-1
 410    continue
        ty(jxy) = arg
        ny = ny+1
      !  restart the computations with the new set of knots.
 420  continue
      !  test whether the least-squares polynomial is a solution of our
      !  approximation problem.
 430  if(ier==(-2)) go to 830
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! part 2: determination of the smoothing spline sp(x,y)                c
      ! *****************************************************                c
      ! we have determined the number of knots and their position. we now    c
      ! compute the b-spline coefficients of the smoothing spline sp(x,y).   c
      ! the observation matrix a is extended by the rows of a matrix,        c
      ! expressing that sp(x,y) must be a polynomial of degree kx in x and   c
      ! ky in y. the corresponding weights of these additional rows are set  c
      ! to 1./p.  iteratively we than have to determine the value of p       c
      ! such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c
      ! we already know that the least-squares polynomial corresponds to     c
      ! p=0  and that the least-squares spline corresponds to p=infinity.    c
      ! the iteration process which is proposed here makes use of rational   c
      ! interpolation. since f(p) is a convex and strictly decreasing        c
      ! function of p, it can be approximated by a rational function r(p)=   c
      ! (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
      ! of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
      ! new value of p such that r(p)=s. convergence is guaranteed by taking c
      ! f1 > 0 and f3 < 0.                                                   c
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      kx2 = kx1+1
      !  test whether there are interior knots in the x-direction.
      if(nk1x==kx1) go to 440
      !  evaluate the discotinuity jumps of the kx-th order derivative of
      !  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
      call fpdisc(tx,nx,kx2,bx,nmax)
 440  ky2 = ky1 + 1
      !  test whether there are interior knots in the y-direction.
      if(nk1y==ky1) go to 450
      !  evaluate the discontinuity jumps of the ky-th order derivative of
      !  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
      call fpdisc(ty,ny,ky2,by,nmax)
      !  initial value for p.
 450  p1 = zero
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = sum(a(:ncof,1))
      rn = ncof
      p = rn/p
      !  find the bandwidth of the extended observation matrix.
      iband3 = kx1*nk1y
      iband4 = iband3 +1
      ich1 = 0
      ich3 = 0
      !  iteration process to find the root of f(p)=s.
      do 770 iter=1,maxit
        pinv = one/p
      !  store the triangularized observation matrix into q.
        ff(1:ncof) = f(1:ncof)
        do i=1,ncof
          q(i,1:iband) = a(i,1:iband)
          ibb = iband+1
          q(i,ibb:iband4) = zero
        end do

        if(nk1y==ky1) go to 560
      !  extend the observation matrix with the rows of a matrix, expressing
      !  that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
        do 550 i=ky2,nk1y
          ii = i-ky1
          do 551 j=1,nk1x
      !  initialize the new row.
            h(1:iband) = zero


      !  fill in the non-zero elements of the row. jrot records the column
      !  number of the first non-zero element in the row.
            h(1:ky2) = by(ii,1:ky2)*pinv
            zi   = zero
            jrot = (j-1)*nk1y+ii
      !  rotate the new row into triangle by givens transformations without
      !  square roots.
            do 540 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if (piv==zero) then
                if (i2<=0) go to 551
                go to 520
              endif
      !  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
      !  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2==0) go to 551
      !  apply that givens transformation to the left hand side.
              do 510 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 510          continue
 520          h(1:i2+1) = [h(2:i2+1),zero]
 540        continue
 551      continue
 550    continue
 560    if(nk1x==kx1) go to 640
      !  extend the observation matrix with the rows of a matrix expressing
      !  that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
        do 630 i=kx2,nk1x
          ii = i-kx1
          do 631 j=1,nk1y
      !  initialize the new row
            h(1:iband4) = zero
      !  fill in the non-zero elements of the row. jrot records the column
      !  number of the first non-zero element in the row.
            j1 = 1
            do 580 l=1,kx2
              h(j1) = bx(ii,l)*pinv
              j1 = j1+nk1y
 580        continue
            zi = zero
            jrot = (i-kx2)*nk1y+j
      !  rotate the new row into triangle by givens transformations .
            do 620 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if (piv==zero) then
                if (i2<=0) go to 631
                go to 600
              endif
      !  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
      !  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2==0) go to 631
      !  apply that givens transformation to the left hand side.
              do 590 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 590          continue
 600          h(1:i2+1) = [h(2:i2+1),zero]
 620        continue
 631      continue
 630    continue
      !  find dmax, the maximum value for the diagonal elements in the
      !  reduced triangle.
 640    dmax = max(zero,maxval(q(1:ncof,1),1))
      !  check whether the matrix is rank deficient.
        sigma = max(eps*dmax,maxval(q(1:ncof,1),1))
      !  backward substitution in case of full rank.
        c(:ncof) = fpback(q,ff,ncof,iband4,nc)
        rank = ncof
        go to 675
      !  in case of rank deficiency, find the minimum norm solution.
        lwest = ncof*iband4+ncof+iband4
        if(lwrk<lwest) go to 780
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,nc,sigma,c,sq,rank,wrk(la), &
         wrk(lf),wrk(lh))
 675    q(1:ncof,1) = q(1:ncof,1)/dmax

        !  compute f(p).
        fp = zero
        do 720 num = 1,nreg
          num1 = num-1
          lx = num1/nyy
          ly = num1-lx*nyy
          jrot = lx*nk1y+ly
          in = index(num)
 690      if(in==0) go to 720
          store = zero
          i1 = jrot
          do 710 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 700 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 700        continue
            i1 = i1+nk1y
 710      continue
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
          go to 690
 720    continue
      !  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms)<=acc) go to 820
      !  test whether the maximum allowable number of iterations has been
      !  reached.
        if(iter==maxit) go to 795
      !  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3/=0) go to 740
        if((f2-f3)>acc) go to 730
      !  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 + p2*con1
        go to 770
 730    if(f2<zero) ich3 = 1
 740    if(ich1/=0) go to 760
        if((f1-f2)>acc) go to 750
      !  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3<zero) go to 770
        if(p>=p3) p = p2*con1 + p3*con9
        go to 770
 750    if(f2>zero) ich1 = 1
      !  test whether the iteration process proceeds as theoretically
      !  expected.
 760    if(f2>=f1 .or. f2<=f3) go to 800
      !  find the new value of p.
        call fprati(p1,f1,p2,f2,p3,f3,p)
 770  continue
      !  error codes and messages.
 780  ier = lwest
      go to 830
 785  ier = 5
      go to 830
 790  ier = 4
      go to 830
 795  ier = 3
      go to 830
 800  ier = 2
      go to 830
 810  ier = 1
      go to 830
 815  ier = -1
      fp = zero
 820  if(ncof/=rank) ier = -rank
      !  test whether x and y are in the original order.
 830  if(ichang<0) go to 930
      !  if not, interchange x and y once more.
      l1 = 1
      do i=1,nk1x
        l2 = i
        do j=1,nk1y
          f(l2) = c(l1)
          l1 = l1+1
          l2 = l2+nk1x
        end do
      end do
      c(1:ncof) = f(1:ncof)
      do i=1,m
        store = x(i)
        x(i) = y(i)
        y(i) = store
      end do
      n = min0(nx,ny)
      do i=1,n
        store = tx(i)
        tx(i) = ty(i)
        ty(i) = store
      end do
      n1 = n+1
      if (nx<ny) go to 880
      if (nx==ny) go to 920
      go to 900
 880  tx(n1:ny) = ty(n1:ny)
      go to 920
 900  ty(n1:nx) = tx(n1:nx)
 920  l = nx ! Swap nx,ny
      nx = ny
      ny = l
 930  if(iopt<0) go to 940
      nx0 = nx
      ny0 = ny
 940  return
      end subroutine fpsurf



      ! subroutine fpsysy solves a linear n x n symmetric system   (a) * (b) = (g), with n<=6
      ! on input, vector g contains the right hand side ; on output it will contain the solution (b).
      pure subroutine fpsysy(a,n,g)

      !  ..scalar arguments..
      integer,     intent(in)    ::  n
      !  ..array arguments..
      real(RKIND), intent(inout) :: a(6,6)
      real(RKIND), intent(inout) :: g(6)

      !  ..local scalars..
      real(RKIND) :: fac
      integer :: i,i1,j,k
      !  ..
      g(1) = g(1)/a(1,1)
      if (n==1) return
      !  decomposition of the symmetric matrix (a) = (l) * (d) *(l)'  with (l) a unit lower triangular
      !  matrix and (d) a diagonal matrix
      a(2:n,1) = a(2:n,1)/a(1,1)

      do i=2,n
         i1 = i-1
         do k=i,n
            fac = a(k,i)
            do j=1,i1
               fac = fac-a(j,j)*a(k,j)*a(i,j)
            end do
            a(k,i) = fac
            if (k>i) a(k,i) = fac/a(i,i)
         end do
      end do
      !  solve the system (l)*(d)*(l)'*(b) = (g).

      !  first step : solve (l)*(d)*(c) = (g).
      do i=2,n
         i1 = i-1
         fac = g(i)
         do j=1,i1
            fac = fac-g(j)*a(j,j)*a(i,j)
         end do
         g(i) = fac/a(i,i)
      end do

      !  second step : solve (l)'*(b) = (c)
      i = n
      do j=2,n
         i1 = i
         i = i-1
         fac = g(i)
         do k=i1,n
            fac = fac-g(k)*a(k,i)
         end do
         g(i) = fac
      end do
      return
      end subroutine fpsysy


      !  subroutine fptrnp reduces the (m+n-7) x (n-4) matrix a to upper triangular form and applies the
      !  same givens transformations to the (m) x (mm) x (idim) matrix z to obtain the (n-4) x (mm) x
      !  (idim) matrix q
      pure subroutine fptrnp(m,mm,idim,n,nr,sp,p,b,z,a,q,right)

      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(in) :: p
      integer,     intent(in) :: m,mm,idim,n
      !  ..array arguments..
      real(RKIND), intent(in)  :: sp(m,4),b(n,5),z(m*mm*idim)
      real(RKIND), intent(out) :: q((n-4)*mm*idim),a(n,5),right(mm*idim)
      integer,     intent(in)  :: nr(m)

      !  ..local scalars..
      real(RKIND) :: cos,pinv,piv,sin
      integer     :: i,iband,irot,it,ii,i2,i3,j,jj,l,mid,nmd,m2,m3,nrold,n4,number,n1
      !  ..local arrays..
      real(RKIND) :: h(SIZ_K+1)
      !  ..subroutine references..
      !    fpgivs,fprota
      !  ..
      pinv = merge(one/p,one,p>zero)
      n4   = n-4
      mid  = mm*idim
      m2   = m*mm
      m3   = n4*mm
      !  reduce the matrix (a) to upper triangular form (r) using givens
      !  rotations. apply the same transformations to the rows of matrix z
      !  to obtain the mm x (n-4) matrix g.
      !  store matrix (r) into (a) and g into q.
      !  initialization.
      nmd         = n4*mid
      q(1:nmd)    = zero
      a(1:n4,1:5) = zero


      !  iband denotes the bandwidth of the matrices (a) and (r).
      iband = 4
      nrold  = 0
      number = nr(1)
      rows: do it=1,m
         number = nr(it)
         increase_old: do ! 150
             if (nrold/=number) then

                 if (p<=zero) then
                    nrold = nrold+1
                    cycle increase_old
                 end if
                 iband = 5
                 ! fetch a new row of matrix (b).
                 n1 = nrold+1
                 h(1:5) = b(n1,1:5)*pinv
                 !  find the appropriate column of q.
                 right(1:mid) = zero
                 irot = nrold

             else
                 ! fetch a new row of matrix (sp).
                 h(iband) = zero
                 h(1:4) = sp(it,1:4)
                 ! find the appropriate column of q.
                 j = 0
                 do ii=1,idim
                   l = (ii-1)*m2+(it-1)*mm
                   do jj=1,mm
                     j = j+1
                     l = l+1
                     right(j) = z(l)
                   end do
                 end do
                 irot = number
             endif

             ! rotate the new row of matrix (a) into triangle.
             rotate: do i=1,iband
                irot = irot+1
                piv  = h(i); if (piv==zero) cycle rotate

                ! calculate the parameters of the givens transformation.
                call fpgivs(piv,a(irot,1),cos,sin)

                ! apply that transformation to the rows of matrix q.
                j = 0
                do ii=1,idim
                   l = (ii-1)*m3+irot
                   do jj=1,mm
                      j = j+1
                      call fprota(cos,sin,right(j),q(l))
                      l = l+n4
                   end do
                end do

                ! apply that transformation to the columns of (a).
                if (i<iband) then
                    i2 = 1
                    i3 = i+1
                    do j=i3,iband
                       i2 = i2+1
                       call fprota(cos,sin,h(j),a(irot,i2))
                    end do
                endif
             end do rotate

             if (nrold==number) cycle rows
             nrold = nrold+1
         end do increase_old
      end do rows
      return
      end subroutine fptrnp


      !  subroutine fptrpe reduces the (m+n-7) x (n-7) cyclic bandmatrix a to upper triangular form and
      !  applies the same givens transformations to the (m) x (mm) x (idim) matrix z to obtain the (n-7) x
      !  (mm) x (idim) matrix q.
      recursive subroutine fptrpe(m,mm,idim,n,nr,sp,p,b,z,a,aa,q,right)


      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(in) :: p
      integer, intent(in) :: m,mm,idim,n
      !  ..array arguments..
      real(RKIND) sp(m,4),b(n,5),z(m*mm*idim),a(n,5),aa(n,4),q((n-7)*mm*idim), right(mm*idim)
      integer nr(m)
      !  ..local scalars..
      real(RKIND) co,pinv,piv,si
      integer i,irot,it,ii,i2,i3,j,jj,l,mid,nmd,m2,m3,nrold,n4,number,n1,n7,n11,m1
      integer i1, ij,j1,jk,jper,l0,l1, ik
      !  ..local arrays..
      real(RKIND) h(5),h1(5),h2(4)
      !  ..subroutine references..
      !    fpgivs,fprota
      !  ..
      pinv = merge(one/p,one,p>zero)
      n4 = n-4
      n7 = n-7
      n11 = n-11
      mid = mm*idim
      m2 = m*mm
      m3 = n7*mm
      m1 = m-1
      !  we determine the matrix (a) and then we reduce her to upper triangular form (r) using givens
      !  rotations. we apply the same transformations to the rows of matrix z to obtain the (mm) x (n-7)
      !  matrix g. we store matrix (r) into a and aa, g into q. the n7 x n7 upper triangular matrix (r)
      !  has the form
      !             | a1 '     |
      !       (r) = |    ' a2  |
      !             |  0 '     |
      !  with (a2) a n7 x 4 matrix and (a1) a n11 x n11 upper triangular matrix of bandwidth 5.

      !  initialization.
      nmd = n7*mid
      q (1:nmd)    = zero
      aa(1:n4,1:4) = zero
      a (1:n4,1:5) = zero
      jper = 0
      nrold = 0
      do 760 it=1,m1
        number = nr(it)
 120    if(nrold==number) go to 180
        if(p<=zero) go to 740
      !  fetch a new row of matrix (b).
        n1 = nrold+1
        h  = b(n1,:)*pinv
      !  find the appropriate row of q.
        right(:mid) = zero
        go to 240
      !  fetch a new row of matrix (sp)
 180    h = [sp(it,:),zero]
      !  find the appropriate row of q.
        j = 0
        do 220 ii=1,idim
          l = (ii-1)*m2+(it-1)*mm
          do jj=1,mm
            j = j+1
            l = l+1
            right(j) = z(l)
          end do
 220    continue
      !  test whether there are non-zero values in the new row of (a)
      !  corresponding to the b-splines n(j,*),j=n7+1,...,n4.
 240     if(nrold<n11) go to 640
         if(jper/=0) go to 320
      !  initialize the matrix (aa).
         jk = n11+1
         do 300 i=1,4
            ik = jk
            do 260 j=1,5
               if(ik<=0) go to 280
               aa(ik,i) = a(ik,j)
               ik = ik-1
 260        continue
 280        jk = jk+1
 300     continue
         jper = 1
      !  if one of the non-zero elements of the new row corresponds to one of
      !  the b-splines n(j;*),j=n7+1,...,n4,we take account of the periodicity
      !  conditions for setting up this row of (a).
 320     h1 = zero
         h2 = zero
         j = nrold-n11
         do 420 i=1,5
            j = j+1
            l0 = j
 360        l1 = l0-4
            if(l1<=0) go to 400
            if(l1<=n11) go to 380
            l0 = l1-n11
            go to 360
 380        h1(l1) = h(i)
            go to 420
 400        h2(l0) = h2(l0) + h(i)
 420     continue
      !  rotate the new row of (a) into triangle.
         if(n11<=0) go to 560
      !  rotations with the rows 1,2,...,n11 of (a).
         do 540 irot=1,n11
            piv = h1(1)
            i2 = min0(n11-irot,4)
            if (piv==zero) go to 500
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
      !  apply that transformation to the columns of matrix q.
            j = 0
            do 440 ii=1,idim
               l = (ii-1)*m3+irot
               do jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
               end do
 440        continue
      !  apply that transformation to the rows of (a) with respect to aa.
            do i=1,4
               call fprota(co,si,h2(i),aa(irot,i))
            end do
      !  apply that transformation to the rows of (a) with respect to a.
            if(i2==0) go to 560
            do 480 i=1,i2
               i1 = i+1
               call fprota(co,si,h1(i1),a(irot,i1))
 480        continue
 500        do 520 i=1,i2
               h1(i) = h1(i+1)
 520        continue
            h1(i2+1) = zero
 540     continue
      !  rotations with the rows n11+1,...,n7 of a.
 560     do 620 irot=1,4
            ij = n11+irot
            if(ij<=0) go to 620
            piv = h2(irot)
            if (piv==zero) go to 620
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,aa(ij,irot),co,si)
      !  apply that transformation to the columns of matrix q.
            j = 0
            do 580 ii=1,idim
               l = (ii-1)*m3+ij
               do jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
               end do
 580        continue
            if(irot==4) go to 620
      !  apply that transformation to the rows of (a) with respect to aa.
            j1 = irot+1
            do 600 i=j1,4
               call fprota(co,si,h2(i),aa(ij,i))
 600        continue
 620     continue
         go to 720
      !  rotation into triangle of the new row of (a), in case the elements
      !  corresponding to the b-splines n(j;*),j=n7+1,...,n4 are all zero.
 640     irot =nrold
         do 700 i=1,5
            irot = irot+1
            piv = h(i)
            if (piv==zero) go to 700
      !  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
      !  apply that transformation to the columns of matrix g.
            j = 0
            do ii=1,idim
               l = (ii-1)*m3+irot
               do jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
               end do
            end do
      !  apply that transformation to the rows of (a).
            if(i==5) go to 700
            i2 = 1
            i3 = i+1
            do 680 j=i3,5
               i2 = i2+1
               call fprota(co,si,h(j),a(irot,i2))
 680        continue
 700     continue
 720     if(nrold==number) go to 760
 740     nrold = nrold+1
         go to 120
 760  continue
      return
      end subroutine fptrpe


      !  subroutine insert inserts a new knot x into a spline function s(x) of degree k and calculates
      !  the b-spline representation of s(x) with respect to the new set of knots. in addition, if
      !  iopt/=0, s(x) will be considered as a periodic spline with period per=t(n-k)-t(k+1) satisfying
      !  the boundary constraints
      !       t(i+n-2*k-1) = t(i)+per  ,i=1,2,...,2*k+1
      !       c(i+n-2*k-1) = c(i)      ,i=1,2,...,k
      !  in that case, the knots and b-spline coefficients returned will also satisfy the periodic BCs, i.e.
      !       tt(i+nn-2*k-1) = tt(i)+per  ,i=1,2,...,2*k+1
      !       cc(i+nn-2*k-1) = cc(i)      ,i=1,2,...,k

      pure subroutine insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)

      !
      !  calling sequence:
      !     call insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)
      !
      !  input parameters:
      !    iopt : integer flag, specifying whether (iopt/=0) or not (iopt=0) the given spline must be
      !           considered as being periodic.
      !    t    : array,length nest, which contains the position of the knots.
      !    n    : integer, giving the total number of knots of s(x).
      !    c    : array,length nest, which contains the b-spline coefficients.
      !    k    : integer, giving the degree of s(x).
      !    x    : real, which gives the location of the knot to be inserted.
      !    nest : integer specifying the dimension of the arrays t,c,tt and cc. nest > n.
      !
      !  output parameters:
      !    tt   : array,length nest, which contains the position of the knots after insertion.
      !    nn   : integer, giving the total number of knots after insertion
      !    cc   : array,length nest, which contains the b-spline coefficients of s(x) with respect to the
      !           new set of knots.
      !    ier  : error flag
      !      ier = 0 : normal return
      !      ier =10 : invalid input data (see restrictions)
      !
      !  restrictions:
      !    nest > n
      !    t(k+1) <= x <= t(n-k)
      !    in case of a periodic spline (iopt/=0) there must be either at least k interior knots t(j)
      !       satisfying t(k+1)<t(j)<=x or at least k interior knots t(j) satisfying x<=t(j)<t(n-k)
      !
      !  other subroutines required: fpinst.
      !
      !  further comments:
      !   subroutine insert may be called as follows
      !        call insert(iopt,t,n,c,k,x,t,n,c,nest,ier)
      !   in which case the new representation will simply replace the old one
      !
      !  references :
      !    boehm w : inserting new knots into b-spline curves. computer aided design 12 (1980) 199-201.
      !   dierckx p. : curve and surface fitting with splines, monographs on numerical analysis, oxford
      !                university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer, intent(in) :: iopt,n,k,nest
      integer, intent(out) :: nn,ier
      real(RKIND), intent(in) :: x
      !  ..array arguments..
      real(RKIND), intent(in) :: t(nest),c(nest)
      real(RKIND), intent(out) :: tt(nest),cc(nest)
      !  ..local scalars..
      integer :: kk,k1,l,nk
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (nest<=n) return
      k1 = k+1
      nk = n-k
      if (x<t(k1) .or. x>t(nk)) return
      !  search for knot interval t(l) <= x < t(l+1).
      l = k1
      do while (l<nk .and. .not.x<t(l+1))
         l = l+1
      end do
      !  if no interval found above, then reverse the search and
      !  look for knot interval t(l) < x <= t(l+1).
      l = nk-1
      do while (l>k .and. .not.x>t(l))
        l = l-1
      end do

      !  no interval found in whole range
      if(t(l)>=t(l+1)) return

      if(iopt/=0) then
         kk = 2*k
         if (l<=kk .and. l>=(n-kk)) return
      endif

      ier = FITPACK_OK
      !  insert the new knot.
      call fpinst(iopt,t,n,c,k,x,l,tt,nn,cc,nest)

      end subroutine insert


      recursive subroutine parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)

      !  given the ordered set of m points x(i) in the idim-dimensional space and given also a corresponding
      !  set of strictly increasing values u(i) and the set of positive numbers w(i),i=1,2,...,m, subroutine
      !  parcur determines a smooth approximating spline curve s(u), i.e.
      !      x1 = s1(u)
      !      x2 = s2(u)       ub <= u <= ue
      !      .........
      !      xidim = sidim(u)
      !  with sj(u),j=1,2,...,idim spline functions of degree k with common knots t(j),j=1,2,...,n.
      !  if ipar=1 the values ub,ue and u(i),i=1,2,...,m must be supplied by the user. if ipar=0 these values
      !  are chosen automatically by parcur as
      !      v(1) = 0
      !      v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
      !      u(i) = v(i)/v(m) ,i=1,2,...,m
      !      ub = u(1) = 0, ue = u(m) = 1.
      !  if iopt=-1 parcur calculates the weighted least-squares spline according to a given set of knots.
      !  if iopt>=0 the number of knots of the splines sj(u) and the position t(j),j=1,2,...,n is chosen
      !  automatically by the routine. the smoothness of s(u) is then achieved by minimalizing the
      !  discontinuity jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,n-k-1. the amount
      !  of smoothness is determined by the condition that f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s,
      !  with s a given non-negative constant, called the smoothing factor. the fit s(u) is given in the
      !  b-spline representation and can be evaluated by means of subroutine curev.
      !
      !  calling sequence:
      !     call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)
      !
      !  parameters:
      !   iopt  : integer flag. on entry iopt must specify whether a weighted least-squares spline curve
      !           (iopt=-1) or a smoothing spline curve (iopt=0 or 1) must be determined.if iopt=0 the routine
      !           will start with an initial set of knots t(i)=ub,t(i+k+1)=ue, i=1,2,...,k+1. if iopt=1
      !           the routine will continue with the knots found at the last call of the routine. attention:
      !           a call with iopt=1 must always be immediately preceded by another call with iopt=1 or iopt=0.
      !           unchanged on exit.
      !   ipar  : integer flag. on entry ipar must specify whether (ipar=1) the user will supply the parameter
      !           values u(i),ub and ue or whether (ipar=0) these values are to be calculated by parcur.
      !           unchanged on exit.
      !   idim  : integer. on entry idim must specify the dimension of the curve. 0 < idim < 11.
      !           unchanged on exit.
      !   m     : integer. on entry m must specify the number of data points. m > k. unchanged on exit.
      !   u     : real array of dimension at least (m). in case ipar=1,before entry, u(i) must be set to the
      !           i-th value of the parameter variable u for i=1,2,...,m. these values must then be supplied
      !           in strictly ascending order and will be unchanged on exit. in case ipar=0, on exit,array
      !           u will contain the values u(i) as determined by parcur.
      !   mx    : integer. on entry mx must specify the actual dimension of the array x as declared in the
      !           calling (sub)program. mx must not be too small (see x). unchanged on exit.
      !   x     : real array of dimension at least idim*m.
      !           before entry, x(idim*(i-1)+j) must contain the j-th coordinate of the i-th data point for
      !           i=1,2,...,m and j=1,2,...,idim. unchanged on exit.
      !   w     : real array of dimension at least (m). before entry, w(i) must be set to the i-th value in
      !           the set of weights. the w(i) must be strictly positive. unchanged on exit. see also further
      !           comments.
      !   ub,ue : real values. on entry (in case ipar=1) ub and ue must contain the lower and upper bound for
      !           the parameter u. ub <=u(1), ue>= u(m). if ipar = 0 these values will automatically be set
      !           to 0 and 1 by parcur.
      !   k     : integer. on entry k must specify the degree of the splines. 1<=k<=5. it is recommended to
      !           use cubic splines (k=3). the user is strongly dissuaded from choosing k even,together
      !           with a small s-value. unchanged on exit.
      !   s     : real.on entry (in case iopt>=0) s must specify the smoothing factor. s >=0. unchanged on exit.
      !           for advice on the choice of s see further comments.
      !   nest  : integer. on entry nest must contain an over-estimate of the total number of knots of the
      !           splines returned, to indicate the storage space available to the routine. nest >=2*k+2.
      !           in most practical situation nest=m/2 will be sufficient. always large enough is nest=m+k+1,
      !           the number of knots needed for interpolation (s=0). unchanged on exit.
      !   n     : integer.
      !           unless ier = 10 (in case iopt >=0), n will contain the total number of knots of the
      !           smoothing spline curve returned if the computation mode iopt=1 is used this value of n
      !           should be left unchanged between subsequent calls. in case iopt=-1, the value of n must be
      !           specified on entry.
      !   t     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the knots of the spline curve,i.e. the position
      !           of the interior knots t(k+2),t(k+3),..,t(n-k-1) as well as the position of the additional
      !           t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for the b-spline representation.
      !           if the computation mode iopt=1 is used, the values of t(1),t(2),...,t(n) should be left
      !           unchanged between subsequent calls. if the computation mode iopt=-1 is used, the values
      !           t(k+2),...,t(n-k-1) must be supplied by the user, before entry. see also the restrictions
      !           (ier=10).
      !   nc    : integer. on entry nc must specify the actual dimension of the array c as declared in the
      !           calling (sub)program. nc must not be too small (see c). unchanged on exit.
      !   c     : real array of dimension at least (nest*idim). on successful exit, this array will contain
      !           the coefficients in the b-spline representation of the spline curve s(u),i.e. the b-spline
      !           coefficients of the spline sj(u) will be given in c(n*(j-1)+i),i=1,2,...,n-k-1 for
      !           j=1,2,...,idim.
      !   fp    : real. unless ier = 10, fp contains the weighted sum of squared residuals of the spline
      !           curve returned.
      !   wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k). used as working space. if the
      !           computation mode iopt=1 is used, the values wrk(1),...,wrk(n) should be left unchanged
      !           between subsequent calls.
      !   lwrk  : integer. on entry,lwrk must specify the actual dimension of the array wrk as declared in
      !           the calling (sub)program. lwrk must not be too small (see wrk). unchanged on exit.
      !   iwrk  : integer array of dimension at least (nest). used as working space. if the computation
      !           mode iopt=1 is used,the values iwrk(1),...,iwrk(n) should be left unchanged between
      !           subsequent calls.
      !   ier   : integer. unless the routine detects an error, ier contains a non-positive value on exit, i.e.
      !    ier=0  : normal return. the curve returned has a residual sum of squares fp such that abs(fp-s)/s
      !             <= tol with tol a relative tolerance set to 0.001 by the program.
      !    ier=-1 : normal return. the curve returned is an interpolating spline curve (fp=0).
      !    ier=-2 : normal return. the curve returned is the weighted least squares polynomial curve of
      !             degree k. in this extreme case fp gives the upper bound fp0 for the smoothing factor s.
      !    ier=1  : error. the required storage space exceeds the available storage space, as specified by the
      !             parameter nest.
      !             likely causes : nest too small. if nest is already large (say nest > m/2), it may also
      !             indicate that s is too small. the approximation returned is the least-squares spline
      !             curve according to the knots t(1),t(2),...,t(n). (n=nest) the parameter fp gives the
      !             corresponding weighted sum of squared residuals (fp>s).
      !    ier=2  : error. a theoretically impossible result was found during the iteration process for
      !             finding a smoothing spline curve with fp = s. probably causes : s too small.
      !             there is an approximation returned but the corresponding weighted sum of squared residuals
      !             does not satisfy the condition abs(fp-s)/s < tol.
      !    ier=3  : error. the maximal number of iterations maxit (set to 20 by the program) allowed for
      !             finding a smoothing curve with fp=s has been reached. probably causes : s too small
      !             there is an approximation returned but the corresponding weighted sum of squared residuals
      !             does not satisfy the condition abs(fp-s)/s < tol.
      !    ier=10 : error. on entry, the input data are controlled on validity the following restrictions
      !             must be satisfied.
      !             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
      !             0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
      !             nc>=nest*idim
      !             if ipar=0: sum j=1,idim (x(idim*i+j)-x(idim*(i-1)+j))**2>0 i=1,2,...,m-1.
      !             if ipar=1: ub<=u(1)<u(2)<...<u(m)<=ue
      !             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
      !                         ub<t(k+2)<t(k+3)<...<t(n-k-1)<ue
      !                            (ub=0 and ue=1 in case ipar=0)
      !                       the schoenberg-whitney conditions, i.e. there must be a subset of data points
      !                       uu(j) such that
      !                         t(j) < uu(j) < t(j+k+1), j=1,2,...,n-k-1
      !             if iopt>=0: s>=0
      !                         if s=0 : nest >= m+k+1
      !             if one of these conditions is found to be violated,control is immediately repassed to the
      !             calling program. in that case there is no approximation returned.
      !
      !  further comments:
      !   by means of the parameter s, the user can control the tradeoff between closeness of fit and
      !   smoothness of fit of the approximation. if s is too large, the curve will be too smooth and signal
      !   will be lost ; if s is too small the curve will pick up too much noise. in the extreme cases the
      !   program will return an interpolating curve if s=0 and the least-squares polynomial curve of degree
      !   k if s is very large. between these extremes, a properly chosen s will result in a good compromise
      !   between closeness of fit and smoothness of fit. to decide whether an approximation, corresponding
      !   to a certain s is satisfactory the user is highly recommended to inspect the fits graphically.
      !   recommended values for s depend on the weights w(i). if these are taken as 1/d(i) with d(i) an
      !   estimate of the standard deviation of x(i), a good s-value should be found in the range
      !   (m-sqrt(2*m),m+sqrt(2*m)). if nothing is known about the statistical error in x(i) each w(i) can be
      !   set equal to one and s determined by trial and error, taking account of the comments above. the best
      !   is then to start with a very large value of s ( to determine the least-squares polynomial curve and
      !   the upper bound fp0 for s) and then to progressively decrease the value of s ( say by a factor 10
      !   in the beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the approximating curve shows
      !   more detail) to obtain closer fits. to economize the search for a good s-value the program provides
      !   with different modes of computation. at the first call of the routine, or whenever he wants to
      !   restart with the initial set of knots the user must set iopt=0. if iopt=1 the program will continue
      !   with the set of knots found at the last call of the routine. this will save a lot of computation
      !   time if parcur is called repeatedly for different values of s. the number of knots of the spline
      !   returned and their location will depend on the value of s and on the complexity of the shape of the
      !   curve underlying the data. but, if the computation mode iopt=1 is used, the knots returned may also
      !   depend on the s-values at previous calls (if these were smaller). therefore, if after a number of
      !   trials with different s-values and iopt=1, the user can finally accept a fit as satisfactory, it
      !   may be worthwhile for him to call parcur once more with the selected value for s but now with iopt=0.
      !   indeed, parcur may then return an approximation of the same quality of fit but with fewer knots and
      !   therefore better if data reduction is also an important objective for the user.
      !
      !   the form of the approximating curve can strongly be affected by the choice of the parameter values
      !   u(i). if there is no physical reason for choosing a particular parameter u, often good results will
      !   be obtained with the choice of parcur (in case ipar=0), i.e.
      !        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
      !   where
      !        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
      !   other possibilities for q(i) are
      !        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
      !        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
      !        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
      !        q(i)= 1
      !
      !  other subroutines required:
      !    fpback,fpbspl,fpchec,fppara,fpdisc,fpgivs,fpknot,fprati,fprota
      !
      !  references:
      !   dierckx p. : algorithms for smoothing data with periodic and parametric splines, computer graphics
      !                and image processing 20 (1982) 171-184.
      !   dierckx p. : algorithms for smoothing data with periodic and parametric splines, report tw55, dept.
      !                computer science, k.u.leuven, 1981.
      !   dierckx p. : curve and surface fitting with splines, monographs on numerical analysis, oxford
      !                university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !
      !  ..scalar arguments..
      real(RKIND), intent(inout) :: ub,ue,s
      real(RKIND), intent(out) :: fp
      integer, intent(in) :: iopt,ipar,idim,m,mx,k,nest,n,lwrk,nc
      integer, intent(out) :: ier
      !  ..array arguments..
      real(RKIND), intent(in) :: x(idim,m)
      real(RKIND), intent(inout) :: u(m),w(m),t(nest),c(nc),wrk(lwrk)
      integer, intent(inout) :: iwrk(nest)
      !  ..local scalars..
      integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,nmin,ncc


      !  we set up the parameters tol and maxit
      integer, parameter :: maxit = 20
      real(RKIND), parameter :: tol = smallnum03

      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (iopt<(-1) .or. iopt>1)       return
      if (ipar<0 .or. ipar>1)          return
      if (idim<=0 .or. idim>SIZ_IDIM)  return
      if (k<=0 .or. k>5)               return

      k1 = k+1
      k2 = k1+1
      nmin = 2*k1
      if (m<k1 .or. nest<nmin)         return
      ncc = nest*idim
      if (mx<m*idim .or. nc<ncc)       return
      lwest = m*k1+nest*(6+idim+3*k)
      if (lwrk<lwest)                  return

      ! Normalize coordinates
      if (ipar==0 .and. iopt<=0) then

          ! Point coordinates are stored in x(:), offset by idim values
          u(1) = zero
          do i=2,m
             u(i) = u(i-1) + norm2(x(:,i)-x(:,i-1))
          end do
          if (u(m)<=zero) return

          u(2:) = u(2:)/u(m)
          ub    = zero
          ue    = one
          u(m)  = ue
      endif

      if (ub>u(1) .or. ue<u(m) .or. w(1)<=zero) return
      if (any(u(2:)>=u(:m-1) .or. w(2:)<=zero)) return
      if (iopt<0) then
          if (n<nmin .or. n>nest) return
          j = n
          do i=1,k1
             t(i) = ub
             t(j) = ue
             j = j-1
          end do
          call fpchec(u,m,t,n,k,ier); if (ier/=FITPACK_OK) return

      else

          if (s<zero) return
          if (s==zero .and. nest<(m+k1)) return

          ier = FITPACK_OK

      endif

      ! we partition the working space and determine the spline curve.
      ifp = 1
      iz = ifp+nest
      ia = iz+ncc
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
         call fppara(iopt,idim,m,u,mx,x,w,ub,ue,k,s,nest,tol,maxit,k1,k2, &
                     n,t,ncc,c,fp,wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
      return
      end subroutine parcur


      recursive subroutine parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx, &
       y,my,z,wrk,lwrk,iwrk,kwrk,ier)

      !  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
      !  ,my the partial derivative ( order nux,nuy) of a bivariate spline
      !  s(x,y) of degrees kx and ky, given in the b-spline representation.
      !
      !  calling sequence:
      !     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
      !    * iwrk,kwrk,ier)
      !
      !  input parameters:
      !   tx    : real array, length nx, which contains the position of the
      !           knots in the x-direction.
      !   nx    : integer, giving the total number of knots in the x-direction
      !   ty    : real array, length ny, which contains the position of the
      !           knots in the y-direction.
      !   ny    : integer, giving the total number of knots in the y-direction
      !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
      !           b-spline coefficients.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   nux   : integer values, specifying the order of the partial
      !   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
      !   x     : real array of dimension (mx).
      !           before entry x(i) must be set to the x co-ordinate of the
      !           i-th grid point along the x-axis.
      !           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
      !   mx    : on entry mx must specify the number of grid points along
      !           the x-axis. mx >=1.
      !   y     : real array of dimension (my).
      !           before entry y(j) must be set to the y co-ordinate of the
      !           j-th grid point along the y-axis.
      !           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
      !   my    : on entry my must specify the number of grid points along
      !           the y-axis. my >=1.
      !   wrk   : real array of dimension lwrk. used as workspace.
      !   lwrk  : integer, specifying the dimension of wrk.
      !           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
      !   iwrk  : integer array of dimension kwrk. used as workspace.
      !   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
      !
      !  output parameters:
      !   z     : real array of dimension (mx*my).
      !           on successful exit z(my*(i-1)+j) contains the value of the
      !           specified partial derivative of s(x,y) at the point
      !           (x(i),y(j)),i=1,...,mx;j=1,...,my.
      !   ier   : integer error flag
      !    ier=0 : normal return
      !    ier=10: invalid input data (see restrictions)
      !
      !  restrictions:
      !   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my
      !   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
      !   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
      !   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
      !
      !  other subroutines required:
      !    fpbisp,fpbspl
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory
      !                6 (1972) 50-62.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1989
      !
      !  ..scalar arguments..
      integer nx,ny,kx,ky,nux,nuy,mx,my,lwrk,kwrk,ier
      !  ..array arguments..
      integer iwrk(kwrk)
      real(RKIND) tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my), &
       wrk(lwrk)
      !  ..local scalars..
      integer i,iwx,iwy,j,kkx,kky,kx1,ky1,lx,ly,lwest,l1,l2,m,m0,m1, &
       nc,nkx1,nky1,nxx,nyy
      real(RKIND) ak,fac
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      if(nux<0 .or. nux>=kx) go to 400
      if(nuy<0 .or. nuy>=ky) go to 400
      lwest = nc +(kx1-nux)*mx+(ky1-nuy)*my
      if(lwrk<lwest) go to 400
      if(kwrk<(mx+my)) go to 400
      if (mx<1) go to 400
      if (mx==1) go to 30
      go to 10
  10  if(any(x(2:mx)<x(1:mx-1))) go to 400
  30  if (my<1) go to 400
      if (my==1) go to 60
      go to 40
  40  if(any(y(2:my)<y(1:my-1))) go to 400
  60  ier = FITPACK_OK
      nxx = nkx1
      nyy = nky1
      kkx = kx
      kky = ky
      !  the partial derivative of order (nux,nuy) of a bivariate spline of
      !  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
      !  we calculate the b-spline coefficients of this spline
      wrk(1:nc) = c(1:nc)
      if(nux==0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = kkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+kkx
          fac = tx(l2)-tx(l1)
          if(fac<=0.) go to 90
          do 80 m=1,nyy
            m1 = m0+nyy
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        kkx = kkx-1
 100  continue
 200  if(nuy==0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = kky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+kky
          fac = ty(l2)-ty(l1)
          if(fac<=0.) go to 220
          m0 = i
          do 210 m=1,nxx
            m1 = m0+1
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        kky = kky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 m=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          wrk(m0) = wrk(m1)
 240    continue
        m1 = m1+nuy
 250  continue
      !  we partition the working space and evaluate the partial derivative
 300  iwx = 1+nxx*nyy
      iwy = iwx+mx*(kx1-nux)
      call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,wrk,kkx,kky, &
                  x,mx,y,my,z,wrk(iwx),wrk(iwy),iwrk(1),iwrk(mx+1))
 400  return
      end subroutine parder



      recursive subroutine pardeu(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,y,z,m, &
       wrk,lwrk,iwrk,kwrk,ier)

      !  subroutine pardeu evaluates on a set of points (x(i),y(i)),i=1,...,m
      !  the partial derivative ( order nux,nuy) of a bivariate spline
      !  s(x,y) of degrees kx and ky, given in the b-spline representation.
      !
      !  calling sequence:
      !     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
      !    * iwrk,kwrk,ier)
      !
      !  input parameters:
      !   tx    : real array, length nx, which contains the position of the
      !           knots in the x-direction.
      !   nx    : integer, giving the total number of knots in the x-direction
      !   ty    : real array, length ny, which contains the position of the
      !           knots in the y-direction.
      !   ny    : integer, giving the total number of knots in the y-direction
      !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
      !           b-spline coefficients.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   nux   : integer values, specifying the order of the partial
      !   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   x     : real array of dimension (mx).
      !   y     : real array of dimension (my).
      !   m     : on entry m must specify the number points. m >= 1.
      !   wrk   : real array of dimension lwrk. used as workspace.
      !   lwrk  : integer, specifying the dimension of wrk.
      !           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
      !   iwrk  : integer array of dimension kwrk. used as workspace.
      !   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
      !
      !  output parameters:
      !   z     : real array of dimension (m).
      !           on successful exit z(i) contains the value of the
      !           specified partial derivative of s(x,y) at the point
      !           (x(i),y(i)),i=1,...,m.
      !   ier   : integer error flag
      !   ier=0 : normal return
      !   ier=10: invalid input data (see restrictions)
      !
      !  restrictions:
      !   lwrk>=m*(kx+1-nux)+m*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
      !
      !  other subroutines required:
      !    fpbisp,fpbspl
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory
      !                6 (1972) 50-62.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1989
      !
      !  ..scalar arguments..
      integer nx,ny,kx,ky,m,lwrk,kwrk,ier,nux,nuy
      !  ..array arguments..
      integer iwrk(kwrk)
      real(RKIND) tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(m),y(m),z(m), &
       wrk(lwrk)
      !  ..local scalars..
      integer i,iwx,iwy,j,kkx,kky,kx1,ky1,lx,ly,lwest,l1,l2,mm,m0,m1, &
       nc,nkx1,nky1,nxx,nyy
      real(RKIND) ak,fac
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = 10
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      if(nux<0 .or. nux>=kx) go to 400
      if(nuy<0 .or. nuy>=ky) go to 400
      lwest = nc +(kx1-nux)*m+(ky1-nuy)*m
      if(lwrk<lwest) go to 400
      if(kwrk<(m+m)) go to 400
      if (m<1) go to 400
      ier = 0
      nxx = nkx1
      nyy = nky1
      kkx = kx
      kky = ky
      !  the partial derivative of order (nux,nuy) of a bivariate spline of
      !  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
      !  we calculate the b-spline coefficients of this spline
      do 70 i=1,nc
        wrk(i) = c(i)
  70  continue
      if(nux==0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = kkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+kkx
          fac = tx(l2)-tx(l1)
          if(fac<=0.) go to 90
          do 80 mm=1,nyy
            m1 = m0+nyy
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        kkx = kkx-1
 100  continue
 200  if(nuy==0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = kky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+kky
          fac = ty(l2)-ty(l1)
          if(fac<=0.) go to 220
          m0 = i
          do 210 mm=1,nxx
            m1 = m0+1
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        kky = kky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 mm=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          wrk(m0) = wrk(m1)
 240    continue
        m1 = m1+nuy
 250  continue
      !  we partition the working space and evaluate the partial derivative
 300  iwx = 1+nxx*nyy
      iwy = iwx+m*(kx1-nux)
      do 390 i=1,m
        call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,wrk,kkx,kky, &
         x(i),1,y(i),1,z(i),wrk(iwx),wrk(iwy),iwrk(1),iwrk(2))
 390  continue
 400  return
      end subroutine pardeu


      recursive subroutine pardtc(tx,nx,ty,ny,c,kx,ky,nux,nuy, &
         newc,ier)

      !  subroutine pardtc takes the knots and coefficients of a bivariate
      !  spline, and returns the coefficients for a new bivariate spline that
      !  evaluates the partial derivative (order nux, nuy) of the original
      !  spline.
      !
      !  calling sequence:
      !     call pardtc(tx,nx,ty,ny,c,kx,ky,nux,nuy,newc,ier)
      !
      !  input parameters:
      !   tx    : real array, length nx, which contains the position of the
      !           knots in the x-direction.
      !   nx    : integer, giving the total number of knots in the x-direction
      !           (hidden)
      !   ty    : real array, length ny, which contains the position of the
      !           knots in the y-direction.
      !   ny    : integer, giving the total number of knots in the y-direction
      !           (hidden)
      !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
      !           b-spline coefficients.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   nux   : integer values, specifying the order of the partial
      !   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
      !
      !  output parameters:
      !   newc  : real array containing the coefficients of the derivative.
      !           the dimension is (nx-nux-kx-1)*(ny-nuy-ky-1).
      !   ier   : integer error flag
      !    ier=0 : normal return
      !    ier=10: invalid input data (see restrictions)
      !
      !  restrictions:
      !   0 <= nux < kx, 0 <= nuy < kyc
      !
      !  other subroutines required:
      !    none
      !
      !  references :
      !   de boor c  : on calculating with b-splines, j. approximation theory
      !                6 (1972) 50-62.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  based on the subroutine "parder" by Paul Dierckx.
      !
      !  author :
      !    Cong Ma
      !    Department of Mathematics and Applied Mathematics, U. of Cape Town
      !    Cross Campus Road, Rondebosch 7700, Cape Town, South Africa.
      !    e-mail : cong.ma@uct.ac.za
      !
      !  latest update : may 2019
      !
      !  ..scalar arguments..
      integer nx,ny,kx,ky,nux,nuy,ier, nc
      !  ..array arguments..
      real(RKIND) tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)), &
       newc((nx-kx-1)*(ny-ky-1))
      !  ..local scalars..
      integer i,j,kx1,ky1,lx,ly,l1,l2,m,m0,m1, &
       nkx1,nky1,nxx,nyy,newkx,newky
      real(RKIND) ak,fac
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nux<0 .or. nux>=kx) go to 400
      if(nuy<0 .or. nuy>=ky) go to 400
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      ier = 0
      nxx = nkx1
      nyy = nky1
      newkx = kx
      newky = ky
      !  the partial derivative of order (nux,nuy) of a bivariate spline of
      !  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
      !  we calculate the b-spline coefficients of this spline
      !  that is to say newkx = kx - nux, newky = ky - nuy
      do 70 i=1,nc
        newc(i) = c(i)
  70  continue
      if(nux==0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = newkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+newkx
          fac = tx(l2)-tx(l1)
          if(fac<=0.) go to 90
          do 80 m=1,nyy
            m1 = m0+nyy
            newc(m0) = (newc(m1)-newc(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        newkx = newkx-1
 100  continue
 200  if(nuy==0) go to 400
      ! orig: if(nuy==0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = newky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+newky
          fac = ty(l2)-ty(l1)
          if(fac<=0.) go to 220
          m0 = i
          do 210 m=1,nxx
            m1 = m0+1
            newc(m0) = (newc(m1)-newc(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        newky = newky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 m=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          newc(m0) = newc(m1)
 240    continue
        m1 = m1+nuy
 250  continue
      !300  iwx = 1+nxx*nyy
      !     iwy = iwx+mx*(kx1-nux)
      !
      ! from parder.f:
      !     call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,newc,newkx,newky,
      !    * x,mx,y,my,z,newc(iwx),newc(iwy),iwrk(1),iwrk(mx+1))
      !
      ! from bispev.f:
      !     call fpbisp(tx,       nx,      ty,       ny,      c,   kx,   ky,
      !    * x,mx,y,my,z,wrk(1),   wrk(iw),  iwrk(1),iwrk(mx+1))
      !
      ! from fpbisp.f:
      !          fpbisp(tx,       nx,      ty,       ny,      c,   kx,   ky,
      !    * x,mx,y,my,z,wx,       wy,       lx,     ly)
 400  return
      end subroutine pardtc



      recursive subroutine parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest, &
       nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)

      !  given the set of ordered points f(i,j) in the idim-dimensional space,
      !  corresponding to grid values (u(i),v(j)) ,i=1,...,mu ; j=1,...,mv,
      !  parsur determines a smooth approximating spline surface s(u,v) , i.e.
      !    f1 = s1(u,v)
      !      ...                u(1) <= u <= u(mu) ; v(1) <= v <= v(mv)
      !    fidim = sidim(u,v)
      !  with sl(u,v), l=1,2,...,idim bicubic spline functions with common
      !  knots tu(i),i=1,...,nu in the u-variable and tv(j),j=1,...,nv in the
      !  v-variable.
      !  in addition, these splines will be periodic in the variable u if
      !  ipar(1) = 1 and periodic in the variable v if ipar(2) = 1.
      !  if iopt=-1, parsur determines the least-squares bicubic spline
      !  surface according to a given set of knots.
      !  if iopt>=0, the number of knots of s(u,v) and their position
      !  is chosen automatically by the routine. the smoothness of s(u,v) is
      !  achieved by minimalizing the discontinuity jumps of the derivatives
      !  of the splines at the knots. the amount of smoothness of s(u,v) is
      !  determined by the condition that
      !  fp=sumi=1,mu(sumj=1,mv(dist(f(i,j)-s(u(i),v(j)))**2))<=s,
      !  with s a given non-negative constant.
      !  the fit s(u,v) is given in its b-spline representation and can be
      !  evaluated by means of routine surev.
      !
      ! calling sequence:
      !     call parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,nu,tu,
      !    *  nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      !
      ! parameters:
      !  iopt  : integer flag. unchanged on exit.
      !          on entry iopt must specify whether a least-squares surface
      !          (iopt=-1) or a smoothing surface (iopt=0 or 1)must be
      !          determined.
      !          if iopt=0 the routine will start with the initial set of
      !          knots needed for determining the least-squares polynomial
      !          surface.
      !          if iopt=1 the routine will continue with the set of knots
      !          found at the last call of the routine.
      !          attention: a call with iopt=1 must always be immediately
      !          preceded by another call with iopt = 1 or iopt = 0.
      !  ipar  : integer array of dimension 2. unchanged on exit.
      !          on entry ipar(1) must specify whether (ipar(1)=1) or not
      !          (ipar(1)=0) the splines must be periodic in the variable u.
      !          on entry ipar(2) must specify whether (ipar(2)=1) or not
      !          (ipar(2)=0) the splines must be periodic in the variable v.
      !  idim  : integer. on entry idim must specify the dimension of the
      !          surface. 1 <= idim <= 3. unchanged on exit.
      !  mu    : integer. on entry mu must specify the number of grid points
      !          along the u-axis. unchanged on exit.
      !          mu >= mumin where mumin=4-2*ipar(1)
      !  u     : real array of dimension at least (mu). before entry, u(i)
      !          must be set to the u-co-ordinate of the i-th grid point
      !          along the u-axis, for i=1,2,...,mu. these values must be
      !          supplied in strictly ascending order. unchanged on exit.
      !  mv    : integer. on entry mv must specify the number of grid points
      !          along the v-axis. unchanged on exit.
      !          mv >= mvmin where mvmin=4-2*ipar(2)
      !  v     : real array of dimension at least (mv). before entry, v(j)
      !          must be set to the v-co-ordinate of the j-th grid point
      !          along the v-axis, for j=1,2,...,mv. these values must be
      !          supplied in strictly ascending order. unchanged on exit.
      !  f     : real array of dimension at least (mu*mv*idim).
      !          before entry, f(mu*mv*(l-1)+mv*(i-1)+j) must be set to the
      !          l-th co-ordinate of the data point corresponding to the
      !          the grid point (u(i),v(j)) for l=1,...,idim ,i=1,...,mu
      !          and j=1,...,mv. unchanged on exit.
      !          if ipar(1)=1 it is expected that f(mu*mv*(l-1)+mv*(mu-1)+j)
      !          = f(mu*mv*(l-1)+j), l=1,...,idim ; j=1,...,mv
      !          if ipar(2)=1 it is expected that f(mu*mv*(l-1)+mv*(i-1)+mv)
      !          = f(mu*mv*(l-1)+mv*(i-1)+1), l=1,...,idim ; i=1,...,mu
      !  s     : real. on entry (if iopt>=0) s must specify the smoothing
      !          factor. s >=0. unchanged on exit.
      !          for advice on the choice of s see further comments
      !  nuest : integer. unchanged on exit.
      !  nvest : integer. unchanged on exit.
      !          on entry, nuest and nvest must specify an upper bound for the
      !          number of knots required in the u- and v-directions respect.
      !          these numbers will also determine the storage space needed by
      !          the routine. nuest >= 8, nvest >= 8.
      !          in most practical situation nuest = mu/2, nvest=mv/2, will
      !          be sufficient. always large enough are nuest=mu+4+2*ipar(1),
      !          nvest = mv+4+2*ipar(2), the number of knots needed for
      !          interpolation (s=0). see also further comments.
      !  nu    : integer.
      !          unless ier=10 (in case iopt>=0), nu will contain the total
      !          number of knots with respect to the u-variable, of the spline
      !          surface returned. if the computation mode iopt=1 is used,
      !          the value of nu should be left unchanged between subsequent
      !          calls. in case iopt=-1, the value of nu should be specified
      !          on entry.
      !  tu    : real array of dimension at least (nuest).
      !          on successful exit, this array will contain the knots of the
      !          splines with respect to the u-variable, i.e. the position of
      !          the interior knots tu(5),...,tu(nu-4) as well as the position
      !          of the additional knots tu(1),...,tu(4) and tu(nu-3),...,
      !          tu(nu) needed for the b-spline representation.
      !          if the computation mode iopt=1 is used,the values of tu(1)
      !          ...,tu(nu) should be left unchanged between subsequent calls.
      !          if the computation mode iopt=-1 is used, the values tu(5),
      !          ...tu(nu-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  nv    : integer.
      !          unless ier=10 (in case iopt>=0), nv will contain the total
      !          number of knots with respect to the v-variable, of the spline
      !          surface returned. if the computation mode iopt=1 is used,
      !          the value of nv should be left unchanged between subsequent
      !          calls. in case iopt=-1, the value of nv should be specified
      !          on entry.
      !  tv    : real array of dimension at least (nvest).
      !          on successful exit, this array will contain the knots of the
      !          splines with respect to the v-variable, i.e. the position of
      !          the interior knots tv(5),...,tv(nv-4) as well as the position
      !          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
      !          tv(nv) needed for the b-spline representation.
      !          if the computation mode iopt=1 is used,the values of tv(1)
      !          ...,tv(nv) should be left unchanged between subsequent calls.
      !          if the computation mode iopt=-1 is used, the values tv(5),
      !          ...tv(nv-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  c     : real array of dimension at least (nuest-4)*(nvest-4)*idim.
      !          on successful exit, c contains the coefficients of the spline
      !          approximation s(u,v)
      !  fp    : real. unless ier=10, fp contains the sum of squared
      !          residuals of the spline surface returned.
      !  wrk   : real array of dimension (lwrk). used as workspace.
      !          if the computation mode iopt=1 is used the values of
      !          wrk(1),...,wrk(4) should be left unchanged between subsequent
      !          calls.
      !  lwrk  : integer. on entry lwrk must specify the actual dimension of
      !          the array wrk as declared in the calling (sub)program.
      !          lwrk must not be too small.
      !           lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+
      !           4*(mu+mv)+q*idim where q is the larger of mv and nuest.
      !  iwrk  : integer array of dimension (kwrk). used as workspace.
      !          if the computation mode iopt=1 is used the values of
      !          iwrk(1),.,iwrk(3) should be left unchanged between subsequent
      !          calls.
      !  kwrk  : integer. on entry kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program.
      !          kwrk >= 3+mu+mv+nuest+nvest.
      !  ier   : integer. unless the routine detects an error, ier contains a
      !          non-positive value on exit, i.e.
      !   ier=0  : normal return. the surface returned has a residual sum of
      !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !            ive tolerance set to 0.001 by the program.
      !   ier=-1 : normal return. the spline surface returned is an
      !            interpolating surface (fp=0).
      !   ier=-2 : normal return. the surface returned is the least-squares
      !            polynomial surface. in this extreme case fp gives the
      !            upper bound for the smoothing factor s.
      !   ier=1  : error. the required storage space exceeds the available
      !            storage space, as specified by the parameters nuest and
      !            nvest.
      !            probably causes : nuest or nvest too small. if these param-
      !            eters are already large, it may also indicate that s is
      !            too small
      !            the approximation returned is the least-squares surface
      !            according to the current set of knots. the parameter fp
      !            gives the corresponding sum of squared residuals (fp>s).
      !   ier=2  : error. a theoretically impossible result was found during
      !            the iteration process for finding a smoothing surface with
      !            fp = s. probably causes : s too small.
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=3  : error. the maximal number of iterations maxit (set to 20
      !            by the program) allowed for finding a smoothing surface
      !            with fp=s has been reached. probably causes : s too small
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=10 : error. on entry, the input data are controlled on validity
      !            the following restrictions must be satisfied.
      !            -1<=iopt<=1, 0<=ipar(1)<=1, 0<=ipar(2)<=1, 1 <=idim<=3
      !            mu >= 4-2*ipar(1),mv >= 4-2*ipar(2), nuest >=8, nvest >= 8,
      !            kwrk>=3+mu+mv+nuest+nvest,
      !            lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))
      !             +4*(mu+mv)+max(nuest,mv)*idim
      !            u(i-1)<u(i),i=2,..,mu, v(i-1)<v(i),i=2,...,mv
      !            if iopt=-1: 8<=nu<=min(nuest,mu+4+2*ipar(1))
      !                        u(1)<tu(5)<tu(6)<...<tu(nu-4)<u(mu)
      !                        8<=nv<=min(nvest,mv+4+2*ipar(2))
      !                        v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(mv)
      !                    the schoenberg-whitney conditions, i.e. there must
      !                    be subset of grid co-ordinates uu(p) and vv(q) such
      !                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
      !                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
      !                     (see fpchec or fpchep)
      !            if iopt>=0: s>=0
      !                       if s=0: nuest>=mu+4+2*ipar(1)
      !                               nvest>=mv+4+2*ipar(2)
      !            if one of these conditions is found to be violated,control
      !            is immediately repassed to the calling program. in that
      !            case there is no approximation returned.
      !
      ! further comments:
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the surface will be too smooth and signal will be
      !   lost ; if s is too small the surface will pick up too much noise. in
      !   the extreme cases the program will return an interpolating surface
      !   if s=0 and the constrained least-squares polynomial surface if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the accuracy of the data values.
      !   if the user has an idea of the statistical errors on the data, he
      !   can also find a proper estimate for s. for, by assuming that, if he
      !   specifies the right s, parsur will return a surface s(u,v) which
      !   exactly reproduces the surface underlying the data he can evaluate
      !   the sum(dist(f(i,j)-s(u(i),v(j)))**2) to find a good estimate for s.
      !   for example, if he knows that the statistical errors on his f(i,j)-
      !   values is not greater than 0.1, he may expect that a good s should
      !   have a value not larger than mu*mv*(0.1)**2.
      !   if nothing is known about the statistical error in f(i,j), s must
      !   be determined by trial and error, taking account of the comments
      !   above. the best is then to start with a very large value of s (to
      !   determine the le-sq polynomial surface and the corresponding upper
      !   bound fp0 for s) and then to progressively decrease the value of s
      !   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
      !   and more carefully as the approximation shows more detail) to
      !   obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=0.
      !   if iopt = 1 the program will continue with the knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if parsur is called repeatedly for different values of s.
      !   the number of knots of the surface returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   surface underlying the data. if the computation mode iopt = 1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt=1,the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   parsur once more with the chosen value for s but now with iopt=0.
      !   indeed, parsur may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !   the number of knots may also depend on the upper bounds nuest and
      !   nvest. indeed, if at a certain stage in parsur the number of knots
      !   in one direction (say nu) has reached the value of its upper bound
      !   (nuest), then from that moment on all subsequent knots are added
      !   in the other (v) direction. this may indicate that the value of
      !   nuest is too small. on the other hand, it gives the user the option
      !   of limiting the number of knots the routine locates in any direction
      !   for example, by setting nuest=8 (the lowest allowable value for
      !   nuest), the user can indicate that he wants an approximation with
      !   splines which are simple cubic polynomials in the variable u.
      !
      !  other subroutines required:
      !    fppasu,fpchec,fpchep,fpknot,fprati,fpgrpa,fptrnp,fpback,
      !    fpbacp,fpbspl,fptrpe,fpdisc,fpgivs,fprota
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1989
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) s,fp
      integer iopt,idim,mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
      !  ..array arguments..
      real(RKIND) u(mu),v(mv),f(mu*mv*idim),tu(nuest),tv(nvest), &
       c((nuest-4)*(nvest-4)*idim),wrk(lwrk)
      integer ipar(2),iwrk(kwrk)
      !  ..local scalars..
      real(RKIND) tol,ub,ue,vb,ve,peru,perv
      integer i,j,jwrk,kndu,kndv,knru,knrv,kwest,l1,l2,l3,l4, &
       lfpu,lfpv,lwest,lww,maxit,nc,mf,mumin,mvmin
      !  ..subroutine references..
      !    fppasu,fpchec,fpchep
      !  ..
      !  we set up the parameters tol and maxit.
      maxit = 20
      tol = smallnum03
      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt<(-1) .or. iopt>1) go to 200
      if(ipar(1)<0 .or. ipar(1)>1) go to 200
      if(ipar(2)<0 .or. ipar(2)>1) go to 200
      if(idim<=0 .or. idim>3) go to 200
      mumin = 4-2*ipar(1)
      if(mu<mumin .or. nuest<8) go to 200
      mvmin = 4-2*ipar(2)
      if(mv<mvmin .or. nvest<8) go to 200
      mf = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+ &
       4*(mu+mv)+max0(nuest,mv)*idim
      kwest = 3+mu+mv+nuest+nvest
      if(lwrk<lwest .or. kwrk<kwest) go to 200
      do 10 i=2,mu
        if(u(i-1)>=u(i)) go to 200
  10  continue
      do 20 i=2,mv
        if(v(i-1)>=v(i)) go to 200
  20  continue
      if(iopt>=0) go to 100
      if(nu<8 .or. nu>nuest) go to 200
      ub = u(1)
      ue = u(mu)
      if (ipar(1)/=0) go to 40
      j = nu
      do 30 i=1,4
        tu(i) = ub
        tu(j) = ue
        j = j-1
  30  continue
      call fpchec(u,mu,tu,nu,3,ier)
      if(ier/=0) go to 200
      go to 60
  40  l1 = 4
      l2 = l1
      l3 = nu-3
      l4 = l3
      peru = ue-ub
      tu(l2) = ub
      tu(l3) = ue
      do 50 j=1,3
        l1 = l1+1
        l2 = l2-1
        l3 = l3+1
        l4 = l4-1
        tu(l2) = tu(l4)-peru
        tu(l3) = tu(l1)+peru
  50  continue
      ier = fpchep(u,mu,tu,nu,3)
      if(ier/=0) go to 200
  60  if(nv<8 .or. nv>nvest) go to 200
      vb = v(1)
      ve = v(mv)
      if (ipar(2)/=0) go to 80
      j = nv
      do 70 i=1,4
        tv(i) = vb
        tv(j) = ve
        j = j-1
  70  continue
      call fpchec(v,mv,tv,nv,3,ier)
      if(ier/=0) go to 200
      go to 150
  80  l1 = 4
      l2 = l1
      l3 = nv-3
      l4 = l3
      perv = ve-vb
      tv(l2) = vb
      tv(l3) = ve
      do 90 j=1,3
        l1 = l1+1
        l2 = l2-1
        l3 = l3+1
        l4 = l4-1
        tv(l2) = tv(l4)-perv
        tv(l3) = tv(l1)+perv
  90  continue
      ier = fpchep(v,mv,tv,nv,3)
      if (ier==0) go to 150
      go to 200
 100  if(s<0.) go to 200
      if(s==zero .and. (nuest<(mu+4+2*ipar(1)) .or. &
       nvest<(mv+4+2*ipar(2))) )go to 200
      ier = 0
      !  we partition the working space and determine the spline approximation
 150  lfpu = 5
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-4-nuest-nvest
      knru = 4
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fppasu(iopt,ipar,idim,u,mu,v,mv,f,mf,s,nuest,nvest, &
       tol,maxit,nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4), &
       wrk(lfpu),wrk(lfpv),iwrk(1),iwrk(2),iwrk(3),iwrk(knru), &
       iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww),jwrk,ier)
 200  return
      end subroutine parsur



      recursive subroutine percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp, &
       wrk,lwrk,iwrk,ier)

      !  given the set of data points (x(i),y(i)) and the set of positive
      !  numbers w(i),i=1,2,...,m-1, subroutine percur determines a smooth
      !  periodic spline approximation of degree k with period per=x(m)-x(1).
      !  if iopt=-1 percur calculates the weighted least-squares periodic
      !  spline according to a given set of knots.
      !  if iopt>=0 the number of knots of the spline s(x) and the position
      !  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
      !  ness of s(x) is then achieved by minimalizing the discontinuity
      !  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
      !  n-k-1. the amount of smoothness is determined by the condition that
      !  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
      !  negative constant, called the smoothing factor.
      !  the fit s(x) is given in the b-spline representation (b-spline coef-
      !  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
      !  subroutine splev.
      !
      !  calling sequence:
      !     call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,
      !    * lwrk,iwrk,ier)
      !
      !  parameters:
      !   iopt  : integer flag. on entry iopt must specify whether a weighted
      !           least-squares spline (iopt=-1) or a smoothing spline (iopt=
      !           0 or 1) must be determined. if iopt=0 the routine will start
      !           with an initial set of knots t(i)=x(1)+(x(m)-x(1))*(i-k-1),
      !           i=1,2,...,2*k+2. if iopt=1 the routine will continue with
      !           the knots found at the last call of the routine.
      !           attention: a call with iopt=1 must always be immediately
      !           preceded by another call with iopt=1 or iopt=0.
      !           unchanged on exit.
      !   m     : integer. on entry m must specify the number of data points.
      !           m > 1. unchanged on exit.
      !   x     : real array of dimension at least (m). before entry, x(i)
      !           must be set to the i-th value of the independent variable x,
      !           for i=1,2,...,m. these values must be supplied in strictly
      !           ascending order. x(m) only indicates the length of the
      !           period of the spline, i.e per=x(m)-x(1).
      !           unchanged on exit.
      !   y     : real array of dimension at least (m). before entry, y(i)
      !           must be set to the i-th value of the dependent variable y,
      !           for i=1,2,...,m-1. the element y(m) is not used.
      !           unchanged on exit.
      !   w     : real array of dimension at least (m). before entry, w(i)
      !           must be set to the i-th value in the set of weights. the
      !           w(i) must be strictly positive. w(m) is not used.
      !           see also further comments. unchanged on exit.
      !   k     : integer. on entry k must specify the degree of the spline.
      !           1<=k<=5. it is recommended to use cubic splines (k=3).
      !           the user is strongly dissuaded from choosing k even,together
      !           with a small s-value. unchanged on exit.
      !   s     : real.on entry (in case iopt>=0) s must specify the smoothing
      !           factor. s >=0. unchanged on exit.
      !           for advice on the choice of s see further comments.
      !   nest  : integer. on entry nest must contain an over-estimate of the
      !           total number of knots of the spline returned, to indicate
      !           the storage space available to the routine. nest >=2*k+2.
      !           in most practical situation nest=m/2 will be sufficient.
      !           always large enough is nest=m+2*k,the number of knots needed
      !           for interpolation (s=0). unchanged on exit.
      !   n     : integer.
      !           unless ier = 10 (in case iopt >=0), n will contain the
      !           total number of knots of the spline approximation returned.
      !           if the computation mode iopt=1 is used this value of n
      !           should be left unchanged between subsequent calls.
      !           in case iopt=-1, the value of n must be specified on entry.
      !   t     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the knots of the
      !           spline,i.e. the position of the interior knots t(k+2),t(k+3)
      !           ...,t(n-k-1) as well as the position of the additional knots
      !           t(1),t(2),...,t(k+1)=x(1) and t(n-k)=x(m),..,t(n) needed for
      !           the b-spline representation.
      !           if the computation mode iopt=1 is used, the values of t(1),
      !           t(2),...,t(n) should be left unchanged between subsequent
      !           calls. if the computation mode iopt=-1 is used, the values
      !           t(k+2),...,t(n-k-1) must be supplied by the user, before
      !           entry. see also the restrictions (ier=10).
      !   c     : real array of dimension at least (nest).
      !           on successful exit, this array will contain the coefficients
      !           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
      !   fp    : real. unless ier = 10, fp contains the weighted sum of
      !           squared residuals of the spline approximation returned.
      !   wrk   : real array of dimension at least (m*(k+1)+nest*(8+5*k)).
      !           used as working space. if the computation mode iopt=1 is
      !           used, the values wrk(1),...,wrk(n) should be left unchanged
      !           between subsequent calls.
      !   lwrk  : integer. on entry,lwrk must specify the actual dimension of
      !           the array wrk as declared in the calling (sub)program. lwrk
      !           must not be too small (see wrk). unchanged on exit.
      !   iwrk  : integer array of dimension at least (nest).
      !           used as working space. if the computation mode iopt=1 is
      !           used,the values iwrk(1),...,iwrk(n) should be left unchanged
      !           between subsequent calls.
      !   ier   : integer. unless the routine detects an error, ier contains a
      !           non-positive value on exit, i.e.
      !    ier=0  : normal return. the spline returned has a residual sum of
      !             squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !             ive tolerance set to 0.001 by the program.
      !    ier=-1 : normal return. the spline returned is an interpolating
      !             periodic spline (fp=0).
      !    ier=-2 : normal return. the spline returned is the weighted least-
      !             squares constant. in this extreme case fp gives the upper
      !             bound fp0 for the smoothing factor s.
      !    ier=1  : error. the required storage space exceeds the available
      !             storage space, as specified by the parameter nest.
      !             probably causes : nest too small. if nest is already
      !             large (say nest > m/2), it may also indicate that s is
      !             too small
      !             the approximation returned is the least-squares periodic
      !             spline according to the knots t(1),t(2),...,t(n). (n=nest)
      !             the parameter fp gives the corresponding weighted sum of
      !             squared residuals (fp>s).
      !    ier=2  : error. a theoretically impossible result was found during
      !             the iteration process for finding a smoothing spline with
      !             fp = s. probably causes : s too small.
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=3  : error. the maximal number of iterations maxit (set to 20
      !             by the program) allowed for finding a smoothing spline
      !             with fp=s has been reached. probably causes : s too small
      !             there is an approximation returned but the corresponding
      !             weighted sum of squared residuals does not satisfy the
      !             condition abs(fp-s)/s < tol.
      !    ier=10 : error. on entry, the input data are controlled on validity
      !             the following restrictions must be satisfied.
      !             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,...,m-1
      !             x(1)<x(2)<...<x(m), lwrk>=(k+1)*m+nest*(8+5*k)
      !             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k)
      !                         x(1)<t(k+2)<t(k+3)<...<t(n-k-1)<x(m)
      !                       the schoenberg-whitney conditions, i.e. there
      !                       must be a subset of data points xx(j) with
      !                       xx(j) = x(i) or x(i)+(x(m)-x(1)) such that
      !                         t(j) < xx(j) < t(j+k+1), j=k+1,...,n-k-1
      !             if iopt>=0: s>=0
      !                         if s=0 : nest >= m+2*k
      !             if one of these conditions is found to be violated,control
      !             is immediately repassed to the calling program. in that
      !             case there is no approximation returned.
      !
      !  further comments:
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating periodic
      !   spline if s=0 and the weighted least-squares constant if s is very
      !   large. between these extremes, a properly chosen s will result in
      !   a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the weights w(i). if these are
      !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
      !   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
      !   sqrt(2*m)). if nothing is known about the statistical error in y(i)
      !   each w(i) can be set equal to one and s determined by trial and
      !   error, taking account of the comments above. the best is then to
      !   start with a very large value of s ( to determine the least-squares
      !   constant and the corresponding upper bound fp0 for s) and then to
      !   progressively decrease the value of s ( say by a factor 10 in the
      !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
      !   approximation shows more detail) to obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=0.
      !   if iopt=1 the program will continue with the set of knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if percur is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. but, if the computation mode iopt=1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt=1, the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   percur once more with the selected value for s but now with iopt=0.
      !   indeed, percur may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !
      !  other subroutines required:
      !    fpbacp,fpbspl,fpchep,fpperi,fpdisc,fpgivs,fpknot,fprati,fprota
      !
      !  references:
      !   dierckx p. : algorithms for smoothing data with periodic and
      !                parametric splines, computer graphics and image
      !                processing 20 (1982) 171-184.
      !   dierckx p. : algorithms for smoothing data with periodic and param-
      !                etric splines, report tw55, dept. computer science,
      !                k.u.leuven, 1981.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !  latest update : march 1987
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) s,fp
      integer iopt,m,k,nest,n,lwrk,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
      !  ..local scalars..
      real(RKIND) per,tol
      integer i,ia1,ia2,ib,ifp,ig1,ig2,iq,iz,i1,i2,j1,j2,k1,k2,lwest, &
       maxit,m1,nmin
      !  ..subroutine references..
      !    perper,pcheck
      !  ..
      !  we set up the parameters tol and maxit
      maxit = 20
      tol = smallnum03
      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k<=0 .or. k>5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt<(-1) .or. iopt>1) go to 50
      nmin = 2*k1
      if(m<2 .or. nest<nmin) go to 50
      lwest = m*k1+nest*(8+5*k)
      if(lwrk<lwest) go to 50
      m1 = m-1
      do 10 i=1,m1
         if(x(i)>=x(i+1) .or. w(i)<=0.) go to 50
  10  continue
      if(iopt>=0) go to 30
      if(n<=nmin .or. n>nest) go to 50
      per = x(m)-x(1)
      j1 = k1
      t(j1) = x(1)
      i1 = n-k
      t(i1) = x(m)
      j2 = j1
      i2 = i1
      do 20 i=1,k
         i1 = i1+1
         i2 = i2-1
         j1 = j1+1
         j2 = j2-1
         t(j2) = t(i2)-per
         t(i1) = t(j1)+per
  20  continue
      ier = fpchep(x,m,t,n,k)
      if (ier==0) go to 40
      go to 50
  30  if(s<0.) go to 50
      if(s==zero .and. nest<(m+2*k)) go to 50
      ier = 0
      ! we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia1 = iz+nest
      ia2 = ia1+nest*k1
      ib = ia2+nest*k
      ig1 = ib+nest*k2
      ig2 = ig1+nest*k2
      iq = ig2+nest*k1
      call fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit,k1,k2,n,t,c,fp, &
       wrk(ifp),wrk(iz),wrk(ia1),wrk(ia2),wrk(ib),wrk(ig1),wrk(ig2), &
       wrk(iq),iwrk,ier)
  50  return
      end subroutine percur


      recursive subroutine pogrid(iopt,ider,mu,u,mv,v,z,z0,r,s, &
       nuest,nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)

      !  subroutine pogrid fits a function f(x,y) to a set of data points
      !  z(i,j) given at the nodes (x,y)=(u(i)*cos(v(j)),u(i)*sin(v(j))),
      !  i=1,...,mu ; j=1,...,mv , of a radius-angle grid over a disc
      !    x ** 2  +  y ** 2  <=  r ** 2 .
      !
      !  this approximation problem is reduced to the determination of a
      !  bicubic spline s(u,v) smoothing the data (u(i),v(j),z(i,j)) on the
      !  rectangle 0<=u<=r, v(1)<=v<=v(1)+2*pi
      !  in order to have continuous partial derivatives
      !              i+j
      !             d   f(0,0)
      !    g(i,j) = ----------
      !                i   j
      !              dx  dy
      !
      !  s(u,v)=f(x,y) must satisfy the following conditions
      !
      !    (1) s(0,v) = g(0,0)   v(1)<=v<= v(1)+2*pi
      !
      !        d s(0,v)
      !    (2) -------- = cos(v)*g(1,0)+sin(v)*g(0,1)  v(1)<=v<= v(1)+2*pi
      !        d u
      !
      !  moreover, s(u,v) must be periodic in the variable v, i.e.
      !
      !         j            j
      !        d s(u,vb)   d s(u,ve)
      !    (3) ---------- = ---------   0 <=u<= r, j=0,1,2 , vb=v(1),
      !           j            j                             ve=vb+2*pi
      !        d v          d v
      !
      !  the number of knots of s(u,v) and their position tu(i),i=1,2,...,nu;
      !  tv(j),j=1,2,...,nv, is chosen automatically by the routine. the
      !  smoothness of s(u,v) is achieved by minimalizing the discontinuity
      !  jumps of the derivatives of the spline at the knots. the amount of
      !  smoothness of s(u,v) is determined by the condition that
      !  fp=sumi=1,mu(sumj=1,mv((z(i,j)-s(u(i),v(j)))**2))+(z0-g(0,0))**2<=s,
      !  with s a given non-negative constant.
      !  the fit s(u,v) is given in its b-spline representation and can be
      !  evaluated by means of routine bispev. f(x,y) = s(u,v) can also be
      !  evaluated by means of function program evapol.
      !
      ! calling sequence:
      !     call pogrid(iopt,ider,mu,u,mv,v,z,z0,r,s,nuest,nvest,nu,tu,
      !    *  ,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      !
      ! parameters:
      !  iopt  : integer array of dimension 3, specifying different options.
      !          unchanged on exit.
      !  iopt(1):on entry iopt(1) must specify whether a least-squares spline
      !          (iopt(1)=-1) or a smoothing spline (iopt(1)=0 or 1) must be
      !          determined.
      !          if iopt(1)=0 the routine will start with an initial set of
      !          knots tu(i)=0,tu(i+4)=r,i=1,...,4;tv(i)=v(1)+(i-4)*2*pi,i=1,.
      !          ...,8.
      !          if iopt(1)=1 the routine will continue with the set of knots
      !          found at the last call of the routine.
      !          attention: a call with iopt(1)=1 must always be immediately
      !          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
      !  iopt(2):on entry iopt(2) must specify the requested order of conti-
      !          nuity for f(x,y) at the origin.
      !          if iopt(2)=0 only condition (1) must be fulfilled and
      !          if iopt(2)=1 conditions (1)+(2) must be fulfilled.
      !  iopt(3):on entry iopt(3) must specify whether (iopt(3)=1) or not
      !          (iopt(3)=0) the approximation f(x,y) must vanish at the
      !          boundary of the approximation domain.
      !  ider  : integer array of dimension 2, specifying different options.
      !          unchanged on exit.
      !  ider(1):on entry ider(1) must specify whether (ider(1)=0 or 1) or not
      !          (ider(1)=-1) there is a data value z0 at the origin.
      !          if ider(1)=1, z0 will be considered to be the right function
      !          value, and it will be fitted exactly (g(0,0)=z0=c(1)).
      !          if ider(1)=0, z0 will be considered to be a data value just
      !          like the other data values z(i,j).
      !  ider(2):on entry ider(2) must specify whether (ider(2)=1) or not
      !          (ider(2)=0) f(x,y) must have vanishing partial derivatives
      !          g(1,0) and g(0,1) at the origin. (in case iopt(2)=1)
      !  mu    : integer. on entry mu must specify the number of grid points
      !          along the u-axis. unchanged on exit.
      !          mu >= mumin where mumin=4-iopt(3)-ider(2) if ider(1)<0
      !                                 =3-iopt(3)-ider(2) if ider(1)>=0
      !  u     : real array of dimension at least (mu). before entry, u(i)
      !          must be set to the u-co-ordinate of the i-th grid point
      !          along the u-axis, for i=1,2,...,mu. these values must be
      !          positive and supplied in strictly ascending order.
      !          unchanged on exit.
      !  mv    : integer. on entry mv must specify the number of grid points
      !          along the v-axis. mv > 3 . unchanged on exit.
      !  v     : real array of dimension at least (mv). before entry, v(j)
      !          must be set to the v-co-ordinate of the j-th grid point
      !          along the v-axis, for j=1,2,...,mv. these values must be
      !          supplied in strictly ascending order. unchanged on exit.
      !          -pi <= v(1) < pi , v(mv) < v(1)+2*pi.
      !  z     : real array of dimension at least (mu*mv).
      !          before entry, z(mv*(i-1)+j) must be set to the data value at
      !          the grid point (u(i),v(j)) for i=1,...,mu and j=1,...,mv.
      !          unchanged on exit.
      !  z0    : real value. on entry (if ider(1) >=0 ) z0 must specify the
      !          data value at the origin. unchanged on exit.
      !  r     : real value. on entry r must specify the radius of the disk.
      !          r>=u(mu) (>u(mu) if iopt(3)=1). unchanged on exit.
      !  s     : real. on entry (if iopt(1)>=0) s must specify the smoothing
      !          factor. s >=0. unchanged on exit.
      !          for advice on the choice of s see further comments
      !  nuest : integer. unchanged on exit.
      !  nvest : integer. unchanged on exit.
      !          on entry, nuest and nvest must specify an upper bound for the
      !          number of knots required in the u- and v-directions respect.
      !          these numbers will also determine the storage space needed by
      !          the routine. nuest >= 8, nvest >= 8.
      !          in most practical situation nuest = mu/2, nvest=mv/2, will
      !          be sufficient. always large enough are nuest=mu+5+iopt(2)+
      !          iopt(3), nvest = mv+7, the number of knots needed for
      !          interpolation (s=0). see also further comments.
      !  nu    : integer.
      !          unless ier=10 (in case iopt(1)>=0), nu will contain the total
      !          number of knots with respect to the u-variable, of the spline
      !          approximation returned. if the computation mode iopt(1)=1 is
      !          used, the value of nu should be left unchanged between sub-
      !          sequent calls. in case iopt(1)=-1, the value of nu should be
      !          specified on entry.
      !  tu    : real array of dimension at least (nuest).
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the u-variable, i.e. the position of
      !          the interior knots tu(5),...,tu(nu-4) as well as the position
      !          of the additional knots tu(1)=...=tu(4)=0 and tu(nu-3)=...=
      !          tu(nu)=r needed for the b-spline representation.
      !          if the computation mode iopt(1)=1 is used,the values of tu(1)
      !          ...,tu(nu) should be left unchanged between subsequent calls.
      !          if the computation mode iopt(1)=-1 is used, the values tu(5),
      !          ...tu(nu-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  nv    : integer.
      !          unless ier=10 (in case iopt(1)>=0), nv will contain the total
      !          number of knots with respect to the v-variable, of the spline
      !          approximation returned. if the computation mode iopt(1)=1 is
      !          used, the value of nv should be left unchanged between sub-
      !          sequent calls. in case iopt(1) = -1, the value of nv should
      !          be specified on entry.
      !  tv    : real array of dimension at least (nvest).
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the v-variable, i.e. the position of
      !          the interior knots tv(5),...,tv(nv-4) as well as the position
      !          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
      !          tv(nv) needed for the b-spline representation.
      !          if the computation mode iopt(1)=1 is used,the values of tv(1)
      !          ...,tv(nv) should be left unchanged between subsequent calls.
      !          if the computation mode iopt(1)=-1 is used, the values tv(5),
      !          ...tv(nv-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  c     : real array of dimension at least (nuest-4)*(nvest-4).
      !          on successful exit, c contains the coefficients of the spline
      !          approximation s(u,v)
      !  fp    : real. unless ier=10, fp contains the sum of squared
      !          residuals of the spline approximation returned.
      !  wrk   : real array of dimension (lwrk). used as workspace.
      !          if the computation mode iopt(1)=1 is used the values of
      !          wrk(1),...,wrk(8) should be left unchanged between subsequent
      !          calls.
      !  lwrk  : integer. on entry lwrk must specify the actual dimension of
      !          the array wrk as declared in the calling (sub)program.
      !          lwrk must not be too small.
      !           lwrk >= 8+nuest*(mv+nvest+3)+nvest*21+4*mu+6*mv+q
      !           where q is the larger of (mv+nvest) and nuest.
      !  iwrk  : integer array of dimension (kwrk). used as workspace.
      !          if the computation mode iopt(1)=1 is used the values of
      !          iwrk(1),.,iwrk(4) should be left unchanged between subsequent
      !          calls.
      !  kwrk  : integer. on entry kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program.
      !          kwrk >= 4+mu+mv+nuest+nvest.
      !  ier   : integer. unless the routine detects an error, ier contains a
      !          non-positive value on exit, i.e.
      !   ier=0  : normal return. the spline returned has a residual sum of
      !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !            ive tolerance set to 0.001 by the program.
      !   ier=-1 : normal return. the spline returned is an interpolating
      !            spline (fp=0).
      !   ier=-2 : normal return. the spline returned is the least-squares
      !            constrained polynomial. in this extreme case fp gives the
      !            upper bound for the smoothing factor s.
      !   ier=1  : error. the required storage space exceeds the available
      !            storage space, as specified by the parameters nuest and
      !            nvest.
      !            probably causes : nuest or nvest too small. if these param-
      !            eters are already large, it may also indicate that s is
      !            too small
      !            the approximation returned is the least-squares spline
      !            according to the current set of knots. the parameter fp
      !            gives the corresponding sum of squared residuals (fp>s).
      !   ier=2  : error. a theoretically impossible result was found during
      !            the iteration process for finding a smoothing spline with
      !            fp = s. probably causes : s too small.
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=3  : error. the maximal number of iterations maxit (set to 20
      !            by the program) allowed for finding a smoothing spline
      !            with fp=s has been reached. probably causes : s too small
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=10 : error. on entry, the input data are controlled on validity
      !            the following restrictions must be satisfied.
      !            -1<=iopt(1)<=1, 0<=iopt(2)<=1, 0<=iopt(3)<=1,
      !            -1<=ider(1)<=1, 0<=ider(2)<=1, ider(2)=0 if iopt(2)=0.
      !            mu >= mumin (see above), mv >= 4, nuest >=8, nvest >= 8,
      !            kwrk>=4+mu+mv+nuest+nvest,
      !            lwrk >= 8+nuest*(mv+nvest+3)+nvest*21+4*mu+6*mv+
      !             max(nuest,mv+nvest)
      !            0< u(i-1)<u(i)<=r,i=2,..,mu, (< r if iopt(3)=1)
      !            -pi<=v(1)< pi, v(1)<v(i-1)<v(i)<v(1)+2*pi, i=3,...,mv
      !            if iopt(1)=-1: 8<=nu<=min(nuest,mu+5+iopt(2)+iopt(3))
      !                           0<tu(5)<tu(6)<...<tu(nu-4)<r
      !                           8<=nv<=min(nvest,mv+7)
      !                           v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(1)+2*pi
      !                    the schoenberg-whitney conditions, i.e. there must
      !                    be subset of grid co-ordinates uu(p) and vv(q) such
      !                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
      !                     (iopt(2)=1 and iopt(3)=1 also count for a uu-value
      !                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
      !                     (vv(q) is either a value v(j) or v(j)+2*pi)
      !            if iopt(1)>=0: s>=0
      !                       if s=0: nuest>=mu+5+iopt(2)+iopt(3), nvest>=mv+7
      !            if one of these conditions is found to be violated,control
      !            is immediately repassed to the calling program. in that
      !            case there is no approximation returned.
      !
      ! further comments:
      !   pogrid does not allow individual weighting of the data-values.
      !   so, if these were determined to widely different accuracies, then
      !   perhaps the general data set routine polar should rather be used
      !   in spite of efficiency.
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating spline if
      !   s=0 and the constrained least-squares polynomial(degrees 3,0)if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the accuracy of the data values.
      !   if the user has an idea of the statistical errors on the data, he
      !   can also find a proper estimate for s. for, by assuming that, if he
      !   specifies the right s, pogrid will return a spline s(u,v) which
      !   exactly reproduces the function underlying the data he can evaluate
      !   the sum((z(i,j)-s(u(i),v(j)))**2) to find a good estimate for this s
      !   for example, if he knows that the statistical errors on his z(i,j)-
      !   values is not greater than 0.1, he may expect that a good s should
      !   have a value not larger than mu*mv*(0.1)**2.
      !   if nothing is known about the statistical error in z(i,j), s must
      !   be determined by trial and error, taking account of the comments
      !   above. the best is then to start with a very large value of s (to
      !   determine the least-squares polynomial and the corresponding upper
      !   bound fp0 for s) and then to progressively decrease the value of s
      !   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
      !   and more carefully as the approximation shows more detail) to
      !   obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt(1)=0.
      !   if iopt(1) = 1 the program will continue with the knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if pogrid is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. if the computation mode iopt(1) = 1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt(1)=1,the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   pogrid once more with the chosen value for s but now with iopt(1)=0.
      !   indeed, pogrid may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !   the number of knots may also depend on the upper bounds nuest and
      !   nvest. indeed, if at a certain stage in pogrid the number of knots
      !   in one direction (say nu) has reached the value of its upper bound
      !   (nuest), then from that moment on all subsequent knots are added
      !   in the other (v) direction. this may indicate that the value of
      !   nuest is too small. on the other hand, it gives the user the option
      !   of limiting the number of knots the routine locates in any direction
      !   for example, by setting nuest=8 (the lowest allowable value for
      !   nuest), the user can indicate that he wants an approximation which
      !   is a simple cubic polynomial in the variable u.
      !
      !  other subroutines required:
      !    fppogr,fpchec,fpchep,fpknot,fpopdi,fprati,fpgrdi,fpsysy,fpback,
      !    fpbacp,fpbspl,fpcyt1,fpcyt2,fpdisc,fpgivs,fprota
      !
      !  references:
      !   dierckx p. : fast algorithms for smoothing data over a disc or a
      !                sphere using tensor product splines, in "algorithms
      !                for approximation", ed. j.c.mason and m.g.cox,
      !                clarendon press oxford, 1987, pp. 51-65
      !   dierckx p. : fast algorithms for smoothing data over a disc or a
      !                sphere using tensor product splines, report tw73, dept.
      !                computer science,k.u.leuven, 1985.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : july 1985
      !  latest update : march 1989
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) z0,r,s,fp
      integer mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
      !  ..array arguments..
      integer iopt(3),ider(2),iwrk(kwrk)
      real(RKIND) u(mu),v(mv),z(mu*mv),c((nuest-4)*(nvest-4)),tu(nuest), &
       tv(nvest),wrk(lwrk)
      !  ..local scalars..
      real(RKIND) per,tol,uu,ve,zmax,zmin,rn,zb
      integer i,i1,i2,j,jwrk,j1,j2,kndu,kndv,knru,knrv,kwest,l, &
       ldz,lfpu,lfpv,lwest,lww,m,maxit,mumin,muu,nc

      !  set constants
      per = pi+pi
      ve = v(1)+per
      !  we set up the parameters tol and maxit.
      maxit = 20
      tol = smallnum03
      !  before starting computations, a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt(1)<(-1) .or. iopt(1)>1) go to 200
      if(iopt(2)<0 .or. iopt(2)>1) go to 200
      if(iopt(3)<0 .or. iopt(3)>1) go to 200
      if(ider(1)<(-1) .or. ider(1)>1) go to 200
      if(ider(2)<0 .or. ider(2)>1) go to 200
      if(ider(2)==1 .and. iopt(2)==0) go to 200
      mumin = 4-iopt(3)-ider(2)
      if(ider(1)>=0) mumin = mumin-1
      if(mu<mumin .or. mv<4) go to 200
      if(nuest<8 .or. nvest<8) go to 200
      m = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 8+nuest*(mv+nvest+3)+21*nvest+4*mu+6*mv+ &
       max0(nuest,mv+nvest)
      kwest = 4+mu+mv+nuest+nvest
      if(lwrk<lwest .or. kwrk<kwest) go to 200
      if(u(1)<=0. .or. u(mu)>r) go to 200
      if(iopt(3)==0) go to 10
      if(u(mu)==r) go to 200
  10  if(mu==1) go to 30
      do 20 i=2,mu
        if(u(i-1)>=u(i)) go to 200
  20  continue
  30  if(v(1)< (-pi) .or. v(1)>=pi ) go to 200
      if(v(mv)>=v(1)+per) go to 200
      do 40 i=2,mv
        if(v(i-1)>=v(i)) go to 200
  40  continue
      if(iopt(1)>0) go to 140
      !  if not given, we compute an estimate for z0.
      if(ider(1)<0) go to 50
      zb = z0
      go to 70
  50  zb = sum(z(1:mv))
      rn = mv
      zb = zb/rn
      !  we determine the range of z-values.
  70  zmin = zb
      zmax = zb
      do 80 i=1,m
         if(z(i)<zmin) zmin = z(i)
         if(z(i)>zmax) zmax = z(i)
  80  continue
      wrk(5) = zb
      wrk(6) = zero
      wrk(7) = zero
      wrk(8) = zmax -zmin
      iwrk(4) = mu
      if(iopt(1)==0) go to 140
      if(nu<8 .or. nu>nuest) go to 200
      if(nv<11 .or. nv>nvest) go to 200
      j = nu
      do 90 i=1,4
        tu(i) = zero
        tu(j) = r
        j = j-1
  90  continue
      l = 9
      wrk(l) = zero
      if(iopt(2)==0) go to 100
      l = l+1
      uu = u(1)
      if(uu>tu(5)) uu = tu(5)
      wrk(l) = uu*half
 100  do 110 i=1,mu
        l = l+1
        wrk(l) = u(i)
 110  continue
      if(iopt(3)==0) go to 120
      l = l+1
      wrk(l) = r
 120  muu = l-8
      call fpchec(wrk(9),muu,tu,nu,3,ier)
      if(ier/=0) go to 200
      j1 = 4
      tv(j1) = v(1)
      i1 = nv-3
      tv(i1) = ve
      j2 = j1
      i2 = i1
      do 130 i=1,3
        i1 = i1+1
        i2 = i2-1
        j1 = j1+1
        j2 = j2-1
        tv(j2) = tv(i2)-per
        tv(i1) = tv(j1)+per
 130  continue
      l = 9
      do 135 i=1,mv
        wrk(l) = v(i)
        l = l+1
 135  continue
      wrk(l) = ve
      ier = fpchep(wrk(9),mv+1,tv,nv,3)
      if (ier==0) go to 150
      go to 200
 140  if(s<zero) go to 200
      if(s==zero .and. (nuest<(mu+5+iopt(2)+iopt(3)) .or. nvest<(mv+7)) ) go to 200
      !  we partition the working space and determine the spline approximation
 150  ldz = 5
      lfpu = 9
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-8-nuest-nvest
      knru = 5
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fppogr(iopt,ider,u,mu,v,mv,z,m,zb,r,s,nuest,nvest,tol,maxit, &
       nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),wrk(lfpu), &
       wrk(lfpv),wrk(ldz),wrk(8),iwrk(1),iwrk(2),iwrk(3),iwrk(4), &
       iwrk(knru),iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww),jwrk,ier)
 200  return
      end subroutine pogrid



      recursive subroutine polar(iopt,m,x,y,z,w,rad,s,nuest,nvest, &
        eps,nu,tu,nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

      !  subroutine polar fits a smooth function f(x,y) to a set of data
      !  points (x(i),y(i),z(i)) scattered arbitrarily over an approximation
      !  domain  x**2+y**2 <= rad(atan(y/x))**2. through the transformation
      !    x = u*rad(v)*cos(v) , y = u*rad(v)*sin(v)
      !  the approximation problem is reduced to the determination of a bi-
      !  cubic spline s(u,v) fitting a corresponding set of data points
      !  (u(i),v(i),z(i)) on the rectangle 0<=u<=1,-pi<=v<=pi.
      !  in order to have continuous partial derivatives
      !              i+j
      !             d   f(0,0)
      !    g(i,j) = ----------
      !                i   j
      !              dx  dy
      !
      !  s(u,v)=f(x,y) must satisfy the following conditions
      !
      !    (1) s(0,v) = g(0,0)   -pi <=v<= pi.
      !
      !        d s(0,v)
      !    (2) -------- = rad(v)*(cos(v)*g(1,0)+sin(v)*g(0,1))
      !        d u
      !                                                    -pi <=v<= pi
      !         2
      !        d s(0,v)         2       2             2
      !    (3) -------- = rad(v)*(cos(v)*g(2,0)+sin(v)*g(0,2)+sin(2*v)*g(1,1))
      !           2
      !        d u                                         -pi <=v<= pi
      !
      !  moreover, s(u,v) must be periodic in the variable v, i.e.
      !
      !         j            j
      !        d s(u,-pi)   d s(u,pi)
      !    (4) ---------- = ---------   0 <=u<= 1, j=0,1,2
      !           j           j
      !        d v         d v
      !
      !  if iopt(1) < 0 circle calculates a weighted least-squares spline
      !  according to a given set of knots in u- and v- direction.
      !  if iopt(1) >=0, the number of knots in each direction and their pos-
      !  ition tu(j),j=1,2,...,nu ; tv(j),j=1,2,...,nv are chosen automatical-
      !  ly by the routine. the smoothness of s(u,v) is then achieved by mini-
      !  malizing the discontinuity jumps of the derivatives of the spline
      !  at the knots. the amount of smoothness of s(u,v) is determined  by
      !  the condition that fp = sum((w(i)*(z(i)-s(u(i),v(i))))**2) be <= s,
      !  with s a given non-negative constant.
      !  the bicubic spline is given in its standard b-spline representation
      !  and the corresponding function f(x,y) can be evaluated by means of
      !  function program evapol.
      !
      ! calling sequence:
      !     call polar(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,
      !    *  nv,tv,u,v,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      !
      ! parameters:
      !  iopt  : integer array of dimension 3, specifying different options.
      !          unchanged on exit.
      !  iopt(1):on entry iopt(1) must specify whether a weighted
      !          least-squares polar spline (iopt(1)=-1) or a smoothing
      !          polar spline (iopt(1)=0 or 1) must be determined.
      !          if iopt(1)=0 the routine will start with an initial set of
      !          knots tu(i)=0,tu(i+4)=1,i=1,...,4;tv(i)=(2*i-9)*pi,i=1,...,8.
      !          if iopt(1)=1 the routine will continue with the set of knots
      !          found at the last call of the routine.
      !          attention: a call with iopt(1)=1 must always be immediately
      !          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
      !  iopt(2):on entry iopt(2) must specify the requested order of conti-
      !          nuity for f(x,y) at the origin.
      !          if iopt(2)=0 only condition (1) must be fulfilled,
      !          if iopt(2)=1 conditions (1)+(2) must be fulfilled and
      !          if iopt(2)=2 conditions (1)+(2)+(3) must be fulfilled.
      !  iopt(3):on entry iopt(3) must specify whether (iopt(3)=1) or not
      !          (iopt(3)=0) the approximation f(x,y) must vanish at the
      !          boundary of the approximation domain.
      !  m     : integer. on entry m must specify the number of data points.
      !          m >= 4-iopt(2)-iopt(3) unchanged on exit.
      !  x     : real array of dimension at least (m).
      !  y     : real array of dimension at least (m).
      !  z     : real array of dimension at least (m).
      !          before entry, x(i),y(i),z(i) must be set to the co-ordinates
      !          of the i-th data point, for i=1,...,m. the order of the data
      !          points is immaterial. unchanged on exit.
      !  w     : real array of dimension at least (m). before entry, w(i) must
      !          be set to the i-th value in the set of weights. the w(i) must
      !          be strictly positive. unchanged on exit.
      !  rad   : real function subprogram defining the boundary of the approx-
      !          imation domain, i.e   x = rad(v)*cos(v) , y = rad(v)*sin(v),
      !          -pi <= v <= pi.
      !          must be declared external in the calling (sub)program.
      !  s     : real. on entry (in case iopt(1) >=0) s must specify the
      !          smoothing factor. s >=0. unchanged on exit.
      !          for advice on the choice of s see further comments
      !  nuest : integer. unchanged on exit.
      !  nvest : integer. unchanged on exit.
      !          on entry, nuest and nvest must specify an upper bound for the
      !          number of knots required in the u- and v-directions resp.
      !          these numbers will also determine the storage space needed by
      !          the routine. nuest >= 8, nvest >= 8.
      !          in most practical situation nuest = nvest = 8+sqrt(m/2) will
      !          be sufficient. see also further comments.
      !  eps   : real.
      !          on entry, eps must specify a threshold for determining the
      !          effective rank of an over-determined linear system of equat-
      !          ions. 0 < eps < 1.  if the number of decimal digits in the
      !          computer representation of a real number is q, then 10**(-q)
      !          is a suitable value for eps in most practical applications.
      !          unchanged on exit.
      !  nu    : integer.
      !          unless ier=10 (in case iopt(1) >=0),nu will contain the total
      !          number of knots with respect to the u-variable, of the spline
      !          approximation returned. if the computation mode iopt(1)=1
      !          is used, the value of nu should be left unchanged between
      !          subsequent calls.
      !          in case iopt(1)=-1,the value of nu must be specified on entry
      !  tu    : real array of dimension at least nuest.
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the u-variable, i.e. the position
      !          of the interior knots tu(5),...,tu(nu-4) as well as the
      !          position of the additional knots tu(1)=...=tu(4)=0 and
      !          tu(nu-3)=...=tu(nu)=1 needed for the b-spline representation
      !          if the computation mode iopt(1)=1 is used,the values of
      !          tu(1),...,tu(nu) should be left unchanged between subsequent
      !          calls. if the computation mode iopt(1)=-1 is used,the values
      !          tu(5),...tu(nu-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  nv    : integer.
      !          unless ier=10 (in case iopt(1)>=0), nv will contain the total
      !          number of knots with respect to the v-variable, of the spline
      !          approximation returned. if the computation mode iopt(1)=1
      !          is used, the value of nv should be left unchanged between
      !          subsequent calls. in case iopt(1)=-1, the value of nv should
      !          be specified on entry.
      !  tv    : real array of dimension at least nvest.
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the v-variable, i.e. the position of
      !          the interior knots tv(5),...,tv(nv-4) as well as the position
      !          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
      !          tv(nv) needed for the b-spline representation.
      !          if the computation mode iopt(1)=1 is used, the values of
      !          tv(1),...,tv(nv) should be left unchanged between subsequent
      !          calls. if the computation mode iopt(1)=-1 is used,the values
      !          tv(5),...tv(nv-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  u     : real array of dimension at least (m).
      !  v     : real array of dimension at least (m).
      !          on successful exit, u(i),v(i) contains the co-ordinates of
      !          the i-th data point with respect to the transformed rectan-
      !          gular approximation domain, for i=1,2,...,m.
      !          if the computation mode iopt(1)=1 is used the values of
      !          u(i),v(i) should be left unchanged between subsequent calls.
      !  c     : real array of dimension at least (nuest-4)*(nvest-4).
      !          on successful exit, c contains the coefficients of the spline
      !          approximation s(u,v).
      !  fp    : real. unless ier=10, fp contains the weighted sum of
      !          squared residuals of the spline approximation returned.
      !  wrk1  : real array of dimension (lwrk1). used as workspace.
      !          if the computation mode iopt(1)=1 is used the value of
      !          wrk1(1) should be left unchanged between subsequent calls.
      !          on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the
      !          values d(i)/max(d(i)),i=1,...,ncof=1+iopt(2)*(iopt(2)+3)/2+
      !          (nv-7)*(nu-5-iopt(2)-iopt(3)) with d(i) the i-th diagonal el-
      !          ement of the triangular matrix for calculating the b-spline
      !          coefficients.it includes those elements whose square is < eps
      !          which are treated as 0 in the case of rank deficiency(ier=-2)
      !  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
      !          the array wrk1 as declared in the calling (sub)program.
      !          lwrk1 must not be too small. let
      !            k = nuest-7, l = nvest-7, p = 1+iopt(2)*(iopt(2)+3)/2,
      !            q = k+2-iopt(2)-iopt(3) then
      !          lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
      !  wrk2  : real array of dimension (lwrk2). used as workspace, but
      !          only in the case a rank deficient system is encountered.
      !  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
      !          the array wrk2 as declared in the calling (sub)program.
      !          lwrk2 > 0 . a save upper bound  for lwrk2 = (p+l*q+1)*(4*l+p)
      !          +p+l*q where p,l,q are as above. if there are enough data
      !          points, scattered uniformly over the approximation domain
      !          and if the smoothing factor s is not too small, there is a
      !          good chance that this extra workspace is not needed. a lot
      !          of memory might therefore be saved by setting lwrk2=1.
      !          (see also ier > 10)
      !  iwrk  : integer array of dimension (kwrk). used as workspace.
      !  kwrk  : integer. on entry kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program.
      !          kwrk >= m+(nuest-7)*(nvest-7).
      !  ier   : integer. unless the routine detects an error, ier contains a
      !          non-positive value on exit, i.e.
      !   ier=0  : normal return. the spline returned has a residual sum of
      !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !            ive tolerance set to 0.001 by the program.
      !   ier=-1 : normal return. the spline returned is an interpolating
      !            spline (fp=0).
      !   ier=-2 : normal return. the spline returned is the weighted least-
      !            squares constrained polynomial . in this extreme case
      !            fp gives the upper bound for the smoothing factor s.
      !   ier<-2 : warning. the coefficients of the spline returned have been
      !            computed as the minimal norm least-squares solution of a
      !            (numerically) rank deficient system. (-ier) gives the rank.
      !            especially if the rank deficiency which can be computed as
      !            1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)-iopt(3))+ier
      !            is large the results may be inaccurate.
      !            they could also seriously depend on the value of eps.
      !   ier=1  : error. the required storage space exceeds the available
      !            storage space, as specified by the parameters nuest and
      !            nvest.
      !            probably causes : nuest or nvest too small. if these param-
      !            eters are already large, it may also indicate that s is
      !            too small
      !            the approximation returned is the weighted least-squares
      !            polar spline according to the current set of knots.
      !            the parameter fp gives the corresponding weighted sum of
      !            squared residuals (fp>s).
      !   ier=2  : error. a theoretically impossible result was found during
      !            the iteration process for finding a smoothing spline with
      !            fp = s. probably causes : s too small or badly chosen eps.
      !            there is an approximation returned but the corresponding
      !            weighted sum of squared residuals does not satisfy the
      !            condition abs(fp-s)/s < tol.
      !   ier=3  : error. the maximal number of iterations maxit (set to 20
      !            by the program) allowed for finding a smoothing spline
      !            with fp=s has been reached. probably causes : s too small
      !            there is an approximation returned but the corresponding
      !            weighted sum of squared residuals does not satisfy the
      !            condition abs(fp-s)/s < tol.
      !   ier=4  : error. no more knots can be added because the dimension
      !            of the spline 1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)
      !            -iopt(3)) already exceeds the number of data points m.
      !            probably causes : either s or m too small.
      !            the approximation returned is the weighted least-squares
      !            polar spline according to the current set of knots.
      !            the parameter fp gives the corresponding weighted sum of
      !            squared residuals (fp>s).
      !   ier=5  : error. no more knots can be added because the additional
      !            knot would (quasi) coincide with an old one.
      !            probably causes : s too small or too large a weight to an
      !            inaccurate data point.
      !            the approximation returned is the weighted least-squares
      !            polar spline according to the current set of knots.
      !            the parameter fp gives the corresponding weighted sum of
      !            squared residuals (fp>s).
      !   ier=10 : error. on entry, the input data are controlled on validity
      !            the following restrictions must be satisfied.
      !            -1<=iopt(1)<=1 , 0<=iopt(2)<=2 , 0<=iopt(3)<=1 ,
      !            m>=4-iopt(2)-iopt(3) , nuest>=8 ,nvest >=8, 0<eps<1,
      !            0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m
      !            lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
      !            kwrk >= m+(nuest-7)*(nvest-7)
      !            if iopt(1)=-1:9<=nu<=nuest,9+iopt(2)*(iopt(2)+1)<=nv<=nvest
      !                          0<tu(5)<tu(6)<...<tu(nu-4)<1
      !                          -pi<tv(5)<tv(6)<...<tv(nv-4)<pi
      !            if iopt(1)>=0: s>=0
      !            if one of these conditions is found to be violated,control
      !            is immediately repassed to the calling program. in that
      !            case there is no approximation returned.
      !   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
      !            space for computing the minimal least-squares solution of
      !            a rank deficient system of linear equations. ier gives the
      !            requested value for lwrk2. there is no approximation re-
      !            turned but, having saved the information contained in nu,
      !            nv,tu,tv,wrk1,u,v and having adjusted the value of lwrk2
      !            and the dimension of the array wrk2 accordingly, the user
      !            can continue at the point the program was left, by calling
      !            polar with iopt(1)=1.
      !
      ! further comments:
      !  by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating spline if
      !   s=0 and the constrained weighted least-squares polynomial if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the weights w(i). if these are
      !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
      !   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
      !   sqrt(2*m)). if nothing is known about the statistical error in z(i)
      !   each w(i) can be set equal to one and s determined by trial and
      !   error, taking account of the comments above. the best is then to
      !   start with a very large value of s ( to determine the least-squares
      !   polynomial and the corresponding upper bound fp0 for s) and then to
      !   progressively decrease the value of s ( say by a factor 10 in the
      !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
      !   approximation shows more detail) to obtain closer fits.
      !   to choose s very small is strongly discouraged. this considerably
      !   increases computation time and memory requirements. it may also
      !   cause rank-deficiency (ier<-2) and endager numerical stability.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt(1)=0.
      !   if iopt(1)=1 the program will continue with the set of knots found
      !   at the last call of the routine. this will save a lot of computation
      !   time if polar is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. if the computation mode iopt(1)=1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt(1)=1,the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   polar once more with the selected value for s but now with iopt(1)=0
      !   indeed, polar may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !   the number of knots may also depend on the upper bounds nuest and
      !   nvest. indeed, if at a certain stage in polar the number of knots
      !   in one direction (say nu) has reached the value of its upper bound
      !   (nuest), then from that moment on all subsequent knots are added
      !   in the other (v) direction. this may indicate that the value of
      !   nuest is too small. on the other hand, it gives the user the option
      !   of limiting the number of knots the routine locates in any direction
      !
      !  other subroutines required:
      !    fpback,fpbspl,fppola,fpdisc,fpgivs,fprank,fprati,fprota,fporde,
      !    fprppo
      !
      !  references:
      !   dierckx p.: an algorithm for fitting data over a circle using tensor
      !               product splines,j.comp.appl.maths 15 (1986) 161-173.
      !   dierckx p.: an algorithm for fitting data on a circle using tensor
      !               product splines, report tw68, dept. computer science,
      !               k.u.leuven, 1984.
      !   dierckx p.: curve and surface fitting with splines, monographs on
      !               numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : june 1984
      !  latest update : march 1989
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(in) :: s,eps,fp
      integer m,nuest,nvest,nu,nv,lwrk1,lwrk2,kwrk,ier
      !  ..array arguments..
      real(RKIND) x(m),y(m),z(m),w(m),tu(nuest),tv(nvest),u(m),v(m), &
       c((nuest-4)*(nvest-4)),wrk1(lwrk1),wrk2(lwrk2)
      integer iopt(3),iwrk(kwrk)
      !  ..user specified function
      real(RKIND) rad
      !  ..local scalars..
      real(RKIND) pi,dist,r
      integer i,ib1,ib3,ki,kn,kwest,la,lbu,lcc,lcs,lro,j, &
       lbv,lco,lf,lff,lfp,lh,lq,lsu,lsv,lwest,ncest,ncc,nuu, &
       nvv,nreg,nrint,nu4,nv4,iopt1,iopt2,iopt3,ipar,nvmin

      !  set up constants
      !  we set up the parameters tol and maxit.
      integer, parameter :: maxit = 20
      real(RKIND), parameter :: tol = smallnum03

      !  before starting computations a data check is made. if the input data
      !  are invalid,control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if(eps<=zero .or. eps>=one) go to 60
      iopt1 = iopt(1)
      if(iopt1<(-1) .or. iopt1>1) go to 60
      iopt2 = iopt(2)
      if(iopt2<0 .or. iopt2>2) go to 60
      iopt3 = iopt(3)
      if(iopt3<0 .or. iopt3>1) go to 60
      if(m<(4-iopt2-iopt3)) go to 60
      if(nuest<8 .or. nvest<8) go to 60
      nu4 = nuest-4
      nv4 = nvest-4
      ncest = nu4*nv4
      nuu = nuest-7
      nvv = nvest-7
      ipar = 1+iopt2*(iopt2+3)/2
      ncc = ipar+nvv*(nuest-5-iopt2-iopt3)
      nrint = nuu+nvv
      nreg = nuu*nvv
      ib1 = 4*nvv
      ib3 = ib1+ipar
      lwest = ncc*(1+ib1+ib3)+2*nrint+ncest+m*8+ib3+5*nuest+12*nvest
      kwest = m+nreg
      if(lwrk1<lwest .or. kwrk<kwest) go to 60
      if(iopt1>0) go to 40
      do 10 i=1,m
        if(w(i)<=0.) go to 60
        dist = x(i)**2+y(i)**2
        u(i) = zero
        v(i) = zero
        if(dist<=zero) go to 10
        v(i) = datan2(y(i),x(i))
        r = rad(v(i))
        if(r<=0.) go to 60
        u(i) = sqrt(dist)/r
        if(u(i)>one) go to 60
  10  continue
      if(iopt1==0) go to 40
      nuu = nu-8
      if(nuu<1 .or. nu>nuest) go to 60
      tu(4) = zero
      do 20 i=1,nuu
         j = i+4
         if(tu(j)<=tu(j-1) .or. tu(j)>=one) go to 60
  20  continue
      nvv = nv-8
      nvmin = 9+iopt2*(iopt2+1)
      if(nv<nvmin .or. nv>nvest) go to 60
      pi = datan2(0d0,-one)
      tv(4) = -pi
      do 30 i=1,nvv
         j = i+4
         if(tv(j)<=tv(j-1) .or. tv(j)>=pi) go to 60
  30  continue
      go to 50
  40  if(s<zero) go to 60
  50  ier = FITPACK_OK
      !  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncc*ib3
      lf = la+ncc*ib1
      lff = lf+ncc
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbu = lh+ib3
      lbv = lbu+5*nuest
      lro = lbv+5*nvest
      lcc = lro+nvest
      lcs = lcc+nvest
      lsu = lcs+nvest*5
      lsv = lsu+m*4
      call fppola(iopt1,iopt2,iopt3,m,u,v,z,w,rad,s,nuest,nvest,eps,tol, &
                  maxit,ib1,ib3,ncest,ncc,nrint,nreg,nu,tu,nv,tv,c,fp,wrk1(1), &
                  wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff),wrk1(lro),wrk1(lcc), &
                  wrk1(lcs),wrk1(la),wrk1(lq),wrk1(lbu),wrk1(lbv),wrk1(lsu), &
                  wrk1(lsv),wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  60  return
      end subroutine polar



      recursive subroutine profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)

      !  if iopt=0 subroutine profil calculates the b-spline coefficients of
      !  the univariate spline f(y) = s(u,y) with s(x,y) a bivariate spline of
      !  degrees kx and ky, given in the b-spline representation.
      !  if iopt = 1 it calculates the b-spline coefficients of the univariate
      !  spline g(x) = s(x,u)
      !
      !  calling sequence:
      !     call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
      !
      !  input parameters:
      !   iopt  : integer flag, specifying whether the profile f(y) (iopt=0)
      !           or the profile g(x) (iopt=1) must be determined.
      !   tx    : real array, length nx, which contains the position of the
      !           knots in the x-direction.
      !   nx    : integer, giving the total number of knots in the x-direction
      !   ty    : real array, length ny, which contains the position of the
      !           knots in the y-direction.
      !   ny    : integer, giving the total number of knots in the y-direction
      !   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
      !           b-spline coefficients.
      !   kx,ky : integer values, giving the degrees of the spline.
      !   u     : real value, specifying the requested profile.
      !           tx(kx+1)<=u<=tx(nx-kx), if iopt=0.
      !           ty(ky+1)<=u<=ty(ny-ky), if iopt=1.
      !   nu    : on entry nu must specify the dimension of the array cu.
      !           nu >= ny if iopt=0, nu >= nx if iopt=1.
      !
      !  output parameters:
      !   cu    : real array of dimension (nu).
      !           on successful exit this array contains the b-spline
      !   ier   : integer error flag
      !    ier=0 : normal return
      !    ier=10: invalid input data (see restrictions)
      !
      !  restrictions:
      !   if iopt=0 : tx(kx+1) <= u <= tx(nx-kx), nu >=ny.
      !   if iopt=1 : ty(ky+1) <= u <= ty(ny-ky), nu >=nx.
      !
      !  other subroutines required:
      !    fpbspl
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1987
      !
      !  ..scalar arguments..
      integer iopt,nx,ny,kx,ky,nu,ier
      real(RKIND) u
      !  ..array arguments..
      real(RKIND) tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),cu(nu)
      !  ..local scalars..
      integer i,j,kx1,ky1,l,l1,m,m0,nkx1,nky1
      real(RKIND) sum
      !  ..local array
      real(RKIND) h(SIZ_K+1)
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      ier = 10

      select case (iopt)

         case (0)

                  if(nu<ny) go to 300
                  if(u<tx(kx1) .or. u>tx(nkx1+1)) go to 300
                  !  the b-splinecoefficients of f(y) = s(u,y).
                  ier = 0
                  l = kx1
                  l1 = l+1
             110  if(u<tx(l1) .or. l==nkx1) go to 120
                  l = l1
                  l1 = l+1
                  go to 110
             120  call fpbspl(tx,nx,kx,u,l,h)
                  m0 = (l-kx1)*nky1+1
                  do 140 i=1,nky1
                    m = m0
                    sum = zero
                    do 130 j=1,kx1
                      sum = sum+h(j)*c(m)
                      m = m+nky1
             130    continue
                    cu(i) = sum
                    m0 = m0+1
             140  continue
                  go to 300

         case default

            if(nu<nx) go to 300
                  if(u<ty(ky1) .or. u>ty(nky1+1)) go to 300
                  !  the b-splinecoefficients of g(x) = s(x,u).
                  ier = 0
                  l = ky1
                  l1 = l+1
             210  if(u<ty(l1) .or. l==nky1) go to 220
                  l = l1
                  l1 = l+1
                  go to 210
             220  call fpbspl(ty,ny,ky,u,l,h)
                  m0 = l-ky
                  do 240 i=1,nkx1
                    m = m0
                    cu(i) = dot_product(h(1:ky1),c(m0:m0+ky))
                    m0 = m0+nky1
             240  continue

                  end select

 300  return
      end subroutine profil



      recursive subroutine regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s, &
       nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

      ! given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
      ! i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
      ! iate spline approximation s(x,y) of degrees kx and ky on the rect-
      ! angle xb <= x <= xe, yb <= y <= ye.
      ! if iopt = -1 regrid calculates the least-squares spline according
      ! to a given set of knots.
      ! if iopt >= 0 the total numbers nx and ny of these knots and their
      ! position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
      ! ally by the routine. the smoothness of s(x,y) is then achieved by
      ! minimalizing the discontinuity jumps in the derivatives of s(x,y)
      ! across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
      ! the amounth of smoothness is determined by the condition that f(p) =
      ! sum ((z(i,j)-s(x(i),y(j))))**2) be <= s, with s a given non-negative
      ! constant, called the smoothing factor.
      ! the fit is given in the b-spline representation (b-spline coefficients
      ! c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
      ! uated by means of subroutine bispev.
      !
      ! calling sequence:
      !     call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
      !    *  nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      !
      ! parameters:
      !  iopt  : integer flag. on entry iopt must specify whether a least-
      !          squares spline (iopt=-1) or a smoothing spline (iopt=0 or 1)
      !          must be determined.
      !          if iopt=0 the routine will start with an initial set of knots
      !          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
      !          1,...,ky+1. if iopt=1 the routine will continue with the set
      !          of knots found at the last call of the routine.
      !          attention: a call with iopt=1 must always be immediately pre-
      !                     ceded by another call with iopt=1 or iopt=0 and
      !                     s/=0.
      !          unchanged on exit.
      !  mx    : integer. on entry mx must specify the number of grid points
      !          along the x-axis. mx > kx . unchanged on exit.
      !  x     : real array of dimension at least (mx). before entry, x(i)
      !          must be set to the x-co-ordinate of the i-th grid point
      !          along the x-axis, for i=1,2,...,mx. these values must be
      !          supplied in strictly ascending order. unchanged on exit.
      !  my    : integer. on entry my must specify the number of grid points
      !          along the y-axis. my > ky . unchanged on exit.
      !  y     : real array of dimension at least (my). before entry, y(j)
      !          must be set to the y-co-ordinate of the j-th grid point
      !          along the y-axis, for j=1,2,...,my. these values must be
      !          supplied in strictly ascending order. unchanged on exit.
      !  z     : real array of dimension at least (mx*my).
      !          before entry, z(my*(i-1)+j) must be set to the data value at
      !          the grid point (x(i),y(j)) for i=1,...,mx and j=1,...,my.
      !          unchanged on exit.
      !  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
      !  yb,ye   aries of the rectangular approximation domain.
      !          xb<=x(i)<=xe,i=1,...,mx; yb<=y(j)<=ye,j=1,...,my.
      !          unchanged on exit.
      !  kx,ky : integer values. on entry kx and ky must specify the degrees
      !          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
      !          (kx=ky=3) splines. unchanged on exit.
      !  s     : real. on entry (in case iopt>=0) s must specify the smoothing
      !          factor. s >=0. unchanged on exit.
      !          for advice on the choice of s see further comments
      !  nxest : integer. unchanged on exit.
      !  nyest : integer. unchanged on exit.
      !          on entry, nxest and nyest must specify an upper bound for the
      !          number of knots required in the x- and y-directions respect.
      !          these numbers will also determine the storage space needed by
      !          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
      !          in most practical situation nxest = mx/2, nyest=my/2, will
      !          be sufficient. always large enough are nxest=mx+kx+1, nyest=
      !          my+ky+1, the number of knots needed for interpolation (s=0).
      !          see also further comments.
      !  nx    : integer.
      !          unless ier=10 (in case iopt >=0), nx will contain the total
      !          number of knots with respect to the x-variable, of the spline
      !          approximation returned. if the computation mode iopt=1 is
      !          used, the value of nx should be left unchanged between sub-
      !          sequent calls.
      !          in case iopt=-1, the value of nx should be specified on entry
      !  tx    : real array of dimension nmax.
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the x-variable, i.e. the position of
      !          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
      !          position of the additional knots tx(1)=...=tx(kx+1)=xb and
      !          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
      !          if the computation mode iopt=1 is used, the values of tx(1),
      !          ...,tx(nx) should be left unchanged between subsequent calls.
      !          if the computation mode iopt=-1 is used, the values tx(kx+2),
      !          ...tx(nx-kx-1) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  ny    : integer.
      !          unless ier=10 (in case iopt >=0), ny will contain the total
      !          number of knots with respect to the y-variable, of the spline
      !          approximation returned. if the computation mode iopt=1 is
      !          used, the value of ny should be left unchanged between sub-
      !          sequent calls.
      !          in case iopt=-1, the value of ny should be specified on entry
      !  ty    : real array of dimension nmax.
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the y-variable, i.e. the position of
      !          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
      !          position of the additional knots ty(1)=...=ty(ky+1)=yb and
      !          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
      !          if the computation mode iopt=1 is used, the values of ty(1),
      !          ...,ty(ny) should be left unchanged between subsequent calls.
      !          if the computation mode iopt=-1 is used, the values ty(ky+2),
      !          ...ty(ny-ky-1) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
      !          on successful exit, c contains the coefficients of the spline
      !          approximation s(x,y)
      !  fp    : real. unless ier=10, fp contains the sum of squared
      !          residuals of the spline approximation returned.
      !  wrk   : real array of dimension (lwrk). used as workspace.
      !          if the computation mode iopt=1 is used the values of wrk(1),
      !          ...,wrk(4) should be left unchanged between subsequent calls.
      !  lwrk  : integer. on entry lwrk must specify the actual dimension of
      !          the array wrk as declared in the calling (sub)program.
      !          lwrk must not be too small.
      !           lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
      !            my*(ky+1) +u
      !           where u is the larger of my and nxest.
      !  iwrk  : integer array of dimension (kwrk). used as workspace.
      !          if the computation mode iopt=1 is used the values of iwrk(1),
      !          ...,iwrk(3) should be left unchanged between subsequent calls
      !  kwrk  : integer. on entry kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program.
      !          kwrk >= 3+mx+my+nxest+nyest.
      !  ier   : integer. unless the routine detects an error, ier contains a
      !          non-positive value on exit, i.e.
      !   ier=0  : normal return. the spline returned has a residual sum of
      !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !            ive tolerance set to 0.001 by the program.
      !   ier=-1 : normal return. the spline returned is an interpolating
      !            spline (fp=0).
      !   ier=-2 : normal return. the spline returned is the least-squares
      !            polynomial of degrees kx and ky. in this extreme case fp
      !            gives the upper bound for the smoothing factor s.
      !   ier=1  : error. the required storage space exceeds the available
      !            storage space, as specified by the parameters nxest and
      !            nyest.
      !            probably causes : nxest or nyest too small. if these param-
      !            eters are already large, it may also indicate that s is
      !            too small
      !            the approximation returned is the least-squares spline
      !            according to the current set of knots. the parameter fp
      !            gives the corresponding sum of squared residuals (fp>s).
      !   ier=2  : error. a theoretically impossible result was found during
      !            the iteration process for finding a smoothing spline with
      !            fp = s. probably causes : s too small.
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=3  : error. the maximal number of iterations maxit (set to 20
      !            by the program) allowed for finding a smoothing spline
      !            with fp=s has been reached. probably causes : s too small
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=10 : error. on entry, the input data are controlled on validity
      !            the following restrictions must be satisfied.
      !            -1<=iopt<=1, 1<=kx,ky<=5, mx>kx, my>ky, nxest>=2*kx+2,
      !            nyest>=2*ky+2, kwrk>=3+mx+my+nxest+nyest,
      !            lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
      !             my*(ky+1) +max(my,nxest),
      !            xb<=x(i-1)<x(i)<=xe,i=2,..,mx,yb<=y(j-1)<y(j)<=ye,j=2,..,my
      !            if iopt=-1: 2*kx+2<=nx<=min(nxest,mx+kx+1)
      !                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
      !                        2*ky+2<=ny<=min(nyest,my+ky+1)
      !                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
      !                    the schoenberg-whitney conditions, i.e. there must
      !                    be subset of grid co-ordinates xx(p) and yy(q) such
      !                    that   tx(p) < xx(p) < tx(p+kx+1) ,p=1,...,nx-kx-1
      !                           ty(q) < yy(q) < ty(q+ky+1) ,q=1,...,ny-ky-1
      !            if iopt>=0: s>=0
      !                        if s=0 : nxest>=mx+kx+1, nyest>=my+ky+1
      !            if one of these conditions is found to be violated,control
      !            is immediately repassed to the calling program. in that
      !            case there is no approximation returned.
      !
      ! further comments:
      !   regrid does not allow individual weighting of the data-values.
      !   so, if these were determined to widely different accuracies, then
      !   perhaps the general data set routine surfit should rather be used
      !   in spite of efficiency.
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating spline if
      !   s=0 and the least-squares polynomial (degrees kx,ky) if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the accuracy of the data values.
      !   if the user has an idea of the statistical errors on the data, he
      !   can also find a proper estimate for s. for, by assuming that, if he
      !   specifies the right s, regrid will return a spline s(x,y) which
      !   exactly reproduces the function underlying the data he can evaluate
      !   the sum((z(i,j)-s(x(i),y(j)))**2) to find a good estimate for this s
      !   for example, if he knows that the statistical errors on his z(i,j)-
      !   values is not greater than 0.1, he may expect that a good s should
      !   have a value not larger than mx*my*(0.1)**2.
      !   if nothing is known about the statistical error in z(i,j), s must
      !   be determined by trial and error, taking account of the comments
      !   above. the best is then to start with a very large value of s (to
      !   determine the least-squares polynomial and the corresponding upper
      !   bound fp0 for s) and then to progressively decrease the value of s
      !   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
      !   and more carefully as the approximation shows more detail) to
      !   obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=0.
      !   if iopt=1 the program will continue with the set of knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if regrid is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. if the computation mode iopt=1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt=1, the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   regrid once more with the selected value for s but now with iopt=0.
      !   indeed, regrid may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !   the number of knots may also depend on the upper bounds nxest and
      !   nyest. indeed, if at a certain stage in regrid the number of knots
      !   in one direction (say nx) has reached the value of its upper bound
      !   (nxest), then from that moment on all subsequent knots are added
      !   in the other (y) direction. this may indicate that the value of
      !   nxest is too small. on the other hand, it gives the user the option
      !   of limiting the number of knots the routine locates in any direction
      !   for example, by setting nxest=2*kx+2 (the lowest allowable value for
      !   nxest), the user can indicate that he wants an approximation which
      !   is a simple polynomial of degree kx in the variable x.
      !
      !  other subroutines required:
      !    fpback,fpbspl,fpregr,fpdisc,fpgivs,fpgrre,fprati,fprota,fpchec,
      !    fpknot
      !
      !  references:
      !   dierckx p. : a fast algorithm for smoothing data on a rectangular
      !                grid while using spline functions, siam j.numer.anal.
      !                19 (1982) 1286-1304.
      !   dierckx p. : a fast algorithm for smoothing data on a rectangular
      !                grid while using spline functions, report tw53, dept.
      !                computer science,k.u.leuven, 1980.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !  latest update : march 1989
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) xb,xe,yb,ye,s,fp
      integer iopt,mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk,ier
      !  ..array arguments..
      real(RKIND) x(mx),y(my),z(mx*my),tx(nxest),ty(nyest), &
       c((nxest-kx-1)*(nyest-ky-1)),wrk(lwrk)
      integer iwrk(kwrk)
      !  ..local scalars..
      real(RKIND) tol
      integer i,j,jwrk,kndx,kndy,knrx,knry,kwest,kx1,kx2,ky1,ky2, &
       lfpx,lfpy,lwest,lww,maxit,nc,nminx,nminy,mz
      !  ..subroutine references..
      !    fpregr,fpchec
      !  ..
      !  we set up the parameters tol and maxit.
      maxit = 20
      tol = smallnum03
      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(kx<=0 .or. kx>5) go to 70
      kx1 = kx+1
      kx2 = kx1+1
      if(ky<=0 .or. ky>5) go to 70
      ky1 = ky+1
      ky2 = ky1+1
      if(iopt<(-1) .or. iopt>1) go to 70
      nminx = 2*kx1
      if(mx<kx1 .or. nxest<nminx) go to 70
      nminy = 2*ky1
      if(my<ky1 .or. nyest<nminy) go to 70
      mz = mx*my
      nc = (nxest-kx1)*(nyest-ky1)
      lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+ &
       my*ky1+max0(nxest,my)
      kwest = 3+mx+my+nxest+nyest
      if(lwrk<lwest .or. kwrk<kwest) go to 70
      if(xb>x(1) .or. xe<x(mx)) go to 70
      do 10 i=2,mx
        if(x(i-1)>=x(i)) go to 70
  10  continue
      if(yb>y(1) .or. ye<y(my)) go to 70
      do 20 i=2,my
        if(y(i-1)>=y(i)) go to 70
  20  continue
      if(iopt>=0) go to 50
      if(nx<nminx .or. nx>nxest) go to 70
      j = nx
      do 30 i=1,kx1
        tx(i) = xb
        tx(j) = xe
        j = j-1
  30  continue
      call fpchec(x,mx,tx,nx,kx,ier)
      if(ier/=0) go to 70
      if(ny<nminy .or. ny>nyest) go to 70
      j = ny
      do 40 i=1,ky1
        ty(i) = yb
        ty(j) = ye
        j = j-1
  40  continue
      call fpchec(y,my,ty,ny,ky,ier)
      if (ier==0) go to 60
      go to 70
  50  if(s<0.) go to 70
      if(s==zero .and. (nxest<(mx+kx1) .or. nyest<(my+ky1)) ) &
       go to 70
      ier = 0
      !  we partition the working space and determine the spline approximation
  60  lfpx = 5
      lfpy = lfpx+nxest
      lww = lfpy+nyest
      jwrk = lwrk-4-nxest-nyest
      knrx = 4
      knry = knrx+mx
      kndx = knry+my
      kndy = kndx+nxest
      call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
       tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4), &
       wrk(lfpx),wrk(lfpy),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx), &
       iwrk(knry),iwrk(kndx),iwrk(kndy),wrk(lww),jwrk,ier)
  70  return
      end subroutine regrid


      !  subroutine spalde evaluates at a point x ALL the derivatives
      !              (j-1)
      !      d(j) = s     (x) , j=1,2,...,k1
      !  of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline representation.
      pure subroutine spalde(t,n,c,k1,x,d,ier)

      !  calling sequence:
      !     call spalde(t,n,c,k1,x,d,ier)
      !
      !  input parameters:
      !    t    : array,length n, which contains the position of the knots.
      !    n    : integer, giving the total number of knots of s(x).
      !    c    : array,length n, which contains the b-spline coefficients.
      !    k1   : integer, giving the order of s(x) (order=degree+1)
      !    x    : real, which contains the point where the derivatives must
      !           be evaluated.
      !
      !  output parameters:
      !    d    : array,length k1, containing the derivative values of s(x).
      !    ier  : error flag
      !      ier = 0 : normal return
      !      ier =10 : invalid input data (see restrictions)
      !
      !  restrictions:
      !    t(k1) <= x <= t(n-k1+1)
      !
      !  further comments:
      !    if x coincides with a knot, right derivatives are computed
      !    ( left derivatives if x = t(n-k1+1) ).
      !
      !  other subroutines required: fpader.
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory
      !                6 (1972) 50-62.
      !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
      !                applics 10 (1972) 134-149.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1987
      !
      !  ..scalar arguments..
      integer, intent(in)      :: n,k1
      integer, intent(out)     :: ier

      real(RKIND), intent(in)  :: x
      !  ..array arguments..
      real(RKIND), intent(in)  :: t(n),c(n)
      real(RKIND), intent(out) :: d(k1)

      !  ..local scalars..
      integer :: l,nk1
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      nk1 = n-k1

      if(x<t(k1) .or. x>t(nk1+1)) return

      !  search for knot interval t(l) <= x < t(l+1)
      l = k1
      do while (.not.(x<t(l+1) .or. l==nk1))
         l = l+1
      end do

      if(t(l)>=t(l+1)) return

      !  calculate the derivatives.
      ier = FITPACK_OK
      call fpader(t,n,c,k1,x,l,d)

      end subroutine spalde


      recursive subroutine spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s, &
       nuest,nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)

      !  given the function values r(i,j) on the latitude-longitude grid
      !  (u(i),v(j)), i=1,...,mu ; j=1,...,mv , spgrid determines a smooth
      !  bicubic spline approximation on the rectangular domain 0<=u<=pi,
      !  vb<=v<=ve (vb = v(1), ve=vb+2*pi).
      !  this approximation s(u,v) will satisfy the properties
      !
      !    (1) s(0,v) = s(0,0) = dr(1)
      !
      !        d s(0,v)           d s(0,0)           d s(0,pi/2)
      !    (2) -------- = cos(v)* -------- + sin(v)* -----------
      !        d u                d u                d u
      !
      !                 = cos(v)*dr(2)+sin(v)*dr(3)
      !                                                     vb <= v <= ve
      !    (3) s(pi,v) = s(pi,0) = dr(4)
      !
      !        d s(pi,v)           d s(pi,0)           d s(pi,pi/2)
      !    (4) -------- = cos(v)*  --------- + sin(v)* ------------
      !        d u                 d u                 d u
      !
      !                 = cos(v)*dr(5)+sin(v)*dr(6)
      !
      !  and will be periodic in the variable v, i.e.
      !
      !         j           j
      !        d s(u,vb)   d s(u,ve)
      !    (5) --------- = ---------   0 <=u<= pi , j=0,1,2
      !           j           j
      !        d v         d v
      !
      !  the number of knots of s(u,v) and their position tu(i),i=1,2,...,nu;
      !  tv(j),j=1,2,...,nv, is chosen automatically by the routine. the
      !  smoothness of s(u,v) is achieved by minimalizing the discontinuity
      !  jumps of the derivatives of the spline at the knots. the amount of
      !  smoothness of s(u,v) is determined by the condition that
      !  fp=sumi=1,mu(sumj=1,mv((r(i,j)-s(u(i),v(j)))**2))+(r0-s(0,v))**2
      !  + (r1-s(pi,v))**2 <= s, with s a given non-negative constant.
      !  the fit s(u,v) is given in its b-spline representation and can be
      !  evaluated by means of routine bispev
      !
      ! calling sequence:
      !     call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,nu,tu,
      !    *  ,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      !
      ! parameters:
      !  iopt  : integer array of dimension 3, specifying different options.
      !          unchanged on exit.
      !  iopt(1):on entry iopt(1) must specify whether a least-squares spline
      !          (iopt(1)=-1) or a smoothing spline (iopt(1)=0 or 1) must be
      !          determined.
      !          if iopt(1)=0 the routine will start with an initial set of
      !          knots tu(i)=0,tu(i+4)=pi,i=1,...,4;tv(i)=v(1)+(i-4)*2*pi,
      !          i=1,...,8.
      !          if iopt(1)=1 the routine will continue with the set of knots
      !          found at the last call of the routine.
      !          attention: a call with iopt(1)=1 must always be immediately
      !          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
      !  iopt(2):on entry iopt(2) must specify the requested order of conti-
      !          nuity at the pole u=0.
      !          if iopt(2)=0 only condition (1) must be fulfilled and
      !          if iopt(2)=1 conditions (1)+(2) must be fulfilled.
      !  iopt(3):on entry iopt(3) must specify the requested order of conti-
      !          nuity at the pole u=pi.
      !          if iopt(3)=0 only condition (3) must be fulfilled and
      !          if iopt(3)=1 conditions (3)+(4) must be fulfilled.
      !  ider  : integer array of dimension 4, specifying different options.
      !          unchanged on exit.
      !  ider(1):on entry ider(1) must specify whether (ider(1)=0 or 1) or not
      !          (ider(1)=-1) there is a data value r0 at the pole u=0.
      !          if ider(1)=1, r0 will be considered to be the right function
      !          value, and it will be fitted exactly (s(0,v)=r0).
      !          if ider(1)=0, r0 will be considered to be a data value just
      !          like the other data values r(i,j).
      !  ider(2):on entry ider(2) must specify whether (ider(2)=1) or not
      !          (ider(2)=0) the approximation has vanishing derivatives
      !          dr(2) and dr(3) at the pole u=0  (in case iopt(2)=1)
      !  ider(3):on entry ider(3) must specify whether (ider(3)=0 or 1) or not
      !          (ider(3)=-1) there is a data value r1 at the pole u=pi.
      !          if ider(3)=1, r1 will be considered to be the right function
      !          value, and it will be fitted exactly (s(pi,v)=r1).
      !          if ider(3)=0, r1 will be considered to be a data value just
      !          like the other data values r(i,j).
      !  ider(4):on entry ider(4) must specify whether (ider(4)=1) or not
      !          (ider(4)=0) the approximation has vanishing derivatives
      !          dr(5) and dr(6) at the pole u=pi (in case iopt(3)=1)
      !  mu    : integer. on entry mu must specify the number of grid points
      !          along the u-axis. unchanged on exit.
      !          mu >= 1, mu >=mumin=4-i0-i1-ider(2)-ider(4) with
      !            i0=min(1,ider(1)+1), i1=min(1,ider(3)+1)
      !  u     : real array of dimension at least (mu). before entry, u(i)
      !          must be set to the u-co-ordinate of the i-th grid point
      !          along the u-axis, for i=1,2,...,mu. these values must be
      !          supplied in strictly ascending order. unchanged on exit.
      !          0 < u(i) < pi.
      !  mv    : integer. on entry mv must specify the number of grid points
      !          along the v-axis. mv > 3 . unchanged on exit.
      !  v     : real array of dimension at least (mv). before entry, v(j)
      !          must be set to the v-co-ordinate of the j-th grid point
      !          along the v-axis, for j=1,2,...,mv. these values must be
      !          supplied in strictly ascending order. unchanged on exit.
      !          -pi <= v(1) < pi , v(mv) < v(1)+2*pi.
      !  r     : real array of dimension at least (mu*mv).
      !          before entry, r(mv*(i-1)+j) must be set to the data value at
      !          the grid point (u(i),v(j)) for i=1,...,mu and j=1,...,mv.
      !          unchanged on exit.
      !  r0    : real value. on entry (if ider(1) >=0 ) r0 must specify the
      !          data value at the pole u=0. unchanged on exit.
      !  r1    : real value. on entry (if ider(1) >=0 ) r1 must specify the
      !          data value at the pole u=pi. unchanged on exit.
      !  s     : real. on entry (if iopt(1)>=0) s must specify the smoothing
      !          factor. s >=0. unchanged on exit.
      !          for advice on the choice of s see further comments
      !  nuest : integer. unchanged on exit.
      !  nvest : integer. unchanged on exit.
      !          on entry, nuest and nvest must specify an upper bound for the
      !          number of knots required in the u- and v-directions respect.
      !          these numbers will also determine the storage space needed by
      !          the routine. nuest >= 8, nvest >= 8.
      !          in most practical situation nuest = mu/2, nvest=mv/2, will
      !          be sufficient. always large enough are nuest=mu+6+iopt(2)+
      !          iopt(3), nvest = mv+7, the number of knots needed for
      !          interpolation (s=0). see also further comments.
      !  nu    : integer.
      !          unless ier=10 (in case iopt(1)>=0), nu will contain the total
      !          number of knots with respect to the u-variable, of the spline
      !          approximation returned. if the computation mode iopt(1)=1 is
      !          used, the value of nu should be left unchanged between sub-
      !          sequent calls. in case iopt(1)=-1, the value of nu should be
      !          specified on entry.
      !  tu    : real array of dimension at least (nuest).
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the u-variable, i.e. the position of
      !          the interior knots tu(5),...,tu(nu-4) as well as the position
      !          of the additional knots tu(1)=...=tu(4)=0 and tu(nu-3)=...=
      !          tu(nu)=pi needed for the b-spline representation.
      !          if the computation mode iopt(1)=1 is used,the values of tu(1)
      !          ...,tu(nu) should be left unchanged between subsequent calls.
      !          if the computation mode iopt(1)=-1 is used, the values tu(5),
      !          ...tu(nu-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  nv    : integer.
      !          unless ier=10 (in case iopt(1)>=0), nv will contain the total
      !          number of knots with respect to the v-variable, of the spline
      !          approximation returned. if the computation mode iopt(1)=1 is
      !          used, the value of nv should be left unchanged between sub-
      !          sequent calls. in case iopt(1) = -1, the value of nv should
      !          be specified on entry.
      !  tv    : real array of dimension at least (nvest).
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the v-variable, i.e. the position of
      !          the interior knots tv(5),...,tv(nv-4) as well as the position
      !          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
      !          tv(nv) needed for the b-spline representation.
      !          if the computation mode iopt(1)=1 is used,the values of tv(1)
      !          ...,tv(nv) should be left unchanged between subsequent calls.
      !          if the computation mode iopt(1)=-1 is used, the values tv(5),
      !          ...tv(nv-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  c     : real array of dimension at least (nuest-4)*(nvest-4).
      !          on successful exit, c contains the coefficients of the spline
      !          approximation s(u,v)
      !  fp    : real. unless ier=10, fp contains the sum of squared
      !          residuals of the spline approximation returned.
      !  wrk   : real array of dimension (lwrk). used as workspace.
      !          if the computation mode iopt(1)=1 is used the values of
      !          wrk(1),..,wrk(12) should be left unchanged between subsequent
      !          calls.
      !  lwrk  : integer. on entry lwrk must specify the actual dimension of
      !          the array wrk as declared in the calling (sub)program.
      !          lwrk must not be too small.
      !           lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+q
      !           where q is the larger of (mv+nvest) and nuest.
      !  iwrk  : integer array of dimension (kwrk). used as workspace.
      !          if the computation mode iopt(1)=1 is used the values of
      !          iwrk(1),.,iwrk(5) should be left unchanged between subsequent
      !          calls.
      !  kwrk  : integer. on entry kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program.
      !          kwrk >= 5+mu+mv+nuest+nvest.
      !  ier   : integer. unless the routine detects an error, ier contains a
      !          non-positive value on exit, i.e.
      !   ier=0  : normal return. the spline returned has a residual sum of
      !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !            ive tolerance set to 0.001 by the program.
      !   ier=-1 : normal return. the spline returned is an interpolating
      !            spline (fp=0).
      !   ier=-2 : normal return. the spline returned is the least-squares
      !            constrained polynomial. in this extreme case fp gives the
      !            upper bound for the smoothing factor s.
      !   ier=1  : error. the required storage space exceeds the available
      !            storage space, as specified by the parameters nuest and
      !            nvest.
      !            probably causes : nuest or nvest too small. if these param-
      !            eters are already large, it may also indicate that s is
      !            too small
      !            the approximation returned is the least-squares spline
      !            according to the current set of knots. the parameter fp
      !            gives the corresponding sum of squared residuals (fp>s).
      !   ier=2  : error. a theoretically impossible result was found during
      !            the iteration process for finding a smoothing spline with
      !            fp = s. probably causes : s too small.
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=3  : error. the maximal number of iterations maxit (set to 20
      !            by the program) allowed for finding a smoothing spline
      !            with fp=s has been reached. probably causes : s too small
      !            there is an approximation returned but the corresponding
      !            sum of squared residuals does not satisfy the condition
      !            abs(fp-s)/s < tol.
      !   ier=10 : error. on entry, the input data are controlled on validity
      !            the following restrictions must be satisfied.
      !            -1<=iopt(1)<=1, 0<=iopt(2)<=1, 0<=iopt(3)<=1,
      !            -1<=ider(1)<=1, 0<=ider(2)<=1, ider(2)=0 if iopt(2)=0.
      !            -1<=ider(3)<=1, 0<=ider(4)<=1, ider(4)=0 if iopt(3)=0.
      !            mu >= mumin (see above), mv >= 4, nuest >=8, nvest >= 8,
      !            kwrk>=5+mu+mv+nuest+nvest,
      !            lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+
      !             max(nuest,mv+nvest)
      !            0< u(i-1)<u(i)< pi,i=2,..,mu,
      !            -pi<=v(1)< pi, v(1)<v(i-1)<v(i)<v(1)+2*pi, i=3,...,mv
      !            if iopt(1)=-1: 8<=nu<=min(nuest,mu+6+iopt(2)+iopt(3))
      !                           0<tu(5)<tu(6)<...<tu(nu-4)< pi
      !                           8<=nv<=min(nvest,mv+7)
      !                           v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(1)+2*pi
      !                    the schoenberg-whitney conditions, i.e. there must
      !                    be subset of grid co-ordinates uu(p) and vv(q) such
      !                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
      !                     (iopt(2)=1 and iopt(3)=1 also count for a uu-value
      !                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
      !                     (vv(q) is either a value v(j) or v(j)+2*pi)
      !            if iopt(1)>=0: s>=0
      !                       if s=0: nuest>=mu+6+iopt(2)+iopt(3), nvest>=mv+7
      !            if one of these conditions is found to be violated,control
      !            is immediately repassed to the calling program. in that
      !            case there is no approximation returned.
      !
      ! further comments:
      !   spgrid does not allow individual weighting of the data-values.
      !   so, if these were determined to widely different accuracies, then
      !   perhaps the general data set routine sphere should rather be used
      !   in spite of efficiency.
      !   by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating spline if
      !   s=0 and the constrained least-squares polynomial(degrees 3,0)if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the accuracy of the data values.
      !   if the user has an idea of the statistical errors on the data, he
      !   can also find a proper estimate for s. for, by assuming that, if he
      !   specifies the right s, spgrid will return a spline s(u,v) which
      !   exactly reproduces the function underlying the data he can evaluate
      !   the sum((r(i,j)-s(u(i),v(j)))**2) to find a good estimate for this s
      !   for example, if he knows that the statistical errors on his r(i,j)-
      !   values is not greater than 0.1, he may expect that a good s should
      !   have a value not larger than mu*mv*(0.1)**2.
      !   if nothing is known about the statistical error in r(i,j), s must
      !   be determined by trial and error, taking account of the comments
      !   above. the best is then to start with a very large value of s (to
      !   determine the least-squares polynomial and the corresponding upper
      !   bound fp0 for s) and then to progressively decrease the value of s
      !   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
      !   and more carefully as the approximation shows more detail) to
      !   obtain closer fits.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt(1)=0.
      !   if iopt(1) = 1 the program will continue with the knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if spgrid is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. if the computation mode iopt(1) = 1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt(1)=1,the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   spgrid once more with the chosen value for s but now with iopt(1)=0.
      !   indeed, spgrid may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !   the number of knots may also depend on the upper bounds nuest and
      !   nvest. indeed, if at a certain stage in spgrid the number of knots
      !   in one direction (say nu) has reached the value of its upper bound
      !   (nuest), then from that moment on all subsequent knots are added
      !   in the other (v) direction. this may indicate that the value of
      !   nuest is too small. on the other hand, it gives the user the option
      !   of limiting the number of knots the routine locates in any direction
      !   for example, by setting nuest=8 (the lowest allowable value for
      !   nuest), the user can indicate that he wants an approximation which
      !   is a simple cubic polynomial in the variable u.
      !
      !  other subroutines required:
      !    fpspgr,fpchec,fpchep,fpknot,fpopsp,fprati,fpgrsp,fpsysy,fpback,
      !    fpbacp,fpbspl,fpcyt1,fpcyt2,fpdisc,fpgivs,fprota
      !
      !  references:
      !   dierckx p. : fast algorithms for smoothing data over a disc or a
      !                sphere using tensor product splines, in "algorithms
      !                for approximation", ed. j.c.mason and m.g.cox,
      !                clarendon press oxford, 1987, pp. 51-65
      !   dierckx p. : fast algorithms for smoothing data over a disc or a
      !                sphere using tensor product splines, report tw73, dept.
      !                computer science,k.u.leuven, 1985.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : july 1985
      !  latest update : march 1989
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) r0,r1,s,fp
      integer mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
      !  ..array arguments..
      integer iopt(3),ider(4),iwrk(kwrk)
      real(RKIND) u(mu),v(mv),r(mu*mv),c((nuest-4)*(nvest-4)),tu(nuest), &
       tv(nvest),wrk(lwrk)
      !  ..local scalars..
      real(RKIND) per,tol,uu,ve,rmax,rmin,rn,rb,re
      integer i,i1,i2,j,jwrk,j1,j2,kndu,kndv,knru,knrv,kwest,l, &
       ldr,lfpu,lfpv,lwest,lww,m,maxit,mumin,muu,nc

      !  set constants
      per = pi+pi
      ve = v(1)+per
      !  we set up the parameters tol and maxit.
      maxit = 20
      tol = smallnum03
      !  before starting computations, a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt(1)<(-1) .or. iopt(1)>1) go to 200
      if(iopt(2)<0 .or. iopt(2)>1) go to 200
      if(iopt(3)<0 .or. iopt(3)>1) go to 200
      if(ider(1)<(-1) .or. ider(1)>1) go to 200
      if(ider(2)<0 .or. ider(2)>1) go to 200
      if(ider(2)==1 .and. iopt(2)==0) go to 200
      if(ider(3)<(-1) .or. ider(3)>1) go to 200
      if(ider(4)<0 .or. ider(4)>1) go to 200
      if(ider(4)==1 .and. iopt(3)==0) go to 200
      mumin = 4
      if(ider(1)>=0) mumin = mumin-1
      if(iopt(2)==1 .and. ider(2)==1) mumin = mumin-1
      if(ider(3)>=0) mumin = mumin-1
      if(iopt(3)==1 .and. ider(4)==1) mumin = mumin-1
      if(mumin==0) mumin = 1
      if(mu<mumin .or. mv<4) go to 200
      if(nuest<8 .or. nvest<8) go to 200
      m = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 12+nuest*(mv+nvest+3)+24*nvest+4*mu+8*mv+ &
       max0(nuest,mv+nvest)
      kwest = 5+mu+mv+nuest+nvest
      if(lwrk<lwest .or. kwrk<kwest) go to 200
      if(u(1)<=0. .or. u(mu)>=pi) go to 200
      if(mu==1) go to 30
      do 20 i=2,mu
        if(u(i-1)>=u(i)) go to 200
  20  continue
  30  if(v(1)< (-pi) .or. v(1)>=pi ) go to 200
      if(v(mv)>=v(1)+per) go to 200
      do 40 i=2,mv
        if(v(i-1)>=v(i)) go to 200
  40  continue
      if(iopt(1)>0) go to 140
      !  if not given, we compute an estimate for r0.
      rn = mv
      if(ider(1)<0) go to 45
      rb = r0
      go to 55
  45  rb = sum(r(1:mv))/rn
      !  if not given, we compute an estimate for r1.
  55  if(ider(3)<0) go to 60
      re = r1
      go to 70
  60  re = zero
      j = m
      do 65 i=1,mv
         re = re+r(j)
         j = j-1
  65  continue
      re = re/rn
      !  we determine the range of r-values.
  70  rmin = rb
      rmax = re
      do 80 i=1,m
         if(r(i)<rmin) rmin = r(i)
         if(r(i)>rmax) rmax = r(i)
  80  continue
      wrk(5) = rb
      wrk(6) = zero
      wrk(7) = zero
      wrk(8) = re
      wrk(9) = zero
      wrk(10) = zero
      wrk(11) = rmax -rmin
      wrk(12) = wrk(11)
      iwrk(4) = mu
      iwrk(5) = mu
      if(iopt(1)==0) go to 140
      if(nu<8 .or. nu>nuest) go to 200
      if(nv<11 .or. nv>nvest) go to 200
      j = nu
      do 90 i=1,4
        tu(i) = zero
        tu(j) = pi
        j = j-1
  90  continue
      l = 13
      wrk(l) = zero
      if(iopt(2)==0) go to 100
      l = l+1
      uu = u(1)
      if(uu>tu(5)) uu = tu(5)
      wrk(l) = uu*half
 100  do 110 i=1,mu
        l = l+1
        wrk(l) = u(i)
 110  continue
      if(iopt(3)==0) go to 120
      l = l+1
      uu = u(mu)
      if(uu<tu(nu-4)) uu = tu(nu-4)
      wrk(l) = uu+(pi-uu)*half
 120  l = l+1
      wrk(l) = pi
      muu = l-12
      call fpchec(wrk(13),muu,tu,nu,3,ier)
      if(ier/=0) go to 200
      j1 = 4
      tv(j1) = v(1)
      i1 = nv-3
      tv(i1) = ve
      j2 = j1
      i2 = i1
      do 130 i=1,3
        i1 = i1+1
        i2 = i2-1
        j1 = j1+1
        j2 = j2-1
        tv(j2) = tv(i2)-per
        tv(i1) = tv(j1)+per
 130  continue
      l = 13
      do 135 i=1,mv
        wrk(l) = v(i)
        l = l+1
 135  continue
      wrk(l) = ve
      ier = fpchep(wrk(13),mv+1,tv,nv,3)
      if (ier==0) go to 150
      go to 200
 140  if(s<0.) go to 200
      if(s==zero .and. (nuest<(mu+6+iopt(2)+iopt(3)) .or. &
       nvest<(mv+7)) ) go to 200
      !  we partition the working space and determine the spline approximation
 150  ldr = 5
      lfpu = 13
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-12-nuest-nvest
      knru = 6
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fpspgr(iopt,ider,u,mu,v,mv,r,m,rb,re,s,nuest,nvest,tol,maxit, &
                  nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),wrk(lfpu), &
                  wrk(lfpv),wrk(ldr),wrk(11),iwrk(1),iwrk(2),iwrk(3),iwrk(4), &
                  iwrk(5),iwrk(knru),iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww), &
                  jwrk,ier)
 200  return
      end subroutine spgrid


      recursive subroutine sphere(iopt,m,teta,phi,r,w,s,ntest,npest, &
        eps,nt,tt,np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

      !  subroutine sphere determines a smooth bicubic spherical spline
      !  approximation s(teta,phi), 0 <= teta <= pi ; 0 <= phi <= 2*pi
      !  to a given set of data points (teta(i),phi(i),r(i)),i=1,2,...,m.
      !  such a spline has the following specific properties
      !
      !    (1) s(0,phi)  = constant   0 <=phi<= 2*pi.
      !
      !    (2) s(pi,phi) = constant   0 <=phi<= 2*pi
      !
      !         j             j
      !        d s(teta,0)   d s(teta,2*pi)
      !    (3) ----------- = ------------   0 <=teta<=pi, j=0,1,2
      !             j             j
      !        d phi         d phi
      !
      !        d s(0,phi)    d s(0,0)             d s(0,pi/2)
      !    (4) ----------  = -------- *cos(phi) + ----------- *sin(phi)
      !        d teta        d teta               d teta
      !
      !        d s(pi,phi)   d s(pi,0)            d s(pi,pi/2)
      !    (5) ----------- = ---------*cos(phi) + ------------*sin(phi)
      !        d teta        d teta               d teta
      !
      !  if iopt =-1 sphere calculates a weighted least-squares spherical
      !  spline according to a given set of knots in teta- and phi- direction.
      !  if iopt >=0, the number of knots in each direction and their position
      !  tt(j),j=1,2,...,nt ; tp(j),j=1,2,...,np are chosen automatically by
      !  the routine. the smoothness of s(teta,phi) is then achieved by mini-
      !  malizing the discontinuity jumps of the derivatives of the spline
      !  at the knots. the amount of smoothness of s(teta,phi) is determined
      !  by the condition that fp = sum((w(i)*(r(i)-s(teta(i),phi(i))))**2)
      !  be <= s, with s a given non-negative constant.
      !  the spherical spline is given in the standard b-spline representation
      !  of bicubic splines and can be evaluated by means of subroutine bispev
      !
      ! calling sequence:
      !     call sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,
      !    *  nt,tt,np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      !
      ! parameters:
      !  iopt  : integer flag. on entry iopt must specify whether a weighted
      !          least-squares spherical spline (iopt=-1) or a smoothing
      !          spherical spline (iopt=0 or 1) must be determined.
      !          if iopt=0 the routine will start with an initial set of knots
      !          tt(i)=0,tt(i+4)=pi,i=1,...,4;tp(i)=0,tp(i+4)=2*pi,i=1,...,4.
      !          if iopt=1 the routine will continue with the set of knots
      !          found at the last call of the routine.
      !          attention: a call with iopt=1 must always be immediately pre-
      !                     ceded by another call with iopt=1 or iopt=0.
      !          unchanged on exit.
      !  m     : integer. on entry m must specify the number of data points.
      !          m >= 2. unchanged on exit.
      !  teta  : real array of dimension at least (m).
      !  phi   : real array of dimension at least (m).
      !  r     : real array of dimension at least (m).
      !          before entry,teta(i),phi(i),r(i) must be set to the spherical
      !          co-ordinates of the i-th data point, for i=1,...,m.the order
      !          of the data points is immaterial. unchanged on exit.
      !  w     : real array of dimension at least (m). before entry, w(i) must
      !          be set to the i-th value in the set of weights. the w(i) must
      !          be strictly positive. unchanged on exit.
      !  s     : real. on entry (in case iopt>=0) s must specify the smoothing
      !          factor. s >=0. unchanged on exit.
      !          for advice on the choice of s see further comments
      !  ntest : integer. unchanged on exit.
      !  npest : integer. unchanged on exit.
      !          on entry, ntest and npest must specify an upper bound for the
      !          number of knots required in the teta- and phi-directions.
      !          these numbers will also determine the storage space needed by
      !          the routine. ntest >= 8, npest >= 8.
      !          in most practical situation ntest = npest = 8+sqrt(m/2) will
      !          be sufficient. see also further comments.
      !  eps   : real.
      !          on entry, eps must specify a threshold for determining the
      !          effective rank of an over-determined linear system of equat-
      !          ions. 0 < eps < 1.  if the number of decimal digits in the
      !          computer representation of a real number is q, then 10**(-q)
      !          is a suitable value for eps in most practical applications.
      !          unchanged on exit.
      !  nt    : integer.
      !          unless ier=10 (in case iopt >=0), nt will contain the total
      !          number of knots with respect to the teta-variable, of the
      !          spline approximation returned. if the computation mode iopt=1
      !          is used, the value of nt should be left unchanged between
      !          subsequent calls.
      !          in case iopt=-1, the value of nt should be specified on entry
      !  tt    : real array of dimension at least ntest.
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the teta-variable, i.e. the position
      !          of the interior knots tt(5),...,tt(nt-4) as well as the
      !          position of the additional knots tt(1)=...=tt(4)=0 and
      !          tt(nt-3)=...=tt(nt)=pi needed for the b-spline representation
      !          if the computation mode iopt=1 is used, the values of tt(1),
      !          ...,tt(nt) should be left unchanged between subsequent calls.
      !          if the computation mode iopt=-1 is used, the values tt(5),
      !          ...tt(nt-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  np    : integer.
      !          unless ier=10 (in case iopt >=0), np will contain the total
      !          number of knots with respect to the phi-variable, of the
      !          spline approximation returned. if the computation mode iopt=1
      !          is used, the value of np should be left unchanged between
      !          subsequent calls.
      !          in case iopt=-1, the value of np (>=9) should be specified
      !          on entry.
      !  tp    : real array of dimension at least npest.
      !          on successful exit, this array will contain the knots of the
      !          spline with respect to the phi-variable, i.e. the position of
      !          the interior knots tp(5),...,tp(np-4) as well as the position
      !          of the additional knots tp(1),...,tp(4) and tp(np-3),...,
      !          tp(np) needed for the b-spline representation.
      !          if the computation mode iopt=1 is used, the values of tp(1),
      !          ...,tp(np) should be left unchanged between subsequent calls.
      !          if the computation mode iopt=-1 is used, the values tp(5),
      !          ...tp(np-4) must be supplied by the user, before entry.
      !          see also the restrictions (ier=10).
      !  c     : real array of dimension at least (ntest-4)*(npest-4).
      !          on successful exit, c contains the coefficients of the spline
      !          approximation s(teta,phi).
      !  fp    : real. unless ier=10, fp contains the weighted sum of
      !          squared residuals of the spline approximation returned.
      !  wrk1  : real array of dimension (lwrk1). used as workspace.
      !          if the computation mode iopt=1 is used the value of wrk1(1)
      !          should be left unchanged between subsequent calls.
      !          on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the
      !          values d(i)/max(d(i)),i=1,...,ncof=6+(np-7)*(nt-8)
      !          with d(i) the i-th diagonal element of the reduced triangular
      !          matrix for calculating the b-spline coefficients. it includes
      !          those elements whose square is less than eps,which are treat-
      !          ed as 0 in the case of presumed rank deficiency (ier<-2).
      !  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
      !          the array wrk1 as declared in the calling (sub)program.
      !          lwrk1 must not be too small. let
      !            u = ntest-7, v = npest-7, then
      !          lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
      !  wrk2  : real array of dimension (lwrk2). used as workspace, but
      !          only in the case a rank deficient system is encountered.
      !  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
      !          the array wrk2 as declared in the calling (sub)program.
      !          lwrk2 > 0 . a save upper bound  for lwrk2 = 48+21*v+7*u*v+
      !          4*(u-1)*v**2 where u,v are as above. if there are enough data
      !          points, scattered uniformly over the approximation domain
      !          and if the smoothing factor s is not too small, there is a
      !          good chance that this extra workspace is not needed. a lot
      !          of memory might therefore be saved by setting lwrk2=1.
      !          (see also ier > 10)
      !  iwrk  : integer array of dimension (kwrk). used as workspace.
      !  kwrk  : integer. on entry kwrk must specify the actual dimension of
      !          the array iwrk as declared in the calling (sub)program.
      !          kwrk >= m+(ntest-7)*(npest-7).
      !  ier   : integer. unless the routine detects an error, ier contains a
      !          non-positive value on exit, i.e.
      !   ier=0  : normal return. the spline returned has a residual sum of
      !            squares fp such that abs(fp-s)/s <= tol with tol a relat-
      !            ive tolerance set to 0.001 by the program.
      !   ier=-1 : normal return. the spline returned is a spherical
      !            interpolating spline (fp=0).
      !   ier=-2 : normal return. the spline returned is the weighted least-
      !            squares constrained polynomial . in this extreme case
      !            fp gives the upper bound for the smoothing factor s.
      !   ier<-2 : warning. the coefficients of the spline returned have been
      !            computed as the minimal norm least-squares solution of a
      !            (numerically) rank deficient system. (-ier) gives the rank.
      !            especially if the rank deficiency which can be computed as
      !            6+(nt-8)*(np-7)+ier, is large the results may be inaccurate
      !            they could also seriously depend on the value of eps.
      !   ier=1  : error. the required storage space exceeds the available
      !            storage space, as specified by the parameters ntest and
      !            npest.
      !            probably causes : ntest or npest too small. if these param-
      !            eters are already large, it may also indicate that s is
      !            too small
      !            the approximation returned is the weighted least-squares
      !            spherical spline according to the current set of knots.
      !            the parameter fp gives the corresponding weighted sum of
      !            squared residuals (fp>s).
      !   ier=2  : error. a theoretically impossible result was found during
      !            the iteration process for finding a smoothing spline with
      !            fp = s. probably causes : s too small or badly chosen eps.
      !            there is an approximation returned but the corresponding
      !            weighted sum of squared residuals does not satisfy the
      !            condition abs(fp-s)/s < tol.
      !   ier=3  : error. the maximal number of iterations maxit (set to 20
      !            by the program) allowed for finding a smoothing spline
      !            with fp=s has been reached. probably causes : s too small
      !            there is an approximation returned but the corresponding
      !            weighted sum of squared residuals does not satisfy the
      !            condition abs(fp-s)/s < tol.
      !   ier=4  : error. no more knots can be added because the dimension
      !            of the spherical spline 6+(nt-8)*(np-7) already exceeds
      !            the number of data points m.
      !            probably causes : either s or m too small.
      !            the approximation returned is the weighted least-squares
      !            spherical spline according to the current set of knots.
      !            the parameter fp gives the corresponding weighted sum of
      !            squared residuals (fp>s).
      !   ier=5  : error. no more knots can be added because the additional
      !            knot would (quasi) coincide with an old one.
      !            probably causes : s too small or too large a weight to an
      !            inaccurate data point.
      !            the approximation returned is the weighted least-squares
      !            spherical spline according to the current set of knots.
      !            the parameter fp gives the corresponding weighted sum of
      !            squared residuals (fp>s).
      !   ier=10 : error. on entry, the input data are controlled on validity
      !            the following restrictions must be satisfied.
      !            -1<=iopt<=1,  m>=2, ntest>=8 ,npest >=8, 0<eps<1,
      !            0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m
      !            lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
      !            kwrk >= m+(ntest-7)*(npest-7)
      !            if iopt=-1: 8<=nt<=ntest , 9<=np<=npest
      !                        0<tt(5)<tt(6)<...<tt(nt-4)<pi
      !                        0<tp(5)<tp(6)<...<tp(np-4)<2*pi
      !            if iopt>=0: s>=0
      !            if one of these conditions is found to be violated,control
      !            is immediately repassed to the calling program. in that
      !            case there is no approximation returned.
      !   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
      !            space for computing the minimal least-squares solution of
      !            a rank deficient system of linear equations. ier gives the
      !            requested value for lwrk2. there is no approximation re-
      !            turned but, having saved the information contained in nt,
      !            np,tt,tp,wrk1, and having adjusted the value of lwrk2 and
      !            the dimension of the array wrk2 accordingly, the user can
      !            continue at the point the program was left, by calling
      !            sphere with iopt=1.
      !
      ! further comments:
      !  by means of the parameter s, the user can control the tradeoff
      !   between closeness of fit and smoothness of fit of the approximation.
      !   if s is too large, the spline will be too smooth and signal will be
      !   lost ; if s is too small the spline will pick up too much noise. in
      !   the extreme cases the program will return an interpolating spline if
      !   s=0 and the constrained weighted least-squares polynomial if s is
      !   very large. between these extremes, a properly chosen s will result
      !   in a good compromise between closeness of fit and smoothness of fit.
      !   to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits
      !   graphically.
      !   recommended values for s depend on the weights w(i). if these are
      !   taken as 1/d(i) with d(i) an estimate of the standard deviation of
      !   r(i), a good s-value should be found in the range (m-sqrt(2*m),m+
      !   sqrt(2*m)). if nothing is known about the statistical error in r(i)
      !   each w(i) can be set equal to one and s determined by trial and
      !   error, taking account of the comments above. the best is then to
      !   start with a very large value of s ( to determine the least-squares
      !   polynomial and the corresponding upper bound fp0 for s) and then to
      !   progressively decrease the value of s ( say by a factor 10 in the
      !   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
      !   approximation shows more detail) to obtain closer fits.
      !   to choose s very small is strongly discouraged. this considerably
      !   increases computation time and memory requirements. it may also
      !   cause rank-deficiency (ier<-2) and endager numerical stability.
      !   to economize the search for a good s-value the program provides with
      !   different modes of computation. at the first call of the routine, or
      !   whenever he wants to restart with the initial set of knots the user
      !   must set iopt=0.
      !   if iopt=1 the program will continue with the set of knots found at
      !   the last call of the routine. this will save a lot of computation
      !   time if sphere is called repeatedly for different values of s.
      !   the number of knots of the spline returned and their location will
      !   depend on the value of s and on the complexity of the shape of the
      !   function underlying the data. if the computation mode iopt=1
      !   is used, the knots returned may also depend on the s-values at
      !   previous calls (if these were smaller). therefore, if after a number
      !   of trials with different s-values and iopt=1, the user can finally
      !   accept a fit as satisfactory, it may be worthwhile for him to call
      !   sphere once more with the selected value for s but now with iopt=0.
      !   indeed, sphere may then return an approximation of the same quality
      !   of fit but with fewer knots and therefore better if data reduction
      !   is also an important objective for the user.
      !   the number of knots may also depend on the upper bounds ntest and
      !   npest. indeed, if at a certain stage in sphere the number of knots
      !   in one direction (say nt) has reached the value of its upper bound
      !   (ntest), then from that moment on all subsequent knots are added
      !   in the other (phi) direction. this may indicate that the value of
      !   ntest is too small. on the other hand, it gives the user the option
      !   of limiting the number of knots the routine locates in any direction
      !   for example, by setting ntest=8 (the lowest allowable value for
      !   ntest), the user can indicate that he wants an approximation which
      !   is a cubic polynomial in the variable teta.
      !
      !  other subroutines required:
      !    fpback,fpbspl,fpsphe,fpdisc,fpgivs,fprank,fprati,fprota,fporde,
      !    fprpsp
      !
      !  references:
      !   dierckx p. : algorithms for smoothing data on the sphere with tensor
      !                product splines, computing 32 (1984) 319-342.
      !   dierckx p. : algorithms for smoothing data on the sphere with tensor
      !                product splines, report tw62, dept. computer science,
      !                k.u.leuven, 1983.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : july 1983
      !  latest update : march 1989
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND) s,eps,fp
      integer iopt,m,ntest,npest,nt,np,lwrk1,lwrk2,kwrk,ier
      !  ..array arguments..
      real(RKIND) teta(m),phi(m),r(m),w(m),tt(ntest),tp(npest), &
       c((ntest-4)*(npest-4)),wrk1(lwrk1),wrk2(lwrk2)
      integer iwrk(kwrk)
      !  ..local scalars..
      real(RKIND) tol
      integer i,ib1,ib3,ki,kn,kwest,la,lbt,lcc,lcs,lro,j, &
       lbp,lco,lf,lff,lfp,lh,lq,lst,lsp,lwest,maxit,ncest,ncc,ntt, &
       npp,nreg,nrint,ncof,nt4,np4

      !  we set up the parameters tol and maxit.
      maxit = 20
      tol = smallnum03
      !  before starting computations a data check is made. if the input data
      !  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps<=0. .or. eps>=1.) go to 80
      if(iopt<(-1) .or. iopt>1) go to 80
      if(m<2) go to 80
      if(ntest<8 .or. npest<8) go to 80
      nt4 = ntest-4
      np4 = npest-4
      ncest = nt4*np4
      ntt = ntest-7
      npp = npest-7
      ncc = 6+npp*(ntt-1)
      nrint = ntt+npp
      nreg = ntt*npp
      ncof = 6+3*npp
      ib1 = 4*npp
      ib3 = ib1+3
      if(ncof>ib1) ib1 = ncof
      if(ncof>ib3) ib3 = ncof
      lwest = 185+52*npp+10*ntt+14*ntt*npp+8*(m+(ntt-1)*npp**2)
      kwest = m+nreg
      if(lwrk1<lwest .or. kwrk<kwest) go to 80
      if(iopt>0) go to 60

      do 20 i=1,m
        if(w(i)<=0.) go to 80
        if(teta(i)<0. .or. teta(i)>pi) go to 80
        if(phi(i) <zero .or. phi(i)>pi2) go to 80
  20  continue
      if(iopt==0) go to 60
      ntt = nt-8
      if(ntt<0 .or. nt>ntest) go to 80
      if(ntt==0) go to 40
      tt(4) = zero
      do 30 i=1,ntt
         j = i+4
         if(tt(j)<=tt(j-1) .or. tt(j)>=pi) go to 80
  30  continue
  40  npp = np-8
      if(npp<1 .or. np>npest) go to 80
      tp(4) = zero
      do 50 i=1,npp
         j = i+4
         if(tp(j)<=tp(j-1) .or. tp(j)>=pi2) go to 80
  50  continue
      go to 70
  60  if(s<zero) go to 80
  70  ier = 0
      !  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncc*ib3
      lf = la+ncc*ib1
      lff = lf+ncc
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbt = lh+ib3
      lbp = lbt+5*ntest
      lro = lbp+5*npest
      lcc = lro+npest
      lcs = lcc+npest
      lst = lcs+npest
      lsp = lst+m*4
      call fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eps,tol,maxit, &
       ib1,ib3,ncest,ncc,nrint,nreg,nt,tt,np,tp,c,fp,wrk1(1),wrk1(lfp), &
       wrk1(lco),wrk1(lf),wrk1(lff),wrk1(lro),wrk1(lcc),wrk1(lcs), &
       wrk1(la),wrk1(lq),wrk1(lbt),wrk1(lbp),wrk1(lst),wrk1(lsp), &
       wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  80  return
      end subroutine sphere



      !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m the derivative of order nu
      !  of a spline s(x) of degree k,given in its b-spline representation.
      pure subroutine splder(t,n,c,k,nu,x,y,m,e,wrk,ier)

      !
      !  calling sequence:
      !     call splder(t,n,c,k,nu,x,y,m,e,wrk,ier)
      !
      !  input parameters:
      !    t    : array,length n, which contains the position of the knots.
      !    n    : integer, giving the total number of knots of s(x).
      !    c    : array,length n, which contains the b-spline coefficients.
      !    k    : integer, giving the degree of s(x).
      !    nu   : integer, specifying the order of the derivative. 0<=nu<=k
      !    x    : array,length m, which contains the points where the derivative of s(x) must be evaluated.
      !    m    : integer, giving the number of points where the derivative of s(x) must be evaluated
      !    e    : integer, if 0 the spline is extrapolated from the end spans for points not in the
      !           support, if 1 the spline evaluates to zero for those points, and if 2 ier is set to
      !           1 and the subroutine returns.
      !    wrk  : real array of dimension n. used as working space.
      !
      !  output parameters:
      !    y    : array,length m, giving the value of the derivative of s(x) at the different points.
      !    ier  : error flag
      !      ier = 0 : normal return
      !      ier = 1 : argument out of bounds and e == 2
      !      ier =10 : invalid input data (see restrictions)
      !
      !  restrictions:
      !    0 <= nu <= k
      !    m >= 1
      !    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
      !
      !  other subroutines required: fpbspl
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory 6 (1972) 50-62.
      !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths applics 10 (1972) 134-149.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : march 1987
      !
      !++ pearu: 13 aug 2003
      !++   - disabled cliping x values to interval [min(t),max(t)]
      !++   - removed the restriction of the orderness of x values
      !++   - fixed initialization of sp to real(RKIND) value
      !
      !  ..scalar arguments..
      integer, intent(in) :: n,k,nu,m,e
      integer, intent(out) :: ier
      !  ..array arguments..
      real(RKIND), intent(in) :: t(n),c(n),x(m)
      real(RKIND), intent(out) :: y(m)
      real(RKIND), intent(inout) :: wrk(n)
      !  ..local scalars..
      integer :: i,j,k1,k2,k3,l,l1,l2,nk1,nk2,kk
      real(RKIND) :: ak,arg,fac,tb,te
      !  ..local arrays ..
      real(RKIND) :: h(SIZ_K+1)
      logical :: nonflat

      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (nu<0 .or. nu>k) return
      if (m<1) return

      kk  = k-nu

      ier = FITPACK_OK

      !  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      k3 = k1+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)

      !  the derivative of order nu of a spline of degree k is a spline of
      !  degree k-nu,the b-spline coefficients wrk(i) of which can be found
      !  using the recurrence scheme of de boor.
      l  = 1
      wrk(1:nk1) = c(1:nk1)
      if (nu/=0) then
          nk2 = nk1
          de_boor: do j=1,nu
             ak  = k1-j
             nk2 = nk2-1
             l1  = l
             do i=1,nk2
                l1 = l1+1
                l2 = l1+k1-j
                fac = t(l2)-t(l1)
                if (fac>zero) wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
             end do
             l = l+1
          end do de_boor
      endif

      nonflat = nu==0 .or. k/=nu

      l  = k1
      l1 = l+1
      k2 = k1-nu
      j  = 1
      !  main loop for the different points.
      user_points: do i=1,m
        ! fetch a new x-value arg.
        arg = x(i)
        ! check if arg is in the support
        if (arg < tb .or. arg > te) then
            select case (e)
               case (OUTSIDE_EXTRAPOLATE)
                  ! continue like any other point
               case (OUTSIDE_ZERO)
                  y(i) = zero
                  cycle user_points
               case (OUTSIDE_NOT_ALLOWED)
                  ier = FITPACK_INSUFFICIENT_STORAGE
                  return
            end select
        endif
        ! search for knot interval t(l) <= arg < t(l+1)
        do while (.not.(arg>=t(l) .or. l1==k3))
           l1 = l
           l  = l-1
           j  = j-1
        end do
        ! ++
        do while (.not.(arg<t(l1) .or. l==nk1))
           l  = l1
           l1 = l+1
           j  = j+1
        end do

        if (nonflat) then
           !  evaluate the non-zero b-splines of degree k-nu at arg.
           call fpbspl(t,n,kk,arg,l,h)
           !  find the value of the derivative at x=arg.
           y(i) = dot_product(h(1:k2),wrk(l-k:l-nu))
        else
           ! if nu=k the derivative is a piecewise constant function
           y(i) = wrk(j)
        endif
      end do user_points

      return
      end subroutine splder

      ! subroutine splev evaluates in a number of points x(i),i=1,2,...,m a spline s(x) of degree k,
      ! given in its b-spline representation.
      pure subroutine splev(t,n,c,k,x,y,m,e,ier)

      !  calling sequence:
      !     call splev(t,n,c,k,x,y,m,e,ier)
      !
      !  input parameters:
      !    t    : array,length n, which contains the position of the knots.
      !    n    : integer, giving the total number of knots of s(x).
      !    c    : array,length n, which contains the b-spline coefficients.
      !    k    : integer, giving the degree of s(x).
      !    x    : array,length m, which contains the points where s(x) must be evaluated.
      !    m    : integer, giving the number of points where s(x) must be evaluated.
      !    e    : integer, boundary condition for points outside the support
      !           0 = the spline is extrapolated from the end spans
      !           1 = the spline evaluates to zero for those points,
      !           2 = extrapolation not allowed, ier is set to 1 and the subroutine returns,
      !           3 = the spline evaluates to the value of the nearest boundary point.
      !
      !  output parameter:
      !    y    : array,length m, giving the value of s(x) at the different points.
      !    ier  : error flag
      !
      !  restrictions:
      !    m >= 1
      !
      !  other subroutines required: fpbspl.
      !
      !  references :
      !    de boor c  : on calculating with b-splines, j. approximation theory 6 (1972) 50-62.
      !    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths applics 10 (1972) 134-149.
      !    dierckx p. : curve and surface fitting with splines, monographs on numerical analysis, oxford
      !                 university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer,     intent(in)  :: n, k, m, e
      integer,     intent(out) :: ier
      !  ..array arguments..
      real(RKIND), intent(in)  :: t(n), c(n), x(m)
      real(RKIND), intent(out) :: y(m)

      !  ..local scalars..
      integer :: i, k1, l, l1, nk1,k2
      real(RKIND) :: arg, tb, te
      !  ..local array..
      real(RKIND) :: h(SIZ_K+1)
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (m<1) return

      ier = FITPACK_OK
      !  fetch tb and te, the boundaries of the approximation interval.
      k1  = k  + 1
      k2  = k1 + 1
      nk1 = n - k1
      tb  = t(k1)
      te  = t(nk1 + 1)
      l   = k1
      l1  = l + 1

      !  main loop for the different points.
      user_points: do i = 1, m

        ! fetch a new x-value arg.
        arg = x(i)

        ! check if arg is in the support
        outside: if (arg<tb .or. arg>te) then
            select case (e)
               case (OUTSIDE_EXTRAPOLATE)
                ! Continue normally
               case (OUTSIDE_ZERO)
                  y(i) = zero
                  cycle user_points
               case (OUTSIDE_NOT_ALLOWED)
                  ier = FITPACK_INVALID_RANGE
                  return
               case (OUTSIDE_NEAREST_BND)
                  arg = max(min(arg,te),tb)
            end select
        endif outside

        ! search for knot interval t(l) <= arg < t(l+1)
        do while (arg<t(l) .and. l1/=k2)
          l1 = l
          l  = l - 1
        end do
        do while (arg>=t(l1) .and. l/=nk1)
          l  = l1
          l1 = l + 1
        end do

        ! evaluate the non-zero b-splines at arg.
        call fpbspl(t, n, k, arg, l, h)

        ! find the value of s(x) at x=arg.
        y(i) = dot_product(c(l-k:l),h(1:k1))
      end do user_points

      end subroutine splev

      !  function splint calculates the integral of a spline function s(x)
      !  of degree k, which is given in its normalized b-spline representation

      real(RKIND) function splint(t,n,c,k,a,b,wrk) result(integral)

      !  calling sequence:
      !     aint = splint(t,n,c,k,a,b,wrk)
      !
      !  input parameters:
      !    t    : array,length n,which contains the position of the knots of s(x).
      !    n    : integer, giving the total number of knots of s(x).
      !    c    : array,length n, containing the b-spline coefficients.
      !    k    : integer, giving the degree of s(x).
      !    a,b  : real values, containing the end points of the integration interval. s(x) is considered to be identically
      !           zero outside the interval (t(k+1),t(n-k)).
      !
      !  output parameter:
      !    aint : real, containing the integral of s(x) between a and b.
      !    wrk  : real array, length n, used as working space. on output, wrk will contain the integrals of the normalized
      !           b-splines defined on the set of knots.
      !
      !  other subroutines required: fpintb.
      !
      !  references :
      !    gaffney p.w. : the calculation of indefinite integrals of b-splines. j. inst. maths applics 17 (1976) 37-41.
      !    dierckx p. : curve and surface fitting with splines, monographs on numerical analysis, oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      real(RKIND), intent(in) :: a,b
      integer    , intent(in) :: n,k
      !  ..array arguments..
      real(RKIND), intent(in)    :: t(n),c(n)
      real(RKIND), intent(inout) :: wrk(n)

      !  ..local scalars..
      integer :: nk1

      nk1 = n-k-1

      !  calculate the integrals wrk(i) of the normalized b-splines ni,k+1(x), i=1,2,...nk1.
      call fpintb(t,n,wrk,nk1,a,b)

      !  calculate the integral of s(x).
      integral = dot_product(c(1:nk1),wrk(1:nk1))
      return
      end function splint

      !  subroutine sproot finds the zeros of a cubic spline s(x),which is given in its normalized
      ! b-spline representation.
      pure subroutine sproot(t,n,c,zeros,mest,m,ier)

      !
      !  input parameters:
      !    t    : real array,length n, containing the knots of s(x).
      !    n    : integer, containing the number of knots.  n>=8
      !    c    : real array,length n, containing the b-spline coefficients.
      !    mest : integer, specifying the dimension of array zero.
      !
      !  output parameters:
      !    zeros : real array,length mest, containing the zeros of s(x).
      !    m     : integer,giving the number of zeros.
      !    ier   : error flag:
      !      ier = 0: normal return.
      !      ier = 1: the number of zeros exceeds mest.
      !      ier =10: invalid input data (see restrictions).
      !
      !  other subroutines required: fpcuro
      !
      !  restrictions:
      !    1) n>= 8.
      !    2) t(4) < t(5) < ... < t(n-4) < t(n-3).
      !       t(1) <= t(2) <= t(3) <= t(4)
      !       t(n-3) <= t(n-2) <= t(n-1) <= t(n)
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  latest update : october 2022
      !
      ! ..
      ! ..scalar arguments..
      integer, intent(in)  :: n,mest
      integer, intent(out) :: m,ier
      !  ..array arguments..
      real(RKIND), intent(in)  :: t(n),c(n)
      real(RKIND), intent(out) :: zeros(mest)
      !  ..local scalars..
      integer :: i,j,j1,l,n4
      real(RKIND) :: ah,a0,a1,a2,a3,bh,b0,b1,c1,c2,c3,c4,c5,d4,d5,h1,h2,t1,t2,t3,t4,t5
      logical :: z0,z1,z2,z3,z4,nz0,nz1,nz2,nz3,nz4
      !  ..local array..
      real(RKIND) :: y(3)
      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid, control is immediately repassed to the calling program.
      n4  = n-4
      ier = FITPACK_INPUT_ERROR
      if(n<8) return
      j = n
      do i=1,3
        if(t(i)>t(i+1)) return
        if(t(j)<t(j-1)) return
        j = j-1
      end do
      if (any(t(4:n4)>=t(5:n4+1))) return

      !  the problem considered reduces to finding the zeros of the cubic polynomials pl(x) which define
      !  the cubic spline in each knot interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
      !  the condition that it belongs to the knot interval. the cubic polynomial pl(x) is determined by
      !  computing s(t(l)), s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute s(t(l+1))
      !  and s'(t(l+1)); because of the continuity conditions of splines and their derivatives, the value
      !  of s(t(l)) and s'(t(l)) is already known from the foregoing knot interval.
      ier = FITPACK_OK

      !  evaluate some constants for the first knot interval
      h1 = t(4)-t(3)
      h2 = t(5)-t(4)
      t1 = t(4)-t(2)
      t2 = t(5)-t(3)
      t3 = t(6)-t(4)
      t4 = t(5)-t(2)
      t5 = t(6)-t(3)
      !  calculate a0 = s(t(4)) and ah = s'(t(4)).
      c1 = c(1)
      c2 = c(2)
      c3 = c(3)
      c4 = (c2-c1)/t4
      c5 = (c3-c2)/t5
      d4 = (h2*c1+t1*c2)/t4
      d5 = (t3*c2+h1*c3)/t5
      a0 = (h2*d4+h1*d5)/t2
      ah = three*(h2*c4+h1*c5)/t2

      z1  = .not.ah<zero
      nz1 = .not.z1

      m = 0
      !  main loop for the different knot intervals.
      knot_intervals: do l=4,n4

      !  evaluate some constants for the knot interval t(l) <= x <= t(l+1).
        h1 = h2
        h2 = t(l+2)-t(l+1)
        t1 = t2
        t2 = t3
        t3 = t(l+3)-t(l+1)
        t4 = t5
        t5 = t(l+3)-t(l)

      !  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
        c1 = c2
        c2 = c3
        c3 = c(l)
        c4 = c5
        c5 = (c3-c2)/t5
        d4 = (h2*c1+t1*c2)/t4
        d5 = (h1*c3+t3*c2)/t5
        b0 = (h2*d4+h1*d5)/t2
        bh = three*(h2*c4+h1*c5)/t2

      !  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
      !  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
        a1 = ah*h1
        b1 = bh*h1
        a2 = three*(b0-a0)-b1-two*a1
        a3 = two*(a0-b0)+b1+a1

      !  test whether or not pl(x) could have a zero in the range
      !  t(l) <= x <= t(l+1).
        z3  = .not.b1<zero
        nz3 = .not.z3
        if (a0*b0<=zero) go to 100

        z0  = .not.a0<zero
        nz0 = .not.z0
        z2  = .not.a2<zero
        nz2 = .not.z2
        z4  = .not.three*a3+a2<zero
        nz4 = .not.z4

        ! find the zeros of ql(y).
        if (.not.((z0.and.(nz1.and.(z3.or.z2.and.nz4).or.nz2.and.z3.and.z4) &
              .or.nz0.and.(z1.and.(nz3.or.nz2.and.z4).or.z2.and.nz3.and.nz4)))) go to 200

 100    call fpcuro(a3,a2,a1,a0,y,j)

        if (j/=0) then
            ! find which zeros of pl(x) are zeros of s(x).
            which_zeros: do i=1,j
              if(y(i)<zero .or. y(i)>one) cycle which_zeros
              ! test whether the number of zeros of s(x) exceeds mest.
              if (m>=mest) then
                 ier = FITPACK_INSUFFICIENT_STORAGE
                 return
              end if
              m = m+1
              zeros(m) = t(l)+h1*y(i)
            end do which_zeros
        endif

 200    a0 = b0
        ah = bh
        z1 = z3
        nz1 = nz3
      end do knot_intervals

      !  the zeros of s(x) are arranged in increasing order.
      if (m<2) return

      ! FP this double loop can be made more efficient
      sort_zeros: do j=1,m
        inner_loop: do j1 = j+1,m
          if (zeros(j1)<zeros(j)) call swap_RKIND(zeros(j),zeros(j1))
        end do inner_loop
      end do sort_zeros

      ! Filter duplicates
      j = m
      m = 1
      filter_duplicates: do i=2,j
        if (zeros(i)==zeros(m)) cycle filter_duplicates
        m = m+1
        zeros(m) = zeros(i)
      end do filter_duplicates
      return

      end subroutine sproot

      !  subroutine surev evaluates on a grid (u(i),v(j)),i=1,...,mu; j=1,...,mv a bicubic spline
      !  surface of dimension idim, given in the b-spline representation.
      pure subroutine surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,iwrk,kwrk,ier)

      !
      !  calling sequence:
      !     call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,iwrk,kwrk,ier)
      !
      !  input parameters:
      !   idim  : integer, specifying the dimension of the spline surface.
      !   tu    : real array, length nu, which contains the position of the knots in the u-direction.
      !   nu    : integer, giving the total number of knots in the u-direction
      !   tv    : real array, length nv, which contains the position of the knots in the v-direction.
      !   nv    : integer, giving the total number of knots in the v-direction
      !   c     : real array, length (nu-4)*(nv-4)*idim, which contains the b-spline coefficients.
      !   u     : real array of dimension (mu).
      !           before entry u(i) must be set to the u co-ordinate of the i-th grid point along the u-axis.
      !           tu(4)<=u(i-1)<=u(i)<=tu(nu-3), i=2,...,mu.
      !   mu    : on entry mu must specify the number of grid points along the u-axis. mu >=1.
      !   v     : real array of dimension (mv).
      !           before entry v(j) must be set to the v co-ordinate of the j-th grid point along the
      !           v-axis.  tv(4)<=v(j-1)<=v(j)<=tv(nv-3), j=2,...,mv.
      !   mv    : on entry mv must specify the number of grid points along the v-axis. mv >=1.
      !   mf    : on entry, mf must specify the dimension of the array f. mf >= mu*mv*idim
      !   wrk   : real array of dimension lwrk. used as workspace.
      !   lwrk  : integer, specifying the dimension of wrk. lwrk >= 4*(mu+mv)
      !   iwrk  : integer array of dimension kwrk. used as workspace.
      !   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mu+mv.
      !
      !  output parameters:
      !   f     : real array of dimension (mf).
      !           on successful exit f(mu*mv*(l-1)+mv*(i-1)+j) contains the l-th co-ordinate of the bicubic
      !           spline surface at the point (u(i),v(j)),l=1,...,idim,i=1,...,mu;j=1,...,mv.
      !   ier   : integer error flag
      !    ier=0 : normal return
      !    ier=10: invalid input data (see restrictions)
      !
      !  restrictions:
      !   mu >=1, mv >=1, lwrk>=4*(mu+mv), kwrk>=mu+mv , mf>=mu*mv*idim
      !   tu(4) <= u(i-1) <= u(i) <= tu(nu-3), i=2,...,mu
      !   tv(4) <= v(j-1) <= v(j) <= tv(nv-3), j=2,...,mv
      !
      !  other subroutines required:
      !    fpsuev,fpbspl
      !
      !  references :
      !    de boor c : on calculating with b-splines, j. approximation theory 6 (1972) 50-62.
      !    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths applics 10 (1972) 134-149.
      !    dierckx p. : curve and surface fitting with splines, monographs on numerical analysis,
      !                 oxford university press, 1993.
      !
      !  author :
      !    p.dierckx
      !    dept. computer science, k.u.leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  ..scalar arguments..
      integer,     intent(in)    :: idim,nu,nv,mu,mv,mf,lwrk,kwrk
      integer,     intent(out)   :: ier
      !  ..array arguments..
      integer,     intent(inout) :: iwrk(kwrk)
      real(RKIND), intent(in)    :: tu(nu),tv(nv),c((nu-4)*(nv-4)*idim),u(mu),v(mv)
      real(RKIND), intent(inout) :: wrk(lwrk)
      real(RKIND), intent(out)   :: f(mf)
      !  ..local scalars..
      integer :: muv

      !  ..
      !  before starting computations a data check is made. if the input data
      !  are invalid control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      muv = mu+mv

      if (mf<mu*mv*idim .or. &
          lwrk<4*muv    .or. &
          kwrk<muv      .or. &
          mu<1          .or. &
          mv<1)         return

      if (mu>1 .and. any(u(2:mu)<u(1:mu-1))) return
      if (mv>1 .and. any(v(2:mv)<v(1:mv-1))) return

      ! All checks passed
      ier = FITPACK_OK
      call fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,wrk(1),wrk(4*mu+1),iwrk(1),iwrk(mu+1))

      end subroutine surev


      recursive subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nmax,eps,nx,tx,ny,ty,&
                                  c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

      ! given the set of data points (x(i),y(i),z(i)) and the set of positive numbers w(i),i=1,...,m,
      ! subroutine surfit determines a smooth bivariate spline approximation s(x,y) of degrees kx and
      ! ky on the rect angle xb <= x <= xe, yb <= y <= ye.
      ! if iopt = -1 surfit calculates the weighted least-squares spline according to a given set of knots.
      ! if iopt >= 0 the total numbers nx and ny of these knots and their position tx(j),j=1,...,nx and
      ! ty(j),j=1,...,ny are chosen automatically by the routine. the smoothness of s(x,y) is then achieved
      ! by minimizing the discontinuity jumps in the derivatives of s(x,y) across the boundaries of the
      ! subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
      ! The amounth of smoothness is determined by the condition that, for a given smoothing factor s>=0,
      ! f(p) = sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s. the fit is given in the b-spline representation
      ! (b-spline coefficients c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be evaluated
      ! by means of subroutine bispev.
      !
      ! calling sequence:
      !     call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
      !    *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      !
      ! parameters:
      !  iopt  : integer flag. on entry iopt must specify whether a weighted least-squares spline (iopt=-1)
      !          or a smoothing spline (iopt=0 or 1) must be determined.
      !          if iopt=0 the routine will start with an initial set of knots
      !          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=1,...,ky+1.
      !          if iopt=1 the routine will continue with the set of knots found at the last call of the
      !          routine. attention: a call with iopt=1 must always be immediately preceded by another call
      !          with iopt=1 or iopt=0. unchanged on exit.
      !  m     : integer. on entry m must specify the number of data points.
      !          m >= (kx+1)*(ky+1). unchanged on exit.
      !  x     : real array of dimension at least (m).
      !  y     : real array of dimension at least (m).
      !  z     : real array of dimension at least (m).
      !          before entry, x(i),y(i),z(i) must be set to the co-ordinates of the i-th data point, for
      !          i=1,...,m. the order of the data points is immaterial. unchanged on exit.
      !  w     : real array of dimension at least (m). before entry, w(i) must be set to the i-th value in
      !          the set of weights. the w(i) must be strictly positive. unchanged on exit.
      !  xb,xe : real values. on entry xb,xe,yb and ye must specify the boundaries of the rectangular
      !  yb,ye   approximation domain. xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
      !  kx,ky : integer values. on entry kx and ky must specify the degrees of the spline. 1<=kx,ky<=5. it
      !          is recommended to use bicubic (kx=ky=3) splines. unchanged on exit.
      !  s     : real. on entry (in case iopt>=0) s must specify the smoothing factor. s>=0. unchanged
      !          on exit. for advice on the choice of s see further comments
      !  nxest : integer. unchanged on exit.
      !  nyest : integer. unchanged on exit.
      !          on entry, nxest and nyest must specify an upper bound for the number of knots required in
      !          the x- and y-directions respect. these numbers will also determine the storage space needed
      !          by the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1). in most practical situation
      !          nxest = kx+1+sqrt(m/2), nyest = ky+1+sqrt(m/2) will be sufficient. see also further comments.
      !  nmax  : integer. on entry nmax must specify the actual dimension of the arrays tx and ty.
      !          nmax >= nxest, nmax >=nyest. unchanged on exit.
      !  eps   : real. on entry, eps must specify a threshold for determining the effective rank of an
      !          over-determined linear system of equations. 0 < eps < 1.  if the number of decimal digits
      !          in the computer representation of a real number is q, then 10**(-q)
      !          is a suitable value for eps in most practical applications. unchanged on exit.
      !  nx    : integer. unless ier=10 (in case iopt >=0), nx will contain the total number of knots with
      !          respect to the x-variable, of the spline approximation returned. if the computation mode
      !          iopt=1 is used, the value of nx should be left unchanged between subsequent calls.
      !          in case iopt=-1, the value of nx should be specified on entry.
      !  tx    : real array of dimension nmax.
      !          on successful exit, this array will contain the knots of the spline with respect to the
      !          x-variable, i.e. the position of the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
      !          position of the additional knots tx(1)=...=tx(kx+1)=xb and tx(nx-kx)=...=tx(nx)=xe needed
      !          for the b-spline representation. if the computation mode iopt=1 is used, the values of tx(1),
      !          ...,tx(nx) should be left unchanged between subsequent calls. if the computation mode
      !          iopt=-1 is used, the values tx(kx+2),...tx(nx-kx-1) must be supplied by the user, before
      !          entry. see also the restrictions (ier=10).
      !  ny    : integer. unless ier=10 (in case iopt >=0), ny will contain the total number of knots with
      !          respect to the y-variable, of the spline approximation returned. if the computation mode
      !          iopt=1 is used, the value of ny should be left unchanged between subsequent calls.
      !          in case iopt=-1, the value of ny should be specified on entry
      !  ty    : real array of dimension nmax. on successful exit, this array will contain the knots of the
      !          spline with respect to the y-variable, i.e. the position of the interior knots ty(ky+2),...,
      !          ty(ny-ky-1) as well as the position of the additional knots ty(1)=...=ty(ky+1)=yb and
      !          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representation. if the computation mode
      !          iopt=1 is used, the values of ty(1),...,ty(ny) should be left unchanged between subsequent
      !          calls. if the computation mode iopt=-1 is used, the values ty(ky+2),...ty(ny-ky-1) must be
      !          supplied by the user, before entry. see also the restrictions (ier=10).
      !  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
      !          on successful exit, c contains the coefficients of the spline approximation s(x,y)
      !  fp    : real. unless ier=10, fp contains the weighted sum of squared residuals of the spline
      !          approximation returned.
      !  wrk1  : real array of dimension (lwrk1). used as workspace.
      !          if the computation mode iopt=1 is used the value of wrk1(1) should be left unchanged between
      !          subsequent calls. on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will contain the
      !          values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1) with d(i) the i-th diagonal element of the
      !          reduced triangular matrix for calculating the b-spline coefficients. it includes those
      !          elements whose square is less than eps,which are treated as 0 in the case of presumed rank
      !          deficiency (ier<-2).
      !  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of the array wrk1 as declared in
      !          the calling (sub)program. lwrk1 must not be too small. let
      !            u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
      !            ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
      !            if(bx<=by) b1 = bx, b2 = b1+v-ky
      !            if(bx >by) b1 = by, b2 = b1+u-kx  then
      !          lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      !  wrk2  : real array of dimension (lwrk2). used as workspace, but only in the case a rank deficient
      !          system is encountered.
      !  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of the array wrk2 as declared in
      !          the calling (sub)program. lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2 where u,v
      !          and b2 are as above. if there are enough data points, scattered uniformly over the
      !          approximation domain and if the smoothing factor s is not too small, there is a good chance
      !          that this extra workspace is not needed. a lot of memory might therefore be saved by setting
      !          lwrk2=1. (see also ier > 10)
      !  iwrk  : integer array of dimension (kwrk). used as workspace.
      !  kwrk  : integer. on entry kwrk must specify the actual dimension of the array iwrk as declared in
      !          the calling (sub)program.  kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
      !  ier   : integer. unless the routine detects an error, ier contains a non-positive value on exit, i.e.
      !   ier=0  : normal return. the spline returned has a residual sum of squares fp such that
      !            abs(fp-s)/s <= tol with tol a relative tolerance set to 0.001 by the program.
      !   ier=-1 : normal return. the spline returned is an interpolating spline (fp=0).
      !   ier=-2 : normal return. the spline returned is the weighted least-squares polynomial of degrees kx
      !            and ky. in this extreme case fp gives the upper bound for the smoothing factor s.
      !   ier<-2 : warning. the coefficients of the spline returned have been computed as the minimal norm
      !            least-squares solution of a (numerically) rank deficient system. (-ier) gives the rank.
      !            especially if the rank deficiency which can be computed as (nx-kx-1)*(ny-ky-1)+ier, is
      !            large the results may be inaccurate. they could also seriously depend on the value of eps.
      !   ier=1  : error. the required storage space exceeds the available storage space, as specified by the
      !            parameters nxest and nyest.
      !            probably causes : nxest or nyest too small. if these parameters are already large, it may
      !            also indicate that s is too small. the approximation returned is the weighted least-squares
      !            spline according to the current set of knots. parameter fp gives the corresponding weighted
      !            sum of squared residuals (fp>s).
      !   ier=2  : error. a theoretically impossible result was found during the iteration process for finding
      !            a smoothing spline with fp = s. probably causes : s too small or badly chosen eps.
      !            there is an approximation returned but the corresponding weighted sum of squared residuals
      !            does not satisfy the condition abs(fp-s)/s < tol.
      !   ier=3  : error. the maximal number of iterations maxit (set to 20 by the program) allowed for
	  !            finding a smoothing spline with fp=s has been reached. probably causes : s too small there
	  !            is an approximation returned but the corresponding weighted sum of squared residuals does
	  !            not satisfy the condition abs(fp-s)/s < tol.
      !   ier=4  : error. no more knots can be added because the number of b-spline coefficients
	  !            (nx-kx-1)*(ny-ky-1) already exceeds the number of data points m. likely causes: either s
	  !            or m too small. the approximation returned is the weighted least-squares spline according
	  !            to the current set of knots. the parameter fp gives the corresponding weighted sum of
      !            squared residuals (fp>s).
      !   ier=5  : error. no more knots can be added because the additional knot would (quasi) coincide with
	  !            an old one. likely causes : s too small or too large a weight to an inaccurate data point.
      !            the approximation returned is the weighted least-squares spline according to the current
	  !            set of knots. the parameter fp gives the corresponding weighted sum of squared residuals (fp>s).
      !   ier=10 : error. on entry, the input data are controlled on validity the following restrictions must
	  !            be satisfied.
      !            -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
      !            nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
      !            xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
      !            lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      !            kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
      !            if iopt=-1: 2*kx+2<=nx<=nxest
      !                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
      !                        2*ky+2<=ny<=nyest
      !                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
      !            if iopt>=0: s>=0
      !            if one of these conditions is found to be violated,control is returned to the calling program.
	  !            in that case there is no approximation returned.
      !   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-space for computing the minimal
	  !            least-squares solution of a rank deficient system of linear equations. ier gives the
	  !            requested value for lwrk2. there is no approximation returned but, having saved the
	  !            information contained in nx,ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and
      !            the dimension of the array wrk2 accordingly, the user cancontinue at the point the program
	  !            was left, by calling surfit with iopt=1.
      !
      ! further comments:
      !  by means of the parameter s, the user can control the tradeoff between closeness of fit and smoothness
      !   of fit of the approximation. if s is too large, the spline will be too smooth and signal will be lost;
      !   if s is too small the spline will pick up too much noise. in the extreme cases the program will return
	  !   an interpolating spline ifs=0 and the weighted least-squares polynomial (degrees kx,ky)if s is very
	  !   large. between these extremes, a properly chosen s will result in a good compromise between closeness
	  !   of fit and smoothness of fit. to decide whether an approximation, corresponding to a certain s is
      !   satisfactory the user is highly recommended to inspect the fits graphically.
      !   recommended values for s depend on the weights w(i). if these are taken as 1/d(i) with d(i) an estimate
	  !   of the standard deviation of z(i), a good s-value should be found in the range (m-sqrt(2*m),m+sqrt(2*m)).
      !   if nothing is known about the statistical error in z(i) each w(i) can be set equal to one and s
	  !   determined by trial and error, taking account of the comments above. the best is then to start with a
	  !   very large value of s ( to determine the least-squares polynomial and the corresponding upper bound fp0
	  !   for s) and then to progressively decrease the value of s ( say by a factor 10 in the beginning, i.e.
	  !   s=fp0/10, fp0/100,...and more carefully as the approximation shows more detail) to obtain closer fits.
      !   to choose s very small is strongly discouraged. this considerably increases computation time and memory
	  !   requirements. it may also cause rank-deficiency (ier<-2) and endager numerical stability. to economize
	  !   the search for a good s-value the program provides with different modes of computation. at the first
	  !   call of the routine, or whenever he wants to restart with the initial set of knots the user must set
	  !   iopt=0.
      !   if iopt=1 the program will continue with the set of knots found at the last call of the routine. this
	  !   will save a lot of computation time if surfit is called repeatedly for different values of s. the number
	  !   of knots of the spline returned and their location will depend on the value of s and on the complexity
	  !   of the shape of the function underlying the data. if the computation mode iopt=1 is used, the knots
	  !   returned may also depend on the s-values at previous calls (if these were smaller). therefore, if after
	  !   a number of trials with different s-values and iopt=1, the user can finally accept a fit as satisfactory,
	  !   it may be worthwhile for him to call surfit once more with the selected value for s but now with iopt=0.
	  !   indeed, surfit may then return an approximation of the same quality of fit but with fewer knots and
	  !   therefore better if data reduction is also an important objective for the user. the number of knots may
	  !   also depend on the upper bounds nxest and nyest. indeed, if at a certain stage in surfit the number of
	  !   knots in one direction (say nx) has reached the value of its upper bound (nxest), then from that moment
	  !   on all subsequent knots are added in the other (y) direction. this may indicate that the value of nxest
	  !   is too small. on the other hand, it gives the user the option of limiting the number of knots the
	  !   routine locates in any direction for example, by setting nxest=2*kx+2 (the lowest allowable value for
      !   nxest), the user can indicate that he wants an approximation which is a simple polynomial of degree kx
	  !   in the variable x.
      !
      !  references:
      !   dierckx p. : an algorithm for surface fitting with spline functions
      !                ima j. numer. anal. 1 (1981) 267-283.
      !   dierckx p. : an algorithm for surface fitting with spline functions
      !                report tw50, dept. computer science,k.u.leuven, 1980.
      !   dierckx p. : curve and surface fitting with splines, monographs on
      !                numerical analysis, oxford university press, 1993.
      !
      !  author:
      !    p.dierckx
      !    dept. computer science, k.u. leuven
      !    celestijnenlaan 200a, b-3001 heverlee, belgium.
      !    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
      !
      !  creation date : may 1979
      !
      !  ..
      !  ..scalar arguments..
      real(RKIND), intent(in) :: xb,xe,yb,ye,s,eps,fp
      integer, intent(in) :: iopt
      integer :: m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
      !  ..array arguments..
      real(RKIND) :: x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
      integer :: iwrk(kwrk)
      !  ..local scalars..
      integer :: ib1,ib3,jb1,ki,kmax,km1,km2,kn,kwest,kx1,ky1,la,lbx,lby,lco,lf,lff,lfp,lh,lq,lsx,lsy, &
                 lwest,ncest,nest,nek,nminx,nminy,nmx,nmy,nreg,nrint,nxk,nyk

      !  we set up the parameters tol and maxit.
      logical, parameter :: verbose = .true.
      integer, parameter :: maxit = 20
      real(RKIND), parameter :: tol = smallnum03

      ! Size parameters
      kx1   = kx+1
      ky1   = ky+1
      kmax  = max(kx,ky)
      km1   = kmax+1
      km2   = km1+1
      nminx = 2*kx1
      nminy = 2*ky1

      !  before starting computations a data check is made. if the input data
      !  are invalid,control is immediately repassed to the calling program.
      ier = FITPACK_INPUT_ERROR
      if (.not.(eps>zero .and. eps<one))         goto 1
      if (.not.(kx>0 .and. kx<=5))               goto 1
      if (.not.(ky>0 .and. ky<=5))               goto 1
      if (.not.(iopt>=(-1) .and. iopt<=1))       goto 1
      if (.not.m>=(kx1*ky1))                     goto 1
      if (.not.(nxest>=nminx .and. nxest<=nmax)) goto 1
      if (.not.(nyest>=nminy .and. nyest<=nmax)) goto 1
      nest = max(nxest,nyest)
      nxk = nxest-kx1
      nyk = nyest-ky1
      ncest = nxk*nyk
      nmx = nxest-nminx+1
      nmy = nyest-nminy+1
      nrint = nmx+nmy
      nreg = nmx*nmy
      jb1 = ky*nxk+kx1
      ib1 = min(kx*nyk+ky1,jb1)
      ib3 = kx1*nyk+1
      lwest = ncest*(2+ib1+ib3)+2*(nrint+nest*km2+m*km1)+ib3
      kwest = m+nreg
      if (.not.(lwrk1>=lwest .and. kwrk>=kwest)) goto 1
      if (.not.(xb<xe .and. yb<ye))              goto 1
      if (any(w<=zero))                          goto 1
      if (any(x<xb .or. x>xe))                   goto 1
      if (any(y<xb .or. y>xe))                   goto 1

      if (iopt>=0) then

          if (s<zero)                            goto 1

      else

          ! Check that the pre-existing x, y knot locations are monotonic
          if (nx<nminx .or. nx>nxest)            goto 1
          nxk       = nx-kx1
          tx(kx1)   = xb
          tx(nxk+1) = xe
          if (any(tx(kx1+1:nxk+1)<=tx(kx1:nxk))) goto 2

          if (ny<nminy .or. ny>nyest)            goto 1
          nyk       = ny-ky1
          ty(ky1)   = yb
          ty(nyk+1) = ye
          if (any(ty(ky1+1:nyk+1)<=ty(ky1:nyk))) goto 3

      endif

      ! All checks passed
      ier = FITPACK_OK

      !  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncest*ib3
      lf = la+ncest*ib1
      lff = lf+ncest
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbx = lh+ib3
      nek = nest*km2
      lby = lbx+nek
      lsx = lby+nek
      lsy = lsx+m*km1
      call fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
                  eps,tol,maxit,nest,km1,km2,ib1,ib3,ncest,nrint,nreg,nx,tx, &
                  ny,ty,c,fp,wrk1(1),wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff), &
                  wrk1(la),wrk1(lq),wrk1(lbx),wrk1(lby),wrk1(lsx),wrk1(lsy), &
                  wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
      return

      ! Input parameter error
  1   if (verbose) then
         print "('[fitpack] surfit input parameter error: ')"
         print "('[fitpack] iopt=',i0,' kx=',i0,' ky=',i0,' m=',i0)", iopt,kx,ky,m
         print "('[fitpack] nxest=',i0,' nyest=',i0,' nmax=',i0)", nxest,nyest,nmax
         print "('[fitpack] lwrk1=',i0,' lwrk2=',i0,' kwrk=',i0)", lwrk1,lwrk2,kwrk
         print "('[fitpack] xbounts=[',g0,',',g0,'] ybounds=[',g0,',',g0,']')",xb,xe,yb,ye
         print "('[fitpack] eps=',g0,' s=',g0)", eps,s
      endif
      return
  2   if (verbose) print "('[fitpack] x knot locations are not monotonic: '/,10x,'tx=',*(1x,g0))", tx
      return
  3   if (verbose) print "('[fitpack] y knot locations are not monotonic: '/,10x,'ty=',*(1x,g0))", ty
      return
      end subroutine surfit

      ! Swap two real numbers
      elemental subroutine swap_RKIND(a,b)
         real(RKIND), intent(inout) :: a,b
         real(RKIND) :: tmp

         tmp = a
         a   = b
         b   = tmp

      end subroutine swap_RKIND

      elemental subroutine swap_int(a,b)
         integer, intent(inout) :: a,b
         integer :: tmp

         tmp = a
         a   = b
         b   = tmp

      end subroutine swap_int


end module fitpack_core
