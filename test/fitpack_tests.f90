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
module fitpack_tests
    use fitpack_core
    use fitpack_test_data
    use iso_fortran_env, only: output_unit
    implicit none
    private

    public :: mnbisp ! test bispev: Evaluation of a bivariate spline function
    public :: mncloc ! test clocur: closed curve
    public :: mncoco ! test concon: smoothing with convexity constraints
    public :: mnconc ! test concur: smoothing with endpoint derivative constraints
    public :: mncosp ! test cocosp: least-squares fitting with convexity constraints
    public :: mncual ! test cualde: derivatives of a closed planar spline curve
    public :: mncurf ! test curfit: General curve fitting
    public :: mnfour ! test fourco: Fourier coefficient calculation
    public :: mnist  ! test insert: knot insertion in periodic and non-periodic splines
    public :: mnpade ! test parder: Partial derivatives of a bivariate spline
    public :: mnparc
    public :: mnperc
    public :: mnpogr
    public :: mnpola
    public :: mnprof
    public :: mnregr
    public :: mnspal
    public :: mnspde
    public :: mnspev
    public :: mnsphe
    public :: mnspin
    public :: mnspro
    public :: mnsuev
    public :: mnsurf
    public :: mncuev
    public :: mndbin
    public :: mnevpo
    public :: mnpasu
    public :: mnspgr


    contains

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c              mnbisp : bispev test program                          cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mnbisp(iunit) result(success)
          integer, optional, intent(in) :: iunit
          integer, parameter :: mx = 6, my = 6
          real(RKIND) :: fac,facx
          integer     :: i,ier,j,kx,kx1,ky,ky1,m0,m1,m2,m3,nc,nkx1,nky1,nx,ny
          real(RKIND) :: tx(15),ty(15),c(100),x(mx),y(my),z(mx*my),wrk(100)
          integer     :: iwrk(20),useUnit

          ! Output unit
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          ! Initialize success
          success = .true.

          ! we set up the grid points for evaluating the tensor product splines.
          x = 0.2_RKIND*[0,1,2,3,4,5]
          y = x

          ! loop for different spline degrees with respect to the x-variable
          spline_x_degree: do kx=1,5,2
              ! the knots in the x-direction
              tx(kx+2) = 0.4
              tx(kx+3) = 0.7
              tx(kx+4) = 0.9
              kx1 = kx+1
              nx  = 3+2*kx1
              tx(1:kx1)    = zero
              tx(nx-kx:nx) = one

              ! loop for different spline degrees with respect to the y-variable
              spline_y_degree: do ky=2,3

                  ! the knots in the y-direction
                  ty(ky+2) = 0.3
                  ty(ky+3) = 0.8
                  ky1 = ky+1
                  ny = 2+2*ky1

                  ty(1:ky1)    = zero
                  ty(ny-ky:ny) = one

                  ! we generate the b-spline coefficients for the test function x*y
                  nkx1 = nx-kx1
                  nky1 = ny-ky1

                  c = zero

                  fac = kx*ky
                  m0 = 1
                  do i=2,nkx1
                    m1 = m0+nky1
                    facx = (tx(i+kx)-tx(i))/fac
                    do j=2,nky1
                      m2 = m0+1
                      m3 = m1+1
                      c(m3) = c(m1)+c(m2)-c(m0)+facx*(ty(j+ky)-ty(j))
                      m0 = m0+1
                      m1 = m1+1
                    end do
                  m0 = m0+1
                  end do

                  ! evaluation of the spline
                  call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,100,iwrk,20,ier)

                  ! printing of the results
                  write(useUnit,900) kx,ky
                  write(useUnit,910)
                  write(useUnit,920) (tx(i),i=1,nx)
                  write(useUnit,930)
                  write(useUnit,920) (ty(i),i=1,ny)
                  nc = nkx1*nky1
                  write(useUnit,940)
                  write(useUnit,950) (c(i),i=1,nc)
                  write(useUnit,960)
                  write(useUnit,970) (y(i),i=1,my)
                  write(useUnit,980)
                  m2 = 0
                  do i=1,mx
                      m1 = m2+1
                      m2 = m2+my
                      write(useUnit,990) x(i),(z(j),j=m1,m2)
                  end do

                  ! return on error
                  if (.not.FITPACK_SUCCESS(ier)) then
                      write(useUnit,1000) kx,ky,FITPACK_MESSAGE(ier)
                      success = .false.
                  end if

              end do spline_y_degree
          end do spline_x_degree

          !  format statements.
          900  format(33h0tensor product spline of degrees,2i3)
          910  format(1x,40hposition of the knots in the x-direction)
          920  format(1x,15f5.1)
          930  format(1x,40hposition of the knots in the y-direction)
          940  format(23h b-spline coefficients )
          950  format(1x,8f9.4)
          960  format(1h0,37hspline values at selected grid points)
          970  format(1h0,8x,1hy,4x,6(4x,f4.1))
          980  format(1h ,7x,1hx)
          990  format(6x,f4.1,5x,6f8.2)
         1000  format(1x,'[mnbisp] b-spline evaluation at x-order=',i0,' y-order=',i0,' failed: ',a)

      end function mnbisp


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                mncloc : clocur test program                        cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mncloc(iunit) result(success)
          integer, intent(in), optional :: iunit

          ! m denotes the number of data points
          integer, parameter :: m = 19
          integer, parameter :: nest = 40
          integer, parameter :: lwrk = 1500
          integer, parameter :: nc = 80
          real(RKIND) :: x(2*m),w(m),u(m),t(nest),c(nc),wrk(lwrk),sp(nest)
          integer     :: iwrk(40)
          real(RKIND) :: al,del,fp,s
          integer     :: i,idim,ier,iopt,ipar,is,i1,i2,j,j1,k,l,l1,mx,n,nk1,useUnit

          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          ! Initialize success
          success = .true.

          ! the data absciss values
          x(1:35:2) = [-4.7,-7.048,-6.894,-3.75,-1.042,0.938,2.5,3.524,4.511,5.0,&
                       4.886,3.524,3.2,1.302,-1.424,-3.0,-3.064,-3.665]

          !  the data ordinate values
          x(2:36:2) = [0.0,2.565,5.785,6.495,5.909,5.318,4.33,2.957,1.642,0.0,-1.779, &
                       -2.957,-5.543,-7.386,-8.075,-5.196,-2.571,-1.334]


          ! the first and last data point coincide
          x(2*m-1:2*m) = x(1:2)

          ! we set up the weights and parameter values of the data points
          w = one
          u = [(20*(i-1),i=1,m)]

          ! we set up the dimension information.
          mx   = 38

          ! we will determine a planar closed curve x=sx(u) , y=sy(u)
          idim = 2

          ! for the first approximations we will use cubic splines
          k = 3

          ! we will also supply the parameter values u(i)
          ipar = 1

          !  loop for the different approximating spline curves
          approximations: do is=1,9

             select case (is)

                case (1)

                   ! we start computing the least-squares point ( s very large)
                   iopt = 0
                   s = 900.0_RKIND

                case (2)

                   ! iopt =  1 from the second call on
                   iopt = 1
                   s = 10.0_RKIND

                case (3)

                   ! a smaller value for s to get a closer approximation
                   s = 0.1_RKIND

                case (4)

                   ! a larger value for s to get a smoother approximation
                   s = 0.5_RKIND

                case (5)

                   ! if a satisfactory fit is obtained we can calculate a curve of equal quality
                   ! of fit (same value for s) but possibly with fewer knots by specifying iopt=0
                   iopt = 0
                   s = 0.5_RKIND

                case (6)

                   ! we determine a spline curve with respect to the same smoothing
                   ! factor s, but now we let the program determine parameter values u(i)
                   ipar = 0
                   iopt = 0
                   s = 0.5_RKIND

                case (7)

                   ! we choose a different degree of spline approximation
                   k = 5
                   iopt = 0
                   s = 0.5_RKIND

                case (8)

                   ! we determine an interpolating curve
                   s = zero

                case (9)

                   ! finally we calculate a least-squares spline curve with specified knots
                   iopt = -1
                   n = 9+2*k
                   j = k+2
                   del = (u(m)-u(1))*0.125_RKIND
                   do l=1,7
                      al = l
                      t(j) = u(1)+al*del
                      j = j+1
                   end do

             end select

             ! determine the approximating closed curve
             call clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)

             ! return on error
             if (.not.FITPACK_SUCCESS(ier)) then
                 write(useUnit,1000) is,FITPACK_MESSAGE(ier)
                 success = .false.
                 stop
             end if

             ! printing of the results.
             if (iopt>=0) then
                 write(useUnit,915) k,ipar
                 write(useUnit,920) s
             else
                 write(useUnit,910) k,ipar
             end if

             write(useUnit,925) fp,ier
             write(useUnit,930) n
             write(useUnit,935)
             if (ipar==1) write(useUnit,940) (t(i),i=1,n)
             if (ipar==0) write(useUnit,950) (t(i),i=1,n)
             nk1 = n-k-1
             write(useUnit,945)
             write(useUnit,950) (c(l),l=1,nk1)
             write(useUnit,955)
             i1 = n+1
             i2 = n+nk1
             write(useUnit,950) (c(l),l=i1,i2)
             write(useUnit,960)

             ! we evaluate the spline curve
             call curev(idim,t,n,c,nc,k,u,m,sp,mx,ier)
             do i=1,9
                 l = (i-1)*4+1
                 l1 = l+1
                 j = l+2
                 j1 = j+1
                 write(useUnit,965) x(l),x(l1),sp(l),sp(l1),x(j),x(j1),sp(j),sp(j1)
             end do

             ! return on error
             if (.not.FITPACK_SUCCESS(ier)) then
                 write(useUnit,1000) is,FITPACK_MESSAGE(ier)
                 success = .false.
             end if

          end do approximations

          ! Format statements
          910  format(38h0least-squares closed curve of degree ,i1,7h  ipar=,i1)
          915  format(34h0smoothing closed curve of degree ,i1,7h  ipar=,i1)
          920  format(20h smoothing factor s=,f7.1)
          925  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
          930  format(1x,24htotal number of knots n=,i3)
          935  format(1x,22hposition of the knots )
          940  format(5x,10f6.0)
          945  format(1x,30hb-spline coefficients of sx(u))
          950  format(5x,8f9.4)
          955  format(1x,30hb-spline coefficients of sy(u))
          960  format(1h0,2(4x,2hxi,7x,2hyi,6x,6hsx(ui),3x,6hsy(ui)))
          965  format(1h ,8f9.4)
         1000  format(1x,'[mncloc] closed-curve test ',i0,' failed: ',a)

      end function mncloc



      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mncoco : concon test program                       cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mncoco(iunit) result(success)
          integer, optional, intent(in) :: iunit

          ! m denotes the number of data points.
          integer, parameter :: m = 16
          integer, parameter :: nest = 20
          integer, parameter :: lwrk = 550
          integer, parameter :: kwrk = 450

          real(RKIND) :: x(m),y(m),w(m),v(m),sx(m),s2(m),t(nest),c(nest),wrk(lwrk),s(3)
          integer     :: iwrk(kwrk)
          logical     :: bind(nest)
          integer     :: i,ier,iopt,is,j,maxbin,maxtr,n,useUnit
          real(RKIND) :: sq

          ! Decide output unit
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          success = .true.

          !  the absciss values of the data points.
          x = [0.1,0.3,0.5,0.7,0.9,1.25,1.75,2.25,2.75,3.5,4.5,5.5,6.5,7.5,8.5,9.5]
          !  the ordinate values of the data points.
          y = [0.124,0.234,0.256,0.277,0.278,0.291,0.308,0.311,0.315,0.322,0.317,0.326,&
               0.323,0.321,0.322,0.328]

          !  we set up the weights of the data points.
          w = one
          w(1) = 10.0_RKIND
          w(2) = 3.0_RKIND
          w(16) = 10.0_RKIND

          !  we will determine concave approximations
          v = one

          !  w set up the dimension information.
          maxtr  = 100
          maxbin = 10

          !  we set up the different s-values.
          s(1) = 0.2
          s(2) = 0.04
          s(3) = 0.0002

          !  loop for the different spline approximations
          approximations: do is=1,3

             ! iopt=0: initialization. iopt=1 from the second call on.
             iopt = merge(1,0,is>1)

             ! we determine the concave spline approximation.
             call concon(iopt,m,x,y,w,v,s(is),nest,maxtr,maxbin,n,t,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier)

             if (.not.FITPACK_SUCCESS(ier)) then
                 success = .false.
                 write(useUnit,1000) 'fit',is,s(is),FITPACK_MESSAGE(ier)
             end if

             ! printing of the results.
             write(useUnit,900) s(is),ier
             write(useUnit,905) sq
             write(useUnit,910) n
             write(useUnit,915)
             write(useUnit,920) (t(i),i=1,n)
             write(useUnit,925)

             do j=1,n-6
                if(bind(j)) write(useUnit,930) t(j+3)
             end do
             write(useUnit,935)

             write(useUnit,940) (c(i),i=1,n-4)

             ! we evaluate the second order derivative of the spline.
             call splder(t,n,c,3,2,x,s2,m,0,wrk,ier)

             if (.not.FITPACK_SUCCESS(ier)) then
                 success = .false.
                 write(useUnit,1000) 'derivative',is,s(is),FITPACK_MESSAGE(ier)
             end if

             write(useUnit,945)
             do i=1,m
                write(useUnit,950) i,x(i),y(i),sx(i),s2(i)
             end do

          end do approximations

          !  format statements
          900 format(48h0upper limit for the sum of squared residuals s=,e8.1,5x,15herror flag ier=,i2)
          905 format(1x,28hsum of squared residuals sq=,e10.3)
          910 format(1x,24htotal number of knots n=,i2)
          915 format(1x,21hposition of the knots)
          920 format(5x,8f7.2)
          925 format(1x,24hthe knots where s''(x)=0)
          930 format(5x,f7.2)
          935 format(1x,21hb-spline coefficients)
          940 format(5x,4f12.6)
          945 format(3h0 i,5x,4hx(i),6x,4hy(i),4x,7hs(x(i)),4x,9hs''(x(i)))
          950 format(1x,i2,5x,f4.2,5x,f5.3,5x,f5.3,5x,f8.4)
         1000 format(1x,'[mncoco] convex smoothing ',a,' test ',i0,' failed (s=',f6.4,'): ',a)

      end function mncoco


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                mnconc : concur test program                        cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mnconc(iunit) result(success)
          integer, optional, intent(in) :: iunit

          !  m denotes the number of data points
          integer, parameter :: m = 31

          !  we have a planar curve  x = sx(u) , y = sy(u)
          integer, parameter :: idim = 2

          !  we set up the dimension information.
          integer, parameter :: mx = idim*m
          integer, parameter :: lwrk = 1400
          integer, parameter :: ndd = 12
          integer, parameter :: np = 24
          integer, parameter :: nb = 6
          integer, parameter :: ne = 6
          integer, parameter :: nest = 50
          integer, parameter :: nc = 100

          real(RKIND) :: x(mx),w(m),u(m),t(nest),c(nc),wrk(lwrk),xx(mx),db(nb),de(ne),cp(np),dd(ndd),sp(mx)
          integer :: iwrk(50)
          real(RKIND) :: ai,del,fp,s,sigma
          integer :: i,ib,ie,ier,iopt,is,i1,i2,j,j1,k,kk,k1,l,l1,l2,n,nk1,useUnit

          ! Initialization
          success = .true.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          ! the data abscissae values
          x(1:61:2) = [ real(RKIND) :: -3.109,-2.188,-1.351,-0.605,0.093,0.451,0.652,0.701,0.518,0.277,0.008,&
                       -0.291,-0.562,-0.679,-0.637,-0.425,-0.049,0.575,1.334,2.167,3.206,4.099,4.872,5.710,&
                       6.330,6.741,6.928,6.965,6.842,6.593,6.269]
          !  the data ordinate values
          x(2:62:2) = [ real(RKIND) :: 3.040,2.876,2.634,2.183,1.586,1.010,0.382,-0.218,-0.632,-0.879,-0.981,&
                       -0.886,-0.642,-0.195,0.373,1.070,1.607,2.165,2.618,2.905,2.991,2.897,2.615,2.164,1.617,&
                       0.977,0.383,-0.194,-0.665,-0.901,-1.010]

          !  set up the parameter values for the data points
          del =pi*0.1
          do i=1,m
             ai = i-11
             u(i) = ai*del
          end do

          ! the weights are taken as 1/sigma with sigma an estimate of the
          ! standard deviation of the data points.
          sigma = 0.04
          w = one/sigma

          !  accordingly, the smoothing factor is chosen as s = m
          s = m

          !  begin point derivatives of the curve
          db(1:6) = [-pi,three,three,zero,zero,-two]

          !  end point derivatives of the curve
          de(1:6) = [pi2,-one,-one,zero,zero,two]

          !  for the first approximations we will use cubic splines.
          k = 3

          !  loop for the different spline curves
          curve_tests: do is=1,7

              select case (is)
                  case (1)

                     !  no derivative constraints
                     iopt = 0
                     ib = 0
                     ie = 0

                  case (2)

                     !  fixed end points
                     iopt = 0
                     ib = 1
                     ie = 1

                  case (3)

                     ! first derivative constraint at the end point
                     iopt = 0
                     ib = 2
                     ie = 1

                  case (4)

                     ! first derivative constraints at begin and end point.
                     iopt = 0
                     ib = 2
                     ie = 2

                  case (5)

                     !  we choose quintic splines with second derivative constraints.
                     iopt = 0
                     k = 5
                     ib = 3
                     ie = 3

                  case (6)

                     !  we choose another s-value and continue with the set of knots found at
                     !  the last call of concur.
                     iopt = 1
                     s = 26.

                  case (7)

                     ! finally we also calculate a least-squares curve with specified knots
                     iopt = -1
                     j = k+2
                     set_knots: do l=1,5
                        ai = l-2
                        t(j) = ai*pi*half
                        j = j+1
                     end do set_knots
                     n = 7+2*k

              end select

              !  determination of the spline curve.
              call concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,nest,n,t,nc,c,np,cp,fp,&
                          wrk,lwrk,iwrk,ier)

              ! printing of the results.
              if (iopt<0) then
                  write(useUnit,910) k
              else
                  write(useUnit,915) k
                  write(useUnit,920) s
              endif
              write(useUnit,925) ib,ie
              write(useUnit,930) fp,ier
              write(useUnit,935) n
              write(useUnit,940)
              write(useUnit,945) (t(i),i=1,n)
              nk1 = n-k-1
              write(useUnit,950)
              write(useUnit,945) (c(l),l=1,nk1)
              write(useUnit,955)
              i1 = n+1
              i2 = n+nk1
              write(useUnit,945) (c(l),l=i1,i2)

              if (.not.FITPACK_SUCCESS(ier)) then
                  success = .false.
                  write(useUnit,1000) is,FITPACK_MESSAGE(ier)
              end if

              ! calculate derivatives at the begin point.
              k1 = k+1
              kk = k1/2
              call cualde(idim,t,n,c,nc,k1,u(1),dd,ndd,ier)
              write(useUnit,960)
              do i=1,kk
                 l = i-1
                 l1 = l*idim+1
                 l2 = l1+1
                 write(useUnit,970) l,dd(l1),dd(l2)
              end do

              ! calculate derivatives at the end point.
              call cualde(idim,t,n,c,nc,k1,u(m),dd,ndd,ier)
              write(useUnit,965)
              do i=1,kk
                 l = i-1
                 l1 = l*idim+1
                 l2 = l1+1
                 write(useUnit,970) l,dd(l1),dd(l2)
              end do

              ! we evaluate the spline curve
              call curev(idim,t,n,c,nc,k,u,m,sp,mx,ier)
              write(useUnit,975)
              do i=1,5
                 l = (i-1)*12+3
                 l1 = l+1
                 j = l+6
                 j1 = j+1
                 write(useUnit,980) x(l),x(l1),sp(l),sp(l1),x(j),x(j1),sp(j),sp(j1)
              end do
          end do curve_tests



         910  format(31h0least-squares curve of degree ,i1)
         915  format(27h0smoothing curve of degree ,i1)
         920  format(20h smoothing factor s=,f5.0)
         925  format(37h number of derivative constraints ib=,i2,5x,3hie=,i2)
         930  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
         935  format(1x,24htotal number of knots n=,i3)
         940  format(1x,22hposition of the knots )
         945  format(5x,8f8.4)
         950  format(1x,30hb-spline coefficients of sx(u))
         955  format(1x,30hb-spline coefficients of sy(u))
         960  format(1x,30hderivatives at the begin point)
         965  format(1x,28hderivatives at the end point)
         970  format(5x,6horder=,i2,2f9.4)
         975  format(1h0,2(4x,2hxi,7x,2hyi,6x,6hsx(ui),3x,6hsy(ui)))
         980  format(1h ,8f9.4)
        1000  format(1x,'[mnconc] smoothing with endpoint derivative constraints test ',i0,' failed: ',a)

      end function mnconc

      !  test function for the polar package
      elemental real(RKIND) function testpo(x,y)
          real(RKIND), intent(in) ::x,y
          testpo=(x**2+y**2)/((x+y)**2+half)
      end function testpo

      pure real(RKIND) function r1(v)
          real(RKIND), intent(in) :: v
          r1 = one
      end function r1

      pure real(RKIND) function r2(v)
          real(RKIND), intent(in) :: v
          r2 = (one+cos(v)**2)*half
      end function r2

      ! calculate the value of a test function for the sphere package.
      elemental real(RKIND) function testsp(v,u)
          real(RKIND), intent(in) :: u,v
          real(RKIND) ::cu,cv,rad1,rad2,rad3,su,sv
          cu = cos(u)
          cv = cos(v)
          su = sin(u)
          sv = sin(v)
          rad1 = (cu*sv*0.2)**2+(su*sv)**2+(cv*0.5)**2
          rad2 = (cu*sv)**2+(su*sv*0.5)**2+(cv*0.2)**2
          rad3 = (cu*sv*0.5)**2+(su*sv*0.2)**2+cv**2
          testsp = 1./sqrt(rad1) + 1./sqrt(rad2) + 1./sqrt(rad3)
          return
      end function testsp

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mncosp : cocosp test program                       cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mncosp(iunit) result(success)
          integer, optional, intent(in) :: iunit

          ! m denotes the number of data points.
          integer, parameter :: m = 10

          !  we set up the dimension information.
          integer, parameter :: maxtr = 100
          integer, parameter :: maxbin = 10
          integer, parameter :: lwrk = 550
          integer, parameter :: kwrk = 450

          real(RKIND) :: x(m),y(m),w(m),sx(m),s2(m),t(20),c(20),e(20),wrk(550)
          integer     :: iwrk(450)
          logical     :: bind(20)
          integer     :: i,ier,is,j,n,n4,n6,useUnit
          real(RKIND) :: sq

          !  initialization.
          success = .true.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          !  the absciss values of the data points.
          x = [ real(RKIND) :: 0.25,0.5,0.75,1.25,1.75,2.25,2.75,3.25,6.25,12.25 ]
          !  the ordinate values of the data points.
          y = [ real(RKIND) :: 17.0,15.2,13.8,12.2,11.0,10.1,9.4,8.6,6.1,3.5 ]

          !  we set up the weights of the data points.
          w = one

          !  we fetch the knots of the cubic spline
          n  = 11
          n4 = n-4
          n6 = n-6

          !  the interior knots
          t(5) = 1.6
          t(6) = 2.5
          t(7) = 6.0

          !  the boundary knots
          t(1:4)  = x(1)
          t(8:11) = x(m)

          !  loop for the different spline approximations
          approximations: do is=1,3

             select case (is)
                 case (1)

                      !  a convex spline approximation
                      write(useUnit,900)
                      e(1:n6) = -one

                 case (2)

                      !  a concave spline approximation (a straight line)
                      write(useUnit,905)
                      e(1:n6) = one

                 case (3)

                      !  no convexity/concavity constraints
                      write(useUnit,910)
                      e(1:n6) = zero

             end select

             !  we determine the spline approximation.
             call cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier)

             if (.not.FITPACK_SUCCESS(ier)) then
                 success = .false.
                 write(useUnit,1000) is,FITPACK_MESSAGE(ier)
             end if

             !  printing of the results.
             write(useUnit,915) ier
             write(useUnit,920) sq
             write(useUnit,925) n
             write(useUnit,930)
             write(useUnit,935) (t(i),i=1,n)
             write(useUnit,940)
             do j=1,n6
                if (bind(j)) write(useUnit,945) t(j+3)
             end do
             write(useUnit,950)
             write(useUnit,955) (c(i),i=1,n4)

             !  we evaluate the second order derivative of the spline.
             call splder(t,n,c,3,2,x,s2,m,0,wrk,ier)
             write(useUnit,960)
             do i=1,m
                write(useUnit,965) i,x(i),y(i),sx(i),s2(i)
             end do

          end do approximations

          !  format statements
          900  format(28h0convex spline approximation)
          905  format(29h0concave spline approximation)
          910  format(35h0unconstrained spline approximation)
          915  format(16h error flag ier=,i2)
          920  format(1x,28hsum of squared residuals sq=,e10.3)
          925  format(1x,24htotal number of knots n=,i2)
          930  format(1x,21hposition of the knots)
          935  format(5x,8f7.2)
          940  format(1x,24hthe knots where s''(x)=0)
          945  format(5x,f7.2)
          950  format(1x,21hb-spline coefficients)
          955  format(5x,4f12.4)
          960  format(3h0 i,6x,4hx(i),5x,4hy(i),4x,7hs(x(i)),3x,9hs''(x(i)))
          965  format(1x,i2,5x,f5.2,5x,f4.1,5x,f5.2,5x,f5.2)
         1000  format(1x,'[mncosp] smoothing with convexity constraints test ',i0,' failed: ',a)
      end function mncosp


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mncual : cualde test program                       cc
      !c         evaluation of a closed planar spline curve                 cc
      !c                    x = sx(u) , y = sy(u)                           cc
      !c            through its polynomial representation                   cc
      !c                    in each knot interval.                          cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mncual(iunit) result(success)
          integer, optional, intent(in) :: iunit

          !  we have a planar curve
          integer, parameter :: idim = 2
          !  set up the dimension information
          integer, parameter :: nc = 40
          integer, parameter :: nd = 12
          integer, parameter :: m = 20

          real(RKIND) :: t(m),c(nc),u(m),sp(nc),d(nd),cof(2,6)
          integer     :: i,ier,ii,ip,i1,i2,j,jj,jn,j1,j2,j3,j4,k,kk,k1,l,l1,n,nk,nk1,useUnit
          real(RKIND) :: aj,arg,fac,per,pol,tt,uu

          !  initialization.
          success = .true.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          !  set up the points where the curve will be evaluated.
          u = 0.05_RKIND*[(i-1,i=1,m)]

          !  main loop for the different spline degrees.
          spline_degree: do k=3,5,2

              ! the order of the spline.
              k1 = k+1

              ! n denotes the total number of knots.
              n = 2*k1+4

              ! set up the knots of the spline
              t(k1:k1+5) = [ real(RKIND) :: zero,0.1,0.3,0.4,0.8,one]

              ! fetch the b-spline coefficients for sx(u)
              c(1:5) = [one,three,four,five,-one]

              ! fetch the b-spline coefficients for sy(u)
              c(n+1:n+5) = [one,two,-three,two,four]

              ! incorporate the boundary conditions for periodic splines
              nk  = n-k
              per = t(nk)-t(k1)

              bc: do j=1,k
                 !  the boundary knots
                 i1 = nk+j
                 i2 = nk-j
                 j1 = k1+j
                 j2 = k1-j
                 t(i1) = t(j1)+per
                 t(j2) = t(i2)-per

                 ! the boundary coefficients
                 jn = j+n
                 c(j+5) = c(j)
                 c(jn+5) = c(jn)
              end do bc

              ! print the data for the spline.
              write(useUnit,900) k
              write(useUnit,905)
              write(useUnit,910) (t(i),i=1,n)
              write(useUnit,915)
              nk1 = n-k1
              write(useUnit,920) (c(i),i=1,nk1)
              write(useUnit,925)
              i1 = n+1
              i2 = n+nk1
              write(useUnit,920) (c(i),i=i1,i2)
              l = k
              l1 = k1
              kk = k1*idim

              ! main loop for the different points of evaluation.
              ip = 0
              evaluate_spline: do i=1,m
                  arg = u(i)

                  ! search for knot interval t(l)<=u(i)<t(l+1).
                  search_knot: do while (arg>=t(l1) .and. l/=nk1)
                      !  a new knot interval.
                      l = l1
                      l1 = l+1
                      if (t(l)==t(l1)) cycle search_knot

                      write(useUnit,930) t(l),t(l1)

                      ! calculate the spline derivatives at the midpoint tt of the interval
                      tt = (t(l)+t(l1))*half
                      call cualde(idim,t,n,c,nc,k1,tt,d,nd,ier)

                      if (.not.FITPACK_SUCCESS(ier)) then
                          success = .false.
                          write(useUnit,1000) k,FITPACK_MESSAGE(ier)
                      end if

                      write(useUnit,935)
                      write(useUnit,940) (d(j),j=1,kk)

                      !  calculate the coefficients cof in the polynomial representation of
                      !  the spline curve in the current knot interval,i.e.
                      !    sx(u) = cof(1,1)+cof(1,2)*(u-tt)+...+cof(1,k1)*(u-tt)**k
                      !    sy(u) = cof(2,1)+cof(2,2)*(u-tt)+...+cof(2,k1)*(u-tt)**k
                      fac = one
                      jj = 0
                      do j=1,k1
                        do ii=1,idim
                          jj = jj+1
                          cof(ii,j) = d(jj)/fac
                        end do
                        aj  = j
                        fac = fac*aj
                      end do

                      write(useUnit,945)
                      write(useUnit,950) (cof(1,j),j=1,k1)
                      write(useUnit,955)
                      write(useUnit,950) (cof(2,j),j=1,k1)
                  end do search_knot

                  ! evaluate the polynomial curve
                  uu = arg-tt
                  eval_curve: do ii=1,idim
                      pol = cof(ii,k1)
                      jj = k1
                      do j=1,k
                         jj = jj-1
                         pol = pol*uu+cof(ii,jj)
                      end do
                      ip = ip+1
                      sp(ip) = pol
                  end do eval_curve
              end do evaluate_spline

              write(useUnit,960)
              i2 = 0
              j4 = 0
              do j=1,10
                 i1 = i2+1
                 i2 = i1+1
                 j1 = j4+1
                 j2 = j1+1
                 j3 = j2+1
                 j4 = j3+1
                 write(useUnit,965) i1,u(i1),sp(j1),sp(j2),i2,u(i2),sp(j3),sp(j4)
              end do
          end do spline_degree

          !  format statements.
          900  format(31h0degree of the spline curve k =,i2)
          905  format(1x,21hposition of the knots)
          910  format(5x,12f6.1)
          915  format(1x,30hb-spline coefficients of sx(u))
          920  format(5x,14f5.0)
          925  format(1x,30hb-spline coefficients of sy(u))
          930  format(16h0knot interval (,f4.1,1h,,f4.1,1h))
          935  format(1x,49hcurve derivatives at the midpoint of the interval)
          940  format(1x,3(1x,2e12.4))
          945  format(1x,50hcoefficients in the polynomial represent. of sx(u))
          950  format(2x,6e13.5)
          955  format(1x,50hcoefficients in the polynomial represent. of sy(u))
          960  format(1x,2(5x,1hi,3x,4hu(i),4x,8hsx(u(i)),4x,8hsy(u(i))))
          965  format(1x,2(i6,f7.2,2f12.5))
         1000  format(1x,'[mncual] closed planar curve with order ',i0,' failed: ',a)
      end function mncual


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mncuev : curev test program                        cc
      !c             evaluation of a closed planar curve                    cc
      !c                    x = sx(u) , y = sy(u)                           cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mncuev
      real(RKIND) ::u(20),t(20),c(40),sp(40)
      integer i,idim,i1,i2,ier,j,jn,j1,j2,j3,j4,k,k1,m,mx,n,nc,nk,nk1
      real(RKIND) ::ai,per
      !  we have a planar curve
      idim = 2
      !  set up the dimension information
      nc = 40
      mx = 40
      !  set up the points where the curve will be evaluated.
      m = 20
      do 10 i=1,m
        ai = i-1
        u(i) = ai*0.5e-01
  10  continue
      !  main loop for the different spline degrees.
      do 50 k=1,5
        k1 = k+1
      !  n denotes the total number of knots.
        n = 2*k1+4
      !  set up the knots of the spline
        t(k1) = 0.
        t(k1+1) = 0.1e0
        t(k1+2) = 0.3e0
        t(k1+3) = 0.4e0
        t(k1+4) = 0.8e0
        t(k1+5) = 0.1e+01
      !  fetch the b-spline coefficients for sx(u)
        c(1) = 0.1e+01
        c(2) = 0.3e+01
        c(3) = 0.4e+01
        c(4) = 0.5e+01
        c(5) = -0.1e+01
      !  fetch the b-spline coefficients for sy(u)
        c(n+1) = 0.1e+01
        c(n+2) = 0.2e+01
        c(n+3) = -0.3e+01
        c(n+4) = 0.2e+01
        c(n+5) = 0.4e+01
      !  incorporate the boundary conditions for periodic splines
        nk = n-k
        per = t(nk)-t(k1)
        do 20 j=1,k
      !  the boundary knots
          i1 = nk+j
          i2 = nk-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
      !  the boundary coefficients
          jn = j+n
          c(j+5) = c(j)
          c(jn+5) = c(jn)
  20    continue
      !  evaluate the spline curve
        call curev(idim,t,n,c,nc,k,u,m,sp,mx,ier)
      !  print the results.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        nk1 = n-k1
        write(6,920) (c(i),i=1,nk1)
        write(6,925)
        i1 = n+1
        i2 = n+nk1
        write(6,920) (c(i),i=i1,i2)
        write(6,930)
        i2 = 0
        j4 = 0
        do 40 j=1,10
          i1 = i2+1
          i2 = i1+1
          j1 = j4+1
          j2 = j1+1
          j3 = j2+1
          j4 = j3+1
          write(6,935) i1,u(i1),sp(j1),sp(j2),i2,u(i2),sp(j3),sp(j4)
  40    continue
  50  continue
      stop
      !  format statements.
 900  format(31h0degree of the spline curve k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,12f6.1)
 915  format(1x,30hb-spline coefficients of sx(u))
 920  format(5x,14f5.0)
 925  format(1x,30hb-spline coefficients of sy(u))
 930  format(1x,2(5x,1hi,3x,4hu(i),4x,8hsx(u(i)),4x,8hsy(u(i))))
 935  format(1x,2(i6,f7.2,2f12.5))
      end subroutine mncuev


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                mncurf : curfit test program                        cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mncurf(iunit) result(success)
          integer, optional, intent(in) :: iunit

          !  m denotes the number of data points
          integer, parameter :: m = 25

          ! we set up the dimension information
          integer :: nest = 35
          integer, parameter :: lwrk = 1000

          real(RKIND) :: t(35),c(35),wrk(lwrk),sp(m)
          real(RKIND) :: ai,fp,s,xb,xe
          integer     :: i,ier,iopt,is,j,k,l,l1,l2,n,nk1,iwrk(35),useUnit

          ! we set up the abscissae, ordinate values, and weights of the data points
          real(RKIND), parameter :: w(m) = one
          real(RKIND), parameter :: x(m) = [(real(i-1,RKIND),i=1,m)]
          real(RKIND), parameter :: y(m) = [ real(RKIND) :: 1.0,1.0,1.4,1.1,1.0,1.0,4.0,9.0,13.0,13.4, &
                                             12.8,13.1,13.0,14.0,13.0,13.5,10.0,2.0,3.0,2.5,2.5,2.5, &
                                             3.0,4.0,3.5]

          ! Initialization.
          success = .true.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          !  we set up the boundaries of the approximation interval
          xb = x(1)
          xe = x(m)

          !  loop for the different spline degrees.
          OUTSIDE_degrees: do k=3,5,2
          !  loop for the different spline approximations of degree k
             test_case: do is=1,7

                select case (is)
                  case (1)
                     !  we start computing the least-squares polynomial (large value for s).
                     iopt = 0
                     s = 1000.0_RKIND
                  case (2)
                     !  iopt=1 from the second call on
                     iopt = 1
                     s = 60.0_RKIND
                  case (3)
                     !  a smaller value for s to get a closer approximation
                     s = 10.0_RKIND
                  case (4)
                     !  a larger value for s to get a smoother approximation
                     s = 30.0_RKIND
                  case (5)
                     !  if a satisfactory fit is obtained  we can calculate a spline of equal quality
                     !  of fit ( same value for s ) but possibly with fewer knots by specifying iopt=0
                     s = 30.0_RKIND
                     iopt = 0
                  case (6)
                     !  we calculate an interpolating spline
                     s = zero

                  case (7)

                     ! finally, we also calculate a least-squares spline function with specified knots
                     iopt = -1
                     j = k+2
                     do l=1,7
                        ai =3*l
                        t(j) = ai
                        j = j+1
                     end do
                     n = 9+2*k

                end select

                ! Call fitting routine
                call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
                if (.not.FITPACK_SUCCESS(ier)) then
                    success = .false.
                    write(useUnit,1000) is,FITPACK_MESSAGE(ier)
                endif

                !  printing of the results.
                if (iopt<0) then
                    write(useUnit,910) k
                else
                    write(useUnit,915) k
                    write(useUnit,920) s
                endif

                write(useUnit,925) fp,ier
                write(useUnit,930) n
                write(useUnit,935)
                write(useUnit,940) (t(i),i=1,n)
                nk1 = n-k-1
                write(useUnit,945)
                write(useUnit,950) (c(i),i=1,nk1)
                write(useUnit,955)

                ! evaluation of the spline approximation
                call splev(t,n,c,k,x,sp,m,0,ier)
                if (.not.FITPACK_SUCCESS(ier)) then
                    success = .false.
                    write(useUnit,1000)is,FITPACK_MESSAGE(ier)
                endif

                do i=1,5
                   l1 = (i-1)*5+1
                   l2 = l1+4
                   write(useUnit,960) (x(l),y(l),sp(l),l=l1,l2)
                end do
             end do test_case
          end do OUTSIDE_degrees

          return

          910 format(32h0least-squares spline of degree ,i1)
          915 format(28h0smoothing spline of degree ,i1)
          920 format(20h smoothing factor s=,f5.0)
          925 format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
          930 format(1x,24htotal number of knots n=,i3)
          935 format(1x,22hposition of the knots )
          940 format(5x,12f6.1)
          945 format(23h0b-spline coefficients )
          950 format(5x,8f9.4)
          955 format(1h0,5(1x,2hxi,3x,2hyi,2x,5hs(xi),1x))
          960 format(1h ,5(f4.1,1x,f4.1,1x,f4.1,2x))
         1000 format(1x,'[mncurf] curve fit test ',i0,' failed: ',a)

      end function mncurf


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c              mndbin : dblint test program                          cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mndbin
      real(RKIND) ::fac,facx,aint,exint,xb,xe,yb,ye
      integer i,j,kx,kx1,ky,ky1,mx,my,m0,m1,m2,m3,nc,nkx1,nky1,nx,ny
      real(RKIND) ::tx(15),ty(15),c(100),x(6),y(6),wrk(50)
      !  we set up the end points of the integration domains.
      mx = 6
      my = 6
      do 10 i=1,6
      x(i) = (i-1)*0.2
      y(i) = x(i)
  10  continue
      !  loop for different spline degrees with respect to the x-variable
      do 300 kx=1,5,2
      !  the knots in the x-direction
        tx(kx+2) = 0.4
        tx(kx+3) = 0.7
        tx(kx+4) = 0.9
        kx1 = kx+1
        nx = 3+2*kx1
        j = nx
        do 20 i=1,kx1
          tx(i) = 0.
          tx(j) = 1.
          j = j-1
  20    continue
      !  loop for different spline degrees with respect to the y-variable
      do 200 ky=2,3
      !  the knots in the y-direction
        ty(ky+2) = 0.3
        ty(ky+3) = 0.8
        ky1 = ky+1
        ny = 2+2*ky1
        j = ny
        do 30 i=1,ky1
          ty(i) = 0.
          ty(j) = 1.
          j = j-1
  30    continue
      !  we generate the b-spline coefficients for the test function x*y
        nkx1 = nx-kx1
        nky1 = ny-ky1
        do 40 i=1,nky1
          c(i) = 0.
  40    continue
        do 50 i=2,nkx1
          c((i-1)*nky1+1) = 0.
  50    continue
        fac = kx*ky
        m0 = 1
        do 70 i=2,nkx1
          m1 = m0+nky1
          facx = (tx(i+kx)-tx(i))/fac
          do 60 j=2,nky1
            m2 = m0+1
            m3 = m1+1
            c(m3) = c(m1)+c(m2)-c(m0)+facx*(ty(j+ky)-ty(j))
            m0 = m0+1
            m1 = m1+1
  60      continue
          m0 = m0+1
  70    continue
      !  printing of the spline information
        write(6,900) kx,ky
        write(6,910)
        write(6,920) (tx(i),i=1,nx)
        write(6,930)
        write(6,920) (ty(i),i=1,ny)
        nc = nkx1*nky1
        write(6,940)
        write(6,950) (c(i),i=1,nc)
      !  integration of the tensor product spline
        write(6,960)
        j = 6
        do 100 i=1,3
          xb = x(i)
          yb = xb
          xe = x(j)
          ye = xe
          j = j-1
          aint = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
          exint = (xe-xb)*(xe+xb)*(ye-yb)*(ye+yb)*0.25
          write(6,970) xb,xe,yb,ye,aint,exint
 100    continue
 200    continue
 300  continue
      stop
      !  format statements.
 900  format(33h0tensor product spline of degrees,2i3)
 910  format(1x,40hposition of the knots in the x-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the y-direction)
 940  format(23h b-spline coefficients )
 950  format(1x,8f9.4)
 960  format(1h0,5x,2hxb,5x,2hxe,5x,2hyb,5x,2hye,6x,6hdblint,7x,5hexint)
 970  format(1x,4f7.1,2f12.5)
      end subroutine mndbin


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c              mnevpo : evapol test program                          cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnevpo
      real(RKIND) ::fac
      integer i,ir,j,m,mx,my,m0,m1,m2,nc,nu4,nv4,nu,nv
      real(RKIND) ::tu(11),tv(10),c(42),x(6),y(6),f(36)
      !  we set up the grid points for evaluating the polar spline.
      mx = 6
      my = 6
      do 10 i=1,6
      x(i) = (2*i-7)*0.1
      y(i) = x(i)
  10  continue
      !  the interior knots with respect to the u-variable.
      tu(5) = 0.4
      tu(6) = 0.7
      tu(7) = 0.9
      nu = 11
      !  the interior knots with respect to the v-variable.
      tv(5) = 0.3
      tv(6) = 0.8
      nv = 10
      !  the boundary knots
      do 20 i=1,4
        tu(i) = 0.
        tv(i) = 0.
        tu(i+7) = 1.
        tv(i+6) = 1.
  20  continue
      !  the number of b-spline coefficients
      nu4 = nu-4
      nv4 = nv-4
      nc = nu4*nv4
      !  we generate the b-spline coefficients for the function s(u,v)=u**2
      m0 = 1
      m1 = m0+nv4
      c(m0) = 0.
      c(m1) = 0.
      do 70 i=3,nu4
        m2 = m1+nv4
        c(m2) = c(m1)+(tu(i+3)-tu(i))*((c(m1)-c(m0))/(tu(i+2)-tu(i-1)) &
          +(tu(i+2)-tu(i))/3.)
        m0 = m1
        m1 = m2
  70  continue
      do 80 i=1,nu4
        m0 = (i-1)*nv4+1
        fac = c(m0)
        do j=2,nv4
          m0 = m0+1
          c(m0) = fac
        end do
  80  continue
      write(6,900)
      write(6,910)
      write(6,920) (tu(i),i=1,nu)
      write(6,930)
      write(6,920) (tv(i),i=1,nv)
      write(6,940)
      write(6,950) (c(j),j=1,nc)
      ! the spline s(u,v) defines a function f(x,y) through the transformation
      !    x = r(v)*u*cos(v)   y = r(v)*u*sin(v)
      ! we consider two different functions r(v)
      do 300 ir=1,2
        go to (110,130),ir
      ! if r(v) =1 and s(u,v) = u**2 then f(x,y) = x**2+y**2
      ! evaluation of f(x,y)
 110    m = 0
        do i=1,my
          do j=1,mx
            m = m+1
            f(m) = evapol(tu,nu,tv,nv,c,r1,x(j),y(i))
          end do
        end do
        write(6,960)
        go to 200
      ! if r(v) = (1+cos(v)**2)/2 and s(u,v) = u**2 then f(x,y) =
      !    4*(x**2+y**2)**3/(4*x**4 +y**4 +4*x**2*y**2)
      ! evaluation of f(x,y)
 130    m = 0
        do i=1,my
          do j=1,mx
            m = m+1
            f(m) = evapol(tu,nu,tv,nv,c,r2,x(j),y(i))
          end do
        end do
        write(6,965)
 200    write(6,970) (x(i),i=1,mx)
        write(6,975)
        m1 = 0
        do 210 j=1,my
          m0 = m1+1
          m1 = m1+mx
          write(6,980) y(j),(f(m),m=m0,m1)
 210    continue
 300  continue
      stop
      !  format statements.
 900  format(24h0polar spline evaluation)
 910  format(1x,40hposition of the knots in the u-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the v-direction)
 940  format(23h b-spline coefficients )
 950  format(5x,8f9.4)
 960  format(1h0,30hf(x,y) corresponding to r(v)=1)
 965  format(1h0,44hf(x,y) corresponding to r(v)=(1+cos(v)**2)/2)
 970  format(1h0,1hx,2x,6f7.1)
 975  format(3x,1hy)
 980  format(1x,f4.1,6f7.3)
      end subroutine mnevpo




      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mnfour : fourco test program                       cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mnfour(iunit) result(success)
          integer, optional, intent(in) :: iunit

          real(RKIND) :: c(20),t(20),wrk1(20),wrk2(20),alfa(10),ress(10),resc(10)
          integer     :: i,ier,j,k,k1,m,n,nk1,useUnit
          real(RKIND) :: ak,rc,rs

          ! Initialization.
          success = .true.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          !  as an example we calculate some integrals of the form
          !          / 1                               / 1
          !         !    x * sin(alfa*x) dx  and      !   x * cos(alfa*x) dx
          !      0 /                               0 /
          !
          !  we will represent y = x as a cubic spline.
          k = 3
          k1 = k+1
          !  we fetch the knots of the cubic spline
          n = 2*k1+4
          !  the boundary knots
          j = n
          t(1:k1)  = zero
          t(n-k:n) = one

          !  the interior knots
          t(5:8) = [real(RKIND) :: 0.1, 0.3, 0.4, 0.8]

          !  find the b-spline representation of y=x
          nk1  = n-k1
          ak   = k
          c(1) = zero
          do i=2,nk1
             j    = i+k
             c(i) = c(i-1)+(t(j)-t(i))/ak
          end do

          !  print the data for the spline.
          write(6,900) k
          write(6,905)
          write(6,910) (t(i),i=1,n)
          write(6,915)
          write(6,920) (c(i),i=1,nk1)

          !  fetch the different values for alfa
          m = 8
          alfa(1) = zero
          alfa(2) = 0.001_RKIND
          do i=3,m
              alfa(i) = -alfa(i-1)*0.1e+02
          end do

          !  calculate the fourier integrals of the cubic spline
          call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)

          if (.not.FITPACK_SUCCESS(ier)) then
             success = .false.
             write(useUnit,1000) FITPACK_MESSAGE(ier)
          end if

          !  print the results
          write(6,925)
          do i=1,m
              ! fetch the exact values of the integrals
              call exfour(alfa(i),rs,rc)
              write(6,930) alfa(i),ress(i),rs,resc(i),rc

              ! Check that the exact values match the numerical ones
              if (.not.abs(ress(i)-rs)<smallnum10+smallnum06*abs(rs)) then
                 write(useUnit,1000)'int{x * sin(alfa*x) dx} values do not match'
                 success = .false.
              end if
              if (.not.abs(resc(i)-rc)<smallnum10+smallnum06*abs(rc)) then
                 write(useUnit,1000)'int{x * cos(alfa*x) dx} values do not match'
                 success = .false.
              end if

          end do

          return

          !  format statements.
          900  format(25h0degree of the spline k =,i2)
          905  format(1x,21hposition of the knots)
          910  format(5x,12f5.1)
          915  format(1x,21hb-spline coefficients)
          920  format(5x,8f9.5)
          925  format(1h0,2x,4halfa,9x,4hress,9x,4hexas,9x,4hresc,9x,4hexac)
          930  format(1x,e8.1,4f13.5)
         1000  format('[mnfour] fourier integral error: ',a)
      end function mnfour

      !  subroutine exfour calculates the exact values of the integrals
      !                 / 1
      !      rs =      !    x*sin(alfa*x) dx    and
      !             0 /
      !                 / 1
      !      rc =      !    x*cos(alfa*x) dx
      !             0 /
      elemental subroutine exfour(alfa,rs,rc)

          real(RKIND), intent(in) :: alfa
          real(RKIND), intent(out) :: rs,rc

          integer :: k,k2
          real(RKIND) :: aa,ak,cc,c1,ss,s1

          if (alfa==zero) then
              rs = zero
              rc = half
          elseif (abs(alfa)>=one) then
              !  integration by parts
              aa = one/alfa
              cc = cos(alfa)
              ss = sin(alfa)
              rs = (ss*aa-cc)*aa
              rc = ((cc-one)*aa+ss)*aa
          else
              !  using the series expansions of sin(alfa*x) and cos(alfa*x)
              rc = half
              rs = alfa/three
              ss = rs
              cc = rc
              aa = -alfa**2
              do k=1,21
                 k2 = 2*(k-1)
                 ak = (k2+2)*(k2+5)
                 ss = ss*aa/ak
                 s1 = rs+ss
                 if (s1==rs) exit
                 rs = s1
              end do
              do k=1,21
                 k2 = 2*(k-1)
                 ak = (k2+1)*(k2+4)
                 cc = cc*aa/ak
                 c1 = rc+cc
                 if (c1==rc) exit
                 rc = c1
              end do
          endif
          return
      end subroutine exfour


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mninst : insert test program                       cc
      !c      application : to find the sum of two periodic splines         cc
      !c                 with different sets of knots                       cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mnist(iunit) result(success)
          integer, optional, intent(in) :: iunit

          integer, parameter :: m = 21
          integer, parameter :: nest = 30

          real(RKIND) :: t1(nest),c1(nest),t2(nest),c2(nest),x(m),y(m),y1(m),y2(m)
          integer :: i,ier,iopt,ip,i1,i2,j,j1,j2,k,k1,nk,n1,n1k1,n2,n2k1,useUnit
          real(RKIND) :: per

          ! Initialization.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          !  set up the points where the splines will be evaluated.
          x = [(0.05_RKIND*(i-1),i=1,m)]

          !  main loop for the different spline degrees.
          spline_degrees: do k=3,5,2
              k1 = k+1

              periodic_or_not: do ip=1,2
                  !  if iopt = 1 the splines will be considered as periodic splines.
                  !  if iopt = 0 they will be considered as ordinary splines.
                  select case (ip)
                     case (1); iopt = 1; write(useUnit,900)
                     case (2); iopt = 0; write(useUnit,905)
                  end select

                  write(useUnit,910) k

                  !  fetch the knots and b-spline coefficients of the first spline s1(x).
                  n1 = 2*k1+5
                  n1k1 = n1-k1

                  t1(k1:k1+6) = [real(RKIND) :: 0.0,0.2,0.3,0.4,0.7,0.9,1.0]
                  c1(1:6)     = [one,two,-one,three,three,-three]

                  !  fetch the knots and b-spline coefficients of the second spline s2(x).
                  n2 = 2*k1+6
                  n2k1 = n2-k1

                  t2(k1:k1+7) = [real(RKIND) :: 0.0,0.1,0.2,0.3,0.4,0.7,0.8,1.0]
                  c2(1:7)     = [two,-two,one,-three,four,four,four]

                  !  incorporate the boundary conditions for periodic splines.
                  per = one
                  nk = n1-k

                  do j=1,k
                      i1 = nk+j
                      i2 = nk-j
                      j1 = k1+j
                      j2 = k1-j
                      ! the boundary knots
                      t1(i1) = t1(j1)+per
                      t1(j2) = t1(i2)-per
                      t2(i1+1) = t2(j1)+per
                      t2(j2) = t2(i2+1)-per
                      ! the boundary coefficients
                      c1(j+6) = c1(j)
                      c2(j+7) = c2(j)
                  end do

                  if (iopt==0) then
                      !  if iopt=0 we insert k knots at the boundaries of the interval to
                      !  find the representation with coincident boundary knots
                      do j=1,k
                         call insert(iopt,t1,n1,c1,k,one,t1,n1,c1,nest,ier)
                         n1 = n1-1

                         if (.not.FITPACK_SUCCESS(ier)) then
                             success = .false.
                             write(useUnit,1000) iopt,j,1,FITPACK_MESSAGE(ier)
                         end if

                         call insert(iopt,t1,n1,c1,k,zero,t1,n1,c1,nest,ier)
                         n1 = n1-1

                         if (.not.FITPACK_SUCCESS(ier)) then
                             success = .false.
                             write(useUnit,1000) iopt,j,1,FITPACK_MESSAGE(ier)
                         end if

                         t1(1:n1) = t1(2:n1+1)
                         c1(1:n1) = c1(2:n1+1)

                      end do

                      do j=1,k
                          call insert(iopt,t2,n2,c2,k,one,t2,n2,c2,nest,ier)
                          n2 = n2-1

                          if (.not.FITPACK_SUCCESS(ier)) then
                              success = .false.
                              write(useUnit,1000) iopt,j,2,FITPACK_MESSAGE(ier)
                          end if

                          call insert(iopt,t2,n2,c2,k,zero,t2,n2,c2,nest,ier)
                          n2 = n2-1

                          if (.not.FITPACK_SUCCESS(ier)) then
                              success = .false.
                              write(useUnit,1000) iopt,j,2,FITPACK_MESSAGE(ier)
                          end if

                          t2(1:n2) = t2(2:n2+1)
                          c2(1:n2) = c2(2:n2+1)
                      end do
                  end if

                  !  print knots and b-spline coefficients of the two splines.
                  write(useUnit,915)
                  write(useUnit,920) (t1(i),i=1,n1)
                  write(useUnit,925)
                  write(useUnit,930) (c1(i),i=1,n1k1)
                  write(useUnit,935)
                  write(useUnit,920) (t2(i),i=1,n2)
                  write(useUnit,940)
                  write(useUnit,930) (c2(i),i=1,n2k1)

                  !  evaluate the two splines
                  call splev(t1,n1,c1,k,x,y1,m,OUTSIDE_EXTRAPOLATE,ier)
                  call splev(t2,n2,c2,k,x,y2,m,OUTSIDE_EXTRAPOLATE,ier)

                  !  insert the knots of the second spline into those of the first one
                  call insert(iopt,t1,n1,c1,k,0.1_RKIND,t1,n1,c1,nest,ier)
                  if (.not.FITPACK_SUCCESS(ier)) then
                      success = .false.
                      write(useUnit,1000) iopt,k+1,1,FITPACK_MESSAGE(ier)
                  end if

                  call insert(iopt,t1,n1,c1,k,0.8_RKIND,t1,n1,c1,nest,ier)
                  if (.not.FITPACK_SUCCESS(ier)) then
                      success = .false.
                      write(useUnit,1000) iopt,k+2,1,FITPACK_MESSAGE(ier)
                  end if

                  !  insert the knots of the first spline into those of the second one
                  call insert(iopt,t2,n2,c2,k,0.9_RKIND,t2,n2,c2,nest,ier)
                  if (.not.FITPACK_SUCCESS(ier)) then
                      success = .false.
                      write(useUnit,1000) iopt,k+1,2,FITPACK_MESSAGE(ier)
                  end if

                  !  print the knots and coefficients of the splines in their new
                  !  representation
                  n1k1 = n1-k1
                  write(useUnit,945)
                  write(useUnit,920) (t1(i),i=1,n1)
                  write(useUnit,925)
                  write(useUnit,930) (c1(i),i=1,n1k1)
                  write(useUnit,940)
                  write(useUnit,930) (c2(i),i=1,n1k1)

                  !  find the coefficients of the sum of the two splines.
                  c1(1:n1k1) = c1(1:n1k1) + c2(1:n1k1)
                  write(useUnit,950)
                  write(useUnit,930) (c1(i),i=1,n1k1)

                  !  evaluate this new spline and compare results
                  call splev(t1,n1,c1,k,x,y,m,0,ier)
                  write(useUnit,955)
                  do i=1,m
                      write(useUnit,960) i,x(i),y1(i),y2(i),y(i)
                  end do
              end do periodic_or_not
          end do spline_degrees

          !  format statements.
          900  format(41h0insertion algorithm for ordinary splines)
          905  format(41h0insertion algorithm for periodic splines)
          910  format(1x,25hdegree of the splines k =,i2)
          915  format(1x,30hposition of the knots of s1(x))
          920  format(5x,15f5.1)
          925  format(1x,30hb-spline coefficients of s1(x))
          930  format(5x,8f9.5)
          935  format(1x,30hposition of the knots of s2(x))
          940  format(1x,30hb-spline coefficients of s2(x))
          945  format(1x,37hposition of the knots after insertion)
          950  format(1x,33hb-spline coefficients of s1+s2(x))
          955  format(3h0 i,6x,1hx,7x,5hs1(x),7x,5hs2(x),6x,8hs1+s2(x))
          960  format(1x,i2,f8.2,3f12.5)
         1000  format('[mnist] with iopt=',i0,' error inserting ',i0,'-th knot on spline ',i0,': ',a)
      end function mnist


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c              mnpade : parder test program                          cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mnpade(iunit) result(success)
          integer, optional, intent(in) :: iunit


          integer, parameter :: mx = 6
          integer, parameter :: my = 6

          real(RKIND) :: fac,facx
          integer :: i,ix,iy,ier,j,kx,kx1,ky,ky1,m0,m1,m2,m3,nc, &
           nkx1,nky1,nux,nuy,nx,ny,useUnit
          real(RKIND) :: tx(15),ty(15),c(100),x(mx),y(my),z(mx*my),wrk(200)
          integer :: iwrk(20)

          ! Initialization.
          success = .true.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          !  we set up the grid points for evaluating the spline derivatives.
          x = [(0.2_RKIND*(i-1),i=1,6)]
          y = x

          !  loop for different spline degrees with respect to the x-variable
          x_spline_degrees: do kx=3,5,2
              !  the knots in the x-direction
              tx(kx+2:kx+4) = [0.4,0.7,0.9]
              kx1 = kx+1
              nx = 3+2*kx1
              tx(1:kx1) = zero
              tx(nx-kx:nx) = one
              ! loop for different spline degrees with respect to the y-variable
              y_spline_degrees: do ky=2,3
                  ! the knots in the y-direction
                  ty(ky+2:ky+3) = [0.3,0.8]
                  ky1 = ky+1
                  ny = 2+2*ky1
                  ty(1:ky1) = zero
                  ty(ny-ky:ny) = one

                  !  we generate the b-spline coefficients for the test function x*y
                  nkx1 = nx-kx1
                  nky1 = ny-ky1
                  c(1:nky1) = zero
                  do i=2,nkx1
                      c((i-1)*nky1+1) = 0.
                  end do
                  fac = kx*ky
                  m0 = 1
                  do i=2,nkx1
                      m1 = m0+nky1
                      facx = (tx(i+kx)-tx(i))/fac
                      do j=2,nky1
                          m2 = m0+1
                          m3 = m1+1
                          c(m3) = c(m1)+c(m2)-c(m0)+facx*(ty(j+ky)-ty(j))
                          m0 = m0+1
                          m1 = m1+1
                      end do
                      m0 = m0+1
                  end do

                  !  printing of the results
                  write(useUnit,900) kx,ky
                  write(useUnit,910)
                  write(useUnit,920) (tx(i),i=1,nx)
                  write(useUnit,930)
                  write(useUnit,920) (ty(i),i=1,ny)
                  nc = nkx1*nky1
                  write(useUnit,940)
                  write(useUnit,950) (c(i),i=1,nc)
                  !  loop for different orders of spline derivatives
                  do ix=1,2
                      nux = ix-1
                      do iy=1,2
                          nuy = iy-1
                          !  evaluation of the spline derivative
                          call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,200,iwrk,20,ier)

                          if (.not.FITPACK_SUCCESS(ier)) then
                              success = .false.
                              write(useUnit,1000) ky,ky,FITPACK_MESSAGE(ier)
                          end if

                          write(useUnit,960) nux,nuy
                          write(useUnit,970) (y(i),i=1,my)
                          write(useUnit,980)
                          m2 = 0
                          do i=1,mx
                             m1 = m2+1
                             m2 = m2+my
                             write(useUnit,990) x(i),(z(j),j=m1,m2)
                          end do
                      end do
                  end do
              end do y_spline_degrees
          end do x_spline_degrees

          !  format statements.
          900  format(33h0tensor product spline of degrees,2i3)
          910  format(1x,40hposition of the knots in the x-direction)
          920  format(1x,15f5.1)
          930  format(1x,40hposition of the knots in the y-direction)
          940  format(23h b-spline coefficients )
          950  format(1x,8f9.4)
          960  format(1h0,26hspline derivative of order,2i4)
          970  format(1h0,8x,1hy,4x,6(4x,f4.1))
          980  format(1h ,7x,1hx)
          990  format(6x,f4.1,5x,6f8.2)
         1000  format('[mnpade] partial derivative for spline w/ orders (',i0,',',i0,' failed: ',a)
      end function mnpade



      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                mnparc : parcur test program                        cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function mnparc(iunit) result(success)
          integer, optional, intent(in) :: iunit

          !  m denotes the number of data points
          integer, parameter :: m = 32

          !  we will determine a planar curve   x=sx(u) , y=sy(u)
          integer, parameter :: idim = 2

          !  we set up the dimension information.
          integer, parameter :: nest = 40
          integer, parameter :: lwrk = 1200
          integer, parameter :: nc = 80
          integer, parameter :: mx = 64

          real(RKIND) :: x(mx),w(m),u(m),t(nest),c(nc),wrk(lwrk),sp(mx)
          integer     :: iwrk(40)
          real(RKIND) :: al,del,fp,s,ub,ue
          integer     :: i,ier,iopt,ipar,is,i1,i2,j,j1,k,l,l1,n,nk1,useUnit

          ! Initialization
          success = .true.
          if (present(iunit)) then
              useUnit = iunit
          else
              useUnit = output_unit
          end if

          !  the data parameter values
          u = [real(RKIND) :: 120.,128.,133.,136.,138.,141.,144.,146.,149.,151.,154.,161.,170.,180.,190.,&
                200.,210.,220.,230.,240.,250.,262.,269.,273.,278.,282.,287.,291.,295.,299.,305.,315.]

          !  the data absciss values
          x(1:63:2) = [-1.5141,-2.0906,-1.9253,-0.8724,-0.3074,-0.5534,0.0192,1.2298,2.5479,2.4710,1.7063,&
                       1.1183,0.5534,0.4727,0.3574,0.1998,0.2882,0.2613,0.2652,0.2805,0.4112,0.9377,1.3527,&
                       1.5564,1.6141,1.6333,1.1567,0.8109,0.2498,-0.2306,-0.7571,-1.1222]
          !  the data ordinate values
          x(2:64:2) = [0.5150,1.3412,2.6094,3.2358,2.7401,2.7823,3.5932,3.8353,2.5863,1.3105,0.6841,0.2575,&
                       0.2460,0.3689,0.2460,0.2998,0.3651,0.3343,0.3881,0.4573,0.5918,0.7110,0.4035,0.0769,&
                       -0.3920,-0.8570,-1.3412,-1.5641,-1.7409,-1.7178,-1.2989,-0.5572]

          !  we set up the weights of the data points
          w = one

          !  for the first approximations we will use cubic splines
          k = 3
          !  we will also supply the parameter values u(i)
          ipar = 1
          ub = 120.0_RKIND
          ue = 320.0_RKIND

          !  loop for the different approximating spline curves
          approximations: do is=1,9

              select case (is)
                 case (1)

                    !  we start computing a polynomial curve ( s very large)
                    iopt = 0
                    s = 100.0_RKIND

                 case (2)

                    !  iopt =  1 from the second call on
                    iopt = 1
                    s = 1.0_RKIND

                 case (3)

                    !  a smaller value for s to get a closer approximation
                    s = 0.05_RKIND

                 case (4)

                    !  a larger value for s to get a smoother approximation
                    s = 0.25_RKIND

                 case (5)

                    !  if a satisfactory fit is obtained we can calculate a curve of equal quality of
                    !  fit (same value for s) but possibly with fewer knots by specifying iopt=0
                    iopt = 0
                    s = 0.25_RKIND

                 case (6)
                    !  we determine a spline curve with respect to the same smoothing
                    !  factor s,  but now we let the program determine parameter values u(i)
                    ipar = 0
                    iopt = 0
                    s = 0.25_RKIND

                 case (7)

                    !  we choose a different degree of spline approximation
                    k = 5
                    iopt = 0
                    s = 0.25_RKIND

                 case (8)

                    !  we determine an interpolating curve
                    s = zero

                 case (9)

                    !  finally we calculate a least-squares spline curve with specified knots
                    iopt =-1
                    n = 9+2*k
                    j = k+2
                    del = (ue-ub)*0.125_RKIND
                    do l=1,7
                        al = l
                        t(j) = ub+al*del
                        j = j+1
                    end do

              end select

              !  determine the approximating curve
              call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)

              if (.not.FITPACK_SUCCESS(ier)) then
                  success = .false.
                  write(useUnit,1000) is,FITPACK_MESSAGE(ier)
                  stop 'here'
              end if

              ! printing of the results.
              if (iopt<0) then
                  write(useUnit,910) k,ipar
              else
                  write(useUnit,915) k,ipar
                  write(useUnit,920) s
              endif
              write(useUnit,925) fp,ier
              write(useUnit,930) n
              write(useUnit,935)

              if (ipar==1) write(useUnit,940) (t(i),i=1,n)
              if (ipar==0) write(useUnit,950) (t(i),i=1,n)
              nk1 = n-k-1
              write(useUnit,945)
              write(useUnit,950) (c(l),l=1,nk1)
              write(useUnit,955)
              i1 = n+1
              i2 = n+nk1
              write(useUnit,950) (c(l),l=i1,i2)
              write(useUnit,960)

              !  we evaluate the spline curve
              call curev(idim,t,n,c,nc,k,u,m,sp,mx,ier)

              if (.not.FITPACK_SUCCESS(ier)) then
                  success = .false.
                  write(useUnit,1000) is,FITPACK_MESSAGE(ier)
              end if

              do i=1,8
                  l = (i-1)*8+3
                  l1 = l+1
                  j = l+4
                  j1 = j+1
                  write(useUnit,965) x(l),x(l1),sp(l),sp(l1),x(j),x(j1),sp(j),sp(j1)
              end do
          end do approximations

         ! Format statements
         910  format(31h0least-squares curve of degree ,i1,7h  ipar=,i1)
         915  format(27h0smoothing curve of degree ,i1,7h  ipar=,i1)
         920  format(20h smoothing factor s=,f7.2)
         925  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
         930  format(1x,24htotal number of knots n=,i3)
         935  format(1x,22hposition of the knots )
         940  format(5x,10f6.0)
         945  format(1x,30hb-spline coefficients of sx(u))
         950  format(5x,8f8.4)
         955  format(1x,30hb-spline coefficients of sy(u))
         960  format(1h0,2(4x,2hxi,7x,2hyi,6x,6hsx(ui),3x,6hsy(ui)))
         965  format(1h ,8f9.4)
        1000  format('[mnparc] parametric curve test ',i0,' failed: ',a)
      end function mnparc


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c               mnpasu : parsur test program                         cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnpasu(datafile)
        real(RKIND), intent(in) :: datafile(:)
        real(RKIND) ::u(21),v(11),f(693),tu(27),tv(17),c(900),wrk(2000),z(693), wk(128)
      integer iwrk(80),iw(32),ipar(2)
      real(RKIND) ::ai,fp,s
      integer kwrk,lwrk,m,mu,mv,j0,j1,j2,j3,nc,nu,nuest,nv,nvest,i,idim,ier,is,iopt,j,l,pos

      ! Store a pointer to the datafile
      pos = 0

      !  we generate the u-coordinates of the grid.
      mu = 21
      do 10 i=1,mu
        u(i) = i-1
  10  continue
      !  we generate the v-coordinates of the grid.
      mv = 11
      do 20 i=1,mv
        v(i) = u(2*i-1)
  20  continue
      !  the dimension of the surface
      idim = 3
      !  we fetch and print the surface co-ordinates at the grid points
      write(6,900)
      write(6,905) (v(i),i=1,mv)
      write(6,910)
      m = mu*mv
      j0 = 0
      do 40 i=1,mu
        write(6,915) u(i)
        j1 = j0
        do 30 l=1,idim
          j2 = j1+1
          j3 = j1+mv

          ! Read from an array
          f(j2:j3) = datafile(pos+1:pos+mv)
          pos = pos+mv

          write(6,925) (f(j),j=j2,j3)
          j1 = j1+m
  30    continue
        j0 = j0+mv
  40  continue
      !  we set up the dimension information
      nuest = 27
      nvest = 17
      lwrk = 2000
      kwrk = 80
      !  main loop for the different spline approximations
      do 300 is=1,4
        go to (110,120,130,140),is
      !  a smoothing surface with no periodicity conditions
 110    iopt = 0
        s = 0.07
        ipar(1) = 0
        ipar(2) = 0
        go to 200
      !  a smoothing surface periodic in the v-variable
 120    ipar(2) = 1
        go to 200
      !  a smoothing surface periodic in both variables
 130    ipar(1) = 1
        go to 200
      !  finally we also calculate a least-squares spline surface
      !  with specified knots.
 140    iopt = -1
        nu = 11
        nv = 11
        j = 5
        do 150 i=1,3
          ai = 5*i
          tu(j) = ai
          tv(j) = tu(j)
          j = j+1
 150    continue
      !  determination of the spline surface.
 200    call parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest, &
         nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      !  printing of the fitting results.
        if(iopt>=0) go to 210
        write(6,935) ipar(1),ipar(2)
        go to 220
 210    write(6,940) ipar(1),ipar(2)
        write(6,945) s
 220    write(6,950) fp,ier
        write(6,955) nu
        write(6,960)
        write(6,965) (tu(i),i=1,nu)
        write(6,970) nv
        write(6,960)
        write(6,965) (tv(i),i=1,nv)
        nc = (nu-4)*(nv-4)
        write(6,975)
        j1 = 0
        do 230 l=1,idim
          j0 = j1+1
          j1 = j1+nc
          write(6,980) (c(j),j=j0,j1)
 230    continue
      !  evaluation of the spline surface.
        call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,z,693, &
         wk,128,iw,32,ier)
        write(6,985)
        write(6,905) (v(i),i=1,mv,2)
        write(6,910)
        j0 = 0
        do 250 i=1,mu,4
          write(6,915) u(i)
          j1 = j0
          do 240 l=1,idim
            j2 = j1+1
            j3 = j1+mv
            write(6,925) (z(j),j=j2,j3,2)
            j1 = j1+m
 240      continue
          j0 = j0+mv*4
 250    continue
 300  continue
      stop
      !  format statements.
 900  format(15h1the input data)
 905  format(1h0,2x,1hv,11(3x,f4.1))
 910  format(1h ,1x,1hu)
 915  format(1h ,f4.1)
 925  format(5x,11f7.3)
 935  format(37h0least-squares surface of periodicity,2i3)
 940  format(33h0smoothing surface of periodicity,2i3)
 945  format(20h smoothing factor s=,f8.2)
 950  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 955  format(1x,42htotal number of knots in the u-direction =,i3)
 960  format(1x,22hposition of the knots )
 965  format(5x,10f6.2)
 970  format(1x,42htotal number of knots in the v-direction =,i3)
 975  format(23h0b-spline coefficients )
 980  format(5x,8f9.4)
 985  format(1h0,37hspline values at selected grid points)
      end subroutine mnpasu


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                mnperc : percur test program                        cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnperc
      real(RKIND) ::x(27),y(27),w(27),t(37),c(37),wrk(1400),sp(27)
      integer iwrk(37)
      real(RKIND) ::al,fp,s
      integer i,ier,iopt,is,j,k,l,lwrk,l1,l2,m,m1,n,nest,nk1
      !  the data absciss values
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11), &
       x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20),x(21), &
       x(22),x(23),x(24),x(25),x(26)/0.0,3.922,7.843,11.765,15.686, &
       19.608,23.509,27.451,31.373,35.294,39.216,43.137,47.059,50.980, &
       54.902,58.824,62.745,66.667,70.588,74.510,78.431,82.353,86.275, &
       90.196,94.118,98.039/
      !  the data ordinate values
      data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11), &
       y(12),y(13),y(14),y(15),y(16),y(17),y(18),y(19),y(20),y(21), &
       y(22),y(23),y(24),y(25),y(26)/10.099,14.835,21.453,25.022,22.427, &
       22.315,22.070,19.673,16.754,13.983,11.973,12.286,16.129,21.560, &
       28.041,39.205,59.489,72.559,75.960,79.137,75.925,68.809,55.758, &
       39.915,22.006,12.076/
      !  m denotes the number of data points
      m = 27
      !  the period of the spline is determined by x(m)
      x(m) = 100.
      y(m) = y(1)
      !  we set up the weights of the data points
      m1 = m-1
      do 10 i=1,m1
         w(i) = 1.0
  10  continue
      !  we set up the dimension information.
      nest = 37
      lwrk = 1400
      !  loop for the different spline degrees.
      do 400 k=3,5,2
      !  loop for the different spline approximations of degree k
         do 300 is=1,7
            go to (110,120,130,140,150,160,170),is
      !  we start computing the least-squares constant (large value for s).
 110        iopt = 0
            s = 65000.
            go to 200
      !  iopt=1 from the second call on
 120        iopt = 1
            s = 500.
            go to 200
      !  a smaller value for s to get a closer approximation
 130        s = 5.
            go to 200
      !  a larger value for s to get a smoother approximation
 140        s = 20.
            go to 200
      !  if a satisfactory fit is obtained  we can calculate a spline of equal
      !  quality of fit ( same value for s ) but possibly with fewer knots by
      !  specifying iopt=0
 150        s = 20.
            iopt = 0
            go to 200
      !  we calculate an interpolating periodic spline.
 160        s = 0.
            go to 200
      !  finally, we also calculate a least-squares periodic spline function
      !  with specified knots.
 170        iopt = -1
            n = 11+2*k
            j = k+2
            do 180 l=1,9
               al = l*10
               t(j) = al
               j = j+1
 180        continue
      !  determine the periodic spline approximation
 200        call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk, &
             iwrk,ier)
      !  printing of the results.
            if(iopt>=0) go to 210
            write(6,910) k
            go to 220
 210        write(6,915) k
            write(6,920) s
 220        write(6,925) fp,ier
            write(6,930) n
            write(6,935)
            write(6,940) (t(i),i=1,n)
            nk1 = n-k-1
            write(6,945)
            write(6,950) (c(i),i=1,nk1)
            write(6,955)
      !  evaluation of the spline approximation
            call splev(t,n,c,k,x,sp,m,0,ier)
            do 230 i=1,9
               l1 = (i-1)*3+1
               l2 = l1+2
               write(6,960) (x(l),y(l),sp(l),l=l1,l2)
 230        continue
 300     continue
 400  continue
      stop
 910  format(41h0least-squares periodic spline of degree ,i1)
 915  format(37h0smoothing periodic spline of degree ,i1)
 920  format(20h smoothing factor s=,f7.0)
 925  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i2)
 930  format(1x,24htotal number of knots n=,i3)
 935  format(1x,22hposition of the knots )
 940  format(5x,8f8.3)
 945  format(23h0b-spline coefficients )
 950  format(5x,8f8.4)
 955  format(1h0,3(3x,2hxi,6x,2hyi,4x,5hs(xi),3x))
 960  format(1h ,3(f7.3,1x,f7.3,1x,f7.3,2x))
      end subroutine mnperc


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                  mnpogr : pogrid test program                      cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnpogr(datafile)
        real(RKIND), intent(in) :: datafile(*)

        !  ..local scalars..
        real(RKIND) ::cv,ermax,er0,exz0,fp,r,sum,sv,x,y,z0,ai,s
        integer i,ier,is,j,k,kwrk,lwrk,m,nc,nuest,nu,nvest,nv,pos

        ! number of u (radius)-values of the grid.
        integer, parameter :: mu = 9

        ! number of v (angle)-values of the grid
        integer, parameter :: mv = 20

        ! ..local arrays..
        integer :: ider(2),iopt(3),iwrk(100),iw(29)
        real(RKIND) :: u(mu),v(mv),z(180),c(300),tu(50),tv(50),f(180),wk(116), &
                       exact(180),err(mu),sp(mu),wrk(1600)

        ! Store a pointer to the data array
        pos = 0

        ! we set up the radius of the disc
        r = one

        ! set up the u-coordinates of the grid.
        forall(i=1:mu) u(i) = 0.1_RKIND*i

        ! we set up the v-coordinates of the grid.
        forall(j=1:mv) v(j) = (j-1)*pi*0.1_RKIND -pi

        ! we fetch the data values at the grid points.
        m = mu*mv
        z(1:m) = datafile(pos+1:pos+m); pos = pos+m

        ! we fetch the data value at the origin.
        z0     = datafile(pos+1)      ; pos = pos+1

      ! we print the data values at the grid points. we also compute and print
      ! the exact value of the test function underlying the data.
      write(6,905)
      write(6,910) (i,i=1,mu)
      write(6,915)
      exz0 = tespog(zero,zero)
      er0 = abs(exz0-z0)
      ermax = er0
      sum = er0
      do 40 j=1,mv
         cv = cos(v(j))
         sv = sin(v(j))
         k = j
         do 30 i=1,mu
            x = u(i)*cv
            y = u(i)*sv
            exact(k) = tespog(x,y)
            err(i) = abs(exact(k)-z(k))
            sum = sum+err(i)
            if(err(i)>ermax) ermax = err(i)
            k = k+mv
  30     continue
         write(6,920) j,(z(k),k=j,m,mv)
         write(6,925) (exact(k),k=j,m,mv)
  40  continue
      ai = m+1
      sum = sum/ai
      write(6,930) z0,exz0
      write(6,935) sum,ermax
      !  we set up the dimension information
      nuest = 16
      nvest = 27
      kwrk = 100
      lwrk = 1600
      ! main loop for the different spline approximations
      do 300 is=1,6
        go to (110,120,130,140,150,160),is
      !  we start computing a set of spline approximations with
      !  only c0-continuity at the origin,
 110    iopt(2) = 0
        ider(2) = 0
      !  non-vanishing at the boundary of the disc,
        iopt(3) = 0
      !  with a data value at the origin.
        ider(1) = 0
      !  initialisation
        iopt(1) = 0
      !  a large value for s for computing the least-squares polynomial
        s = 5.
        go to 200
      !  iopt(1) = 1 from the second call on
 120    s = 0.1
        iopt(1) = 1
        go to 200
      !  an interpolating spline
 130    s = 0.
        go to 200
      !  a second set of approximations with c1-continuity at the origin
 140    iopt(2) = 1
      !  vanishing at the boundary of the disc.
        iopt(3) = 1
      !  exact value at the origin.
        ider(1) = 1
        z0 = exz0
      ! reinitialization
        iopt(1) = 0
        s = 0.1
        go to 200
      !  no data value at the origin
 150    ider(1) = -1
      !  vanishing partial derivatives at the origin
        ider(2) = 1
      ! reinitialization
        iopt(1) = 0
        go to 200
      ! finally we calculate the least-squares spline according to the current
      !  set of knots
 160    iopt(1) = -1
 200    call pogrid(iopt,ider,mu,u,mv,v,z,z0,r,s,nuest,nvest, &
          nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      ! printing of the fitting results.
        if(iopt(1)>=0) go to 210
        write(6,940)
        go to 220
 210    write(6,945) s
 220    write(6,950) iopt(2)
        if(ider(2)==1) write(6,955)
        if(iopt(3)==1) write(6,960)
        write(6,965) fp,ier
        write(6,970) nu
        write(6,975)
        write(6,980) (tu(i),i=1,nu)
        write(6,985) nv
        write(6,975)
        write(6,980) (tv(i),i=1,nv)
        nc = (nu-4)*(nv-4)
        write(6,990)
        write(6,980) (c(i),i=1,nc)
      !  evaluation of the spline approximation
        call bispev(tu,nu,tv,nv,c,3,3,u,mu,v,mv,f,wk,116,iw,29,ier)
        write(6,995)
        write(6,910) (i,i=1,mu,2)
        write(6,915)
        er0 = abs(exz0-c(1))
        ermax = er0
        sum = er0
        do 240 j=1,mv
          k = j
          do 230 i=1,mu
            sp(i) = f(k)
            err(i) = abs(exact(k)-f(k))
            sum = sum+err(i)
            if(err(i)>ermax) ermax = err(i)
            k = k+mv
 230      continue
          if( (j/3)*3 .ne.j ) go to 240
          write(6,920) j,(sp(i),i=1,mu,2)
          write(6,925) (err(i),i=1,mu,2)
 240    continue
        sum = sum/ai
        write(6,1000) c(1),er0
        write(6,935) sum,ermax
 300  continue
      stop
 905  format(49h1data value (exact function value) at grid points)
 910  format(8h u(i),i=,3x,9(i1,7x))
 915  format(8h v(j),j=)
 920  format(1x,i5,9(2x,f6.3))
 925  format(7x,9(2h (,f5.3,1h)))
 930  format(23h0data value at (0,0) = ,f7.3,5x,14hexact value = ,f7.3)
 935  format(19h0mean abs. error = ,f9.3,5x,18hmax. abs. error = ,f9.3)
 940  format(21h0least-squares spline)
 945  format(25h0smoothing spline with s=,f7.2)
 950  format(1x,35horder of continuity at the origin =,i3)
 955  format(1x,43hvanishing partial derivatives at the origin)
 960  format(1x,37hvanishing at the boundary of the disc)
 965  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 970  format(1x,42htotal number of knots in the u-direction =,i3)
 975  format(1x,22hposition of the knots )
 980  format(5x,8f9.4)
 985  format(1x,42htotal number of knots in the v-direction =,i3)
 990  format(23h0b-spline coefficients )
 995  format(50h0spline value (approximation error) at grid points)
1000  format(25h0spline value at (0,0) = ,f7.3,5x,8herror = ,f7.3)
      end subroutine mnpogr
      !
      real(RKIND) function tespog(x,y)
      ! function program tespog calculates the value of the test function
      ! underlying the data.
      !  ..
      !  ..scalar arguments..
      real(RKIND) ::x,y,f
      !  ..
      f = 1.-((3.*x-1.)**2+(3.*y-1.)**2)/(11.-6.*(x+y))
      tespog = f-(1.-x**2-y**2)*(x+y)*54./121.
      return
      end function tespog


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mnpola : polar test program                        cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnpola(datafile)
        integer, parameter :: m1 = 200
        integer, parameter :: m2 = 90

        real(RKIND), intent(in) :: datafile(:)

        real(RKIND), dimension(m1) :: x,y,z,w,u,v,exact,f
        real(RKIND) :: tu(30),tv(30),c(300),s,fp,eps,sum,ermax,error,ai
        real(RKIND), allocatable :: wrk1(:),wrk2(:)
        integer :: iopt(3),iwrk(500)
        integer :: i,is,ier,kwrk,l,lwrk1,lwrk2,l1,l2,m,nc,nu,nv,nuest,nvest,pos

        allocate(wrk1(15000),wrk2(5700))
        pos = 0

      !  we fetch and print the coordinates and function values of the data.
      write(6,900)
      write(6,905)
      l2 = 0
      do 10 i=1,50
         l1 = l2+1
         l2 = l2+4

         do l=l1,l2
            pos = pos+1; x(l) = datafile(pos)
            pos = pos+1; y(l) = datafile(pos)
            pos = pos+1; z(l) = datafile(pos)
         end do

         !read(5,910) (x(l),y(l),z(l),l=l1,l2)
  10  continue
      write(6,915)(x(l),y(l),z(l),l=1,m1)
      !  we calculate the exact function values and set up the weights w(i)=
      !  (0.01)**(-1) (0.01 is an estimate for the standard deviation of the
      !  error in z(i)). at the same time we calculate the mean and maximum
      !  errors for the data values.
      sum = 0.
      ermax = 0.
      do 20 i=1,m1
         w(i) = 0.1e03
         exact(i) = testpo(x(i),y(i))
         error = abs(z(i)-exact(i))
         sum = sum+error
         if(error>ermax) ermax = error
  20  continue
      ai = m1
      sum = sum/ai
      write(6,920) sum,ermax
      !  we set up the dimension information
      nuest = 15
      nvest = 19
      lwrk1 = 15000
      lwrk2 = 5700
      kwrk = 500
      !  we choose a value for eps
      eps = 0.1e-05
      !  main loop for the different spline approximations
      do 400 is=1,5
        go to (110,120,130,140,160),is
      !  we determine a number of smoothing spline approximations on the unit
      !  disk x**2+y**2 <= 1.
      !  all the data points are considered.
 110    m = m1
      !  we set up the smoothing factor.
        s = 1500.
      !  the approximations are not restricted at the boundaries of the disk
        iopt(3) = 0
      !  we request c2-continuity at the origin.
        iopt(2) = 2
      !  at the first call of polar iopt(1) must be zero.
        iopt(1) = 0
        go to 200
      !  iopt(1) = 1 from the second call on
 120    iopt(1) = 1
        s = 200.
        go to 200
 130    s = 170.
        go to 200
      !  we determine a smoothing spline approximation on the ellips
      !  3*x**2+3*y**2-4*x*y<=1.
      !  we only consider the data points inside this domain.
 140    m = m2
        ai = m
      !  the given function has a constant value 0.4 at the boundary of the
      !  ellips. we calculate new data values by substracting this constant
      !  from the old ones.
        do 150 i=1,m
          z(i) = z(i)-0.4
 150    continue
      !  given these data we will then determine approximations which are
      !  identically zero at the boundary of the ellips.
        iopt(3) = 1
      !  we still request c2-continuity at the origin.
        iopt(2) = 2
      !  reinitialization for the knots.
        iopt(1) = 0
      !  we set up the smoothing factor.
        s = 90.
        go to 250
      !  at the last call we will determine the least-squares spline
      !  approximation corresponding to the current set of knots
 160    iopt(1) = -1
        go to 250
      !  determination of the spline approximation on the disk
 200    call polar(iopt,m,x,y,z,w,rad1,s,nuest,nvest,eps,nu,tu, &
         nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
        nc = (nu-4)*(nv-4)
      !  we calculate the function values at the different points.
        do 220 i=1,m
            f(i) = evapol(tu,nu,tv,nv,c,rad1,x(i),y(i))
 220    continue
        write(6,925) s
        go to 300
      !  determination of the spline approximation on the ellips.
 250    call polar(iopt,m,x,y,z,w,rad2,s,nuest,nvest,eps,nu,tu, &
         nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      !  we determine the b-spline coefficients for the spline approximations
      !  of the given function.
        nc = (nu-4)*(nv-4)
        do 260 i=1,nc
            c(i) = c(i)+0.4
 260    continue
      !  we calculate the function values at the different points.
        do 270 i=1,m
            f(i) = evapol(tu,nu,tv,nv,c,rad2,x(i),y(i))
 270    continue
        if(iopt(1)<0) go to 280
        write(6,930) s
        go to 300
 280    write(6,935)
 300    write(6,940) fp,ier
        write(6,945) nu
        write(6,950)
        write(6,955) (tu(i),i=1,nu)
        write(6,960) nv
        write(6,950)
        write(6,955) (tv(i),i=1,nv)
        write(6,965)
        write(6,970) (c(i),i=1,nc)
      !  we determine mean and maximum errors.
        sum = 0.
        ermax = 0.
        do 350 i=1,m
          error = abs(f(i)-exact(i))
          sum = sum+error
          if(error>ermax) ermax = error
 350    continue
        sum = sum/ai
        write(6,975)
        write(6,980)
        write(6,915)(x(l),y(l),f(l),l=2,m,3)
        write(6,920) sum,ermax
 400  continue
      stop
      !  format statements
 900  format(15h1the input data)
 905  format(1h0,3(3x,1hx,6x,1hy,6x,1hz,5x))
 915  format(1h ,3(3f7.3,2x))
 920  format(14h0mean error = ,f7.4,5x,13hmax. error = ,f7.4)
 925  format(38h0smoothing spline on the disk with s =,f5.0)
 930  format(40h0smoothing spline on the ellips with s =,f5.0)
 935  format(35h0least-squares spline on the ellips)
 940  format(27h0sum of squared residuals =,e15.6,5x,12herror flag =,i5)
 945  format(1x,42htotal number of knots in the u-direction =,i3)
 950  format(1x,22hposition of the knots )
 955  format(5x,8f8.4)
 960  format(1x,42htotal number of knots in the v-direction =,i3)
 965  format(23h0b-spline coefficients )
 970  format(5x,8f9.4)
 975  format(33h0spline values at selected points)
 980  format(1h0,3(3x,1hx,6x,1hy,6x,1hf,5x))
      end subroutine mnpola


      !  the boundary of the approximation domain  x**2+y**2<=1. in polar coordinates
      pure real(RKIND) function rad1(v)
         real(RKIND), intent(in) :: v
         rad1 = one
         return
      end function rad1

      ! the boundary of the approximation domain  3*x**2+3*y**2-4*x*y<=1. in polar coordinates
      pure real(RKIND) function rad2(v)
         real(RKIND), intent(in) :: v
         rad2 = one/sqrt(three-two*sin(2*v))
         return
      end function rad2


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c              mnprof : profil test program                          cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnprof
      real(RKIND) ::fac,facx,u
      integer i,ier,iopt,j,kx,kx1,ky,ky1,m,mx,my,m0,m1,m2,m3,nc, &
       nkx1,nky1,nx,ny
      real(RKIND) ::tx(15),ty(15),c(100),x(6),y(6),z(36),cc(15)
      !  we set up the grid points for evaluating the tensor product splines.
      mx = 6
      my = 6
      m = mx*my
      do 10 i=1,6
      x(i) = (i-1)*0.2
      y(i) = x(i)
  10  continue
      !  loop for different spline degrees with respect to the x-variable
      do 300 kx=3,5,2
      !  the knots in the x-direction
        tx(kx+2) = 0.4
        tx(kx+3) = 0.7
        tx(kx+4) = 0.9
        kx1 = kx+1
        nx = 3+2*kx1
        j = nx
        do 20 i=1,kx1
          tx(i) = 0.
          tx(j) = 1.
          j = j-1
  20    continue
      !  loop for different spline degrees with respect to the y-variable
      do 200 ky=2,3
      !  the knots in the y-direction
        ty(ky+2) = 0.3
        ty(ky+3) = 0.8
        ky1 = ky+1
        ny = 2+2*ky1
        j = ny
        do 30 i=1,ky1
          ty(i) = 0.
          ty(j) = 1.
          j = j-1
  30    continue
      !  we generate the b-spline coefficients for the test function x*y
        nkx1 = nx-kx1
        nky1 = ny-ky1
        do 40 i=1,nky1
          c(i) = 0.
  40    continue
        do 50 i=2,nkx1
          c((i-1)*nky1+1) = 0.
  50    continue
        fac = kx*ky
        m0 = 1
        do 70 i=2,nkx1
          m1 = m0+nky1
          facx = (tx(i+kx)-tx(i))/fac
          do 60 j=2,nky1
            m2 = m0+1
            m3 = m1+1
            c(m3) = c(m1)+c(m2)-c(m0)+facx*(ty(j+ky)-ty(j))
            m0 = m0+1
            m1 = m1+1
  60      continue
          m0 = m0+1
  70    continue
      !  printing of the spline information
        write(6,900) kx,ky
        write(6,910)
        write(6,920) (tx(i),i=1,nx)
        write(6,930)
        write(6,920) (ty(i),i=1,ny)
        nc = nkx1*nky1
        write(6,940)
        write(6,950) (c(i),i=1,nc)
      !  we calculate a number of profiles f(y)=s(u,y)
        iopt = 0
        m0 = 1
        do 80 i=1,mx
          u = x(i)
          call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,15,cc,ier)
          write(6,955) u
          write(6,950) (cc(j),j=1,nky1)
      !  evaluation of the one-dimensional spline f(y)
          call splev(ty,ny,cc,ky,y,z(m0),my,OUTSIDE_EXTRAPOLATE,ier)
          m0 = m0+my
  80    continue
        write(6,960)
        write(6,970) (y(i),i=1,my)
        write(6,980)
        m2 = 0
        do 100 i=1,mx
          m1 = m2+1
          m2 = m2+my
          write(6,990) x(i),(z(j),j=m1,m2)
 100    continue
      !  we calculate a number of profiles g(x)=s(x,u)
        iopt = 1
        m0 = 1
        do 120 i=1,my
          u = y(i)
          call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,15,cc,ier)
          write(6,995) u
          write(6,950) (cc(j),j=1,nkx1)
      !  evaluation of the one-dimensional spline g(x)
          call splev(tx,nx,cc,kx,x,z(m0),mx,OUTSIDE_EXTRAPOLATE,ier)
          m0 = m0+mx
 120    continue
        write(6,960)
        write(6,970) (y(i),i=1,my)
        write(6,980)
        do 140 i=1,mx
          write(6,990) x(i),(z(j),j=i,m,mx)
 140    continue
 200    continue
 300  continue
      stop
      !  format statements.
 900  format(33h0tensor product spline of degrees,2i3)
 910  format(1x,40hposition of the knots in the x-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the y-direction)
 940  format(23h b-spline coefficients )
 950  format(1x,8f9.4)
 955  format(45h0b-spline coefficients of the profile f(y)=s(,f3.1, &
       3h,y))
 960  format(1h0,37hspline values at selected grid points)
 970  format(1h0,8x,1hy,4x,6(4x,f4.1))
 980  format(1h ,7x,1hx)
 990  format(6x,f4.1,5x,6f8.2)
 995  format(47h0b-spline coefficients of the profile g(x)=s(x,,f3.1,1h))
      end subroutine mnprof


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c               mnregr : regrid test program                         cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnregr(x,y,z)
        real(RKIND), intent(in) :: x(:),y(:),z(size(x)*size(y))
        real(RKIND) :: tx(17),ty(17),c(300),wrk(850),f(121), wk(132)
        integer iwrk(60),iw(22)
        real(RKIND) ::ai,fp,s,xb,xe,yb,ye
        integer kx,ky,kwrk,lwrk,m,mx,my,m1,m2,nc,nx,nxest,ny,nyest,i,ier,is,iopt,j

        ! fetch the grid size.
        mx = size(x)
        my = size(y)
        m  = mx*my

        !  printing of the input data.
        write(6,915)
        write(6,920) (y(i),i=1,6)
        write(6,925)
        m1 = 1
        do i=1,mx
          m2 = m1+5
          write(6,930) x(i),(z(j),j=m1,m2)
          m1 = m1+my
        end do

        write(6,920) (y(i),i=7,my)
        write(6,925)
      m1 = 7
      do 20 i=1,mx
        m2 = m1+4
        write(6,930) x(i),(z(j),j=m1,m2)
        m1 = m1+my
  20  continue
      !  we set up the boundaries of the approximation domain.
      xb = x(1)
      yb = y(1)
      xe = x(mx)
      ye = y(my)
      !  we set up the dimension information
      nxest = 17
      nyest = 17
      lwrk = 850
      kwrk = 60
      !  main loop for the different spline approximations
      do 300 is=1,6
        go to (110,120,130,140,150,160),is
      !  we start computing the least-squares bicubic polynomial
 110    iopt = 0
        kx = 3
        ky = 3
        s = 10.
        go to 200
      !  iopt=1 from the second call on
 120    iopt = 1
        s = 0.22
        go to 200
      !  overfitting (s too small)
 130    s = 0.1
        go to 200
      !  an interpolating spline
 140    s = 0.
        go to 200
      !  we change the degrees of the spline
 150    kx = 5
        ky = 5
        s = 0.2
        iopt = 0
        go to 200
      !  finally we also calculate a least-squares spline approximation
      !  with specified knots.
 160    iopt = -1
        kx = 3
        ky = 3
        nx = 11
        ny = 11
        j = kx+2
        do 170 i=1,3
          ai = i-2
          tx(j) = ai*0.5
          ty(j) = tx(j)
          j = j+1
 170    continue
      !  determination of the spline approximation.
 200    call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest, &
         nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      !  printing of the fitting results.
        if(iopt>=0) go to 210
        write(6,935) kx,ky
        go to 220
 210    write(6,940) kx,ky
        write(6,945) s
 220    write(6,950) fp,ier
        write(6,955) nx
        write(6,960)
        write(6,965) (tx(i),i=1,nx)
        write(6,970) ny
        write(6,960)
        write(6,965) (ty(i),i=1,ny)
        nc = (nx-kx-1)*(ny-ky-1)
        write(6,975)
        write(6,980) (c(i),i=1,nc)
      !  evaluation of the spline approximation.
        call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,f, &
         wk,132,iw,22,ier)
        write(6,985)
        write(6,920) (y(i),i=1,my,2)
        write(6,925)
        m1 = 1
        do 230 i=1,mx,2
          m2 = m1+my-1
          write(6,930) x(i),(f(j),j=m1,m2,2)
          m1 = m1+2*my
 230    continue
 300  continue
      stop
      !  format statements.
 915  format(15h1the input data)
 920  format(1h0,8x,1hy,4x,6(4x,f4.1))
 925  format(1h ,7x,1hx)
 930  format(6x,f4.1,5x,6f8.4)
 935  format(32h0least-squares spline of degrees,2i3)
 940  format(28h0smoothing spline of degrees,2i3)
 945  format(20h smoothing factor s=,f7.2)
 950  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 955  format(1x,42htotal number of knots in the x-direction =,i3)
 960  format(1x,22hposition of the knots )
 965  format(5x,10f6.2)
 970  format(1x,42htotal number of knots in the y-direction =,i3)
 975  format(23h0b-spline coefficients )
 980  format(5x,8f9.4)
 985  format(1h0,37hspline values at selected grid points)
      end subroutine mnregr



      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mnspal : spalde test program                       cc
      !c    evaluation of a spline function through its polynomial          cc
      !c            representation in each knot interval.                   cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnspal
      real(RKIND) ::x(21),y(21),t(20),c(20),d(6),cof(6)
      integer i,i1,i2,ier,j,jj,k,k1,l,l1,m,n,nk1
      real(RKIND) ::ai,aj,arg,fac,pol,tt,xx
      !  set up the points where the splines will be evaluated.
      m = 21
      do 10 i=1,m
        ai = i-1
        x(i) = ai*0.5e-01
  10  continue
      !  main loop for the different spline degrees.
      do 100 k=3,5,2
        k1 = k+1
      !  n denotes the total number of knots.
        n = 2*k1+4
      !  set up the knots of the spline
        j = n
      !  the boundary knots
        do 20 i=1,k1
          t(i) = 0.
          t(j) = 0.1e+01
          j = j-1
  20    continue
      !  the interior knots
        t(k1+1) = 0.1e+0
        t(k1+2) = 0.3e+0
        t(k1+3) = 0.4e+0
        t(k1+4) = 0.8e+0
      !  generate the b-spline coefficients.
        nk1 = n-k1
        do 30 i=1,nk1
          ai = i
          c(i) = 0.1e-01*ai*(ai-0.5e01)
  30    continue
      !  print the data for the spline.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
        l = k
        l1 = k1
      !  main loop for the different points of evaluation.
        do 80 i=1,m
          arg = x(i)
      !  search for knot interval t(l)<=x(i)<t(l+1).
  40      if(arg<t(l1) .or. l==nk1) go to 60
      !  a new knot interval.
          l = l1
          l1 = l+1
          if(t(l)==t(l1)) go to 40
          write(6,925) t(l),t(l1)
      !  calculate the spline derivatives at the midpoint tt of the interval
          tt = (t(l)+t(l1))*0.5e0
          call spalde(t,n,c,k1,tt,d,ier)
          write(6,930)
          write(6,935) (d(j),j=1,k1)
      !  calculate the coefficients cof in the polynomial representation of
      !  the spline in the current knot interval,i.e.
      !    s(x) = cof(1)+cof(2)*(x-tt)+...+cof(k1)*(x-tt)**k
          fac = 0.1e01
          do 50 j=1,k1
            cof(j) = d(j)/fac
            aj = j
            fac = fac*aj
  50      continue
          write(6,940)
          write(6,935) (cof(j),j=1,k1)
          go to 40
      !  evaluate the polynomial
  60      xx = arg-tt
          pol = cof(k1)
          jj = k1
          do 70 j=1,k
            jj = jj-1
            pol = pol*xx+cof(jj)
  70      continue
          y(i) = pol
  80    continue
        write(6,945)
        i2 = 0
        do 90 j=1,7
          i1 = i2+1
          i2 = i1+2
          write(6,950) (i,x(i),y(i),i=i1,i2)
  90    continue
 100  continue
      stop
      !  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(16h0knot interval (,f4.1,1h,,f4.1,1h))
 930  format(1x,49hderivative values at the midpoint of the interval)
 935  format(2x,6e13.5)
 940  format(1x,45hcoefficients in the polynomial representation)
 945  format(1h0,3(7x,1hi,3x,4hx(i),4x,7hs(x(i))))
 950  format(1x,3(i8,f7.2,f11.5))
      end subroutine mnspal


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mnspde : splder test program                       cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnspde
      real(RKIND) ::x(7),y(42),t(20),c(20),wrk(20),d(6)
      integer i,ier,j,k,k1,l,m,n,nk1,nu
      real(RKIND) ::ai
      !  set up the points where the spline derivatives will be evaluated.
      m = 7
      x(1) = 0.
      ai = 0.5e-01
      do 10 i=2,m
        x(i) = x(i-1)+ai
        ai = ai+0.5e-01
  10  continue
      x(m) = 0.1e+01
      !  main loop for the different spline degrees.
      do 70 k=1,5
        k1 = k+1
      !  n denotes the total number of knots.
        n = 2*k1+4
      !  set up the knots of the spline
        j = n
      !  the boundary knots
        do 20 i=1,k1
          t(i) = 0.
          t(j) = 0.1e+01
          j = j-1
  20    continue
      !  the interior knots
        t(k1+1) = 0.1e+0
        t(k1+2) = 0.3e+0
        t(k1+3) = 0.4e+0
        t(k1+4) = 0.8e+0
      !  generate the b-spline coefficients.
        nk1 = n-k1
        do 30 i=1,nk1
          ai = i
          c(i) = 0.1e-01*ai*(ai-0.5e01)
  30    continue
      !  evaluate the spline derivatives.
        j = 1
        do 40 i=1,k1
      !  nu denotes the order of the derivative
          nu = i-1
          call splder(t,n,c,k,nu,x,y(j),m,OUTSIDE_EXTRAPOLATE,wrk,ier)
          j = j+m
  40    continue
      !  print the results.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
        write(6,925) (i,i=1,5)
        write(6,930)
        do 60 i=1,m
          j = i
          do 50 l=1,k1
            d(l) = y(j)
            j = j+m
  50      continue
          write(6,935) i,x(i),(d(l),l=1,k1)
  60    continue
  70  continue
      stop
      !  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(21x,5(i4,8x))
 930  format(1x,1hi,6h  x(i),3x,8hs(x(i)) ,5(4x,8hs (x(i))))
 935  format(1x,i1,f6.2,6e12.4)
      end subroutine mnspde


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mnspev : splev test program                        cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnspev
      real(RKIND) ::x(21),y(21),t(20),c(20)
      integer i,i1,i2,ier,j,k,k1,m,n,nk1
      real(RKIND) ::ai
      !  set up the points where the splines will be evaluated.
      m = 21
      do 10 i=1,m
        ai = i-1
        x(i) = ai*0.5e-01
  10  continue
      !  main loop for the different spline degrees.
      do 50 k=1,5
        k1 = k+1
      !  n denotes the total number of knots.
        n = 2*k1+4
      !  set up the knots of the spline
        j = n
      !  the boundary knots
        do 20 i=1,k1
          t(i) = 0.
          t(j) = 0.1e+01
          j = j-1
  20    continue
      !  the interior knots
        t(k1+1) = 0.1e+0
        t(k1+2) = 0.3e+0
        t(k1+3) = 0.4e+0
        t(k1+4) = 0.8e+0
      !  generate the b-spline coefficients.
        nk1 = n-k1
        do 30 i=1,nk1
          ai = i
          c(i) = 0.1e-01*ai*(ai-0.5e01)
  30    continue
      !  evaluate the spline.
        call splev(t,n,c,k,x,y,m,OUTSIDE_EXTRAPOLATE,ier)
      !  print the results.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
        write(6,925) ier
        write(6,930)
        i2 = 0
        do 40 j=1,7
          i1 = i2+1
          i2 = i1+2
          write(6,935) (i,x(i),y(i),i=i1,i2)
  40    continue
  50  continue
      stop
      !  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(1x,16herror flag ier =,i3)
 930  format(1x,3(7x,1hi,3x,4hx(i),4x,7hs(x(i))))
 935  format(1x,3(i8,f7.2,f11.5))
      end subroutine mnspev


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                  mnspgr : spgrid test program                      cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnspgr(u,v,r)
        real(RKIND), intent(in) :: u(:),v(:),r(size(u)*size(v))

        !  ..local scalars..
        real(RKIND) :: ermax,erf,exr0,exr1,fp,sum,r0,r1,ai,s
        integer i,ier,is,j,k,kwrk,l,lwrk,m,mu,mv,nc,nuest,nu,nvest,nv
        !  ..local arrays..
        integer ider(4),iopt(3),iwrk(70)
        real(RKIND) :: c(300),tu(25),tv(25),f(154),wk(100),exact(154),err(14),sp(14),wrk(1500)


        ! set up the number of u (latitude)-values of the grid.
        mu = size(u)
        ! set up the number of v (longitude)-values of the grid
        mv = size(v)

        ! we fetch the data values at the grid points.
        m  = mu*mv

        ! we print the data values at the grid points. we also compute and print
        ! the exact value of the test function underlying the data.
        write(6,910)
        write(6,915) (j,j=1,mv,2)
        write(6,920)
        exr0 = tesspg(zero,zero)
        exr1 = tesspg(pi,zero)
        ermax = zero
        sum = zero
      l = 0
      do 40 i=1,mu
        l = (i-1)*mv+1
        k = 1
        do 30 j=1,7
          exact(l) = tesspg(u(i),v(k))
          erf = abs(exact(l)-r(l))
          sum = sum+erf
          if(erf>ermax) ermax = erf
          sp(j) = exact(l)
          err(j) = r(l)
          l = l+2
          k = k+2
  30    continue
        write(6,925) i,(err(j),j=1,7)
        write(6,930) (sp(j),j=1,7)
  40  continue
      write(6,915) (j,j=2,mv,2)
      write(6,920)
      do 60 i=1,mu
        l = (i-1)*mv+2
        k = 2
        do 50 j=1,7
          exact(l) = tesspg(u(i),v(k))
          erf = abs(exact(l)-r(l))
          sum = sum+erf
          if(erf>ermax) ermax = erf
          sp(j) = exact(l)
          err(j) = r(l)
          l = l+2
          k = k+2
  50    continue
        write(6,925) i,(err(j),j=1,7)
        write(6,930) (sp(j),j=1,7)
  60  continue
      ai = m
      sum = sum/ai
      write(6,935) sum,ermax
      write(6,940) exr0,exr1
      !  we set up the dimension information
      nuest = 19
      nvest = 21
      kwrk = 70
      lwrk = 1500
      ! main loop for the different spline approximations
      do 300 is=1,6
        go to (110,120,130,140,150,160),is
      !  we start computing a set of spline approximations with
      !  only c0-continuity at the poles,
 110    iopt(2) = 0
        ider(2) = 0
        iopt(3) = 0
        ider(4) = 0
      !  with no data values at the poles.
        ider(1) = -1
        ider(3) = -1
      !  initialisation
        iopt(1) = 0
      !  a large value for s for computing the least-squares polynomial
        s = 60.
        go to 200
      !  iopt(1) = 1 from the second call on
 120    s = 0.05
        iopt(1) = 1
        go to 200
      !  an interpolating spline
 130    s = zero
        go to 200
      !  a second set of approximations with c1-continuity at the poles
 140    iopt(2) = 1
        iopt(3) = 1
      !  exact values at the poles.
        ider(1) = 1
        ider(3) = 1
        r0 = exr0
        r1 = exr1
      ! reinitialization
        iopt(1) = 0
        s = 0.05
        go to 200
      !  vanishing derivatives at the poles
 150    ider(2) = 1
        ider(4) = 1
      ! reinitialization
        iopt(1) = 0
        go to 200
      ! finally we calculate the least-squares spline according to the current
      !  set of knots
 160    iopt(1) = -1
 200    call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest, &
          nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      ! printing of the fitting results.
        if(iopt(1)>=0) go to 210
        write(6,945)
        go to 220
 210    write(6,950) s
 220    write(6,955) iopt(2),iopt(3)
        if(ider(2)==1) write(6,960)
        if(ider(4)==1) write(6,965)
        write(6,970) fp,ier
        write(6,975) nu
        write(6,980)
        write(6,985) (tu(i),i=1,nu)
        write(6,990) nv
        write(6,980)
        write(6,985) (tv(i),i=1,nv)
        nc = (nu-4)*(nv-4)
        write(6,995)
        write(6,985) (c(i),i=1,nc)
      !  evaluation of the spline approximation
        call bispev(tu,nu,tv,nv,c,3,3,u,mu,v,mv,f,wk,100,iwrk,70,ier)
        write(6,1000)
        write(6,915) (j,j=1,mv,2)
        write(6,920)
        ermax = zero
        sum = zero
        do 240 i=1,mu
          k = i
          do 230 j=1,mv
            sp(j) = f(k)
            err(j) = abs(exact(k)-f(k))
            sum = sum+err(j)
            if(err(j)>ermax) ermax = err(j)
            k = k+1
 230      continue
          if( (i/2)*2 .ne.i ) go to 240
          write(6,925) i,(sp(j),j=1,mv,2)
          write(6,930) (err(j),j=1,mv,2)
 240    continue
        sum = sum/ai
        write(6,935) sum,ermax
        write(6,1005) c(1),c(nc)
 300  continue
      stop
 910  format(49h1data value (exact function value) at grid points)
 915  format(8h0v(j),j=,3x,7(i2,6x))
 920  format(8h u(i),i=)
 925  format(1x,i5,7(2x,f6.3))
 930  format(7x,7(2h (,f5.3,1h)))
 935  format(19h0mean abs. error = ,f9.3,5x,18hmax. abs. error = ,f9.3)
 940  format(30h function values at the poles ,f7.3,5x,f7.3)
 945  format(21h0least-squares spline)
 950  format(25h0smoothing spline with s=,f7.2)
 955  format(1x,34horder of continuity at the poles =,2i5)
 960  format(1x,37hvanishing derivatives at the pole u=0)
 965  format(1x,38hvanishing derivatives at the pole u=pi)
 970  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 975  format(1x,42htotal number of knots in the u-direction =,i3)
 980  format(1x,22hposition of the knots )
 985  format(5x,8f9.4)
 990  format(1x,42htotal number of knots in the v-direction =,i3)
 995  format(23h0b-spline coefficients )
1000  format(50h0spline value (approximation error) at grid points)
1005  format(28h spline values at the poles ,f7.3,5x,f7.3)
      end subroutine mnspgr

      !
      real(RKIND) function tesspg(u,v)
      ! function program tesspg calculates the value of the test function
      ! underlying the data.
      real(RKIND) ::u,v,sin,cos
      tesspg = 2./(4.1+cos(3.*u)+3.*cos(v+v+u*0.25)*sin(u)**2)
      return
      end function tesspg



      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c             mnsphe : sphere test program                           cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnsphe(datafile)
        real(RKIND), intent(in) :: datafile(:)
        real(RKIND) :: teta(192),phi(192),r(192),w(192),tp(30),tt(30),c(300), &
                       p(9),t(9),f(81),eps,fp,s,scale,scp,sct,ai
        real(RKIND), allocatable :: wrk1(:),wrk2(:)
        integer :: i,ier,iopt,j,kwrk,lwrk1,lwrk2,l1,l2,l,m,np,npest,nt,ntest, &
                   is,i1,i2,nc,ntt,npp,iwrk(300),pos

      allocate(wrk1(12000),wrk2(72))



      !  set constants
      scale = pi4/0.45e+02
      !  we fetch the number of data points.
      m = 192
      !  we fetch and print the latitude - longitude coordinates of the data
      !  points (in degrees).
      write(6,900)
      write(6,905)
      l2 = 0
      pos = 0
      do 10 i=1,48
         l1 = l2+1
         l2 = l2+4

         do l=l1,l2
            teta(l) = datafile(pos+1); pos = pos+1
            phi (l) = datafile(pos+1); pos = pos+1
         end do

         write(6,915)(teta(l),phi(l),l=l1,l2)
  10  continue
      !  we set up the weights, scale into radians the latitude-longitude
      !  coordinates and calculate the function values.
      do 20 i=1,m
         w(i) = 0.1e+01
         teta(i) = teta(i)*scale
         phi(i) = phi(i)*scale
         if(teta(i)>pi) teta(i) = pi
         if(phi(i)>pi2) phi(i) = pi2
         r(i) = testsp(teta(i),phi(i))
  20  continue
      !  we set up the coordinates of the grid points for the evaluation of
      !  the spline approximations.
      sct = pi/8
      scp = pi2/8
      do 30 i=1,8
         ai = i-1
         t(i) = ai*sct
         p(i) = ai*scp
  30  continue
      t(9) = pi
      p(9) = pi2
      ! we set up the dimension information
      ntest = 15
      npest = 19
      lwrk1 = 12000
      lwrk2 = 72
      kwrk = 300
      !  we choose a value for eps
      eps = 0.1e-05
      !  main loop for the different spline approximations
      do 300 is=1,4
        go to (110,120,130,140),is
      !  we start computing the least-squares constrained polynomial (large s)
 110    iopt = 0
        s = 500.
        go to 200
      !  iopt = 1 from the second call on.
 120    iopt = 1
        s = 135.
        go to 200
 130    s = 15.
        go to 200
      !  a least-squares spherical spline with specified knots.
 140    iopt = -1
      !  we set up the number of knots.
        nt = 11
        np = 15
      !  we set up the position of the interior knots of the spline.
        ntt = nt-8
        do 150 i=1,ntt
          ai = i
          j = i+4
          tt(j) = ai*pi4
 150    continue
        npp = np-8
        do 160 i=1,npp
          ai = i
          j = i+4
          tp(j) = ai*pi4
 160    continue
      !  determination of the spline approximation.
 200    call sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,nt,tt, &
         np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      !  printing of the fitting results.
        if(iopt>=0) go to 210
        write(6,920)
        go to 220
 210    write(6,925)
        write(6,930) s
 220    write(6,935) fp,ier
        write(6,940) nt
        write(6,945)
        write(6,950) (tt(i),i=1,nt)
        write(6,955) np
        write(6,945)
        write(6,950) (tp(i),i=1,np)
        nc = (nt-4)*(np-4)
        write(6,960)
        write(6,965) (c(i),i=1,nc)
      !  evaluation of the spline approximation.
        call bispev(tt,nt,tp,np,c,3,3,t,9,p,9,f,wrk2,lwrk2, &
         iwrk,kwrk,ier)
        write(6,970) (p(i),i=1,9)
        write(6,975)
        i2 = 0
        do 230 i=1,9
          i1 = i2+1
          i2 = i2+9
          write(6,980) t(i),(f(j),j=i1,i2)
 230    continue
 300  continue
      stop
      !  format statements.
 900  format(55h1latitude-longitude values of the data points (degrees))
 905  format(1h0,4(3x,10hteta   phi,3x))
 915  format(1h ,4(3x,f4.0,2x,f4.0,3x))
 920  format(50h0least-squares spline approximation on the sphere.)
 925  format(32h0smoothing spline on the sphere.)
 930  format(20h smoothing factor s=,f9.0)
 935  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
 940  format(1x,45htotal number of knots in the teta-direction =,i3)
 945  format(1x,22hposition of the knots )
 950  format(5x,8f8.4)
 955  format(1x,44htotal number of knots in the phi-direction =,i3)
 960  format(23h0b-spline coefficients )
 965  format(5x,8f9.4)
 970  format(9h      phi,9f7.3)
 975  format(6h  teta)
 980  format(1h ,f6.3,2x,9f7.3)

      end subroutine mnsphe

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mnspin : splint test program                       cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnspin
      real(RKIND) ::t(20),c(20),wrk(20)
      integer i,j,k,k1,n,nk1,ier
      real(RKIND) ::a,aint,ak,b,exint
      !  as an example we calculate some integrals of the form
      !          / b
      !         !     (1-x)**k  dx
      !      a /
      !
      !  main loop for the different spline degrees.
      do 30 k=1,5
        k1 = k+1
        ak = k1
      !  find the b-spline representation of the polynomial (1-x)**k.
        n = 2*k1
        j = n
        do 10 i=1,k1
          c(i) = 0.
          t(i) = 0.
          t(j) = 0.1e+01
          j = j-1
  10    continue
        c(1) = 0.1e+01
      !  insert a number of knots
        call insert(0,t,n,c,k,0.8d0,t,n,c,20,ier)
        call insert(0,t,n,c,k,0.4d0,t,n,c,20,ier)
        call insert(0,t,n,c,k,0.3d0,t,n,c,20,ier)
        call insert(0,t,n,c,k,0.1d0,t,n,c,20,ier)
      !  print the data for the spline.
        write(6,900) k
        write(6,905)
        write(6,910) (t(i),i=1,n)
        nk1 = n-k1
        write(6,915)
        write(6,920) (c(i),i=1,nk1)
      !  loop for the different integration limits a and b.
        a = 0.
        b = 0.1e+01
        write(6,925)
        do 20 j=1,4
      !  calculate the value of the spline integral
          aint = splint(t,n,c,k,a,b,wrk)
      !  calculate the exact value of the integral
          exint = ((0.1e01-a)**k1-(0.1e01-b)**k1)/ak
          write(6,930) a,b,aint,exint
          a = a+0.1e0
          b = b-0.3e0
  20    continue
  30  continue
      stop
      !  format statements.
 900  format(25h0degree of the spline k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,15f5.1)
 915  format(1x,21hb-spline coefficients)
 920  format(5x,8f9.5)
 925  format(1h0,5x,1ha,6x,1hb,6x,6hsplint,7x,5hexint)
 930  format(1x,2f7.1,2f12.5)
      end subroutine mnspin


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c                 mnspro : sproot test program                       cc
      !c      application : to find the intersection of a planar            cc
      !c      cubic spline curve   x = sx(u)   y = sy(u)   with             cc
      !c      a straight line   alfa*x + beta*y = gamma                     cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnspro
      real(RKIND) ::t(13),c(26),zero(20),sp(40),cc(13)
      integer i,idim,ier,is,i1,i2,j,k,k1,l1,l2,m,mest,n,nc,nk1
      real(RKIND) ::alfa,beta,gamma,per
      !  we have a planar curve
      idim = 2
      !  we have a cubic spline curve.
      k = 3
      k1 = k+1
      !  set up the dimension information
      nc = 26
      mest = 20
      !  n denotes the total number of knots.
      n = 13
      !  set up the knots of the spline curve
      t(4) = 0.
      t(5) = 0.2e0
      t(6) = 0.3e0
      t(7) = 0.5e0
      t(8) = 0.6e0
      t(9) = 0.7e0
      t(10) = 0.1e+01
      !  fetch the b-spline coefficients for sx(u)
      c(1) = 0.1e+01
      c(2) = 0.3e+01
      c(3) = 0.4e+01
      c(4) = 0.5e+01
      c(5) = 0.3e+01
      c(6) = -0.1e+01
      !  fetch the b-spline coefficients for sy(u)
      c(14) = 0.1e+01
      c(15) = 0.2e+01
      c(16) = -0.3e+01
      c(17) = 0.2e+01
      c(18) = 0.1e+01
      c(19) = 0.4e+01
      !  we have a closed curve.
      !  incorporate the boundary conditions for periodic splines
      per = t(10)-t(4)
      do 10 i=1,3
      !  the boundary knots
        t(i) = t(i+6)-per
        t(i+10) = t(i+4)+per
      !  the boundary coefficients
        c(i+6) = c(i)
        c(i+19) = c(i+13)
  10  continue
      !  print the data of the spline curve.
      write(6,900) k
      write(6,905)
      write(6,910) (t(i),i=1,n)
      write(6,915)
      nk1 = n-k1
      write(6,920) (c(i),i=1,nk1)
      write(6,925)
      i1 = n+1
      i2 = n+nk1
      write(6,920) (c(i),i=i1,i2)
      !  loop for the different lines.
      do 200 is=1,5
        go to (110,120,130,140,150),is
      !  fetch the parameters of the straight line.
 110    alfa = 0.
        beta = 0.1e+01
        gamma = 0.
        go to 160
 120    alfa = 0.1e+01
        beta = 0.
        go to 160
 130    beta = -0.1e+01
        go to 160
 140    alfa = 0.4e0
        beta = 0.3e0
        gamma = 0.12e+01
        go to 160
 150    beta = 0.4e0
        gamma = 0.
      !  print the parameters of the straight line.
 160    write(6,930) alfa,beta,gamma
      !  calculate the coefficients of s(u) = sx(u)*alfa + sy(u)*beta - gamma
        do 170 i=1,nk1
          j = i+n
          cc(i) = alfa*c(i)+beta*c(j)-gamma
 170    continue
      !  find the zeros of s(u)
        call sproot(t,n,cc,zero,mest,m,ier)
        write(6,935) m
        if(m==0) go to 200
      !  find the intersection points
        call curev(idim,t,n,c,nc,k,zero,m,sp,nc,ier)
      !  print the intersection points
        write(6,940)
        l2 = 0
        do 180 i=1,m
          l1 = l2+1
          l2 = l1+1
          write(6,945) i,zero(i),sp(l1),sp(l2)
 180    continue
 200  continue
      stop
      !  format statements.
 900  format(31h0degree of the spline curve k =,i2)
 905  format(1x,21hposition of the knots)
 910  format(5x,7f6.1)
 915  format(1x,30hb-spline coefficients of sx(u))
 920  format(5x,14f5.0)
 925  format(1x,30hb-spline coefficients of sy(u))
 930  format(18h0intersection with,f6.1,5h *x +,f5.1,5h *y =,f5.1)
 935  format(1x,33hnumber of intersection points m =,i3)
 940  format(6x,1hi,7x,4hu(i),5x,8hsx(u(i)),4x,8hsy(u(i)))
 945  format(1x,i6,3f12.5)
      end subroutine mnspro


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c              mnsuev : surev test program                           cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnsuev
      real(RKIND) ::fac
      integer i,idim,ier,j,m,mu,mv,m0,m1,m2,m3,nc,nu4,nv4,nu,nv,l
      real(RKIND) ::tu(11),tv(10),c(126),u(6),v(6),f(108),wrk(48)
      integer iwrk(12)
      !  we set up the grid points for evaluating the spline surface.
      mu = 6
      mv = 6
      do 10 i=1,6
      u(i) = (i-1)*0.2
      v(i) = u(i)
  10  continue
      !  the interior knots with respect to the u-variable.
      tu(5) = 0.4
      tu(6) = 0.7
      tu(7) = 0.9
      nu = 11
      !  the interior knots with respect to the v-variable.
      tv(5) = 0.3
      tv(6) = 0.8
      nv = 10
      !  the boundary knots
      do 20 i=1,4
        tu(i) = 0.
        tv(i) = 0.
        tu(i+7) = 1.
        tv(i+6) = 1.
  20  continue
      !  we generate the b-spline coefficients for the test surface
      !        x = u*v    y = v**2    z = u+v     0 <= u,v <= 1
      !  the dimension of the surface
      idim = 3
      !  the number of b-spline coefficients for each co-ordinate
      nu4 = nu-4
      nv4 = nv-4
      nc = nu4*nv4
      !  the coefficients for x = u*v
      do 30 i=1,nv4
        c(i) = 0.
  30  continue
      do 40 i=2,nu4
        c((i-1)*nv4+1) = 0.
  40  continue
      m0 = 1
      do 60 i=2,nu4
        m1 = m0+nv4
        fac = (tu(i+3)-tu(i))/9.
        do 50 j=2,nv4
          m2 = m0+1
          m3 = m1+1
          c(m3) = c(m1)+c(m2)-c(m0)+fac*(tv(j+3)-tv(j))
          m0 = m0+1
          m1 = m1+1
  50    continue
        m0 = m0+1
  60  continue
      !  the coefficients for y = v**2.
      l = nc
      m0 = l+1
      m1 = m0+1
      c(m0) = 0.
      c(m1) = 0.
      do 70 i=3,nv4
        c(m1+1) = c(m1)+(tv(i+3)-tv(i))*((c(m1)-c(m0))/(tv(i+2)-tv(i-1)) &
          +(tv(i+2)-tv(i))/3.)
        m0 = m1
        m1 = m0+1
  70  continue
      do i=1,nv4
        m0 = l+i
        fac = c(m0)
        do j=1,nu4
          m0 = m0+nv4
          c(m0) = fac
        end do
      end do
      !  the coefficients for z = u+v
      l = l+nc
      m0 = l+1
      c(m0) = 0.
      do 90 i=2,nv4
        m1 = m0+1
        c(m1) = c(m0)+(tv(i+3)-tv(i))/3.
        m0 = m1
  90  continue
      do i=1,nv4
        m0 = l+i
        do j=2,nu4
          m1 = m0+nv4
          c(m1) = c(m0)+(tu(j+3)-tu(j))/3.
          m0 = m1
        end do
      end do
      !  evaluation of the spline surface
      call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,108,wrk,48,iwrk,12,ier)
      !  printing of the results
      write(6,900)
      write(6,910)
      write(6,920) (tu(i),i=1,nu)
      write(6,930)
      write(6,920) (tv(i),i=1,nv)
      write(6,940)
      m1 = 0
      do 110 l=1,idim
        m0 = m1+1
        m1 = m1+nc
        write(6,950) (c(j),j=m0,m1)
 110  continue
      write(6,960)
      write(6,970) (v(i),i=1,mv)
      write(6,980)
      m = mu*mv
      m0 = 0
      do 130 i=1,mu
        write(6,990) u(i)
        m1 = m0
        do 120 l=1,idim
          m2 = m1+1
          m3 = m1+mv
          write(6,995) (f(j),j=m2,m3)
          m1 = m1+m
 120    continue
        m0 = m0+mv
 130  continue
      stop
      !  format statements.
 900  format(23h0bicubic spline surface)
 910  format(1x,40hposition of the knots in the u-direction)
 920  format(1x,15f5.1)
 930  format(1x,40hposition of the knots in the v-direction)
 940  format(23h b-spline coefficients )
 950  format(5x,8f9.4)
 960  format(1h0,37hspline values at selected grid points)
 970  format(1h0,2x,1hv,6(3x,f4.1))
 980  format(1h ,1x,1hu)
 990  format(1h ,f4.1)
 995  format(5x,6f7.3)
      end subroutine mnsuev


      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !c                                                                    cc
      !c        mnsurf : surfit test program                                cc
      !c                                                                    cc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnsurf(xyz,delta,iunit)

         ! Surface points fetched from an external file
         real(RKIND), intent(in) :: xyz(:,:), delta
         integer, optional, intent(in) :: iunit

         ! Local variables
         real(RKIND), allocatable :: wrk1(:),wrk2(:)
         real(RKIND), dimension(size(xyz,2)) :: x,y,z,w
         real(RKIND) :: tx(15),ty(15),c(200)
         real(RKIND) :: ai,eps,fp,s,xb,xe,yb,ye,xx(11),yy(11),zz(121)
         integer :: iwrk(300),i,ier,iopt,is,j,kwrk,kx,ky,lwrk1,lwrk2,m,mx,my,nc, &
                    nmax,nx,nxest,ny,nyest,useUnit

         ! Output unit
         if (present(iunit)) then
             useUnit = iunit
         else
             useUnit = output_unit
         end if

         !  we fetch the number of data points
         m = size(xyz,2)
         write(useUnit,905) m

         !  we fetch the co-ordinate and function values of each data point.
         x = xyz(1,:)
         y = xyz(2,:)
         z = xyz(3,:)
         write(useUnit,910)

         ! Only print half of the points
         do i=1,m
            if (mod(i,2)/=0) cycle
            j = i-1
            write(useUnit,920) j,x(j),y(j),z(j),i,x(i),y(i),z(i)
         end do

         !  we fetch an estimate of the standard deviation of the data values.
         write(useUnit,930) delta

         !  the weights are set equal to delta**(-1)
         w = one/delta

         !  we set up the boundaries of the approximation domain.
         xb = -two
         xe = +two
         yb = -two
         ye = +two

         ! we generate a rectangular grid for evaluating the splines.
         mx = 11
         my = 11
         do i=1,11
           ai    = i-6
           xx(i) = ai*0.4
           yy(i) = xx(i)
         end do

         ! we set up the dimension information
         nxest = 15
         nyest = 15
         nmax  = 15
         kwrk  = 300
         lwrk1 = 12000
         lwrk2 = 6000
         allocate(wrk1(lwrk1),wrk2(lwrk2))

         ! we choose a value for eps
         eps=0.1e-05

         ! main loop for the different spline approximations.
         all_tests: do is=1,6

            select case (is)

               case (1) !  we start computing the least-squares bicubic polynomial (large s)
                  iopt = 0
                  kx = 3
                  ky = 3
                  s = 900000.

               case (2) !  iopt=1 from the second call on.
                  iopt = 1
                  s = 200.

               case (3) !  a value for s within its confidence interval
                  s = m

               case (4) !  overfitting (s too small)
                  s = 20.

               case (5) !  we change the degrees of the spline
                  iopt = 0
                  kx   = 5
                  ky   = 5
                  s    = m

               case (6) !  calculate a least-squares spline approximation with specified knots.
                  iopt = -1
                  kx = 3
                  ky = 3
                  nx = 11
                  ny = 11
                  j = kx+2
                  do i=1,3
                     ai = i-2
                     tx(j) = ai
                     ty(j) = ai
                     j = j+1
                  end do

            end select

            ! determination of the spline approximation.
            call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
                        nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
            ! printing of the fitting results.
            if (iopt>=0) then
               write(useUnit,940) kx,ky
               write(useUnit,945) s
            else
               write(useUnit,935) kx,ky
            endif

            write(useUnit,950) fp,ier
            write(useUnit,955) nx
            write(useUnit,960)
            write(useUnit,965) (tx(i),i=1,nx)
            write(useUnit,970) ny
            write(useUnit,960)
            write(useUnit,965) (ty(i),i=1,ny)
            nc = (nx-kx-1)*(ny-ky-1)
            write(useUnit,975)
            write(useUnit,980) (c(i),i=1,nc)

            ! evaluation of the spline approximation.
            call bispev(tx,nx,ty,ny,c,kx,ky,xx,mx,yy,my,zz,wrk2,lwrk2,iwrk,kwrk,ier)

            write(useUnit,1000)
            write(useUnit,985) (xx(i),i=1,mx)
            write(useUnit,990)
            do j=1,my
               write(useUnit,995) yy(j),(zz(i),i=j,121,11)
            end do

         end do all_tests

         ! format statements.
         905  format(1h1,i3,12h data points)
         910  format(1h0,2(2x,1hi,5x,4hx(i),6x,4hy(i),6x,4hz(i),6x))
         920  format(1x,2(i3,3f10.4,5x))
         930  format(1x,40hestimate of standard deviation of z(i) =,e15.6)
         935  format(32h0least-squares spline of degrees,2i3)
         940  format(28h0smoothing spline of degrees,2i3)
         945  format(20h smoothing factor s=,f9.0)
         950  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
         955  format(1x,42htotal number of knots in the x-direction =,i3)
         960  format(1x,22hposition of the knots )
         965  format(5x,10f7.3)
         970  format(1x,42htotal number of knots in the y-direction =,i3)
         975  format(23h0b-spline coefficients )
         980  format(5x,8f9.4)
         985  format(1h0,1hx,2x,11f7.1)
         990  format(3x,1hy)
         995  format(1x,f4.1,11f7.3)
         1000 format(1h0,33hspline evaluation on a given grid)

      end subroutine mnsurf




end module fitpack_tests
