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
module fitpack_curve_tests
    use fitpack
    use fitpack_core
    use fitpack_test_data, only: dapola,dasphe,daregr_x,daregr_y,daregr_z,dapogr, &
                                 daspgr_r,daspgr_u,daspgr_v,dapasu
    use iso_fortran_env, only: output_unit


    implicit none
    private

    public :: test_sine_fit
    public :: test_zeros
    public :: test_periodic_fit
    public :: test_parametric_fit
    public :: test_closed_fit
    public :: test_polar_fit
    public :: test_polar_boundary
    public :: test_sphere_fit
    public :: test_gridded_fit
    public :: test_constrained_curve
    public :: test_gridded_polar
    public :: test_gridded_sphere
    public :: test_parametric_surface
    public :: test_fpknot_crash
    public :: test_curve_comm_roundtrip
    public :: test_comm_roundtrips
    public :: test_umbrella_api
    public :: test_convex_fit
    public :: test_convex_least_squares
    public :: test_insert_knot
    public :: test_insert_knot_periodic
    public :: test_surface_integral
    public :: test_cross_section
    public :: test_derivative_spline


    contains

    ! This test reproduces https://github.com/scipy/scipy/issues/3691
    logical function test_fpknot_crash() result(success)
       
       integer(FP_SIZE) :: j,ierr
       type(fitpack_curve) :: spline
       real(FP_REAL), parameter :: x(*) = [(j, j=0,108)]
       real(FP_REAL), parameter :: y(*) = [ &
            0., 0., 0., 0., 0., 10.9, 0., 11., 0., &
            0., 0., 10.9, 0., 0., 0., 0., 0., 0., &
            10.9, 0., 0., 0., 11., 0., 0., 0., 10.9, &
            0., 0., 0., 10.5, 0., 0., 0., 10.7, 0., &
            0., 0., 11., 0., 0., 0., 0., 0., 0., &
            10.9, 0., 0., 10.7, 0., 0., 0., 10.6, 0., &
            0., 0., 10.5, 0., 0., 10.7, 0., 0., 10.5, &
            0., 0., 11.5, 0., 0., 0., 10.7, 0., 0., &
            10.7, 0., 0., 10.9, 0., 0., 10.8, 0., 0., &
            0., 10.7, 0., 0., 10.6, 0., 0., 0., 10.4, &
            0., 0., 10.6, 0., 0., 10.5, 0., 0., 0., &
            10.7, 0., 0., 0., 10.4, 0., 0., 0., 10.8, 0. ]
  
       
       ierr = spline%new_fit(x,y,order=1)
       
       success = FITPACK_SUCCESS(ierr)
    
    end function test_fpknot_crash

    ! Test parametric curve with constraints
    logical function test_constrained_curve(iunit) result(success)
       integer, optional, intent(in) :: iunit

       integer :: useUnit,ierr,loop,i
       type(fitpack_constrained_curve) :: curve
       integer, parameter :: m    = 31
       integer, parameter :: idim = 2
       real(FP_REAL) :: x(idim,m),ddx_begin(2,0:2),ddx_end(2,0:2)
       real(FP_REAL), allocatable :: w(:),ybegin(:,:),yend(:,:),u(:)
       character(64) :: whereAt

       success = .true.
       if (present(iunit)) then
          useUnit = iunit
       else
          useUnit = output_unit
       end if

       ! Curve coordinates (2D)
       x(1,:) = [-3.109,-2.188,-1.351,-0.605,0.093,0.451,0.652,0.701,0.518,0.277,0.008,-0.291,-0.562,-0.679,&
                 -0.637,-0.425,-0.049,0.575,1.334,2.167,3.206,4.099,4.872,5.710,6.330,6.741,6.928,6.965,6.842,&
                 6.593,6.269]
       x(2,:) = [3.040,2.876,2.634,2.183,1.586,1.010,0.382,-0.218,-0.632,-0.879,-0.981,-0.886,-0.642,-0.195,&
                 0.373,1.070,1.607,2.165,2.618,2.905,2.991,2.897,2.615,2.164,1.617,0.977,0.383,-0.194,-0.665,&
                 -0.901,-1.010]

       ! Data and derivatives at the extremes (point, 1st, 2nd derivative)
       ddx_begin = reshape([-pi,three,  three,zero,  zero,-two],[2,3])
       ddx_end   = reshape([pi2,-one,  -one,zero,  zero,two],[2,3])

       ! weights: 1/sigma with sigma an estimate of the standard deviation of the data points.
       allocate(w(m),source=one/0.04_FP_REAL)

       ! Set up the parameter values for the data points
       u = linspace(zero,3*pi,m)

       do loop=1,6

          select case (loop)


             case (1) ! the smoothing factor is chosen as s = m
                ierr = curve%new_fit(x,u=u,w=w,smoothing=real(m,FP_REAL))
                whereAt = 'no constraints'

             case (2) ! Fix end points

                call curve%set_constraints(ddx_begin(:,0:0),ddx_end(:,0:0),ierr)
                ierr = curve%fit(smoothing=real(m,FP_REAL))
                whereAt = 'fix points only'

             case (3) ! Fix points, 1st derivative (begin)

                call curve%set_constraints(ddx_begin(:,0:1),ddx_end(:,0:0),ierr)
                ierr = curve%fit(smoothing=real(m,FP_REAL))
                whereAt = 'fix points + 1st derivative (left)'

             case (4) ! Fix both first derivatives

                call curve%set_constraints(ddx_begin(:,0:1),ddx_end(:,0:1),ierr)
                ierr = curve%fit(smoothing=real(m,FP_REAL))
                whereAt = 'fix points + 1st derivative (both)'

             case (5)

                whereAt = 'quintic spline, 2st derivative constraints (both)'
                call curve%set_constraints(ddx_begin(:,0:2),ddx_end(:,0:2),ierr)
                ierr = curve%fit(smoothing=real(m,FP_REAL),order=5)

             case (6) ! Reduce smoothing

                whereAt = 'quintic spline, 2st derivative constraints (both), s=26'
                ierr = curve%fit(smoothing=26.0_FP_REAL,order=5)

             case (7)

                whereAt = ''
!                ! finally we also calculate a least-squares curve with specified knots
!                iopt = -1
!                j = k+2
!                set_knots: do l=1,5
!                   ai = l-2
!                   t(j) = ai*pi*half
!                   j = j+1
!                end do set_knots
!                n = 7+2*k

          end select

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
              exit
          end if

          ! Calculate derivatives at the begin point.
          ybegin = curve%dfdx_all(curve%ubegin,ierr)
          write(useUnit,960) trim(whereAt)
          do i=0,curve%order
             write(useUnit,970) i,ybegin(:,i+1)
          end do

          ! Calculate derivatives at the end point
          yend = curve%dfdx_all(curve%uend,ierr)
          write(useUnit,965) trim(whereAt)
          do i=0,curve%order
             write(useUnit,970) i,yend(:,i+1)
          end do

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
              exit
          end if

          ! Verify the constraints
          select case (loop)
             case (2)
               if (.not.maxval(abs(ybegin(:,1)-ddx_begin(:,0)))<smallnum06) success = .false.
               if (.not.maxval(abs(  yend(:,1)-ddx_end  (:,0)))<smallnum06) success = .false.
             case (3)
               if (.not.maxval(abs(ybegin(:,1:2)-ddx_begin(:,0:1)))<smallnum06) success = .false.
               if (.not.maxval(abs(  yend(:,1)-ddx_end  (:,0)))<smallnum06)     success = .false.
             case (4)
               if (.not.maxval(abs(ybegin(:,1:2)-ddx_begin(:,0:1)))<smallnum06) success = .false.
               if (.not.maxval(abs(  yend(:,1:2)-ddx_end  (:,0:1)))<smallnum06) success = .false.
             case (5:)
               if (.not.maxval(abs(ybegin(:,1:3)-ddx_begin(:,0:2)))<smallnum06) success = .false.
               if (.not.maxval(abs(  yend(:,1:3)-ddx_end  (:,0:2)))<smallnum06) success = .false.
          end select

          if (.not.success) then
              write(useUnit,1000) loop,'Constraints failed accuracy test'
              exit
          end if

       end do

        960 format('[test_constrained_curve] begin point (',a,'): ')
        965 format('[test_constrained_curve] endpoint (',a,'): ')
        970 format(5x,'order=',i2,2f9.4)
       1000 format('[test_constrained_curve] test ',i0,' failed: ',a)

    end function test_constrained_curve

    ! Test closed parametric curve
    logical function test_closed_fit(iunit) result(success)
       integer, optional, intent(in) :: iunit

       real(FP_REAL), allocatable :: x(:,:),w(:),y(:,:),u(:)
       type(fitpack_closed_curve) :: curve
       integer :: ierr,loop,useUnit,i

       success = .true.

       if (present(iunit)) then
           useUnit = iunit
       else
           useUnit = output_unit
       end if

       allocate(x(2,19))

       x(1,1:18) = [-4.7,-7.048,-6.894,-3.75,-1.042,0.938,2.5,3.524,4.511,5.0,4.886,3.524,3.2,1.302,-1.424,&
                    -3.0,-3.064,-3.665]
       x(2,1:18) = [0.0,2.565,5.785,6.495,5.909,5.318,4.33,2.957,1.642,0.0,-1.779,-2.957,-5.543,-7.386,-8.075,&
                    -5.196,-2.571,-1.334]

       ! Set closed curve
       x(:,19) = x(:,1)

       ! Use flat weights, do not supply parameter values
       allocate(w(size(x,2)),source=one)

       ! Supply parameter values (also used as evaluation points)
       u = [(20*i,i=0,size(x,2)-1)]


       do loop=1,8

          select case (loop)

             ! start with a least-squares point (s is very large)
             case (1); ierr = curve%new_fit(x,w=w,smoothing=900.0_FP_REAL)

             ! Smaller values of s to get a tighter approximation
             case (2); ierr = curve%fit(smoothing=10.0_FP_REAL)
             case (3); ierr = curve%fit(smoothing=0.1_FP_REAL)

             ! Larger values to get a smoother approximation
             case (4); ierr = curve%fit(smoothing=0.5_FP_REAL)

             ! New fit to get possibly fewer knots
             case (5); ierr = curve%new_fit(x,w=w,smoothing=0.5_FP_REAL)

             ! Let the program determine parameter values u(i)
             case (6); ierr = curve%new_fit(x,w=w,smoothing=0.5_FP_REAL)

             ! Quintic spline approximation
             case (7); ierr = curve%new_fit(x,w=w,smoothing=0.5_FP_REAL,order=5)

             ! Interpolating curve
             case (8); ierr = curve%interpolate()

!                 case (9)
!
!                   ! finally we calculate a least-squares spline curve with specified knots
!                   iopt = -1
!                   n = 9+2*k
!                   j = k+2
!                   del = (u(m)-u(1))*0.125_FP_REAL
!                   do l=1,7
!                      al = l
!                      t(j) = u(1)+al*del
!                      j = j+1
!                   end do

          end select

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
              exit
          end if

          ! Evaluate the spline curve
          y = curve%eval_many(u,ierr)

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
              exit
          end if

       end do

       1000 format('[test_closed_fit] parametric curve test ',i0,' failed: ',a)

    end function test_closed_fit

    ! Test parametric curve
    logical function test_parametric_fit() result(success)

       real(FP_REAL), allocatable :: x(:,:),u(:),w(:),y(:,:)
       type(fitpack_parametric_curve) :: curve
       integer :: ierr,loop,useUnit

       success = .true.

       useUnit = output_unit

       ! Set data points
       u = [real(FP_REAL) :: 120,128,133,136,138,141,144,146,149,151,154,161,170,180,190,&
                           200,210,220,230,240,250,262,269,273,278,282,287,291,295,299,305,315]

       allocate(x(2,32))

       x(1,:) = [-1.5141,-2.0906,-1.9253,-0.8724,-0.3074,-0.5534,0.0192,1.2298,2.5479,2.4710,1.7063,&
                 1.1183,0.5534,0.4727,0.3574,0.1998,0.2882,0.2613,0.2652,0.2805,0.4112,0.9377,1.3527,&
                 1.5564,1.6141,1.6333,1.1567,0.8109,0.2498,-0.2306,-0.7571,-1.1222]
       x(2,:) = [0.5150,1.3412,2.6094,3.2358,2.7401,2.7823,3.5932,3.8353,2.5863,1.3105,0.6841,0.2575,&
                 0.2460,0.3689,0.2460,0.2998,0.3651,0.3343,0.3881,0.4573,0.5918,0.7110,0.4035,0.0769,&
                 -0.3920,-0.8570,-1.3412,-1.5641,-1.7409,-1.7178,-1.2989,-0.5572]

       ! Use flat weights
       allocate(w(size(x,2)),source=one)

       do loop=1,8

          select case (loop)

             ! start with a polynomial curve ( s very large)
             case (1); ierr = curve%new_fit(x,u=u,w=w,smoothing=100.0_FP_REAL)

             ! Smaller values of s to get a tighter approximation
             case (2); ierr = curve%fit(smoothing=1.0_FP_REAL)
             case (3); ierr = curve%fit(smoothing=0.05_FP_REAL)

             ! Larger values to get a smoother approximation
             case (4); ierr = curve%fit(smoothing=0.25_FP_REAL)

             ! New fit to get possibly fewer knots
             case (5); ierr = curve%new_fit(x,u=u,w=w,smoothing=0.25_FP_REAL)

             ! Let the program determine parameter values u(i)
             case (6); ierr = curve%new_fit(x,w=w,smoothing=0.25_FP_REAL)

             ! Quintic spline approximation
             case (7); ierr = curve%new_fit(x,w=w,smoothing=0.25_FP_REAL,order=5)

             ! Interpolating curve
             case (8); ierr = curve%interpolate()

!                 case (9)
!
!                    !  finally we calculate a least-squares spline curve with specified knots
!                    iopt =-1
!                    n = 9+2*k
!                    j = k+2
!                    del = (ue-ub)*0.125_FP_REAL
!                    do l=1,7
!                        al = l
!                        t(j) = ub+al*del
!                        j = j+1
!                    end do

          end select

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
              exit
          end if

          ! Evaluate the spline curve
          y = curve%eval(u,ierr)

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
              exit
          end if

       end do

       1000 format('[test_parametric_fit] parametric curve test ',i0,' failed: ',a)

    end function test_parametric_fit

    ! Dumb test: fit a sine function
    logical function test_sine_fit() result(success)

       integer, parameter     :: N = 200
       real(FP_REAL), parameter :: RTOL = 1.0e-1_FP_REAL
       real(FP_REAL), parameter :: ATOL = 1.0e-2_FP_REAL
       type(fitpack_curve) :: curve
       real(FP_REAL) :: x(N),y(N),xrand(N),yeval,yprime,dfdx(0:3),fint,eint
       integer :: ierr,i,order

       success = .false.

       ! Generate a sine function with 200 points
       x = linspace(zero,pi2,N)
       y = sin(x)

       ! Create INTERPOLATING
       ierr = curve%new_fit(x,y,smoothing=zero)

       ! Failed to create
       if (.not.FITPACK_SUCCESS(ierr)) then
          print *, '[sine_fit] error generating sine function '
          return
       end if

       ! Create 200 points in between the range
       xrand(1:N-1) = half*(x(1:N-1)+x(2:N))
       xrand(N)     = x(N)

       do i=1,n

          ! Evaluate curve
          yeval = curve%eval(xrand(i),ierr)
          if (.not.FITPACK_SUCCESS(ierr)) then
             print *, '[sine_fit] error evaluating sine function '
             return
          end if

          ! Get analytical function and derivatives
          dfdx(0) = sin(xrand(i)) ! the function
          dfdx(1) = cos(xrand(i))
          dfdx(2) = -dfdx(0)
          dfdx(3) = -dfdx(1)

          ! error
          if (abs(yeval-dfdx(0))*rewt(RTOL,ATOL,dfdx(0))>one) then
             print 1, xrand(i),yeval,dfdx(0)
             return
          end if

          ! Evaluate first derivative
          do order = 1,3
             yprime = curve%dfdx(xrand(i),order=order,ierr=ierr)

             ! Check evaluation
             if (.not.FITPACK_SUCCESS(ierr)) then
               print 2, order,xrand(i),FITPACK_MESSAGE(ierr)
               return
             end if

             ! Check error
             if (abs(yprime-dfdx(order))*rewt(RTOL,ATOL,dfdx(order))>one) then
                print 3, order,xrand(i),yprime,dfdx(order)
                return
             end if
          end do

       end do

       ! Evaluate integral, compare with exact
       fint = curve%integral(from=pi/3,to=two*pi)
       eint = cos(pi/3)-cos(two*pi)

       ! Check error
       if (abs(eint-fint)*rewt(RTOL,ATOL,eint)>one) then
          print 4, pi/3,two*pi,fint,eint
          success = .false.
          return
       end if

       ! All checks passed: success!
       success = .true.

       1 format('[sine_fit] sine function error is too large: x=',f6.2,' yspline=',f6.2,' analytical=',f6.2)
       2 format('[sine_fit] cannot evaluate ',i0,'-th derivative at ',f6.2,': ',a)
       3 format('[sine_fit] ',i0,'-th derivative error is too large: x=',f6.2,' yp(spline)=',f6.2,' analytical=',f6.2)
       4 format('[sine_fit] integral error is too large: [a=',f6.2,',b=',f6.2,'] int(spline)=',f6.2,' analytical=',f6.2)

    end function test_sine_fit

    ! Test zeros of a cubic spline
    logical function test_zeros() result(success)

       type(fitpack_curve) :: curve
       real(FP_REAL), allocatable :: x(:),y(:),x0(:)
       integer :: ierr

       success = .false.

       ! Try f(x) = x3 � 3x2 + 2x = x(x � 1)(x � 2), with real roots x = 0, x = 1, and x = 2.
       x = linspace(-10.0_FP_REAL,10.0_FP_REAL,20)
       y = x**3-3*x**2+2*x

       ! Get the interpolating curve
       ierr = curve%new_fit(x,y,smoothing=zero)

       ! Get the zeros
       x0 = curve%zeros(ierr)

       success = FITPACK_SUCCESS(ierr) .and. size(x0)==3 &
                 .and. all(abs(x0(1:3)-[zero,one,two])<sqrt(epsilon(zero)))

    end function test_zeros

    ! Periodic test: fit a cosine function
    logical function test_periodic_fit() result(success)

       integer, parameter     :: N = 200
       real(FP_REAL), parameter :: RTOL = 1.0e-2_FP_REAL
       real(FP_REAL), parameter :: ATOL = 1.0e-4_FP_REAL
       type(fitpack_periodic_curve) :: curve
       real(FP_REAL) :: x(N),y(N),xrand(N),yeval,yprime,dfdx(0:3),fint,eint
       integer :: ierr,i,order

       success = .false.

       ! Generate a periodic function with 200 points
       x = linspace(zero,pi2,N)
       y = fun(x)

       ! Create INTERPOLATING
       ierr = curve%new_fit(x,y,smoothing=zero)

       ! Failed to create
       if (.not.FITPACK_SUCCESS(ierr)) then
          print *, '[periodic_fit] error generating sine function '
          return
       end if

       ! Create 200 points in between the range. Include both extremes
       xrand(1)     = x(1)
       xrand(2:N-1) = half*(x(1:N-2)+x(2:N-1))
       xrand(N)     = x(n)

       do i=1,n

          ! Evaluate curve
          yeval = curve%eval(xrand(i),ierr);
          if (.not.FITPACK_SUCCESS(ierr)) then
             print *, '[periodic_fit] error evaluating sine function '
             return
          end if

          ! Get analytical function and derivatives
          dfdx(0) =  cos(xrand(i)) +   sin(2*xrand(i)) ! the function
          dfdx(1) = -sin(xrand(i)) + 2*cos(2*xrand(i)) ! derivative
          dfdx(2) = -cos(xrand(i)) - 4*sin(2*xrand(i))
          dfdx(3) = +sin(xrand(i)) - 8*cos(2*xrand(i))

          ! error
          if (abs(yeval-dfdx(0))*rewt(RTOL,ATOL,dfdx(0))>one) then
             print 1, xrand(i),yeval,dfdx(0)
             return
          end if

          ! Evaluate first derivative
          do order = 1,3
             yprime = curve%dfdx(xrand(i),order=order,ierr=ierr)

             ! Check evaluation
             if (.not.FITPACK_SUCCESS(ierr)) then
               print 2, order,xrand(i),FITPACK_MESSAGE(ierr)
               return
             end if

             ! Check error
             if (abs(yprime-dfdx(order))*rewt(RTOL,ATOL,dfdx(order))>one) then
                print 3, order,xrand(i),yprime,dfdx(order)
                return
             end if
          end do

       end do

       ! Test integral calculation *across the periodic boundary*
       fint = curve%integral(from=half*pi,to=3*pi)
       eint = intgl(3*pi)-intgl(half*pi)

       ! Check error
       if (abs(eint-fint)*rewt(RTOL,ATOL,eint)>one) then
          print 4, pi/3,two*pi,fint,eint
          return
       end if

       ! All checks passed: success!
       success = .true.

       1 format('[periodic_fit] sine function error is too large: x=',f6.2,' yspline=',f6.2,' analytical=',f6.2)
       2 format('[periodic_fit] cannot evaluate ',i0,'-th derivative at ',f6.2,': ',a)
       3 format('[periodic_fit] ',i0,'-th derivative error is too large: x=',f6.2,' yp(spline)=',f6.2,' analytical=',f6.2)
       4 format('[periodic_fit] integral error is too large: [a=',f6.2,',b=',f6.2,'] int(spline)=',f6.2,' analytical=',f6.2)

       contains

         elemental real(FP_REAL) function fun(x) result(y)
            real(FP_REAL), intent(in) :: x
            y = cos(x) + sin(2*x)
         end function fun

         elemental real(FP_REAL) function intgl(x) result(y)
            real(FP_REAL), intent(in) :: x
            y = sin(x) - half*cos(2*x)
         end function intgl

    end function test_periodic_fit

    logical function test_polar_fit(iunit) result(success)
       integer, optional, intent(in) :: iunit

       type(fitpack_polar) :: polar
       integer :: m,useUnit,loop,ierr
       character(64) :: domain
       real(FP_REAL), allocatable :: x(:),y(:),z(:),w(:),f(:),exact(:)
       real(FP_REAL) :: avg,ermax

       ! Initialization.
       success = .true.
       if (present(iunit)) then
           useUnit = iunit
       else
           useUnit = output_unit
       end if

       ! Get datapoints from array "dapola"
       x = dapola(1:600-2:3)
       y = dapola(2:600-1:3)
       z = dapola(3:600  :3)

       ! And the exact function value
       exact = testpo(x,y)

       ! Number of points
       m = size(x)

       ! Set uniform weights: w(i)=(0.01)**(-1), where 0.01 is an estimate for the standard deviation
       ! of the error in z(i)).
       allocate(w(m),f(m),source=100.0_FP_REAL)

       ! we determine a number of smoothing spline approximations on the unit disk: x**2+y**2 <= 1.
       approximations: do loop=1,6

           ! Compute/update fit
           select case (loop)
              case (1)
                  domain = 'unity disk'
                  ierr = polar%new_fit(x,y,z,unit_disk,w,smoothing=1500.0_FP_REAL)
              case (2)
                  ierr = polar%fit(smoothing=200.0_FP_REAL)
              case (3)
                  ierr = polar%fit(smoothing=170.0_FP_REAL)
              case (4)
                  ! Ellipsoid: Only consider datapoints inside this domain
                  domain = ' ellipsoid'
                  x = x(1:90)
                  y = y(1:90)
                  z = z(1:90)
                  w = w(1:90)

                  ! Impose z=0 at the boundary by shifting all previous z's
                  z     = z - 0.4_FP_REAL
                  exact = testpo(x,y) - 0.4_FP_REAL
                  ierr  = polar%new_fit(x,y,z,rad2_ellipsoid,w,smoothing=90.0_FP_REAL)
              case (5)

                  ! Determine least-squares approximation with current knots
                  ierr = polar%least_squares()

              case (6)

                  ! Determine interpolating spline
                  ierr = polar%interpolate()

           end select

           if (.not.FITPACK_SUCCESS(ierr)) then
               success = .false.
               write(useUnit,1000) loop,trim(domain),FITPACK_MESSAGE(ierr)
               exit
           end if

           ! Evaluate fit at the initial points
           f = polar%eval(x,y,ierr)

           if (.not.FITPACK_SUCCESS(ierr)) then
               success = .false.
               write(useUnit,1000) loop,trim(domain),FITPACK_MESSAGE(ierr)
               exit
           end if

           ! Determine mean and maximum errors
           avg   = sum(abs(f-exact))/m
           ermax = maxval(abs(f-exact),1)

           write(useUnit,920) trim(domain),polar%smoothing,avg,ermax

       end do approximations

        920 format('[test_polar_fit] ',a,': s=',f7.1,', mean error = ',f7.4,5x,'max. error = ',f7.4)
       1000 format('[test_polar_fit] polar test ',i0,' (',a,' domain) failed: ',a)

       contains

          !  test function for the polar package
          elemental real(FP_REAL) function testpo(x,y)
              real(FP_REAL), intent(in) ::x,y
              testpo=(x**2+y**2)/((x+y)**2+half)
          end function testpo

          ! CIRCLE: the boundary of the approximation domain  x**2+y**2<=1. in polar coordinates
          pure real(FP_REAL) function unit_disk(v)
             real(FP_REAL), intent(in) :: v
             unit_disk = one
             return
          end function unit_disk

          ! ELLIPSOID: the boundary of the approximation domain  3*x**2+3*y**2-4*x*y<=1. in polar coordinates
          pure real(FP_REAL) function rad2_ellipsoid(v)
             real(FP_REAL), intent(in) :: v
             rad2_ellipsoid = one/sqrt(three-two*sin(2*v))
             return
          end function rad2_ellipsoid

    end function test_polar_fit

    ! Test polar fitting with various boundary shapes
    ! This specifically tests the abstract interface fitpack_polar_boundary
    logical function test_polar_boundary(iunit) result(success)
       integer, optional, intent(in) :: iunit

       type(fitpack_polar) :: polar
       integer :: useUnit, ierr, m
       real(FP_REAL), allocatable :: theta(:), r(:), x(:), y(:), z(:), w(:), f(:)
       real(FP_REAL) :: avg, ermax

       success = .true.
       useUnit = output_unit; if (present(iunit)) useUnit = iunit

       ! Generate test points inside a cardioid: r = 0.5*(1 + cos(theta))
       m = 100
       theta = linspace(zero, pi2, m)
       r     = 0.4_FP_REAL * (one + cos(theta)) * linspace(0.3_FP_REAL, one, m)
       x     = r * cos(theta)
       y     = r * sin(theta)
       z     = test_func(x, y)
       w     = [(one, ierr=1,m)]

       ! Test 1: Cardioid boundary
       ierr = polar%new_fit(x, y, z, cardioid, w, smoothing=0.1_FP_REAL)
       if (.not.FITPACK_SUCCESS(ierr)) then
           success = .false.
           write(useUnit, 1000) 'cardioid', FITPACK_MESSAGE(ierr)
       else
           f = polar%eval(x, y, ierr)
           if (FITPACK_SUCCESS(ierr)) then
               avg   = sum(abs(f - z)) / m
               ermax = maxval(abs(f - z))
               write(useUnit, 920) 'cardioid', polar%smoothing, avg, ermax
           end if
       end if

       ! Generate test points inside a 3-petal rose: r = |cos(3*theta)|
       theta = linspace(zero, pi, m)
       r     = abs(cos(three * theta)) * linspace(0.3_FP_REAL, 0.9_FP_REAL, m)
       x     = r * cos(theta)
       y     = r * sin(theta)
       z     = test_func(x, y)

       ! Test 2: Rose curve boundary (use higher smoothing to avoid too many knots)
       ierr = polar%new_fit(x, y, z, rose_3petal, w, smoothing=1.0_FP_REAL)
       if (.not.FITPACK_SUCCESS(ierr)) then
           success = .false.
           write(useUnit, 1000) 'rose', FITPACK_MESSAGE(ierr)
       else
           f = polar%eval(x, y, ierr)
           if (FITPACK_SUCCESS(ierr)) then
               avg   = sum(abs(f - z)) / m
               ermax = maxval(abs(f - z))
               write(useUnit, 920) 'rose', polar%smoothing, avg, ermax
           end if
       end if

        920 format('[test_polar_boundary] ',a,': s=',f7.3,', mean error = ',f7.4,5x,'max. error = ',f7.4)
       1000 format('[test_polar_boundary] ',a,' domain failed: ',a)

    contains

       ! Test function
       elemental real(FP_REAL) function test_func(x, y)
           real(FP_REAL), intent(in) :: x, y
           test_func = sin(pi * x) * cos(pi * y)
       end function test_func

       ! Cardioid: r = 0.5*(1 + cos(theta))
       pure real(FP_REAL) function cardioid(theta)
           real(FP_REAL), intent(in) :: theta
           cardioid = half * (one + cos(theta))
       end function cardioid

       ! 3-petal rose: r = |cos(3*theta)|
       pure real(FP_REAL) function rose_3petal(theta)
           real(FP_REAL), intent(in) :: theta
           rose_3petal = max(abs(cos(three * theta)), 0.1_FP_REAL)
       end function rose_3petal

    end function test_polar_boundary

    logical function test_sphere_fit(iunit) result(success)
       integer, optional, intent(in) :: iunit

       type(fitpack_sphere) :: sphere
       integer :: useUnit,loop,ierr,i,j,m
       real(FP_REAL), allocatable, dimension(:) :: theta,phi,r,w,teval,phieval
       real(FP_REAL), allocatable, dimension(:,:) :: exact,f
       real(FP_REAL) :: avg,ermax
       character(64) :: loop_name

       ! Initialization.
       success = .true.
       if (present(iunit)) then
           useUnit = iunit
       else
           useUnit = output_unit
       end if

       ! Get latitude and longitude coordinates from array "dasphe"
       theta = min(deg2rad*dasphe(1:384-1:2),pi)
       phi   = min(deg2rad*dasphe(2:384  :2),pi2)
       r     = testsp(theta,phi) ! calculate the function values.

       ! Set up an evaluation grid at fixed observation points
       teval   = [((pi/32)*i,i=0,32)]
       phieval = [((pi2/32)*i,i=0,32)]
       allocate(exact(33,33)); forall(i=1:33,j=1:33) exact(j,i) = testsp(theta(i),phi(j))

       call write_grid_to_file('exact_grid',teval,phieval,exact)

       ! Number of points
       m = size(theta)

       ! Set uniform weights: w(i)=(0.01)**(-1), where 0.01 is an estimate for the standard deviation
       ! of the error in z(i)).
       allocate(w(m),source=1.0_FP_REAL)

       ! we determine a number of smoothing spline approximations on the unit disk: x**2+y**2 <= 1.
       approximations: do loop=1,3

           ! Compute/update fit
           select case (loop)
              case (1)
                ierr = sphere%new_fit(theta,phi,r,smoothing=500.0_FP_REAL)
              case (2)
                ierr = sphere%fit(smoothing=135.0_FP_REAL)
              case (3)
                ierr = sphere%fit(smoothing=15.0_FP_REAL)
           end select

           if (.not.FITPACK_SUCCESS(ierr)) then
               success = .false.
               write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
               exit
           end if

           ! Evaluate fit at the initial points
           f = sphere%eval(teval,phieval,ierr)

           write(loop_name,10) loop
           call write_grid_to_file(trim(loop_name),teval,phieval,f)

           if (.not.FITPACK_SUCCESS(ierr)) then
               success = .false.
               write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
               exit
           end if

           ! Determine mean and maximum errors
           avg   = sum(abs(f-exact))/size(f)
           ermax = maxval(abs(f-exact))

           !write(useUnit,920) sphere%smoothing,avg,ermax

       end do approximations

         10 format('loop_',i0)
        920 format('[test_polar_fit] s=',f7.1,', mean error = ',f7.4,5x,'max. error = ',f7.4)
       1000 format('[test_polar_fit] polar test ',i0,' failed: ',a)

       contains

          ! Write sphere data to file
          subroutine write_grid_to_file(fileName,theta,phi,data)
              character(*), intent(in) :: fileName
              real(FP_REAL), intent(in) :: theta(:),phi(:),data(:,:) ! [size(phi),size(theta)]

              integer :: iunit,i,j

              open(newunit=iunit,file=fileName,form='formatted',action='write')

              write(iunit,1) zero,theta

              do i=1,size(phi)
                 write(iunit,1) phi(i),(data(i,j),j=1,size(data,2))
              end do

              close(iunit)

              1 format(*(1pe13.5e3))

          end subroutine write_grid_to_file

          ! calculate the value of a test function for the sphere package.
          elemental real(FP_REAL) function testsp(v,u)
              real(FP_REAL), intent(in) :: u,v
              real(FP_REAL) :: cu,cv,rad1,rad2,rad3,su,sv
              cu = cos(u)
              cv = cos(v)
              su = sin(u)
              sv = sin(v)
              rad1 = (cu*sv*0.2_FP_REAL)**2+(su*sv)**2+(cv*0.5_FP_REAL)**2
              rad2 = (cu*sv)**2+(su*sv*0.5_FP_REAL)**2+(cv*0.2_FP_REAL)**2
              rad3 = (cu*sv*0.5_FP_REAL)**2+(su*sv*0.2_FP_REAL)**2+cv**2
              testsp = one/sqrt(rad1) + one/sqrt(rad2) + one/sqrt(rad3)
              return
          end function testsp

    end function test_sphere_fit

    !> Test surface fit from gridded data
    logical function test_gridded_fit(iunit) result(success)
        integer, optional, intent(in) :: iunit
        type(fitpack_grid_surface) :: surf

        integer :: ierr,loop,useUnit
        real(FP_REAL), allocatable :: fit_z(:,:)
        real(FP_REAL), parameter :: RTOL = 1.0e-4_FP_REAL
        real(FP_REAL), parameter :: ATOL = 1.0e-4_FP_REAL

        success = .true.
        if (present(iunit)) then
            useUnit = iunit
        else
            useUnit = output_unit
        end if

        !> Use tabulated data
        do loop = 1,5

            select case (loop)
               case (1)
                  ! we start computing the least-squares bicubic polynomial
                  ierr = surf%new_fit(daregr_x,daregr_y,daregr_z,smoothing=10.0_FP_REAL)
               case (2)
                  ierr = surf%fit(smoothing=0.22_FP_REAL)
               case (3)
                  ! Overfitting (s is too small
                  ierr = surf%fit(smoothing=0.1_FP_REAL)
               case (4)
                  ! Interpolating spline
                  ierr = surf%interpolate()
               case (5)
                  ! Quintic spline
                  ierr = surf%fit(smoothing=0.2_FP_REAL,order=5)
               case (6)
                  ! Finally, calculate a least-squares spline approximation with specified knots
!                  iopt = -1
!                  kx = 3
!                  ky = 3
!                  nx = 11
!                  ny = 11
!                  j  = kx+2
!
!                  tx(kx+2:kx+4) = [(half*(i-2),i=1,3)]
!                  ty(kx+2:kx+4) = tx(kx+2:kx+4)

            end select

            if (.not.FITPACK_SUCCESS(ierr)) then
                success = .false.
                write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
                exit
            end if

            ! Evaluation of the spline approximation
            fit_z = surf%eval(daregr_x,daregr_y,ierr)

            if (.not.FITPACK_SUCCESS(ierr)) then
                success = .false.
                write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
                exit
            end if

            if (loop==4) then
                ! Check that the evaluations match the original data
                success = all(abs(daregr_z-fit_z)*rewt(RTOL,ATOL,daregr_z)<one)
                if (.not.success) then
                    write(useUnit,1100) loop
                    print 900, 'daregr=',daregr_z
                    print 900, 'fit   =',fit_z
                    exit
                end if
            endif

        end do
        
         900 format(a,*(1x,1pe12.3))
        1000 format('[test_gridded_surface] test ',i0,' failed: ',a)
        1100 format('[test_gridded_surface] test ',i0,' evaluation failed ')


    end function test_gridded_fit


    !> Test surface fit from gridded data
    logical function test_gridded_polar(iunit) result(success)
        integer, optional, intent(in) :: iunit
        type(fitpack_grid_polar) :: polr

        integer :: ierr,loop,useUnit,i
        real(FP_REAL) :: z0,r
        real(FP_REAL), allocatable :: u(:),v(:),z(:,:),fit_z(:,:),exact_z(:,:),err(:,:)

        success = .true.
        if (present(iunit)) then
            useUnit = iunit
        else
            useUnit = output_unit
        end if

        !> Setup grid
        u = linspace(0.1_FP_REAL,0.9_FP_REAL,9)      ! u in [0,r]
        v = linspace(-pi,0.9_FP_REAL*pi,19)  ! v in [-pi,pi]
        r = one

        !> Get tabulated data from external datafile
        z  = reshape(dapogr(1:9*19),[19,9])
        z0 = dapogr(9*19+1)

        !> Get analytical solution on the grid points
        exact_z = polar_fun(spread(u,1,19),spread(v,2,9))

        !> Use tabulated data
        do loop = 1,6

            select case (loop)
               case (1)
                  ! Initialize domain
                  call polr%new_points(u,v,r,z)
                  call polr%write('polar_grid.dat')

                  ! Set c0-continuity at the origin
                  call polr%set_origin_BC(z0,exact=.false.,differentiable=.false.)
                  ! Non-vanishing at the boundary
                  polr%bc_boundary = OUTSIDE_EXTRAPOLATE

                  ! Large value of s for computing least-squares polynomial
                  ierr = polr%fit(smoothing=5.0_FP_REAL)
               case (2)
                  ierr = polr%fit(smoothing=0.1_FP_REAL)
               case (3)
                  ierr = polr%interpolate()
               case (4)

                  ! Change BCs
                  polr%bc_boundary = OUTSIDE_ZERO
                  call polr%set_origin_BC(polar_fun(zero,zero),exact=.true.,differentiable=.true.)
                  ierr = polr%fit(smoothing=0.1_FP_REAL)

               case (5)

                  call polr%set_origin_BC(differentiable=.true.)

                  ierr = polr%fit(smoothing=0.1_FP_REAL)

               case (6)

                  ierr = polr%least_squares()

            end select

            if (.not.FITPACK_SUCCESS(ierr)) then
                success = .false.
                write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
                exit
            end if

            ! Evaluation of the spline approximation
            fit_z = polr%eval(u,v,ierr)

            if (.not.FITPACK_SUCCESS(ierr)) then
                success = .false.
                write(useUnit,1100) loop,FITPACK_MESSAGE(ierr)
                exit
            end if

            ! Check error vs analytical solution
            err = abs(exact_z-fit_z)
            write(useUnit,800) polr%c(1),abs(polr%c(1)-polar_fun(zero,zero))
            write(useUnit,935) sum(err)/size(err),maxval(err)

        end do

         800 format('[test_gridded_polar] spline value at (0,0) = ',f7.3,5x,'error = ',f7.3)

         935 format('[test_gridded_polar] error: mean = ',f9.3,'     max. abs. = ',f9.3)

        1000 format('[test_gridded_polar] test ',i0,' failed: ',a)
        1100 format('[test_gridded_polar] test ',i0,' evaluation failed ')

        contains

        ! function program tespog calculates the value of the test function underlying the data.
        elemental real(FP_REAL) function polar_fun(u,v)
           real(FP_REAL), intent(in) :: u,v
           real(FP_REAL) :: x,y,f
           x = u*cos(v)
           y = u*sin(v)
           f = one-((3*x-one)**2+(3*y-one)**2)/(11.-6*(x+y))
           polar_fun = f-(one-x**2-y**2)*(x+y)*54.0_FP_REAL/121.0_FP_REAL
        end function polar_fun

    end function test_gridded_polar

    !> Test surface fit from gridded sphere data
    logical function test_gridded_sphere(iunit) result(success)
        integer, optional, intent(in) :: iunit
        type(fitpack_grid_sphere) :: spgr

        integer :: ierr,loop,useUnit
        real(FP_REAL), allocatable :: u(:),v(:),z(:,:),fit_z(:,:),exact_z(:,:),err(:,:)

        success = .true.
        if (present(iunit)) then
            useUnit = iunit
        else
            useUnit = output_unit
        end if

        !> Load grid and data
        u = daspgr_u
        v = daspgr_v
        z = reshape(daspgr_r,[size(v),size(u)])

        !> Get analytical solution on the grid points
        exact_z = tesspg(spread(u,1,size(v)),spread(v,2,size(u)))

        ! Check error of input data vs analytical solution
        err = abs(z-exact_z)
        write(useUnit,935) sum(err)/size(err),maxval(err)

        !> Use tabulated data
        do loop = 1,6

            select case (loop)
               case (1)

                  ! Initialize data
                  call spgr%new_points(u,v,z)

                  ! Set c0-continuity at the poles, no data values
                  call spgr%BC_south_pole(differentiable=.false.)
                  call spgr%BC_north_pole(differentiable=.false.)

                  ! A large value for s for computing the least-squares polynomial
                  ierr = spgr%fit(smoothing=60.0_FP_REAL)

               case (2)

                  ierr = spgr%fit(smoothing=0.05_FP_REAL)
               case (3)
                  ierr = spgr%interpolate()
               case (4)

                  ! a second set of approximations with c1-continuity at the poles and the exact values
                  call spgr%BC_north_pole(z0=tesspg(zero,zero),exact=.true.,differentiable=.true.,zero_grad=.false.)
                  call spgr%BC_south_pole(z0=tesspg(pi,zero),exact=.true.,differentiable=.true.,zero_grad=.false.)
                  ierr = spgr%fit(smoothing=0.05_FP_REAL)

               case (5)

                  ! Vanishing derivatives at the poles
                  call spgr%BC_north_pole(z0=tesspg(zero,zero),exact=.true.,differentiable=.true.,zero_grad=.true.)
                  call spgr%BC_south_pole(z0=tesspg(pi,zero),exact=.true.,differentiable=.true.,zero_grad=.true.)
                  ierr = spgr%fit(smoothing=0.05_FP_REAL)

               case (6)

                  ! finally we calculate the least-squares spline according to the current set of knots
                  ierr = spgr%least_squares()

            end select

            if (.not.FITPACK_SUCCESS(ierr)) then
                success = .false.
                write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
                exit
            end if

            ! Evaluation of the spline approximation
            fit_z = spgr%eval(u,v,ierr)

            if (.not.FITPACK_SUCCESS(ierr)) then
                success = .false.
                write(useUnit,1100) loop,FITPACK_MESSAGE(ierr)
                exit
            end if

            ! Check error vs analytical solution
            err = abs(exact_z-fit_z)
            write(useUnit,800) spgr%c(1),spgr%c(size(spgr%c))
            write(useUnit,800) tesspg(zero,zero),tesspg(pi,zero)
            write(useUnit,935) sum(err)/size(err),maxval(err)

        end do

         800 format('[test_gridded_sphere] spline value at the poles =',2(1x,f7.3,5x))

         935 format('[test_gridded_sphere] error: mean = ',f9.3,'     max. abs. = ',f9.3)

        1000 format('[test_gridded_sphere] test ',i0,' failed: ',a)
        1100 format('[test_gridded_sphere] test ',i0,' evaluation failed ')

        contains

          ! function program tesspg calculates the value of the test function
          ! underlying the data.
          elemental real(FP_REAL) function tesspg(u,v)
              real(FP_REAL), intent(in) :: u,v
              tesspg = 2.0_FP_REAL/(4.1_FP_REAL+cos(3*u)+3*cos(v+v+u*0.25_FP_REAL)*sin(u)**2)
              return
          end function tesspg

    end function test_gridded_sphere

    logical function test_parametric_surface(iunit) result(success)
        integer, optional, intent(in) :: iunit

        integer, parameter :: mu   = 21
        integer, parameter :: mv   = 11
        integer, parameter :: idim = 3
        integer :: useUnit,i,loop,ier
        type(fitpack_parametric_surface) :: surf
        real(FP_REAL), allocatable :: u(:),v(:),f1(:,:,:),f(:,:,:)


        success = .true.
        if (present(iunit)) then
            useUnit = iunit
        else
            useUnit = output_unit
        end if

        !  we generate the u-coordinates of the grid and load surface data
        u = [(i-1,i=1,mu)]
        v = [(u(2*i-1),i=1,mv)]

        ! dapasu is provided as [v,idim,u], f is requested as [v,u,idim] -> swap dimensions
        f  = reshape(dapasu,[mv,mu,idim],order=[1,3,2])

        !  main loop for the different spline approximations
        approximations: do loop=1,4
           select case (loop)
              case (1)

                 ! Smoothing surface, no periodicity constraint
                 ier = surf%new_fit(u,v,f,smoothing=0.07_FP_REAL,periodic_BC=[.false.,.false.])

              case (2)

                 ! A smoothign surface periodic in the v-variable
                 ier = surf%fit(periodic=[.false.,.true.])

              case (3)

                 ! Smoothing surface with both periodic variables
                 ier = surf%fit(periodic=[.true.,.true.])

              case (4)

                 ! Least-squares with specified knots
                 ier = surf%least_squares(u_knots=[(5.0_FP_REAL*i,i=1,3)],&
                                          v_knots=[(5.0_FP_REAL*i,i=1,3)])

           end select

           if (.not.FITPACK_SUCCESS(ier)) then
               success = .false.
               write(useUnit,1000) loop,FITPACK_MESSAGE(ier)
           end if

           !  printing of the fitting results.
           if (surf%iopt>=0) then
               write(useUnit,940) merge(1,0,surf%periodic_dim)
               write(useUnit,945) surf%smoothing
           else
               write(useUnit,935) merge(1,0,surf%periodic_dim)
           endif

           ! Evaluate the spline surface at the same nodes
           f1 = surf%eval(u,v,ier)

           if (.not.FITPACK_SUCCESS(ier)) then
               success = .false.
               write(useUnit,1000) loop,FITPACK_MESSAGE(ier)
               exit
           end if

        end do approximations

        935  format('0least-squares surface of periodicity',2i3)
        940  format('0smoothing surface of periodicity',2i3)
        945  format(' smoothing factor s=',f8.2)
       1000  format('[test_parametric_surface] case ',i0,' failed with message: ',a)

    end function test_parametric_surface

    !> Test that pack/expand round-trip preserves curve data
    logical function test_curve_comm_roundtrip() result(success)

        integer, parameter :: N = 50
        type(fitpack_curve) :: curve1, curve2
        real(FP_REAL) :: x(N), y(N), xtest(10), y1(10), y2(10)
        real(FP_COMM), allocatable :: buffer(:)
        integer(FP_SIZE) :: buf_size
        integer :: ierr, i

        success = .false.

        ! Generate test data: a simple polynomial
        x = linspace(zero, pi2, N)
        y = sin(x) + half * cos(2 * x)

        ! Create an interpolating curve
        ierr = curve1%new_fit(x, y, smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_curve_comm_roundtrip] error creating curve: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Get buffer size and allocate
        buf_size = curve1%comm_size()
        allocate(buffer(buf_size))

        ! Pack curve into buffer
        call curve1%comm_pack(buffer)

        ! Expand buffer into new curve
        call curve2%comm_expand(buffer)

        ! Generate test points (different from fitting points)
        xtest = linspace(0.1_FP_REAL, pi2 - 0.1_FP_REAL, 10)

        ! Evaluate both curves
        y1 = curve1%eval(xtest)
        y2 = curve2%eval(xtest)

        ! Check that evaluations match
        if (maxval(abs(y1 - y2)) > epsilon(one)) then
            print *, '[test_curve_comm_roundtrip] evaluation mismatch after round-trip'
            print *, '  max difference: ', maxval(abs(y1 - y2))
            return
        end if

        ! Check scalar members
        if (curve1%m /= curve2%m) then
            print *, '[test_curve_comm_roundtrip] m mismatch: ', curve1%m, ' vs ', curve2%m
            return
        end if
        if (curve1%order /= curve2%order) then
            print *, '[test_curve_comm_roundtrip] order mismatch'
            return
        end if
        if (curve1%knots /= curve2%knots) then
            print *, '[test_curve_comm_roundtrip] knots mismatch'
            return
        end if
        if (abs(curve1%smoothing - curve2%smoothing) > epsilon(one)) then
            print *, '[test_curve_comm_roundtrip] smoothing mismatch'
            return
        end if
        if (abs(curve1%fp - curve2%fp) > epsilon(one)) then
            print *, '[test_curve_comm_roundtrip] fp mismatch'
            return
        end if

        ! Check array sizes match
        if (size(curve1%t) /= size(curve2%t)) then
            print *, '[test_curve_comm_roundtrip] t array size mismatch'
            return
        end if
        if (size(curve1%c) /= size(curve2%c)) then
            print *, '[test_curve_comm_roundtrip] c array size mismatch'
            return
        end if

        ! Check knot and coefficient values
        if (maxval(abs(curve1%t - curve2%t)) > epsilon(one)) then
            print *, '[test_curve_comm_roundtrip] t values mismatch'
            return
        end if
        if (maxval(abs(curve1%c - curve2%c)) > epsilon(one)) then
            print *, '[test_curve_comm_roundtrip] c values mismatch'
            return
        end if

        success = .true.

    end function test_curve_comm_roundtrip

    !> Test pack/expand round-trip for all fitter types.
    !> Uses the double-roundtrip pattern: pack obj1 -> expand into obj2 -> pack obj2 -> compare buffers.
    logical function test_comm_roundtrips() result(success)

        character(*), parameter :: PREFIX = '[test_comm_roundtrips] '
        integer(FP_FLAG) :: ierr
        integer :: i, j

        success = .true.

        ! --- fitpack_surface (scattered) ---
        scattered_surface: block
            integer, parameter :: NX = 5, NY = 5, NS = NX*NY
            type(fitpack_surface) :: s1, s2
            real(FP_REAL) :: xs(NS), ys(NS), zs(NS), xg(NX), yg(NY)

            xg = linspace(-one, one, NX)
            yg = linspace(-one, one, NY)
            call meshgrid(xg, yg, xs, ys)
            zs = xs**2 + ys**2

            ierr = s1%new_fit(xs, ys, zs)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'surface fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(s1, s2, 'fitpack_surface')) then
                success = .false.; return
            end if
        end block scattered_surface

        ! --- fitpack_grid_surface ---
        grid_surface: block
            type(fitpack_grid_surface) :: g1, g2

            ierr = g1%new_fit(daregr_x, daregr_y, daregr_z, smoothing=0.5_FP_REAL)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'grid surface fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(g1, g2, 'fitpack_grid_surface')) then
                success = .false.; return
            end if
        end block grid_surface

        ! --- fitpack_polar (scattered) ---
        scattered_polar: block
            type(fitpack_polar) :: p1, p2
            real(FP_REAL), allocatable :: x(:), y(:), z(:), w(:)

            x = dapola(1:600-2:3)
            y = dapola(2:600-1:3)
            z = dapola(3:600  :3)
            allocate(w(size(x)), source=100.0_FP_REAL)

            ierr = p1%new_fit(x, y, z, unit_disk_boundary, w)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'polar fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(p1, p2, 'fitpack_polar')) then
                success = .false.; return
            end if
        end block scattered_polar

        ! --- fitpack_grid_polar ---
        gridded_polar: block
            type(fitpack_grid_polar) :: p1, p2
            real(FP_REAL), allocatable :: u(:), v(:), z(:,:)
            real(FP_REAL) :: z0, r

            u = linspace(0.1_FP_REAL, 0.9_FP_REAL, 9)
            v = linspace(-pi, 0.9_FP_REAL*pi, 19)
            r = one
            z  = reshape(dapogr(1:9*19), [19,9])
            z0 = dapogr(9*19+1)

            call p1%new_points(u, v, r, z, z0)
            ierr = p1%fit(smoothing=0.05_FP_REAL)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'grid polar fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(p1, p2, 'fitpack_grid_polar')) then
                success = .false.; return
            end if
        end block gridded_polar

        ! --- fitpack_sphere (scattered) ---
        scattered_sphere: block
            type(fitpack_sphere) :: s1, s2
            real(FP_REAL), allocatable :: theta(:), phi(:), r(:)

            theta = min(deg2rad*dasphe(1:384-1:2), pi)
            phi   = min(deg2rad*dasphe(2:384  :2), pi2)
            r     = theta + phi ! simple function of coordinates

            ierr = s1%new_fit(theta, phi, r)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'sphere fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(s1, s2, 'fitpack_sphere')) then
                success = .false.; return
            end if
        end block scattered_sphere

        ! --- fitpack_grid_sphere ---
        gridded_sphere: block
            type(fitpack_grid_sphere) :: s1, s2
            real(FP_REAL), allocatable :: z(:,:)

            z = reshape(daspgr_r, [size(daspgr_v), size(daspgr_u)])

            ierr = s1%new_fit(daspgr_u, daspgr_v, z, smoothing=0.05_FP_REAL)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'grid sphere fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(s1, s2, 'fitpack_grid_sphere')) then
                success = .false.; return
            end if
        end block gridded_sphere

        ! --- fitpack_parametric_curve ---
        parametric_curve: block
            type(fitpack_parametric_curve) :: c1, c2
            real(FP_REAL), allocatable :: x(:,:), w(:), u(:)

            u = [real(FP_REAL) :: 120,128,133,136,138,141,144,146,149,151,154,161,170,180,190,&
                                200,210,220,230,240,250,262,269,273,278,282,287,291,295,299,305,315]
            allocate(x(2,32))
            x(1,:) = [-1.5141,-2.0906,-1.9253,-0.8724,-0.3074,-0.5534,0.0192,1.2298,2.5479,2.4710,1.7063,&
                       1.1183,0.5534,0.4727,0.3574,0.1998,0.2882,0.2613,0.2652,0.2805,0.4112,0.9377,1.3527,&
                       1.5564,1.6141,1.6333,1.1567,0.8109,0.2498,-0.2306,-0.7571,-1.1222]
            x(2,:) = [0.5150,1.3412,2.6094,3.2358,2.7401,2.7823,3.5932,3.8353,2.5863,1.3105,0.6841,0.2575,&
                       0.2460,0.3689,0.2460,0.2998,0.3651,0.3343,0.3881,0.4573,0.5918,0.7110,0.4035,0.0769,&
                      -0.3920,-0.8570,-1.3412,-1.5641,-1.7409,-1.7178,-1.2989,-0.5572]
            allocate(w(32), source=one)

            ierr = c1%new_fit(x, w=w, u=u, smoothing=0.2_FP_REAL)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'parametric curve fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(c1, c2, 'fitpack_parametric_curve')) then
                success = .false.; return
            end if
        end block parametric_curve

        ! --- fitpack_closed_curve ---
        closed_curve: block
            type(fitpack_closed_curve) :: c1, c2
            real(FP_REAL) :: x(2,19), w(19), u(19)

            x(1,1:18) = [-4.7,-7.048,-6.894,-3.75,-1.042,0.938,2.5,3.524,4.511,5.0,4.886,3.524,3.2,&
                          1.302,-1.424,-3.0,-3.064,-3.665]
            x(2,1:18) = [0.0,2.565,5.785,6.495,5.909,5.318,4.33,2.957,1.642,0.0,-1.779,-2.957,-5.543,&
                         -7.386,-8.075,-5.196,-2.571,-1.334]
            x(:,19) = x(:,1)
            w = one
            u = [(20*i, i=0,18)]

            ierr = c1%new_fit(x, w=w, u=u, smoothing=100.0_FP_REAL)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'closed curve fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(c1, c2, 'fitpack_closed_curve')) then
                success = .false.; return
            end if
        end block closed_curve

        ! --- fitpack_constrained_curve ---
        constrained_curve: block
            type(fitpack_constrained_curve) :: c1, c2
            real(FP_REAL) :: x(2,31)
            real(FP_REAL), allocatable :: w(:)

            x(1,:) = [-3.109,-2.188,-1.351,-0.605,0.093,0.451,0.652,0.701,0.518,0.277,0.008,-0.291,&
                      -0.562,-0.679,-0.637,-0.425,-0.049,0.575,1.334,2.167,3.206,4.099,4.872,5.710,&
                       6.330,6.741,6.928,6.965,6.842,6.593,6.269]
            x(2,:) = [3.040,2.876,2.634,2.183,1.586,1.010,0.382,-0.218,-0.632,-0.879,-0.981,-0.886,&
                      -0.642,-0.195,0.373,1.070,1.607,2.165,2.618,2.905,2.991,2.897,2.615,2.164,&
                       1.617,0.977,0.383,-0.194,-0.665,-0.901,-1.010]
            allocate(w(31), source=one)

            ierr = c1%new_fit(x, w=w, smoothing=60.0_FP_REAL)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'constrained curve fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(c1, c2, 'fitpack_constrained_curve')) then
                success = .false.; return
            end if
        end block constrained_curve

        ! --- fitpack_parametric_surface ---
        parametric_surface: block
            integer, parameter :: mu = 21, mv = 11, idim = 3
            type(fitpack_parametric_surface) :: s1, s2
            real(FP_REAL), allocatable :: u(:), v(:), f(:,:,:)

            u = [(i-1, i=1,mu)]
            v = [(u(2*i-1), i=1,mv)]
            f = reshape(dapasu, [mv,mu,idim], order=[1,3,2])

            ierr = s1%new_fit(u, v, f)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, PREFIX//'parametric surface fit failed: ', FITPACK_MESSAGE(ierr)
                success = .false.; return
            end if

            if (.not.check_roundtrip(s1, s2, 'fitpack_parametric_surface')) then
                success = .false.; return
            end if
        end block parametric_surface

    contains

        ! Double roundtrip helper: pack obj1 -> expand into obj2 -> pack obj2 -> compare buffers
        logical function check_roundtrip(obj1, obj2, name) result(ok)
            class(fitpack_fitter), intent(in)    :: obj1
            class(fitpack_fitter), intent(inout) :: obj2
            character(*), intent(in) :: name

            real(FP_COMM), allocatable :: buf1(:), buf2(:)
            integer(FP_SIZE) :: n1, n2

            ok = .false.

            ! Pack original
            n1 = obj1%comm_size()
            allocate(buf1(n1))
            call obj1%comm_pack(buf1)

            ! Expand into copy
            call obj2%comm_expand(buf1)

            ! Re-pack the copy
            n2 = obj2%comm_size()
            allocate(buf2(n2))
            call obj2%comm_pack(buf2)

            ! Check buffer sizes
            if (n1 /= n2) then
                print *, PREFIX//name//': buffer size mismatch: ', n1, ' vs ', n2
                return
            end if

            ! Check buffer contents (use epsilon tolerance for workspace padding)
            if (maxval(abs(buf1 - buf2)) > epsilon(1.0_FP_COMM)) then
                print *, PREFIX//name//': buffer contents mismatch, max diff = ', maxval(abs(buf1-buf2))
                return
            end if

            ! Check base fields
            if (abs(obj1%smoothing - obj2%smoothing) > epsilon(one)) then
                print *, PREFIX//name//': smoothing mismatch'
                return
            end if
            if (abs(obj1%fp - obj2%fp) > epsilon(one)) then
                print *, PREFIX//name//': fp mismatch'
                return
            end if
            if (obj1%iopt /= obj2%iopt) then
                print *, PREFIX//name//': iopt mismatch'
                return
            end if
            if (size(obj1%c) /= size(obj2%c)) then
                print *, PREFIX//name//': c array size mismatch'
                return
            end if
            if (maxval(abs(obj1%c - obj2%c)) > epsilon(one)) then
                print *, PREFIX//name//': c values mismatch'
                return
            end if

            ok = .true.

        end function check_roundtrip

        ! Unit disk boundary for polar test
        pure real(FP_REAL) function unit_disk_boundary(v)
            real(FP_REAL), intent(in) :: v
            unit_disk_boundary = one
        end function unit_disk_boundary

    end function test_comm_roundtrips

    ! Test that `use fitpack` alone provides access to constants, error flags,
    ! and new methods (eval, mse, interpolate, least_squares on surface/curve)
    logical function test_umbrella_api() result(success)

       type(fitpack_curve) :: curve
       type(fitpack_parametric_curve) :: pcurve
       type(fitpack_surface) :: surf
       type(fitpack_grid_surface) :: gsurf
       integer(FP_FLAG) :: ierr
       real(FP_REAL) :: mse_val
       real(FP_REAL), allocatable :: fscatter(:),fgrid(:,:)

       ! Curve test data
       integer, parameter :: NC = 50
       real(FP_REAL) :: xc(NC),yc(NC)
       real(FP_REAL) :: pcxy(2,NC)

       ! Surface test data: 5x5 scattered grid
       integer, parameter :: NX = 5, NY = 5, NS = NX*NY
       real(FP_REAL) :: xs(NS),ys(NS),zs(NS),xg(NX),yg(NY),zg(NY,NX)

       integer :: i,j

       success = .false.

       ! ---- Verify constants are accessible from `use fitpack` ----
       ! (These would not compile if not re-exported)
       if (FITPACK_OK /= 0) return
       if (IOPT_NEW_SMOOTHING /= 0) return
       if (OUTSIDE_NEAREST_BND /= 3) return
       if (abs(pi - atan2(zero,-one)) > epsilon(one)) return

       ! ---- Test curve least_squares ----
       xc = linspace(zero,pi,NC)
       yc = sin(xc)
       ierr = curve%new_fit(xc,yc)
       if (.not.FITPACK_SUCCESS(ierr)) return

       mse_val = curve%mse()
       ierr = curve%least_squares()
       if (.not.FITPACK_SUCCESS(ierr)) return

       ! ---- Test parametric curve interpolate/least_squares ----
       pcxy(1,:) = cos(xc)
       pcxy(2,:) = sin(xc)

       ! Smoothing fit, then interpolation: fp must be zero
       ierr = pcurve%new_fit(pcxy)
       if (.not.FITPACK_SUCCESS(ierr)) return
       ierr = pcurve%interpolate()
       if (.not.FITPACK_SUCCESS(ierr)) return
       if (pcurve%mse() > zero) return

       ! Least squares with knot reset: large smoothing -> fewer knots -> fp > 0
       ierr = pcurve%least_squares(smoothing=real(NC,FP_REAL),reset_knots=.true.)
       if (.not.FITPACK_SUCCESS(ierr)) return
       if (.not. pcurve%mse() > zero) return

       ! ---- Test fitpack_surface eval, eval_ongrid, mse ----
       xg = linspace(-one,one,NX)
       yg = linspace(-one,one,NY)

       ! Generate scattered points from the grid
       call meshgrid(xg, yg, xs, ys)
       zs = xs**2 + ys**2

       ierr = surf%new_fit(xs,ys,zs)
       if (.not.FITPACK_SUCCESS(ierr)) return

       ! Evaluate at all scattered points
       fscatter = surf%eval(xs,ys,ierr)
       if (.not.FITPACK_SUCCESS(ierr)) return

       ! Evaluate on a grid
       fgrid = surf%eval_ongrid(xg,yg,ierr)
       if (.not.FITPACK_SUCCESS(ierr)) return
       if (size(fgrid,1)/=NY .or. size(fgrid,2)/=NX) return

       ! MSE
       mse_val = surf%mse()

       ! Interpolate
       ierr = surf%interpolate()
       if (.not.FITPACK_SUCCESS(ierr)) return

       ! Least squares (after interpolation, knots are set)
       ierr = surf%least_squares()
       if (.not.FITPACK_SUCCESS(ierr)) return

       ! ---- Test grid_surface mse ----
       do j = 1,NX
          zg(:,j) = xg(j)**2 + yg**2
       end do

       ierr = gsurf%new_fit(xg,yg,zg)
       if (.not.FITPACK_SUCCESS(ierr)) return

       mse_val = gsurf%mse()

       success = .true.

    end function test_umbrella_api

    !> Test convex curve fitting (concon) using mncoco test data
    logical function test_convex_fit() result(success)

        integer(FP_SIZE), parameter :: m = 16
        type(fitpack_convex_curve) :: curve
        real(FP_REAL) :: x(m),y(m),w(m),v(m),s2(m)
        real(FP_REAL) :: s(3)
        integer(FP_FLAG) :: ierr
        integer(FP_SIZE) :: is

        success = .false.

        ! mncoco test data
        x = [0.1_FP_REAL,0.3_FP_REAL,0.5_FP_REAL,0.7_FP_REAL,0.9_FP_REAL, &
             1.25_FP_REAL,1.75_FP_REAL,2.25_FP_REAL,2.75_FP_REAL,3.5_FP_REAL, &
             4.5_FP_REAL,5.5_FP_REAL,6.5_FP_REAL,7.5_FP_REAL,8.5_FP_REAL,9.5_FP_REAL]
        y = [0.124_FP_REAL,0.234_FP_REAL,0.256_FP_REAL,0.277_FP_REAL,0.278_FP_REAL, &
             0.291_FP_REAL,0.308_FP_REAL,0.311_FP_REAL,0.315_FP_REAL,0.322_FP_REAL, &
             0.317_FP_REAL,0.326_FP_REAL,0.323_FP_REAL,0.321_FP_REAL,0.322_FP_REAL, &
             0.328_FP_REAL]

        w = one
        w(1) = 10.0_FP_REAL
        w(2) = 3.0_FP_REAL
        w(16) = 10.0_FP_REAL

        ! Concave constraints
        v = one

        s = [0.2_FP_REAL, 0.04_FP_REAL, 0.0002_FP_REAL]

        ! First fit
        call curve%new_points(x,y,w)
        ierr = curve%set_convexity(v)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_convex_fit] set_convexity failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        do is=1,3

            ierr = curve%fit(smoothing=s(is),keep_knots=(is>1))

            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, '[test_convex_fit] fit failed at s=',s(is),': ', FITPACK_MESSAGE(ierr)
                return
            end if

            ! Evaluate second derivative at data points
            s2 = curve%dfdx(x,order=2,ierr=ierr)
            if (.not.FITPACK_SUCCESS(ierr)) then
                print *, '[test_convex_fit] derivative failed: ', FITPACK_MESSAGE(ierr)
                return
            end if

            ! Check concavity: s''(x) <= 0 at all data points (with small tolerance)
            if (any(s2 > 1.0e-6_FP_REAL)) then
                print *, '[test_convex_fit] concavity violated at s=',s(is)
                print *, '  max s''(x) = ', maxval(s2)
                return
            end if

        end do

        success = .true.

    end function test_convex_fit

    !> Test convex least-squares fitting (cocosp) using mncosp test data
    logical function test_convex_least_squares() result(success)

        integer(FP_SIZE), parameter :: m = 10
        type(fitpack_convex_curve) :: curve
        real(FP_REAL) :: x(m),y(m),w(m),v_convex(m),v_concave(m),s2(m)
        integer(FP_FLAG) :: ierr

        success = .false.

        ! mncosp test data
        x = [0.25_FP_REAL,0.5_FP_REAL,0.75_FP_REAL,1.25_FP_REAL,1.75_FP_REAL, &
             2.25_FP_REAL,2.75_FP_REAL,3.25_FP_REAL,6.25_FP_REAL,12.25_FP_REAL]
        y = [17.0_FP_REAL,15.2_FP_REAL,13.8_FP_REAL,12.2_FP_REAL,11.0_FP_REAL, &
             10.1_FP_REAL,9.4_FP_REAL,8.6_FP_REAL,6.1_FP_REAL,3.5_FP_REAL]
        w = one

        ! Convex constraints (v=-1 => s'' >= 0)
        v_convex  = -one
        ! Concave constraints (v=1 => s'' <= 0)
        v_concave = one

        ! First: fit with concon to establish knots, then least-squares with cocosp
        call curve%new_points(x,y,w)
        ierr = curve%set_convexity(v_convex)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_convex_least_squares] set_convexity failed'
            return
        end if

        ! Use reset_knots to first fit via concon, then do least-squares via cocosp
        ierr = curve%least_squares(smoothing=1.0e+10_FP_REAL,reset_knots=.true.)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_convex_least_squares] convex LS failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Check convexity: s''(x) >= 0 at data points
        s2 = curve%dfdx(x,order=2,ierr=ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_convex_least_squares] derivative failed'
            return
        end if
        if (any(s2 < -1.0e-6_FP_REAL)) then
            print *, '[test_convex_least_squares] convexity violated'
            print *, '  min s''(x) = ', minval(s2)
            return
        end if

        ! Now test concave
        call curve%new_points(x,y,w)
        ierr = curve%set_convexity(v_concave)
        if (.not.FITPACK_SUCCESS(ierr)) return

        ierr = curve%least_squares(smoothing=1.0e+10_FP_REAL,reset_knots=.true.)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_convex_least_squares] concave LS failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        s2 = curve%dfdx(x,order=2,ierr=ierr)
        if (.not.FITPACK_SUCCESS(ierr)) return
        if (any(s2 > 1.0e-6_FP_REAL)) then
            print *, '[test_convex_least_squares] concavity violated'
            print *, '  max s''(x) = ', maxval(s2)
            return
        end if

        success = .true.

    end function test_convex_least_squares

    !> Test knot insertion for a standard curve
    logical function test_insert_knot() result(success)

        integer, parameter :: N = 50
        type(fitpack_curve) :: curve
        real(FP_REAL) :: x(N), y(N), xtest(20), y_before(20), y_after(20)
        integer(FP_SIZE) :: knots_before
        integer(FP_FLAG) :: ierr
        real(FP_REAL) :: midpt

        success = .false.

        ! Generate sine curve data
        x = linspace(zero, pi2, N)
        y = sin(x)

        ierr = curve%new_fit(x,y,smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_insert_knot] fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Evaluate at test points before insertion
        xtest = linspace(0.1_FP_REAL, pi2 - 0.1_FP_REAL, 20)
        y_before = curve%eval(xtest)
        knots_before = curve%knots

        ! Insert a knot at the midpoint
        midpt = half*(curve%xleft + curve%xright)
        call curve%insert_knot(midpt,ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_insert_knot] insert_knot failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Verify knot count increased by 1
        if (curve%knots /= knots_before + 1) then
            print *, '[test_insert_knot] expected knots=',knots_before+1,' got ',curve%knots
            return
        end if

        ! Evaluate after insertion — values must match to machine precision
        y_after = curve%eval(xtest)
        if (maxval(abs(y_before - y_after)) > 100*epsilon(one)) then
            print *, '[test_insert_knot] values changed after insertion!'
            print *, '  max diff = ', maxval(abs(y_before - y_after))
            return
        end if

        success = .true.

    end function test_insert_knot

    !> Test knot insertion for a periodic curve
    logical function test_insert_knot_periodic() result(success)

        integer, parameter :: N = 200
        type(fitpack_periodic_curve) :: curve
        real(FP_REAL) :: x(N), y(N), xtest(20), y_before(20), y_after(20)
        integer(FP_SIZE) :: knots_before
        integer(FP_FLAG) :: ierr
        real(FP_REAL) :: midpt

        success = .false.

        ! Generate periodic function data (cos+sin is periodic on [0,2pi])
        x = linspace(zero, pi2, N)
        y = cos(x) + sin(2*x)

        ierr = curve%new_fit(x,y,smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_insert_knot_periodic] fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Evaluate at test points before insertion
        xtest = linspace(0.1_FP_REAL, pi2 - 0.1_FP_REAL, 20)
        y_before = curve%eval(xtest)
        knots_before = curve%knots

        ! Insert a knot at the midpoint
        midpt = half*(curve%xleft + curve%xright)
        call curve%insert_knot(midpt,ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_insert_knot_periodic] insert_knot failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Verify knot count increased by 1
        if (curve%knots /= knots_before + 1) then
            print *, '[test_insert_knot_periodic] expected knots=',knots_before+1,' got ',curve%knots
            return
        end if

        ! Evaluate after insertion — values must match to machine precision
        y_after = curve%eval(xtest)
        if (maxval(abs(y_before - y_after)) > 100*epsilon(one)) then
            print *, '[test_insert_knot_periodic] values changed after insertion!'
            print *, '  max diff = ', maxval(abs(y_before - y_after))
            return
        end if

        success = .true.

    end function test_insert_knot_periodic

    !> Test double integration of surface over rectangular domain
    logical function test_surface_integral() result(success)

        type(fitpack_surface) :: surf
        type(fitpack_grid_surface) :: gsurf
        integer(FP_FLAG) :: ierr
        real(FP_REAL) :: val, exact
        real(FP_REAL), parameter :: TOL = 5.0e-2_FP_REAL

        integer, parameter :: NX = 10, NY = 10, NS = NX*NY
        real(FP_REAL) :: xg(NX), yg(NY), zg(NY,NX)
        real(FP_REAL) :: xs(NS), ys(NS), zs(NS)
        integer :: j

        success = .false.

        ! Grid data for z = x*y on [0,1] x [0,1]
        xg = linspace(zero, one, NX)
        yg = linspace(zero, one, NY)
        do j = 1, NX
            zg(:,j) = xg(j) * yg
        end do

        ! Scattered data from the same grid
        call meshgrid(xg, yg, xs, ys)
        zs = xs * ys

        ! Fit scattered surface
        ierr = surf%new_fit(xs, ys, zs, smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_surface_integral] scattered fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Integral over full domain: int_0^1 int_0^1 x*y dx dy = 0.25
        exact = 0.25_FP_REAL
        val = surf%integral([zero,zero], [one,one])
        if (abs(val - exact) > TOL) then
            print *, '[test_surface_integral] scattered full domain: expected', exact, ' got', val
            return
        end if

        ! Integral over sub-domain: int_0^0.5 int_0^0.5 x*y dx dy = 0.015625
        exact = 0.015625_FP_REAL
        val = surf%integral([zero,zero], [half,half])
        if (abs(val - exact) > TOL) then
            print *, '[test_surface_integral] scattered sub-domain: expected', exact, ' got', val
            return
        end if

        ! Fit gridded surface
        ierr = gsurf%new_fit(xg, yg, zg, smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_surface_integral] gridded fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Integral over full domain
        exact = 0.25_FP_REAL
        val = gsurf%integral([zero,zero], [one,one])
        if (abs(val - exact) > TOL) then
            print *, '[test_surface_integral] gridded full domain: expected', exact, ' got', val
            return
        end if

        ! Integral over sub-domain: int_0^0.5 int_0^0.5 x*y dx dy = 0.015625
        exact = 0.015625_FP_REAL
        val = gsurf%integral([zero,zero], [half,half])
        if (abs(val - exact) > TOL) then
            print *, '[test_surface_integral] gridded sub-domain: expected', exact, ' got', val
            return
        end if

        success = .true.

    end function test_surface_integral

    !> Test cross-section extraction from a bivariate surface
    logical function test_cross_section() result(success)

        type(fitpack_surface) :: surf
        type(fitpack_grid_surface) :: gsurf
        type(fitpack_curve) :: cs
        integer(FP_FLAG) :: ierr
        real(FP_REAL) :: yeval(5), feval(5), exact(5)
        real(FP_REAL), parameter :: TOL = 5.0e-2_FP_REAL

        integer, parameter :: NX = 10, NY = 10, NS = NX*NY
        real(FP_REAL) :: xg(NX), yg(NY), zg(NY,NX)
        real(FP_REAL) :: xs(NS), ys(NS), zs(NS)
        integer :: j

        success = .false.

        ! Grid data for z = x*y on [0,1] x [0,1]
        xg = linspace(zero, one, NX)
        yg = linspace(zero, one, NY)
        do j = 1, NX
            zg(:,j) = xg(j) * yg
        end do

        ! Scattered data
        call meshgrid(xg, yg, xs, ys)
        zs = xs * ys

        yeval = linspace(0.1_FP_REAL, 0.9_FP_REAL, 5)

        ! --- Scattered surface ---
        ierr = surf%new_fit(xs, ys, zs, smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_cross_section] scattered fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! y-profile at u=0.5: f(y) = 0.5*y
        cs = surf%cross_section(half, along_y=.true., ierr=ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_cross_section] scattered y-profile failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Validate curve order and knots match the profiled direction
        if (cs%order /= surf%order(2)) then
            print *, '[test_cross_section] scattered y-profile order mismatch: expected', surf%order(2), ' got', cs%order
            return
        end if
        if (cs%knots /= surf%knots(2)) then
            print *, '[test_cross_section] scattered y-profile knots mismatch: expected', surf%knots(2), ' got', cs%knots
            return
        end if

        feval = cs%eval(yeval)
        exact = half * yeval
        if (maxval(abs(feval - exact)) > TOL) then
            print *, '[test_cross_section] scattered y-profile error:', maxval(abs(feval - exact))
            return
        end if

        ! x-profile at u=0.3: g(x) = 0.3*x
        cs = surf%cross_section(0.3_FP_REAL, along_y=.false., ierr=ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_cross_section] scattered x-profile failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Validate curve order and knots match the profiled direction
        if (cs%order /= surf%order(1)) then
            print *, '[test_cross_section] scattered x-profile order mismatch: expected', surf%order(1), ' got', cs%order
            return
        end if
        if (cs%knots /= surf%knots(1)) then
            print *, '[test_cross_section] scattered x-profile knots mismatch: expected', surf%knots(1), ' got', cs%knots
            return
        end if

        feval = cs%eval(yeval)
        exact = 0.3_FP_REAL * yeval
        if (maxval(abs(feval - exact)) > TOL) then
            print *, '[test_cross_section] scattered x-profile error:', maxval(abs(feval - exact))
            return
        end if

        ! --- Gridded surface ---
        ierr = gsurf%new_fit(xg, yg, zg, smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_cross_section] gridded fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! y-profile at u=0.5
        cs = gsurf%cross_section(half, along_y=.true., ierr=ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_cross_section] gridded y-profile failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        if (cs%order /= gsurf%order(2)) then
            print *, '[test_cross_section] gridded y-profile order mismatch: expected', gsurf%order(2), ' got', cs%order
            return
        end if

        feval = cs%eval(yeval)
        exact = half * yeval
        if (maxval(abs(feval - exact)) > TOL) then
            print *, '[test_cross_section] gridded y-profile error:', maxval(abs(feval - exact))
            return
        end if

        ! x-profile at u=0.3: g(x) = 0.3*x
        cs = gsurf%cross_section(0.3_FP_REAL, along_y=.false., ierr=ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_cross_section] gridded x-profile failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        if (cs%order /= gsurf%order(1)) then
            print *, '[test_cross_section] gridded x-profile order mismatch: expected', gsurf%order(1), ' got', cs%order
            return
        end if

        feval = cs%eval(yeval)
        exact = 0.3_FP_REAL * yeval
        if (maxval(abs(feval - exact)) > TOL) then
            print *, '[test_cross_section] gridded x-profile error:', maxval(abs(feval - exact))
            return
        end if

        success = .true.

    end function test_cross_section

    !> Test derivative spline computation from a bivariate surface
    logical function test_derivative_spline() result(success)

        type(fitpack_grid_surface) :: gsurf, dsurf
        type(fitpack_surface) :: surf, dsurf_s
        integer(FP_FLAG) :: ierr
        real(FP_REAL), parameter :: TOL = 5.0e-2_FP_REAL

        integer, parameter :: NX = 10, NY = 10, NS = NX*NY
        real(FP_REAL) :: xg(NX), yg(NY), zg(NY,NX)
        real(FP_REAL) :: xs(NS), ys(NS), zs(NS)
        real(FP_REAL) :: xeval(3), yeval(3), fgrid(3,3)
        integer :: i, j

        success = .false.

        ! Grid data for z = x*y on [0,1] x [0,1]
        xg = linspace(zero, one, NX)
        yg = linspace(zero, one, NY)
        do j = 1, NX
            zg(:,j) = xg(j) * yg
        end do

        ! Evaluation grid
        xeval = linspace(0.2_FP_REAL, 0.8_FP_REAL, 3)
        yeval = linspace(0.2_FP_REAL, 0.8_FP_REAL, 3)

        ! --- Gridded surface ---
        ierr = gsurf%new_fit(xg, yg, zg, smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] gridded fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! ds/dx = y
        dsurf = gsurf%derivative_spline(1_FP_SIZE, 0_FP_SIZE, ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] gridded ds/dx failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Verify reduced order: kx drops by 1, ky unchanged
        if (dsurf%order(1) /= gsurf%order(1) - 1) then
            print *, '[test_derivative_spline] gridded ds/dx order(1) mismatch: expected', &
                     gsurf%order(1)-1, ' got', dsurf%order(1)
            return
        end if
        if (dsurf%order(2) /= gsurf%order(2)) then
            print *, '[test_derivative_spline] gridded ds/dx order(2) mismatch: expected', &
                     gsurf%order(2), ' got', dsurf%order(2)
            return
        end if

        ! Verify reduced knots: nx drops by 2, ny unchanged
        if (dsurf%knots(1) /= gsurf%knots(1) - 2) then
            print *, '[test_derivative_spline] gridded ds/dx knots(1) mismatch: expected', &
                     gsurf%knots(1)-2, ' got', dsurf%knots(1)
            return
        end if
        if (dsurf%knots(2) /= gsurf%knots(2)) then
            print *, '[test_derivative_spline] gridded ds/dx knots(2) mismatch: expected', &
                     gsurf%knots(2), ' got', dsurf%knots(2)
            return
        end if

        fgrid = dsurf%eval(xeval, yeval, ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] gridded ds/dx eval failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! fgrid(j,i) should be approximately yeval(j)
        do i = 1, 3
            do j = 1, 3
                if (abs(fgrid(j,i) - yeval(j)) > TOL) then
                    print *, '[test_derivative_spline] gridded ds/dx mismatch at', i, j, &
                             ': expected', yeval(j), ' got', fgrid(j,i)
                    return
                end if
            end do
        end do

        ! ds/dy = x
        dsurf = gsurf%derivative_spline(0_FP_SIZE, 1_FP_SIZE, ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] gridded ds/dy failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Verify reduced order: kx unchanged, ky drops by 1
        if (dsurf%order(1) /= gsurf%order(1)) then
            print *, '[test_derivative_spline] gridded ds/dy order(1) mismatch: expected', &
                     gsurf%order(1), ' got', dsurf%order(1)
            return
        end if
        if (dsurf%order(2) /= gsurf%order(2) - 1) then
            print *, '[test_derivative_spline] gridded ds/dy order(2) mismatch: expected', &
                     gsurf%order(2)-1, ' got', dsurf%order(2)
            return
        end if

        ! Verify reduced knots: nx unchanged, ny drops by 2
        if (dsurf%knots(1) /= gsurf%knots(1)) then
            print *, '[test_derivative_spline] gridded ds/dy knots(1) mismatch: expected', &
                     gsurf%knots(1), ' got', dsurf%knots(1)
            return
        end if
        if (dsurf%knots(2) /= gsurf%knots(2) - 2) then
            print *, '[test_derivative_spline] gridded ds/dy knots(2) mismatch: expected', &
                     gsurf%knots(2)-2, ' got', dsurf%knots(2)
            return
        end if

        fgrid = dsurf%eval(xeval, yeval, ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] gridded ds/dy eval failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! fgrid(j,i) should be approximately xeval(i)
        do i = 1, 3
            do j = 1, 3
                if (abs(fgrid(j,i) - xeval(i)) > TOL) then
                    print *, '[test_derivative_spline] gridded ds/dy mismatch at', i, j, &
                             ': expected', xeval(i), ' got', fgrid(j,i)
                    return
                end if
            end do
        end do

        ! --- Scattered surface ---
        call meshgrid(xg, yg, xs, ys)
        zs = xs * ys

        ierr = surf%new_fit(xs, ys, zs, smoothing=zero)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] scattered fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! ds/dx = y (scattered)
        dsurf_s = surf%derivative_spline(1_FP_SIZE, 0_FP_SIZE, ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] scattered ds/dx failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Verify reduced order
        if (dsurf_s%order(1) /= surf%order(1) - 1) then
            print *, '[test_derivative_spline] scattered ds/dx order(1) mismatch: expected', &
                     surf%order(1)-1, ' got', dsurf_s%order(1)
            return
        end if

        ! Evaluate scattered derivative on grid and compare with analytical ds/dx = y
        fgrid = dsurf_s%eval_ongrid(xeval, yeval, ierr)
        if (.not.FITPACK_SUCCESS(ierr)) then
            print *, '[test_derivative_spline] scattered ds/dx eval failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        do i = 1, 3
            do j = 1, 3
                if (abs(fgrid(j,i) - yeval(j)) > TOL) then
                    print *, '[test_derivative_spline] scattered ds/dx mismatch at', i, j, &
                             ': expected', yeval(j), ' got', fgrid(j,i)
                    return
                end if
            end do
        end do

        success = .true.

    end function test_derivative_spline

    ! ODE-style reciprocal error weight
    elemental real(FP_REAL) function rewt(RTOL,ATOL,x)
       real(FP_REAL), intent(in) :: RTOL,ATOL,x
       rewt = one/(RTOL*abs(x)+ATOL)
    end function rewt

    ! Simple linspace function
    pure function linspace(x1,x2,n)
       real(FP_REAL), intent(in) :: x1,x2
       integer,       intent(in) :: n
       real(FP_REAL), dimension(max(2,n)) :: linspace

       integer :: nx,i
      real(FP_REAL) :: dx

       nx = max(n,2)
       dx = (x2-x1)/(nx-1)
       forall(i=1:nx) linspace(i) = x1+dx*(i-1)

    end function linspace

    ! Flatten 1D grid vectors into scattered (x,y) coordinate arrays
    pure subroutine meshgrid(xg, yg, xs, ys)
       real(FP_REAL), intent(in)  :: xg(:), yg(:)
       real(FP_REAL), intent(out) :: xs(size(xg)*size(yg)), ys(size(xg)*size(yg))
       integer :: i, j, nx
       nx = size(xg)
       forall(j=1:size(yg), i=1:nx)
          xs((j-1)*nx+i) = xg(i)
          ys((j-1)*nx+i) = yg(j)
       end forall
    end subroutine meshgrid


end module fitpack_curve_tests
