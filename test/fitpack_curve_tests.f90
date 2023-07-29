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
    use iso_fortran_env, only: output_unit


    implicit none
    private

    public :: test_sine_fit
    public :: test_periodic_fit
    public :: test_parametric_fit


    contains

    ! Test parametric curve
    logical function test_parametric_fit() result(success)

       real(RKIND), allocatable :: x(:,:),u(:),w(:),y(:,:)
       type(fitpack_parametric_curve) :: curve
       real(RKIND) :: s
       integer :: ierr,loop,useUnit

       useUnit = output_unit

       ! Set data points
       u = [real(RKIND) :: 120,128,133,136,138,141,144,146,149,151,154,161,170,180,190,&
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
             case (1); ierr = curve%new_fit(x,u=u,w=w,smoothing=100.0_RKIND)

             ! Smaller values of s to get a tighter approximation
             case (2); ierr = curve%fit(smoothing=1.0_RKIND)
             case (3); ierr = curve%fit(smoothing=0.05_RKIND)

             ! Larger values to get a smoother approximation
             case (4); ierr = curve%fit(smoothing=0.25_RKIND)

             ! New fit to get possibly fewer knots
             case (5); ierr = curve%new_fit(x,u=u,w=w,smoothing=0.25_RKIND)

             ! Let the program determine parameter values u(i)
             case (6); ierr = curve%new_fit(x,w=w,smoothing=0.25_RKIND)

             ! Quintic spline approximation
             case (7); ierr = curve%new_fit(x,w=w,smoothing=0.25_RKIND,order=5)

             ! Interpolating curve
             case (8); ierr = curve%interpolate()

!                 case (9)
!
!                    !  finally we calculate a least-squares spline curve with specified knots
!                    iopt =-1
!                    n = 9+2*k
!                    j = k+2
!                    del = (ue-ub)*0.125_RKIND
!                    do l=1,7
!                        al = l
!                        t(j) = ub+al*del
!                        j = j+1
!                    end do

          end select

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
          end if

          ! Evaluate the spline curve
          y = curve%eval(u,ierr)

          if (.not.FITPACK_SUCCESS(ierr)) then
              success = .false.
              write(useUnit,1000) loop,FITPACK_MESSAGE(ierr)
          end if



       end do

       1000  format('[test_parametric_fit] parametric curve test ',i0,' failed: ',a)

    end function test_parametric_fit

    ! Dumb test: fit a sine function
    logical function test_sine_fit() result(success)

       integer, parameter     :: N = 200
       real(RKIND), parameter :: RTOL = 1.0e-1_RKIND
       real(RKIND), parameter :: ATOL = 1.0e-2_RKIND
       type(fitpack_curve) :: curve
       real(RKIND) :: x(N),y(N),xrand(N),yeval,yprime,dfdx(0:3)
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
          yeval = curve%eval(xrand(i),ierr);
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

       ! All checks passed: success!
       success = .true.

       1 format('[sine_fit] sine function error is too large: x=',f6.2,' yspline=',f6.2,' analytical=',f6.2)
       2 format('[sine_fit] cannot evaluate ',i0,'-th derivative at ',f6.2,': ',a)
       3 format('[sine_fit] ',i0,'-th derivative error is too large: x=',f6.2,' yp(spline)=',f6.2,' analytical=',f6.2)

    end function test_sine_fit

    ! Periodic test: fit a cosine function
    logical function test_periodic_fit() result(success)

       integer, parameter     :: N = 200
       real(RKIND), parameter :: RTOL = 1.0e-2_RKIND
       real(RKIND), parameter :: ATOL = 1.0e-4_RKIND
       type(fitpack_periodic_curve) :: curve
       real(RKIND) :: x(N),y(N),xrand(N),yeval,yprime,dfdx(0:3)
       integer :: ierr,i,order

       success = .false.

       ! Generate a periodic function with 200 points
       x = linspace(zero,pi2,N)
       y = cos(x) + sin(2*x)

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

       ! All checks passed: success!
       success = .true.

       1 format('[periodic_fit] sine function error is too large: x=',f6.2,' yspline=',f6.2,' analytical=',f6.2)
       2 format('[periodic_fit] cannot evaluate ',i0,'-th derivative at ',f6.2,': ',a)
       3 format('[periodic_fit] ',i0,'-th derivative error is too large: x=',f6.2,' yp(spline)=',f6.2,' analytical=',f6.2)

    end function test_periodic_fit

    ! ODE-style reciprocal error weight
    elemental real(RKIND) function rewt(RTOL,ATOL,x)
       real(RKIND), intent(in) :: RTOL,ATOL,x
       rewt = one/(RTOL*abs(x)+ATOL)
    end function rewt

    ! Simple linspace function
    pure function linspace(x1,x2,n)
       real(RKIND), intent(in) :: x1,x2
       integer,     intent(in) :: n
       real(RKIND), dimension(max(2,n)) :: linspace

       integer :: nx,i
       real(RKIND) :: dx

       nx = max(n,2)
       dx = (x2-x1)/(nx-1)
       forall(i=1:nx) linspace(i) = x1+dx*(i-1)

    end function linspace


end module fitpack_curve_tests
