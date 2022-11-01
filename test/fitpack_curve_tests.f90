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


    implicit none
    private

    public :: test_sine_fit


    contains

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
