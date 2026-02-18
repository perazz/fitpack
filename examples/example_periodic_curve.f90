! Example: Periodic spline fitting with fitpack_periodic_curve
!
! Fits a periodic spline to f(x) = cos(x) + sin(2x) on [0, 2*pi],
! demonstrating that derivatives match across the period boundary.
program example_periodic_curve
    use fitpack, only: fitpack_periodic_curve, FP_REAL, FP_FLAG, &
                       FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half, pi
    implicit none

    integer, parameter :: m = 200
    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: pi2 = two * pi

    type(fitpack_periodic_curve) :: curve
    integer(FP_FLAG) :: ierr
    real(FP_REAL) :: x(m), y(m)
    real(FP_REAL) :: y_left, y_right, dy_left, dy_right, area, exact_area, dx_eval
    integer :: i

    ! Generate periodic function data on [0, 2*pi]
    ! Include endpoint: percur requires full period coverage
    do i = 1, m
        x(i) = (i - 1) * pi2 / (m - 1)
        y(i) = fun(x(i))
    end do

    ! ---- Smoothing fit ----
    print '(a)', '=== Smoothing periodic fit ==='
    ierr = curve%new_fit(x, y, smoothing=one)
    if (.not. FITPACK_SUCCESS(ierr)) then
        print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
    else
        print '(a,i0)',    '  Knots:    ', curve%knots
        print '(a,es10.3)', '  Residual: ', curve%mse()
    end if

    ! Evaluate at a few points and compare
    print '(/,a)', '=== Evaluation ==='
    print '(a)', '     x         s(x)        exact'
    do i = 0, 8
        dx_eval = i * pi / 4
        print '(f8.4,f12.6,f12.6)', dx_eval, curve%eval(dx_eval, ierr), fun(dx_eval)
    end do

    ! Check periodicity: values and derivatives at endpoints
    y_left   = curve%eval(x(1), ierr)
    y_right  = curve%eval(x(m), ierr)
    dy_left  = curve%dfdx(x(1), order=1, ierr=ierr)
    dy_right = curve%dfdx(x(m), order=1, ierr=ierr)

    print '(/,a)', '=== Periodicity check ==='
    print '(a,f12.8)', '  s(x_1)    = ', y_left
    print '(a,f12.8)', '  s(x_m)    = ', y_right
    print '(a,f12.8)', '  s''(x_1)   = ', dy_left
    print '(a,f12.8)', '  s''(x_m)   = ', dy_right

    ! ---- Integration ----
    print '(/,a)', '=== Integration ==='
    area = curve%integral(half * pi, two * pi)
    ! Exact: integral of cos(x)+sin(2x) from pi/2 to 2*pi
    ! = [sin(x) - cos(2x)/2] from pi/2 to 2*pi
    exact_area = (sin(two * pi) - half * cos(4.0_FP_REAL * pi)) &
               - (sin(half * pi) - half * cos(pi))
    print '(a,f12.8)', '  integral [pi/2, 2*pi] = ', area
    print '(a,f12.8)', '  exact                 = ', exact_area

    ! ---- Interpolating fit ----
    print '(/,a)', '=== Interpolating periodic spline ==='
    ierr = curve%interpolate()
    if (.not. FITPACK_SUCCESS(ierr)) then
        print '(a,a)', '  Interpolation failed: ', FITPACK_MESSAGE(ierr)
    else
        print '(a,i0)',    '  Knots:    ', curve%knots
        print '(a,es10.3)', '  Residual: ', curve%mse()

        y_left  = curve%eval(x(1), ierr)
        y_right = curve%eval(x(m), ierr)
        print '(a,f12.8)', '  s(0)    = ', y_left
        print '(a,f12.8)', '  s(x_m)  = ', y_right
    end if

    ! ---- Smoothing sweep ----
    print '(/,a)', '=== Smoothing sweep ==='
    print '(a)', '     S       knots   residual'
    do i = 1, 4
        ierr = curve%new_fit(x, y, smoothing=10.0_FP_REAL**(2 - i))
        if (FITPACK_SUCCESS(ierr)) then
            print '(es10.1,i8,es12.3)', curve%smoothing, curve%knots, curve%mse()
        end if
    end do

    call curve%destroy()
    print '(/,a)', 'Done.'

contains

    elemental real(FP_REAL) function fun(x) result(f)
        real(FP_REAL), intent(in) :: x
        f = cos(x) + sin(two * x)
    end function

end program example_periodic_curve
