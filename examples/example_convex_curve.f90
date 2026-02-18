! Example: Convexity-constrained curve fitting with fitpack_convex_curve
!
! Fits a concave spline to data inspired by moisture-content measurements:
! monotonically increasing, concave (s''(x) <= 0) function.
program example_convex_curve
    use fitpack, only: fitpack_convex_curve, FP_REAL, FP_FLAG, &
                       FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one
    implicit none

    integer, parameter :: m = 16
    type(fitpack_convex_curve) :: curve
    integer(FP_FLAG) :: ierr
    real(FP_REAL) :: x(m), y(m), w(m), v(m), s2(m)
    real(FP_REAL) :: smoothing_values(3)
    real(FP_REAL) :: y_eval, dx
    integer :: i, is

    ! Data: moisture content vs depth (inspired by Dierckx Ch. 7)
    ! Concave, increasing trend
    x = [0.1_FP_REAL, 0.3_FP_REAL, 0.5_FP_REAL, 0.7_FP_REAL, 0.9_FP_REAL, &
         1.25_FP_REAL, 1.75_FP_REAL, 2.25_FP_REAL, 2.75_FP_REAL, 3.5_FP_REAL, &
         4.5_FP_REAL, 5.5_FP_REAL, 6.5_FP_REAL, 7.5_FP_REAL, 8.5_FP_REAL, 9.5_FP_REAL]

    y = [0.124_FP_REAL, 0.234_FP_REAL, 0.256_FP_REAL, 0.277_FP_REAL, 0.278_FP_REAL, &
         0.291_FP_REAL, 0.308_FP_REAL, 0.311_FP_REAL, 0.315_FP_REAL, 0.322_FP_REAL, &
         0.317_FP_REAL, 0.326_FP_REAL, 0.323_FP_REAL, 0.321_FP_REAL, 0.322_FP_REAL, &
         0.328_FP_REAL]

    ! Weights: emphasize endpoints
    w = one
    w(1) = 10.0_FP_REAL
    w(2) = 3.0_FP_REAL
    w(m) = 10.0_FP_REAL

    ! Concavity constraint: v(i) = 1 means concave (s''(x) <= 0)
    ! Use v = -1 for convex (s''(x) >= 0), v = 0 for unconstrained
    v = one

    ! Load data and set convexity flags
    call curve%new_points(x, y, w)
    ierr = curve%set_convexity(v)
    if (.not. FITPACK_SUCCESS(ierr)) then
        print '(a,a)', 'set_convexity failed: ', FITPACK_MESSAGE(ierr)
        error stop
    end if

    ! Fit with decreasing smoothing
    smoothing_values = [0.2_FP_REAL, 0.04_FP_REAL, 0.0002_FP_REAL]

    print '(a)', '=== Concavity-constrained fits ==='
    print '(a)', '     S       knots   residual   max s''''(x)'

    do is = 1, 3
        ierr = curve%fit(smoothing=smoothing_values(is), keep_knots=(is > 1))
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,es8.1,a,a)', 'Fit at s=', smoothing_values(is), ' failed: ', FITPACK_MESSAGE(ierr)
            cycle
        end if

        ! Evaluate second derivative at data points
        s2 = curve%dfdx(x, order=2, ierr=ierr)

        print '(es10.1, i8, es12.3, es12.3)', &
            smoothing_values(is), curve%knots, curve%mse(), maxval(s2)
    end do

    ! Evaluate the tightest fit at several points
    print '(/,a)', '=== Evaluation (tightest fit) ==='
    print '(a)', '     x       y_data    s(x)      s''''(x)'
    do i = 1, m
        y_eval = curve%eval(x(i), ierr)
        print '(f8.3, f10.4, f10.4, es12.3)', x(i), y(i), y_eval, s2(i)
    end do

    ! Evaluate at intermediate points
    print '(/,a)', '=== Interpolated values ==='
    print '(a)', '     x        s(x)      s''(x)'
    do i = 1, 10
        dx = i * one
        y_eval = curve%eval(dx, ierr)
        print '(f8.3, f10.5, es12.3)', dx, y_eval, curve%dfdx(dx, order=2, ierr=ierr)
    end do

    call curve%destroy()
    print '(/,a)', 'Done.'

end program example_convex_curve
