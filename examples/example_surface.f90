! Example: Scattered surface fitting with fitpack_surface
!
! Fits a bivariate spline to scattered data sampled from Franke's test function,
! demonstrating smoothing, evaluation, partial derivatives, and integration.
program example_surface
    use fitpack, only: fitpack_surface, fitpack_curve, FP_REAL, FP_FLAG, &
                       FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half, pi
    implicit none

    integer, parameter :: m = 200
    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: nine = 9.0_FP_REAL

    real(FP_REAL) :: x(m), y(m), z(m)
    type(fitpack_surface) :: surf
    type(fitpack_curve) :: cross
    integer(FP_FLAG) :: ierr
    real(FP_REAL) :: z_val, dzdx, dzdy, vol, exact_eval
    real(FP_REAL) :: x_eval(5), y_eval(5), z_eval(5)
    integer :: i

    ! Generate scattered data from a smooth test function
    ! f(x,y) = sin(pi*x)*cos(pi*y) on [0,1]x[0,1]
    ! Using a simple quasi-random distribution
    do i = 1, m
        x(i) = mod(i * 0.6180339887_FP_REAL, one)        ! golden ratio mod 1
        y(i) = mod(i * 0.3247179572_FP_REAL, one)        ! sqrt(2)/4 mod 1
        z(i) = test_func(x(i), y(i))
    end do

    ! ---- Smoothing fit ----
    print '(a)', '=== Scattered surface: smoothing fit ==='
    ierr = surf%new_fit(x, y, z, smoothing=real(m, FP_REAL))
    if (.not. FITPACK_SUCCESS(ierr)) then
        print '(a,a)', 'Fit failed: ', FITPACK_MESSAGE(ierr)
        error stop
    end if
    print '(a,i0,a,i0)', '  Knots: ', surf%knots(1), ' x ', surf%knots(2)
    print '(a,es10.3)',   '  Residual: ', surf%mse()

    ! Evaluate at a few test points
    print '(/,a)', '=== Point evaluation ==='
    print '(a)', '     x       y       s(x,y)    exact'
    x_eval = [0.25_FP_REAL, half, 0.75_FP_REAL, 0.1_FP_REAL, 0.9_FP_REAL]
    y_eval = [0.25_FP_REAL, half, 0.75_FP_REAL, 0.9_FP_REAL, 0.1_FP_REAL]
    do i = 1, 5
        z_val = surf%eval(x_eval(i), y_eval(i), ierr)
        exact_eval = test_func(x_eval(i), y_eval(i))
        print '(2f8.3, 2f12.6)', x_eval(i), y_eval(i), z_val, exact_eval
    end do

    ! ---- Partial derivatives ----
    print '(/,a)', '=== Partial derivatives at (0.5, 0.5) ==='
    dzdx = surf%dfdx(half, half, 1, 0, ierr)
    dzdy = surf%dfdx(half, half, 0, 1, ierr)
    print '(a,f12.6)', '  ds/dx  = ', dzdx
    print '(a,f12.6)', '  ds/dy  = ', dzdy
    print '(a,f12.6)', '  exact ds/dx = ', pi * cos(pi * half) * cos(pi * half)
    print '(a,f12.6)', '  exact ds/dy = ', -pi * sin(pi * half) * sin(pi * half)

    ! ---- Integration ----
    print '(/,a)', '=== Integration over [0,1] x [0,1] ==='
    vol = surf%integral([zero, zero], [one, one])
    print '(a,f12.6)', '  integral = ', vol
    ! Exact: integral of sin(pi*x)*cos(pi*y) dxdy over [0,1]^2
    ! = (2/pi) * (0) = 0 (since integral of cos(pi*y) from 0 to 1 = 0)
    print '(a,f12.6)', '  exact    = ', zero

    ! ---- Cross-section ----
    print '(/,a)', '=== Cross-section at y = 0.5 ==='
    cross = surf%cross_section(half, .true., ierr)
    if (FITPACK_SUCCESS(ierr)) then
        print '(a,i0)', '  Cross-section knots: ', cross%knots
        print '(a)', '     x       profile   exact'
        do i = 1, 5
            z_val = cross%eval(x_eval(i), ierr)
            exact_eval = test_func(x_eval(i), half)
            print '(f8.3, 2f12.6)', x_eval(i), z_val, exact_eval
        end do
    end if

    ! ---- Tighter fit ----
    print '(/,a)', '=== Smoothing sweep ==='
    print '(a)', '     S       knots_x  knots_y  residual'
    do i = 1, 4
        ierr = surf%new_fit(x, y, z, smoothing=10.0_FP_REAL**(3 - i))
        if (FITPACK_SUCCESS(ierr)) then
            print '(es10.1, 2i9, es12.3)', surf%smoothing, surf%knots(1), surf%knots(2), surf%mse()
        end if
    end do

    call cross%destroy()
    call surf%destroy()
    print '(/,a)', 'Done.'

contains

    elemental real(FP_REAL) function test_func(x, y) result(f)
        real(FP_REAL), intent(in) :: x, y
        f = sin(pi * x) * cos(pi * y)
    end function

end program example_surface
