! Example: Univariate spline curve fitting with fitpack_curve
!
! Demonstrates smoothing, interpolation, evaluation, derivatives,
! integration, and root-finding on a noisy sine wave.
program example_curve
    use fitpack, only: fitpack_curve, FP_REAL, FP_FLAG, FITPACK_SUCCESS, FITPACK_MESSAGE, &
                       zero, one, half, pi
    implicit none

    integer, parameter :: m = 50
    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: pi2 = two * pi

    real(FP_REAL) :: x(m), y(m), w(m), dx, noise
    type(fitpack_curve) :: curve
    integer(FP_FLAG) :: ierr
    real(FP_REAL) :: y_eval, dydx, area, exact_area
    real(FP_REAL), allocatable :: roots(:)
    integer :: i

    ! Generate noisy sine data on [0, 2*pi]
    dx = pi2 / (m - 1)
    do i = 1, m
        x(i) = (i - 1) * dx
        ! Simple deterministic "noise" for reproducibility
        noise = 0.05_FP_REAL * sin(17.0_FP_REAL * x(i)) * cos(31.0_FP_REAL * x(i))
        y(i) = sin(x(i)) + noise
    end do

    ! Unit weights
    w = one

    ! ---- Smoothing spline (automatic knots) ----
    print '(a)', '=== Smoothing spline ==='
    ierr = curve%new_fit(x, y, w=w, smoothing=real(m, FP_REAL))
    if (.not. FITPACK_SUCCESS(ierr)) then
        print '(a,a)', 'Fit failed: ', FITPACK_MESSAGE(ierr)
        error stop
    end if
    print '(a,i0)',    '  Knots:    ', curve%knots
    print '(a,es10.3)', '  Residual: ', curve%mse()

    ! Evaluate at a single point
    y_eval = curve%eval(pi, ierr)
    print '(a,f8.5,a,f8.5)', '  s(pi) = ', y_eval, '  (exact sin(pi) = ', sin(pi), ')'

    ! ---- Interpolation ----
    print '(/,a)', '=== Interpolating spline ==='
    ierr = curve%new_fit(x, y, smoothing=zero)
    if (.not. FITPACK_SUCCESS(ierr)) then
        print '(a,a)', 'Interpolation failed: ', FITPACK_MESSAGE(ierr)
        error stop
    end if
    print '(a,i0)',    '  Knots:    ', curve%knots
    print '(a,es10.3)', '  Residual: ', curve%mse()

    ! ---- Derivatives ----
    print '(/,a)', '=== Derivatives at x = pi/4 ==='
    y_eval = curve%eval(pi / 4, ierr)
    dydx   = curve%dfdx(pi / 4, order=1, ierr=ierr)
    print '(a,f10.6)', '  s(pi/4)  = ', y_eval
    print '(a,f10.6)', '  s''(pi/4) = ', dydx

    ! ---- Integration ----
    print '(/,a)', '=== Integration ==='
    area = curve%integral(zero, pi)
    exact_area = two  ! integral of sin(x) from 0 to pi
    print '(a,f10.6)', '  integral [0, pi]       = ', area
    print '(a,f10.6)', '  exact (noise-free)     = ', exact_area

    ! ---- Zeros ----
    print '(/,a)', '=== Zeros ==='
    roots = curve%zeros(ierr)
    if (FITPACK_SUCCESS(ierr)) then
        print '(a,i0)', '  Number of zeros: ', size(roots)
        do i = 1, size(roots)
            print '(a,i0,a,f10.6)', '  root(', i, ') = ', roots(i)
        end do
    end if

    ! ---- Smoothing parameter sweep ----
    print '(/,a)', '=== Smoothing sweep ==='
    print '(a)', '     S       knots   residual'
    do i = 1, 5
        ierr = curve%new_fit(x, y, smoothing=10.0_FP_REAL**(3 - i))
        if (FITPACK_SUCCESS(ierr)) then
            print '(es10.1,i8,es12.3)', curve%smoothing, curve%knots, curve%mse()
        end if
    end do

    call curve%destroy()
    print '(/,a)', 'Done.'

end program example_curve
