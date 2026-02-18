! Example: Polar domain fitting with fitpack_polar and fitpack_grid_polar
!
! Fits splines on disc-shaped domains: scattered data on a unit disk
! and gridded data on a polar grid.
program example_polar
    use fitpack, only: fitpack_polar, fitpack_grid_polar, &
                       FP_REAL, FP_FLAG, FITPACK_SUCCESS, FITPACK_MESSAGE, &
                       zero, one, half, pi
    implicit none

    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: pi2 = two * pi

    ! Part 1: Scattered polar
    call demo_scattered_polar()

    ! Part 2: Gridded polar
    call demo_gridded_polar()

    print '(/,a)', 'Done.'

contains

    ! Boundary function: unit disk (radius = 1 for all angles)
    pure real(FP_REAL) function unit_disk(theta)
        real(FP_REAL), intent(in) :: theta
        unit_disk = one
    end function

    ! Test function on the disk
    elemental real(FP_REAL) function test_polar(x, y) result(f)
        real(FP_REAL), intent(in) :: x, y
        f = (x**2 + y**2) / ((x + y)**2 + half)
    end function

    subroutine demo_scattered_polar()
        integer, parameter :: m = 100
        real(FP_REAL) :: x(m), y(m), z(m), w(m)
        real(FP_REAL) :: r, theta
        type(fitpack_polar) :: polar
        integer(FP_FLAG) :: ierr
        real(FP_REAL) :: z_val, exact
        real(FP_REAL) :: x_test(5), y_test(5)
        integer :: i

        print '(a)', '=== Scattered polar: unit disk ==='

        ! Generate scattered points inside the unit disk
        do i = 1, m
            ! Simple quasi-random on disk via rejection-free polar coords
            r = sqrt(mod(i * 0.6180339887_FP_REAL, one))
            theta = pi2 * mod(i * 0.3247179572_FP_REAL, one) - pi
            x(i) = r * cos(theta)
            y(i) = r * sin(theta)
            z(i) = test_polar(x(i), y(i))
        end do

        ! Unit weights
        w = one

        ! Smoothing fit
        ierr = polar%new_fit(x, y, z, unit_disk, w, smoothing=real(m, FP_REAL))
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0,a,i0)', '  Knots: ', polar%knots(1), ' x ', polar%knots(2)
        print '(a,es10.3)',   '  Residual: ', polar%mse()

        ! Evaluate at test points
        print '(/,a)', '  Point evaluation:'
        print '(a)', '     x       y       s(x,y)    exact'
        x_test = [zero, half, -half, 0.3_FP_REAL, -0.2_FP_REAL]
        y_test = [zero, zero,  half, 0.4_FP_REAL, -0.5_FP_REAL]
        do i = 1, 5
            z_val = polar%eval(x_test(i), y_test(i), ierr)
            exact = test_polar(x_test(i), y_test(i))
            print '(2f8.3, 2f12.6)', x_test(i), y_test(i), z_val, exact
        end do

        ! Tighter fit
        ierr = polar%fit(smoothing=10.0_FP_REAL)
        if (FITPACK_SUCCESS(ierr)) then
            print '(/,a,i0,a,i0)', '  Tighter fit knots: ', polar%knots(1), ' x ', polar%knots(2)
            print '(a,es10.3)', '  Residual: ', polar%mse()
        end if

        call polar%destroy()
    end subroutine

    subroutine demo_gridded_polar()
        integer, parameter :: nu = 10, nv = 20
        real(FP_REAL) :: u(nu), v(nv), zg(nv, nu)
        real(FP_REAL) :: du, dv, x_pt, y_pt, z0, z_val
        type(fitpack_grid_polar) :: gpol
        integer(FP_FLAG) :: ierr
        integer :: i, j

        print '(/,a)', '=== Gridded polar: unit disk ==='

        ! Create polar grid: u in (0,1], v in [-pi, pi)
        du = one / nu
        dv = pi2 / nv
        do i = 1, nu
            u(i) = i * du
        end do
        do j = 1, nv
            v(j) = -pi + (j - 1) * dv
        end do

        ! Evaluate test function on grid
        do i = 1, nu
            do j = 1, nv
                x_pt = u(i) * cos(v(j))
                y_pt = u(i) * sin(v(j))
                zg(j, i) = test_polar(x_pt, y_pt)
            end do
        end do

        ! Value at origin
        z0 = test_polar(zero, zero)

        ! Fit with prescribed origin value
        ierr = gpol%new_fit(u, v, one, zg, z0=z0)
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0,a,i0)', '  Knots: ', gpol%knots(1), ' x ', gpol%knots(2)
        print '(a,es10.3)',   '  Residual: ', gpol%mse()

        ! Set origin BC: exact value, differentiable
        call gpol%set_origin_BC(z0=z0, exact=.true., differentiable=.true.)
        ierr = gpol%fit(smoothing=real(nu * nv, FP_REAL))
        if (FITPACK_SUCCESS(ierr)) then
            print '(a,es10.3)', '  Residual (with origin BC): ', gpol%mse()
        end if

        ! Evaluate at origin (u=0 maps to any v)
        z_val = gpol%eval(zero, zero, ierr)
        print '(/,a,f12.6)', '  s(origin) = ', z_val
        print '(a,f12.6)',   '  exact     = ', z0

        call gpol%destroy()
    end subroutine

end program example_polar
