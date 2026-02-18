! Example: Grid surface fitting and parametric surfaces
!
! Demonstrates fitpack_grid_surface on a rectangular grid (peaks-like function)
! and fitpack_parametric_surface for a torus.
program example_grid_surface
    use fitpack, only: fitpack_grid_surface, fitpack_parametric_surface, fitpack_curve, &
                       FP_REAL, FP_FLAG, FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half, pi
    implicit none

    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: pi2 = two * pi

    ! Part 1: Grid surface
    call demo_grid_surface()

    ! Part 2: Parametric surface
    call demo_parametric_surface()

    print '(/,a)', 'Done.'

contains

    subroutine demo_grid_surface()
        integer, parameter :: nx = 20, ny = 25
        real(FP_REAL) :: xg(nx), yg(ny), zg(ny, nx)
        type(fitpack_grid_surface) :: gsurf
        type(fitpack_curve) :: cross
        integer(FP_FLAG) :: ierr
        real(FP_REAL) :: z_val, vol
        real(FP_REAL) :: x_test(3), y_test(3)
        integer :: i, j

        print '(a)', '=== Grid surface: smooth test function ==='

        ! Create grid: x in [-2,2], y in [-2,2]
        do i = 1, nx
            xg(i) = -two + 4.0_FP_REAL * (i - 1) / (nx - 1)
        end do
        do j = 1, ny
            yg(j) = -two + 4.0_FP_REAL * (j - 1) / (ny - 1)
        end do

        ! f(x,y) = exp(-(x^2+y^2)/2) * cos(x) * sin(y)
        do i = 1, nx
            do j = 1, ny
                zg(j, i) = exp(-half * (xg(i)**2 + yg(j)**2)) * cos(xg(i)) * sin(yg(j))
            end do
        end do

        ! Smoothing fit
        ierr = gsurf%new_fit(xg, yg, zg, smoothing=real(nx * ny, FP_REAL))
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0,a,i0)', '  Knots: ', gsurf%knots(1), ' x ', gsurf%knots(2)
        print '(a,es10.3)',   '  Residual: ', gsurf%mse()

        ! Interpolation
        ierr = gsurf%new_fit(xg, yg, zg, smoothing=zero)
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Interpolation failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0,a,i0)', '  Interp knots: ', gsurf%knots(1), ' x ', gsurf%knots(2)

        ! Evaluate at test points
        print '(/,a)', '  Point evaluation:'
        print '(a)', '     x       y       s(x,y)'
        x_test = [zero, one, -one]
        y_test = [zero, one, -one]
        do i = 1, 3
            z_val = gsurf%eval(x_test(i), y_test(i), ierr)
            print '(2f8.3, f12.6)', x_test(i), y_test(i), z_val
        end do

        ! Integration
        vol = gsurf%integral([-two, -two], [two, two])
        print '(/,a,f12.6)', '  Integral [-2,2]x[-2,2] = ', vol

        ! Cross-section at y = 0
        print '(/,a)', '  Cross-section at y = 0:'
        cross = gsurf%cross_section(zero, .true., ierr)
        if (FITPACK_SUCCESS(ierr)) then
            print '(a,i0)', '  Profile knots: ', cross%knots
            do i = 1, 5
                z_val = cross%eval(-two + (i - 1) * one, ierr)
                print '(a,f6.2,a,f12.6)', '    f(', -two + (i - 1) * one, ', 0) = ', z_val
            end do
        end if

        ! Derivative spline
        block
            type(fitpack_grid_surface) :: dsurf
            dsurf = gsurf%derivative_spline(1, 0, ierr)
            if (FITPACK_SUCCESS(ierr)) then
                z_val = dsurf%eval(zero, zero, ierr)
                print '(/,a,f12.6)', '  ds/dx at (0,0) = ', z_val
            end if
            call dsurf%destroy()
        end block

        call cross%destroy()
        call gsurf%destroy()
    end subroutine

    subroutine demo_parametric_surface()
        integer, parameter :: nu = 20, nv = 30, idim = 3
        real(FP_REAL) :: u(nu), v(nv), z(nv, nu, idim)
        real(FP_REAL) :: du, dv, cu, su, cv, sv
        real(FP_REAL), parameter :: R = 2.0_FP_REAL, r_tube = 0.5_FP_REAL
        type(fitpack_parametric_surface) :: psurf
        integer(FP_FLAG) :: ierr
        real(FP_REAL) :: pt(idim)
        integer :: i, j

        print '(/,a)', '=== Parametric surface: torus ==='

        ! Create parameter grid
        du = pi2 / nu
        dv = pi2 / nv
        do i = 1, nu
            u(i) = (i - 1) * du
        end do
        do j = 1, nv
            v(j) = (j - 1) * dv
        end do

        ! Torus: x=(R+r*cos(v))*cos(u), y=(R+r*cos(v))*sin(u), z=r*sin(v)
        do i = 1, nu
            cu = cos(u(i)); su = sin(u(i))
            do j = 1, nv
                cv = cos(v(j)); sv = sin(v(j))
                z(j, i, 1) = (R + r_tube * cv) * cu
                z(j, i, 2) = (R + r_tube * cv) * su
                z(j, i, 3) = r_tube * sv
            end do
        end do

        ! Fit with periodicity in both directions
        ierr = psurf%new_fit(u, v, z, smoothing=zero, periodic_BC=[.true., .true.])
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0,a,i0)', '  Knots: ', psurf%knots(1), ' x ', psurf%knots(2)
        print '(a,es10.3)',   '  Residual: ', psurf%mse()

        ! Evaluate at a few points
        print '(/,a)', '     u       v       x         y         z'
        pt = psurf%eval(zero, zero, ierr)
        print '(2f8.3, 3f10.4)', zero, zero, pt
        pt = psurf%eval(pi, zero, ierr)
        print '(2f8.3, 3f10.4)', pi, zero, pt
        pt = psurf%eval(zero, pi, ierr)
        print '(2f8.3, 3f10.4)', zero, pi, pt

        call psurf%destroy()
    end subroutine

end program example_grid_surface
