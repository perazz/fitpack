! Example: Spherical domain fitting with fitpack_sphere and fitpack_grid_sphere
!
! Fits splines on the unit sphere: scattered data and latitude-longitude grids,
! demonstrating pole boundary conditions.
program example_sphere
    use fitpack, only: fitpack_sphere, fitpack_grid_sphere, &
                       FP_REAL, FP_FLAG, FITPACK_SUCCESS, FITPACK_MESSAGE, &
                       zero, one, half, pi
    implicit none

    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: pi2 = two * pi

    ! Part 1: Scattered sphere
    call demo_scattered_sphere()

    ! Part 2: Gridded sphere
    call demo_gridded_sphere()

    print '(/,a)', 'Done.'

contains

    ! Smooth test function on the sphere
    elemental real(FP_REAL) function test_sphere(theta, phi) result(f)
        real(FP_REAL), intent(in) :: theta, phi
        ! f = 1 + cos(theta) + 0.5*sin(theta)^2*cos(2*phi)
        ! This is a linear combination of spherical harmonics (Y_0^0, Y_1^0, Y_2^2)
        f = one + cos(theta) + half * sin(theta)**2 * cos(two * phi)
    end function

    subroutine demo_scattered_sphere()
        integer, parameter :: m = 200
        real(FP_REAL) :: theta(m), phi(m), r(m), w(m)
        type(fitpack_sphere) :: sph
        integer(FP_FLAG) :: ierr
        real(FP_REAL), allocatable :: f(:,:)
        real(FP_REAL) :: t_eval(5), p_eval(5)
        real(FP_REAL) :: avg_err, max_err
        integer :: i

        print '(a)', '=== Scattered sphere fit ==='

        ! Generate quasi-random points on the sphere
        do i = 1, m
            ! Uniform in cos(theta) for uniform area coverage
            theta(i) = acos(one - two * mod(i * 0.6180339887_FP_REAL, one))
            phi(i) = pi2 * mod(i * 0.3247179572_FP_REAL, one)
            r(i) = test_sphere(theta(i), phi(i))
        end do

        ! Unit weights
        w = one

        ! Smoothing fit
        ierr = sph%new_fit(theta, phi, r, smoothing=real(m, FP_REAL))
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0,a,i0)', '  Knots: ', sph%knots(1), ' x ', sph%knots(2)
        print '(a,es10.3)',   '  Residual: ', sph%mse()

        ! Evaluate on a small grid (eval_many returns grid output)
        t_eval = [0.1_FP_REAL, half, one, two, pi - 0.1_FP_REAL]
        p_eval = [0.1_FP_REAL, one, two, pi, 5.0_FP_REAL]
        f = sph%eval(t_eval, p_eval, ierr)
        ! f has shape (size(p_eval), size(t_eval))
        if (FITPACK_SUCCESS(ierr)) then
            print '(/,a)', '  Grid evaluation sample (theta=0.5, phi=1.0):'
            print '(a,f12.6)', '    s(0.5, 1.0) = ', f(2, 2)
            print '(a,f12.6)', '    exact       = ', test_sphere(t_eval(2), p_eval(2))
        end if

        ! Tighter fit
        ierr = sph%fit(smoothing=10.0_FP_REAL)
        if (FITPACK_SUCCESS(ierr)) then
            print '(/,a)',      '  Tighter fit:'
            print '(a,i0,a,i0)', '    Knots: ', sph%knots(1), ' x ', sph%knots(2)
            print '(a,es10.3)',   '    Residual: ', sph%mse()

            f = sph%eval(t_eval, p_eval, ierr)
            if (FITPACK_SUCCESS(ierr)) then
                print '(a,f12.6)', '    s(0.5, 1.0) = ', f(2, 2)
            end if
        end if

        call sph%destroy()
    end subroutine

    subroutine demo_gridded_sphere()
        integer, parameter :: nu = 15, nv = 30
        real(FP_REAL) :: u(nu), v(nv), zg(nv, nu)
        real(FP_REAL) :: du, dv, z_north, z_south, z_val
        type(fitpack_grid_sphere) :: gsph
        integer(FP_FLAG) :: ierr
        real(FP_REAL), allocatable :: f(:,:)
        real(FP_REAL) :: u_eval(3), v_eval(3)
        integer :: i, j

        print '(/,a)', '=== Gridded sphere fit ==='

        ! Create colatitude-longitude grid
        ! u (colatitude) in (0, pi), v (longitude) in [0, 2*pi)
        du = pi / (nu + 1)
        dv = pi2 / nv
        do i = 1, nu
            u(i) = i * du
        end do
        do j = 1, nv
            v(j) = (j - 1) * dv
        end do

        ! Evaluate function on grid
        do i = 1, nu
            do j = 1, nv
                zg(j, i) = test_sphere(u(i), v(j))
            end do
        end do

        ! Pole values
        z_north = test_sphere(zero, zero)
        z_south = test_sphere(pi, zero)

        ! Fit
        ierr = gsph%new_fit(u, v, zg, smoothing=zero)
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0,a,i0)', '  Knots: ', gsph%knots(1), ' x ', gsph%knots(2)
        print '(a,es10.3)',   '  Residual: ', gsph%mse()

        ! Set pole boundary conditions
        call gsph%BC_north_pole(z0=z_north, exact=.true., differentiable=.true.)
        call gsph%BC_south_pole(z0=z_south, exact=.true., differentiable=.true.)

        ! Re-fit with pole BCs
        ierr = gsph%fit(smoothing=zero)
        if (FITPACK_SUCCESS(ierr)) then
            print '(a,es10.3)', '  Residual (with pole BCs): ', gsph%mse()
        end if

        ! Evaluate
        u_eval = [0.1_FP_REAL, half * pi, pi - 0.1_FP_REAL]
        v_eval = [zero, pi, pi2 - 0.1_FP_REAL]
        f = gsph%eval(u_eval, v_eval, ierr)
        if (FITPACK_SUCCESS(ierr)) then
            print '(/,a)', '  Evaluation at select points:'
            print '(a)', '     theta    phi      s(t,p)    exact'
            do i = 1, 3
                print '(2f9.4, 2f12.6)', u_eval(i), v_eval(i), &
                    f(i, i), test_sphere(u_eval(i), v_eval(i))
            end do
        end if

        ! Check pole values
        z_val = gsph%eval(0.001_FP_REAL, zero, ierr)
        print '(/,a,f12.6)', '  Near north pole: s(0.001, 0) = ', z_val
        print '(a,f12.6)',   '  Exact at pole:               = ', z_north

        call gsph%destroy()
    end subroutine

end program example_sphere
