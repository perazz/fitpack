! Example: Parametric, closed, and constrained curves
!
! Demonstrates fitpack_parametric_curve (Lissajous figure),
! fitpack_closed_curve (ellipse), and fitpack_constrained_curve
! (curve with prescribed endpoint derivatives).
program example_parametric_curves
    use fitpack, only: fitpack_parametric_curve, fitpack_closed_curve, &
                       fitpack_constrained_curve, FP_REAL, FP_FLAG, &
                       FITPACK_SUCCESS, FITPACK_MESSAGE, zero, one, half, pi
    implicit none

    real(FP_REAL), parameter :: two = 2.0_FP_REAL
    real(FP_REAL), parameter :: three = 3.0_FP_REAL
    real(FP_REAL), parameter :: pi2 = two * pi

    ! ---- Part 1: Parametric curve (Lissajous figure) ----
    call demo_parametric()

    ! ---- Part 2: Closed curve (ellipse) ----
    call demo_closed()

    ! ---- Part 3: Constrained curve ----
    call demo_constrained()

    print '(/,a)', 'Done.'

contains

    subroutine demo_parametric()
        integer, parameter :: m = 60
        real(FP_REAL) :: pts(2, m), u(m), du, eval_u(6)
        real(FP_REAL) :: y(2), dy(2)
        type(fitpack_parametric_curve) :: curve
        integer(FP_FLAG) :: ierr
        integer :: i

        print '(a)', '=== Parametric curve: Lissajous figure ==='

        ! Generate Lissajous: x = sin(3t), y = sin(2t)
        du = pi2 / (m - 1)
        do i = 1, m
            u(i) = (i - 1) * du
            pts(1, i) = sin(three * u(i))
            pts(2, i) = sin(two * u(i))
        end do

        ! Fit with automatic parameter values
        ierr = curve%new_fit(pts, smoothing=zero)
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0)', '  Dimensions: ', curve%idim
        print '(a,i0)', '  Knots:      ', curve%knots

        ! Evaluate at a few parameter values
        eval_u = [zero, one, two, three, pi2 - 0.01_FP_REAL, pi2]
        print '(/,a)', '     u       x(u)      y(u)'
        do i = 1, 6
            y = curve%eval_one(eval_u(i), ierr)
            print '(f8.3, 2f10.5)', eval_u(i), y
        end do

        ! Derivative at u = 0
        dy = curve%dfdx(zero, order=1, ierr=ierr)
        print '(/,a,2f10.5)', '  dx/du, dy/du at u=0: ', dy

        call curve%destroy()
    end subroutine

    subroutine demo_closed()
        integer, parameter :: m = 40
        real(FP_REAL) :: pts(2, m), du, y(2)
        type(fitpack_closed_curve) :: curve
        integer(FP_FLAG) :: ierr
        integer :: i

        print '(/,a)', '=== Closed curve: Ellipse ==='

        ! Generate ellipse: x = 2*cos(t), y = sin(t)
        du = pi2 / m
        do i = 1, m
            pts(1, i) = two * cos((i - 1) * du)
            pts(2, i) = sin((i - 1) * du)
        end do

        ! Fit closed curve
        ierr = curve%new_fit(pts, smoothing=zero)
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if
        print '(a,i0)', '  Knots: ', curve%knots

        ! Check closure: eval at start and end of parameter range
        y = curve%eval_one(curve%ubegin, ierr)
        print '(a,2f10.5)', '  s(u_begin) = ', y
        y = curve%eval_one(curve%uend, ierr)
        print '(a,2f10.5)', '  s(u_end)   = ', y

        call curve%destroy()
    end subroutine

    subroutine demo_constrained()
        integer, parameter :: m = 25
        integer, parameter :: idim = 2
        real(FP_REAL) :: pts(idim, m), u(m), du
        real(FP_REAL) :: ddx_begin(idim, 0:1), ddx_end(idim, 0:1)
        real(FP_REAL) :: y_begin(idim), y_end(idim)
        real(FP_REAL), allocatable :: yall(:,:)
        type(fitpack_constrained_curve) :: curve
        integer(FP_FLAG) :: ierr
        integer :: i

        print '(/,a)', '=== Constrained curve: prescribed endpoints ==='

        ! Generate a 2D curve
        du = pi / (m - 1)
        do i = 1, m
            u(i) = (i - 1) * du
            pts(1, i) = u(i) * cos(u(i))
            pts(2, i) = u(i) * sin(u(i))
        end do

        ! Prescribe position and 1st derivative at both endpoints
        ! At u=0: position = (0, 0), derivative = (1, 0)
        ddx_begin(:, 0) = [zero, zero]
        ddx_begin(:, 1) = [one, zero]

        ! At u=pi: position = (-pi, 0), derivative = (-1, -pi)
        ddx_end(:, 0) = [-pi, zero]
        ddx_end(:, 1) = [-one, -pi]

        ! Fit without constraints first
        ierr = curve%new_fit(pts, u=u, smoothing=real(m, FP_REAL))
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Unconstrained fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        y_begin = curve%eval_one(curve%ubegin, ierr)
        y_end   = curve%eval_one(curve%uend, ierr)
        print '(a)',        '  Unconstrained:'
        print '(a,2f10.5)', '    s(0)  = ', y_begin
        print '(a,2f10.5)', '    s(pi) = ', y_end

        ! Apply endpoint constraints and re-fit
        call curve%set_constraints(ddx_begin, ddx_end, ierr)
        ierr = curve%fit(smoothing=real(m, FP_REAL))
        if (.not. FITPACK_SUCCESS(ierr)) then
            print '(a,a)', '  Constrained fit failed: ', FITPACK_MESSAGE(ierr)
            return
        end if

        ! Check: position and derivatives at endpoints
        yall = curve%dfdx_all(curve%ubegin, ierr)
        print '(/,a)',      '  Constrained:'
        print '(a,2f10.5)', '    s(0)    = ', yall(:, 1)
        print '(a,2f10.5)', '    s''(0)   = ', yall(:, 2)

        yall = curve%dfdx_all(curve%uend, ierr)
        print '(a,2f10.5)', '    s(pi)   = ', yall(:, 1)
        print '(a,2f10.5)', '    s''(pi)  = ', yall(:, 2)

        call curve%destroy()
    end subroutine

end program example_parametric_curves
