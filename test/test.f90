program test
    use fitpack_tests
    use fitpack_test_data
    use fitpack_curve_tests
    use iso_fortran_env, only: output_unit
    implicit none

    integer :: passed,failed

    passed = 0
    failed = 0

    ! Run object-oriented tests
    call run_interface_tests

    if (failed>0) then
       print 1, passed, failed, 'interface'
       stop -1
    endif

    ! Run legacy tests
    call run_legacy_tests

    if (failed>0) then
       print 1, passed, failed, 'legacy'
       stop -1
    else
       print 2, passed
       stop 0
    endif


    1 format('[fitpack] there are ',i0,' passed, ',i0,' failed ',a,' tests.')
    2 format(//'[fitpack] SUCCESS! ',i0,' tests passed, none failed.')


    contains


    subroutine run_interface_tests()

        ! Sine function interpolant: test f(x) and df/dx
        call add_test(test_sine_fit())
        call add_test(test_periodic_fit())
        call add_test(test_parametric_fit())
        call add_test(test_closed_fit())
        call add_test(test_polar_fit())
        call add_test(test_sphere_fit())

    end subroutine run_interface_tests


    ! Todo: replace interactive mode
    subroutine run_legacy_tests

        integer :: itest

        ! Perform all legacy tests
        do itest = 1,29
           call add_test(perform_legacy_test(itest))
        end do

    end subroutine run_legacy_tests

    ! Perform FITPACK's i-th legacy test case
    logical function perform_legacy_test(itest,iunit) result(success)
        integer, intent(in) :: itest
        integer, optional, intent(in) :: iunit

        integer :: useUnit

        if (present(iunit)) then
            useUnit = iunit
        else
            useUnit = output_unit
        end if

        success = .true.

        print "(///)"

        select case (itest)
            case (1);  success = mnbisp(iunit)
            case (2);  success = mncloc(iunit)
            case (3);  success = mncoco(iunit)
            case (4);  success = mnconc(iunit)
            case (5);  success = mncosp(iunit)
            case (6);  success = mncual(iunit)
            case (7);  success = mncurf(iunit)
            case (8);  success = mnfour(iunit)
            case (9);  success = mnist (iunit)
            case (10); success = mnpade(iunit)
            case (11); success = mnparc(iunit)
            case (12); success = mnperc(iunit)
            case (13); success = mnpogr(dapogr,iunit)
            case (14); success = mnpola(dapola,iunit)
            case (15); success = mnprof(iunit)
            case (16); success = mnregr(daregr_x,daregr_y,daregr_z,iunit)
            case (17); success = mnspal(iunit)
            case (18); success = mnspde(iunit)
            case (19); success = mnspev(iunit)
            case (20); success = mnsphe(dasphe,iunit)
            case (21); success = mnspin(iunit)
            case (22); success = mnspro(iunit)
            case (23); success = mnsuev(iunit)
            case (24); success = mnsurf(dasurf_xyz,dasurf_delta)
            case (25); success = mncuev(iunit)
            case (26); success = mndbin(iunit)
            case (27); success = mnevpo(iunit)
            case (28); success = mnpasu(dapasu,iunit)
            case (29); success = mnspgr(daspgr_u,daspgr_v,daspgr_r,iunit)
            case default;
                 write(useUnit,*) '[fitpack] invalid test ID: try 1-29'
                 success = .false.
        end select

    end function perform_legacy_test

    subroutine add_test(success)
        logical, intent(in) :: success
        if (success) then
            passed = passed+1
        else
            failed = failed+1
        end if
    end subroutine

end program test
