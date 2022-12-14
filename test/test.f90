program test
    use fitpack_tests
    use fitpack_test_data
    use fitpack_curve_tests
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
    2 format(//'[fitpack] SUCCESS! ',i0,' test passed, none failed.')


    contains


    subroutine run_interface_tests()


        ! Sine function interpolant: test f(x) and df/dx
        call add_test(test_sine_fit())

    end subroutine run_interface_tests


    ! Todo: replace interactive mode
    subroutine run_legacy_tests

        integer :: itest

        ! Perform all le
        do itest = 24,24
           call add_test(perform_legacy_test(itest))
        end do



    end subroutine run_legacy_tests

    ! Perform FITPACK's i-th legacy test case
    logical function perform_legacy_test(itest) result(success)
        integer, intent(in) :: itest

        success = .true.

        select case (itest)
            case (1);  call mnbisp
            case (2);  call mncloc
            case (3);  call mncoco
            case (4);  call mnconc
            case (5);  call mncosp
            case (6);  call mncual
            case (7);  success = mncurf()
            case (8);  call mnfour
            case (9);  call mnist
            case (10); call mnpade
            case (11); call mnparc
            case (12); call mnperc
            case (13); call mnpogr(dapogr)
            case (14); call mnpola(dapola)
            case (15); call mnprof
            case (16); call mnregr(daregr_x,daregr_y,daregr_z)
            case (17); call mnspal
            case (18); call mnspde
            case (19); call mnspev
            case (20); call mnsphe(dasphe)
            case (21); call mnspin
            case (22); call mnspro
            case (23); call mnsuev
            case (24); call mnsurf(dasurf_xyz,dasurf_delta)
            case (25); call mncuev
            case (26); call mndbin
            case (27); call mnevpo
            case (28); call mnpasu(dapasu)
            case (29); call mnspgr(daspgr_u,daspgr_v,daspgr_r)
            case default;
                 print *, 'invalid test ID'
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
