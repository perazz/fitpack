program test
    use fitpack_tests

    ! Get test ID from the useri
    write(*,'(A)',advance='no') 'Enter fitpack test ID [1:29] > '
    read *, itest

    select case (itest)
        case (1); call mnbisp
        case (2); call mncloc
        case (3); call mncoco
        case (4); call mnconc
        case (5); call mncosp
        case (6); call mncual
        case (7); call mncurf
        case (8); call mnfour
        case (9); call mnist
        case (10); call mnpade
        case (11); call mnparc
        case (12); call mnperc
        case (13); call mnpogr
        case (14); call mnpola
        case (15); call mnprof
        case (16); call mnregr
        case (17); call mnspal
        case (18); call mnspde
        case (19); call mnspev
        case (20); call mnsphe
        case (21); call mnspin
        case (22); call mnspro
        case (23); call mnsuev
        case (24); call mnsurf
        case (25); call mncuev
        case (26); call mndbin
        case (27); call mnevpo
        case (28); call mnpasu
        case (29); call mnspgr
        case default; stop 'invalid test ID'
    end select



end program test
