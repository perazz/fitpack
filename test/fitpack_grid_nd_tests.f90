! **************************************************************************************************
!   FITPACK N-D gridded core — bit-for-bit equivalence gate (slice 1)
!
!   Asserts that the dimensionalized *_nd gridded routines reproduce their 2-D originals EXACTLY
!   (bit-for-bit) when run at dims=2, over the standard 11x11 `daregr` battery across all iopt
!   modes and spline orders. This is the oracle that licenses extending the *_nd path to dims>2.
!
!   Gates:
!     A. fpbisp_nd  vs bispev  (gridded evaluation)
!     D. regrid_nd  vs regrid  (gridded fit; transitively covers fpregr_nd and fpgrre_nd, whose
!                               outputs tx,ty,c,fp,nx,ny are exactly what regrid_nd returns)
!     D'. regrid_nd vs regrid on a NON-SQUARE grid (mx /= my) — guards the z = (ny,nx) y-fast data
!                               convention so an x<->y axis/size swap cannot hide behind mx==my
!
!   See todo/fitpack_nd_grids.md (slice 1) and Dierckx, Ch. 5 §5.4 (pp. 98-103).
! **************************************************************************************************
module fitpack_grid_nd_tests
    use fitpack_core
    use fitpack_test_data, only: daregr_x,daregr_y,daregr_z
    use iso_fortran_env, only: output_unit
    implicit none
    private

    public :: test_nd_grid_equivalence

    !> regrid workspace sizing (generous; mirrors legacy mnregr which uses nxest=nyest=17)
    integer(FP_SIZE), parameter :: NXEST = 17, NYEST = 17
    integer(FP_SIZE), parameter :: LWRK  = 2000, KWRK = 200, NCMAX = 400

    !> A fitting case: spline order and smoothing for one regrid call
    type :: grid_case
        integer(FP_SIZE) :: kx,ky
        real(FP_REAL)    :: s
        character(len=24):: label
    end type

    contains

    !> @brief Bit-for-bit equivalence gate for the N-D gridded core at dims=2.
    logical function test_nd_grid_equivalence(iunit) result(success)
        integer, optional, intent(in) :: iunit
        integer :: useUnit

        useUnit = output_unit
        if (present(iunit)) useUnit = iunit

        success = .true.

        ! Gate D: regrid_nd vs regrid (backbone; transitively covers fpregr_nd + fpgrre_nd)
        if (success) success = gate_regrid_nd(useUnit)

        ! Gate D': same, on a NON-SQUARE grid (mx /= my) so an x<->y data/size swap cannot hide
        if (success) success = gate_regrid_nd_nonsquare(useUnit)

        ! Gate A: fpbisp_nd vs bispev
        if (success) success = gate_fpbisp_nd(useUnit)

        if (success) write(useUnit,'(a)') '[test_nd_grid_equivalence] all dims=2 gates bit-for-bit OK'

    end function test_nd_grid_equivalence

    !> @brief The 2-D fitting battery (orders + smoothing) shared by all gates.
    pure function cases() result(cs)
        type(grid_case) :: cs(3)
        cs(1) = grid_case(3,3,0.22_FP_REAL,'cubic smoothing s=0.22')
        cs(2) = grid_case(3,3,0.0_FP_REAL ,'cubic interpolation s=0')
        cs(3) = grid_case(5,5,0.2_FP_REAL ,'quintic smoothing s=0.2')
    end function cases

    !> @brief Run regrid for one case, returning the fitted knots/coefficients.
    subroutine fit_regrid(gc,nx,tx,ny,ty,c,fp,ier)
        type(grid_case),  intent(in)  :: gc
        integer(FP_SIZE), intent(out) :: nx,ny
        real(FP_REAL),    intent(out) :: tx(NXEST),ty(NYEST),c(NCMAX),fp
        integer(FP_FLAG), intent(out) :: ier

        integer(FP_SIZE) :: mx,my
        real(FP_REAL)    :: xb,xe,yb,ye,wrk(LWRK)
        integer(FP_SIZE) :: iwrk(KWRK)
        real(FP_REAL)    :: zin(size(daregr_x)*size(daregr_y))

        mx = size(daregr_x,kind=FP_SIZE)
        my = size(daregr_y,kind=FP_SIZE)
        xb = daregr_x(1); xe = daregr_x(mx)
        yb = daregr_y(1); ye = daregr_y(my)
        zin = reshape(daregr_z,[mx*my])   ! same sequence association legacy mnregr relies on

        call regrid(IOPT_NEW_SMOOTHING,mx,daregr_x,my,daregr_y,zin,xb,xe,yb,ye, &
                    gc%kx,gc%ky,gc%s,NXEST,NYEST,nx,tx,ny,ty,c,fp, &
                    wrk,LWRK,iwrk,KWRK,ier)
    end subroutine fit_regrid

    !> @brief Gate A — fpbisp_nd(dims=2) must equal bispev's gridded evaluation bit-for-bit.
    logical function gate_fpbisp_nd(useUnit) result(ok)
        integer, intent(in) :: useUnit

        type(grid_case)  :: cs(3),gc
        integer(FP_SIZE) :: mx,my,nx,ny,nc,kx,ky,mz,maxn,maxm,maxk1
        integer(FP_FLAG) :: ier
        real(FP_REAL)    :: tx(NXEST),ty(NYEST),c(NCMAX),fp
        ! reference evaluation (bispev)
        real(FP_REAL)    :: zref(size(daregr_x)*size(daregr_y)),ewrk(NCMAX)
        integer(FP_SIZE) :: eiwrk(KWRK)
        ! fpbisp_nd marshalling
        real(FP_REAL)    :: ztest(size(daregr_x)*size(daregr_y))
        real(FP_REAL),    allocatable :: t(:,:),xg(:,:),w(:,:,:)
        integer(FP_SIZE), allocatable :: lidx(:,:)
        integer(FP_SIZE) :: n2(2),k2(2),m2(2)
        integer          :: icase
        real(FP_REAL)    :: maxdiff

        ok = .true.
        cs = cases()
        mx = size(daregr_x,kind=FP_SIZE)
        my = size(daregr_y,kind=FP_SIZE)
        mz = mx*my

        do icase=1,size(cs)
            gc = cs(icase)
            kx = gc%kx; ky = gc%ky

            call fit_regrid(gc,nx,tx,ny,ty,c,fp,ier)
            if (.not.FITPACK_SUCCESS(ier)) then
                ok = .false.
                write(useUnit,1000) 'fit',trim(gc%label),FITPACK_MESSAGE(ier)
                return
            end if

            nc = (nx-kx-1)*(ny-ky-1)

            ! reference: bispev (which internally calls fpbisp)
            call bispev(tx,nx,ty,ny,c,kx,ky,daregr_x,mx,daregr_y,my,zref, &
                        ewrk,size(ewrk,kind=FP_SIZE),eiwrk,size(eiwrk,kind=FP_SIZE),ier)
            if (.not.FITPACK_SUCCESS(ier)) then
                ok = .false.
                write(useUnit,1000) 'bispev',trim(gc%label),FITPACK_MESSAGE(ier)
                return
            end if

            ! test: fpbisp_nd at dims=2
            maxn  = max(nx,ny); maxm = max(mx,my); maxk1 = max(kx,ky)+1
            n2 = [nx,ny]; k2 = [kx,ky]; m2 = [mx,my]
            if (allocated(t)) deallocate(t,xg,w,lidx)
            allocate(t(maxn,2),xg(maxm,2),w(maxm,maxk1,2),lidx(maxm,2))
            t = zero
            t(1:nx,1) = tx(1:nx);      t(1:ny,2) = ty(1:ny)
            xg(1:mx,1) = daregr_x;     xg(1:my,2) = daregr_y

            call fpbisp_nd(2_FP_DIM,t,n2,c(1:nc),k2,xg,m2,ztest(1:mz),w,lidx)

            maxdiff = maxval(abs(zref(1:mz)-ztest(1:mz)))
            if (.not.all(zref(1:mz)==ztest(1:mz))) then
                ok = .false.
                write(useUnit,1100) trim(gc%label),maxdiff
                return
            end if
        end do

        write(useUnit,'(a)') '[gate A] fpbisp_nd == bispev (bit-for-bit, all cases)'

        1000 format('[gate A] ',a,' failed for ',a,': ',a)
        1100 format('[gate A] fpbisp_nd /= bispev for ',a,'  (max |diff| = ',1pe12.3,')')
    end function gate_fpbisp_nd

    !> @brief Gate D — regrid_nd(dims=2) must equal regrid bit-for-bit across all iopt modes.
    !!
    !! Runs regrid (reference) and regrid_nd (test) through the same iopt sequence and asserts
    !! nx,ny,tx,ty,c,fp and ier match to the last bit. The two are independent invocations with
    !! their own preserved workspace; iopt=1 cases continue from the previous knot set on both
    !! sides. A green result here certifies fpregr_nd and fpgrre_nd as well (their outputs ARE
    !! what regrid_nd returns).
    logical function gate_regrid_nd(useUnit) result(ok)
        integer, intent(in) :: useUnit

        ! reference (regrid) state
        integer(FP_SIZE) :: mx,my,nxr,nyr
        real(FP_REAL)    :: txr(NXEST),tyr(NYEST),cr(NCMAX),fpr,wrkr(LWRK)
        integer(FP_SIZE) :: iwrkr(KWRK)
        integer(FP_FLAG) :: ier_r
        ! test (regrid_nd) state
        integer(FP_SIZE) :: ntst(2),kk(2),m2(2),nest2(2)
        real(FP_REAL)    :: ttst(NXEST,2),ctst(NCMAX),fptst,wrkt(LWRK),xgt(size(daregr_x),2)
        integer(FP_SIZE) :: iwrkt(KWRK)
        integer(FP_FLAG) :: ier_t
        ! shared
        real(FP_REAL)    :: zin(size(daregr_x)*size(daregr_y)),xb,xe,yb,ye,s,lo(2),hi(2)
        integer(FP_SIZE) :: kx,ky,iopt,icase,i
        character(len=40):: lbl

        ok = .true.
        mx = size(daregr_x,kind=FP_SIZE); my = size(daregr_y,kind=FP_SIZE)
        xb = daregr_x(1); xe = daregr_x(mx); yb = daregr_y(1); ye = daregr_y(my)
        zin = reshape(daregr_z,[mx*my])

        ! regrid_nd marshalling that is constant across calls
        m2 = [mx,my]; nest2 = [NXEST,NYEST]; lo = [xb,yb]; hi = [xe,ye]
        xgt = zero; xgt(1:mx,1) = daregr_x; xgt(1:my,2) = daregr_y
        wrkr = zero; iwrkr = 0; wrkt = zero; iwrkt = 0

        ! ---- Phase A: cubic, iopt=0 then iopt=1 continuation chain (incl. s=0 interpolation) ----
        kx = 3; ky = 3; kk = [kx,ky]
        do icase=1,4
            select case (icase)
               case (1); iopt = 0; s = 10.0_FP_REAL
               case (2); iopt = 1; s = 0.22_FP_REAL
               case (3); iopt = 1; s = 0.1_FP_REAL
               case (4); iopt = 1; s = 0.0_FP_REAL
            end select
            call regrid(iopt,mx,daregr_x,my,daregr_y,zin,xb,xe,yb,ye,kx,ky,s,NXEST,NYEST, &
                        nxr,txr,nyr,tyr,cr,fpr,wrkr,LWRK,iwrkr,KWRK,ier_r)
            call regrid_nd(iopt,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                           ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
            write(lbl,'(a,i0)') 'cubic iopt-chain #',icase
            if (.not.pair_ok(useUnit,trim(lbl),kk,nxr,nyr,txr,tyr,cr,fpr,ier_r, &
                             ntst,ttst,ctst,fptst,ier_t)) then
                ok = .false.; return
            end if
        end do

        ! ---- Phase B: quintic, fresh iopt=0 ----
        kx = 5; ky = 5; kk = [kx,ky]; iopt = 0; s = 0.2_FP_REAL
        call regrid(iopt,mx,daregr_x,my,daregr_y,zin,xb,xe,yb,ye,kx,ky,s,NXEST,NYEST, &
                    nxr,txr,nyr,tyr,cr,fpr,wrkr,LWRK,iwrkr,KWRK,ier_r)
        call regrid_nd(iopt,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                       ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
        if (.not.pair_ok(useUnit,'quintic fresh iopt=0',kk,nxr,nyr,txr,tyr,cr,fpr,ier_r, &
                         ntst,ttst,ctst,fptst,ier_t)) then
            ok = .false.; return
        end if

        ! ---- Phase C: least-squares on given knots, iopt=-1 ----
        kx = 3; ky = 3; kk = [kx,ky]; s = 0.0_FP_REAL
        nxr = 11; nyr = 11; ntst = [11,11]
        txr = zero; txr(kx+2:kx+4) = [(half*(i-2),i=1,3)]; tyr = zero; tyr(ky+2:ky+4) = txr(kx+2:kx+4)
        ttst = zero
        ttst(kx+2:kx+4,1) = txr(kx+2:kx+4); ttst(ky+2:ky+4,2) = tyr(ky+2:ky+4)
        call regrid(-1_FP_FLAG,mx,daregr_x,my,daregr_y,zin,xb,xe,yb,ye,kx,ky,s,NXEST,NYEST, &
                    nxr,txr,nyr,tyr,cr,fpr,wrkr,LWRK,iwrkr,KWRK,ier_r)
        call regrid_nd(-1_FP_FLAG,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                       ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
        if (.not.pair_ok(useUnit,'lsq given knots iopt=-1',kk,nxr,nyr,txr,tyr,cr,fpr,ier_r, &
                         ntst,ttst,ctst,fptst,ier_t)) then
            ok = .false.; return
        end if

        write(useUnit,'(a)') '[gate D] regrid_nd == regrid (bit-for-bit; iopt=-1,0,1; cubic+quintic)'
    end function gate_regrid_nd

    !> @brief Gate D' — regrid_nd vs regrid on a NON-SQUARE grid (mx /= my).
    !!
    !! The square 11x11 daregr battery cannot expose an x<->y swap in the data tensor or in the
    !! per-axis sizing: with mx==my both indexings address the same elements. A non-square grid
    !! makes any such swap a hard bit-for-bit failure, so this directly guards the z = (ny,nx)
    !! (y-fast) storage convention as the *_nd path heads toward dims>2. The synthetic data is
    !! deliberately NON-symmetric under x<->y (different centres/weights) so even a pure transpose
    !! on a hypothetical square grid would be caught.
    logical function gate_regrid_nd_nonsquare(useUnit) result(ok)
        integer, intent(in) :: useUnit

        integer(FP_SIZE), parameter :: MX = 9, MY = 13
        real(FP_REAL)    :: xg1(MX),yg1(MY),zin(MX*MY)
        ! reference (regrid)
        integer(FP_SIZE) :: nxr,nyr
        real(FP_REAL)    :: txr(NXEST),tyr(NYEST),cr(NCMAX),fpr,wrkr(LWRK)
        integer(FP_SIZE) :: iwrkr(KWRK)
        integer(FP_FLAG) :: ier_r
        ! test (regrid_nd)
        integer(FP_SIZE) :: ntst(2),kk(2),m2(2),nest2(2)
        real(FP_REAL)    :: ttst(NXEST,2),ctst(NCMAX),fptst,wrkt(LWRK),xgt(MY,2)
        integer(FP_SIZE) :: iwrkt(KWRK)
        integer(FP_FLAG) :: ier_t
        ! shared
        real(FP_REAL)    :: xb,xe,yb,ye,s,lo(2),hi(2)
        integer(FP_SIZE) :: i,j,iopt,icase
        character(len=40):: lbl

        ok = .true.

        ! strictly increasing grids on [0,1]; distinct sizes mx=9, my=13
        do i=1,MX; xg1(i) = real(i-1,FP_REAL)/real(MX-1,FP_REAL); end do
        do j=1,MY; yg1(j) = real(j-1,FP_REAL)/real(MY-1,FP_REAL); end do

        ! smooth, NON-symmetric data, stored y-fast: zin((i-1)*MY+j) = f(x_i,y_j)
        do i=1,MX
           do j=1,MY
              zin((i-1)*MY+j) = (xg1(i)-0.5_FP_REAL)**2 + 2.0_FP_REAL*(yg1(j)-0.3_FP_REAL)**2 &
                              + xg1(i)*yg1(j)
           end do
        end do

        xb = xg1(1); xe = xg1(MX); yb = yg1(1); ye = yg1(MY)
        m2 = [MX,MY]; nest2 = [NXEST,NYEST]; lo = [xb,yb]; hi = [xe,ye]
        xgt = zero; xgt(1:MX,1) = xg1; xgt(1:MY,2) = yg1
        wrkr = zero; iwrkr = 0; wrkt = zero; iwrkt = 0
        kk = [3_FP_SIZE,3_FP_SIZE]

        do icase=1,2
            select case (icase)
               case (1); iopt = 0; s = 0.05_FP_REAL   ! smoothing
               case (2); iopt = 0; s = 0.0_FP_REAL    ! interpolation
            end select
            call regrid(iopt,MX,xg1,MY,yg1,zin,xb,xe,yb,ye,kk(1),kk(2),s,NXEST,NYEST, &
                        nxr,txr,nyr,tyr,cr,fpr,wrkr,LWRK,iwrkr,KWRK,ier_r)
            call regrid_nd(iopt,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                           ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
            write(lbl,'(a,i0)') 'non-square 9x13 #',icase
            if (.not.pair_ok(useUnit,trim(lbl),kk,nxr,nyr,txr,tyr,cr,fpr,ier_r, &
                             ntst,ttst,ctst,fptst,ier_t)) then
                ok = .false.; return
            end if
        end do

        write(useUnit,'(a)') '[gate D''] regrid_nd == regrid (bit-for-bit; non-square 9x13)'
    end function gate_regrid_nd_nonsquare

    !> @brief Compare a regrid vs regrid_nd result pair bit-for-bit; report on mismatch.
    logical function pair_ok(useUnit,label,k,nxr,nyr,txr,tyr,cr,fpr,ier_r, &
                             ntst,ttst,ctst,fptst,ier_t) result(ok)
        integer,          intent(in) :: useUnit
        character(*),     intent(in) :: label
        integer(FP_SIZE), intent(in) :: k(2),nxr,nyr,ntst(2)
        real(FP_REAL),    intent(in) :: txr(:),tyr(:),cr(:),fpr,ttst(:,:),ctst(:),fptst
        integer(FP_FLAG), intent(in) :: ier_r,ier_t
        integer(FP_SIZE) :: nc

        ok = .true.

        if (ier_r/=ier_t) then
            write(useUnit,2000) trim(label),'ier',ier_r,ier_t; ok = .false.; return
        end if
        if (nxr/=ntst(1) .or. nyr/=ntst(2)) then
            write(useUnit,2100) trim(label),nxr,ntst(1),nyr,ntst(2); ok = .false.; return
        end if

        nc = (nxr-k(1)-1)*(nyr-k(2)-1)

        if (.not.all(txr(1:nxr)==ttst(1:ntst(1),1))) ok = .false.
        if (.not.all(tyr(1:nyr)==ttst(1:ntst(2),2))) ok = .false.
        if (.not.all(cr(1:nc)==ctst(1:nc)))          ok = .false.
        if (fpr/=fptst)                              ok = .false.

        if (.not.ok) then
            write(useUnit,2200) trim(label), &
                maxval(abs(txr(1:nxr)-ttst(1:ntst(1),1))), &
                maxval(abs(tyr(1:nyr)-ttst(1:ntst(2),2))), &
                maxval(abs(cr(1:nc)-ctst(1:nc))), abs(fpr-fptst)
        end if

        2000 format('[gate D] ',a,': ',a,' mismatch (ref=',i0,' nd=',i0,')')
        2100 format('[gate D] ',a,': knot count mismatch nx ',i0,'/',i0,' ny ',i0,'/',i0)
        2200 format('[gate D] ',a,' NOT bit-for-bit: max|dtx|=',1pe10.2,' max|dty|=',e10.2, &
                    ' max|dc|=',e10.2,' |dfp|=',e10.2)
    end function pair_ok

end module fitpack_grid_nd_tests
