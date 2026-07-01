! **************************************************************************************************
!   FITPACK N-D gridded core — equivalence gate (slice 1) + first dims>2 path (slice 2)
!
!   Gates A/D/D' assert that the dimensionalized *_nd gridded routines reproduce their 2-D originals
!   EXACTLY (bit-for-bit) when run at dims=2, over the standard 11x11 `daregr` battery across all iopt
!   modes and spline orders — the oracle that licenses extending the *_nd path to dims>2. Gate E then
!   exercises the first genuine dims>2 path (the eval contraction at dims=3) against independent refs.
!
!   Gates:
!     A. fpndsp  vs bispev  (gridded evaluation, dims=2 bit-for-bit)
!     D. regrid_nd  vs regrid  (gridded fit; transitively covers fpregr_nd and fpgrre_nd, whose
!                               outputs tx,ty,c,fp,nx,ny are exactly what regrid_nd returns)
!     D'. regrid_nd vs regrid on a NON-SQUARE grid (mx /= my) — guards the z = (ny,nx) y-fast data
!                               convention so an x<->y axis/size swap cannot hide behind mx==my
!     E. fpndsp(dims=3) vs independent references (slice 2) — separable coeffs vs a product of 1-D
!                               splev, and non-separable coeffs vs the dims=2 path combined along z
!     F. regrid_nd(dims=3) least-squares FIT on given knots (slice 3) — the generalized fpgrre_nd solve
!                               recovers an in-space polynomial; checked vs the closed-form (fpndsp)
!
!   See todo/fitpack_nd_grids.md (slices 1–3) and Dierckx, Ch. 5 §5.4 (pp. 98-103).
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

        ! Gate A: fpndsp vs bispev
        if (success) success = gate_fpndsp(useUnit)

        ! Gate E: fpndsp at dims=3 vs two independent references (the first genuine dims>2 path)
        if (success) success = gate_fpndsp_3d(useUnit)

        ! Gate F: regrid_nd at dims=3 — the first genuine dims>2 FIT (generalized fpgrre_nd solve)
        if (success) success = gate_regrid_nd_3d(useUnit)

        ! Gate G: regrid_nd at dims=3 SMOOTHING (s>0) — the generalized knot-direction arbiter
        if (success) success = gate_regrid_nd_3d_smoothing(useUnit)

        if (success) write(useUnit,'(a)') '[test_nd_grid_equivalence] all N-D gridded gates OK'

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

    !> @brief Gate A — the three 2-D entry points (fpndsp kernel, ndspev wrapper, bispev adapter)
    !!        must agree bit-for-bit, and all must match an independent splev separation-of-variables
    !!        oracle. Since slice 6 every 2-D path routes through fpndsp, so the oracle (not the mutual
    !!        identity) is what guards correctness.
    logical function gate_fpndsp(useUnit) result(ok)
        integer, intent(in) :: useUnit

        type(grid_case)  :: cs(3),gc
        integer(FP_SIZE) :: mx,my,nx,ny,nc,kx,ky,mz,maxn,maxm,maxk1
        integer(FP_FLAG) :: ier
        real(FP_REAL)    :: tx(NXEST),ty(NYEST),c(NCMAX),fp
        ! reference evaluation (bispev)
        real(FP_REAL)    :: zref(size(daregr_x)*size(daregr_y)),ewrk(NCMAX)
        integer(FP_SIZE) :: eiwrk(KWRK)
        ! fpndsp marshalling
        real(FP_REAL)    :: ztest(size(daregr_x)*size(daregr_y))
        real(FP_REAL)    :: ztest2(size(daregr_x)*size(daregr_y))
        real(FP_REAL),    allocatable :: t(:,:),xg(:,:),w(:,:,:),nwrk(:)
        integer(FP_SIZE), allocatable :: lidx(:,:),niwrk(:)
        integer(FP_SIZE) :: n2(2),k2(2),m2(2)
        integer          :: icase
        real(FP_REAL)    :: maxdiff
        ! independent oracle (separation of variables via the 1-D splev)
        real(FP_REAL)    :: zoracle(size(daregr_x)*size(daregr_y))
        real(FP_REAL)    :: crow(NYEST),ccol(NXEST),drow(size(daregr_y)),scol(size(daregr_x))
        real(FP_REAL)    :: dmat(NXEST,size(daregr_y))
        integer(FP_SIZE) :: a,ii,jj,nkx1,nky1

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

            ! test: fpndsp at dims=2
            maxn  = max(nx,ny); maxm = max(mx,my); maxk1 = max(kx,ky)+1
            n2 = [nx,ny]; k2 = [kx,ky]; m2 = [mx,my]
            if (allocated(t)) deallocate(t,xg,w,lidx)
            allocate(t(maxn,2),xg(maxm,2),w(maxm,maxk1,2),lidx(maxm,2))
            t = zero
            t(1:nx,1) = tx(1:nx);      t(1:ny,2) = ty(1:ny)
            xg(1:mx,1) = daregr_x;     xg(1:my,2) = daregr_y

            call fpndsp(2_FP_DIM,t,n2,c(1:nc),k2,xg,m2,ztest(1:mz),w,lidx)

            maxdiff = maxval(abs(zref(1:mz)-ztest(1:mz)))
            if (.not.all(zref(1:mz)==ztest(1:mz))) then
                ok = .false.
                write(useUnit,1100) trim(gc%label),maxdiff
                return
            end if

            ! test: public ndspev wrapper (flat workspace, pointer-carved) -- must match bispev too
            if (allocated(nwrk)) deallocate(nwrk,niwrk)
            allocate(nwrk(maxm*maxk1*2),niwrk(maxm*2))
            call ndspev(2_FP_DIM,t,n2,c(1:nc),k2,xg,m2,ztest2(1:mz), &
                        nwrk,size(nwrk,kind=FP_SIZE),niwrk,size(niwrk,kind=FP_SIZE),ier)
            if (.not.FITPACK_SUCCESS(ier) .or. .not.all(zref(1:mz)==ztest2(1:mz))) then
                ok = .false.
                write(useUnit,1200) trim(gc%label),maxval(abs(zref(1:mz)-ztest2(1:mz)))
                return
            end if

            ! independent oracle: separation of variables through the 1-D splev (now that every
            !   2-D path routes through fpndsp, this is the non-tautological correctness check).
            !   d(a,j)     = sum_b c_{a,b} M_b(y_j)   (ny-spline per x-coefficient-row a)
            !   s(x_i,y_j) = sum_a d(a,j) N_a(x_i)    (nx-spline per y-column j)
            nkx1 = nx-kx-1; nky1 = ny-ky-1
            do a=1,nkx1
                crow = zero
                crow(1:nky1) = c((a-1)*nky1+1:(a-1)*nky1+nky1)
                call splev(ty(1:ny),ny,crow(1:ny),ky,daregr_y,drow,my,OUTSIDE_NEAREST_BND,ier)
                dmat(a,1:my) = drow
            end do
            do jj=1,my
                ccol = zero
                ccol(1:nkx1) = dmat(1:nkx1,jj)
                call splev(tx(1:nx),nx,ccol(1:nx),kx,daregr_x,scol,mx,OUTSIDE_NEAREST_BND,ier)
                do ii=1,mx
                    zoracle((ii-1)*my+jj) = scol(ii)
                end do
            end do
            maxdiff = maxval(abs(ztest(1:mz)-zoracle(1:mz)))
            if (maxdiff > 1.0e-10_FP_REAL*max(one,maxval(abs(zoracle(1:mz))))) then
                ok = .false.
                write(useUnit,1300) trim(gc%label),maxdiff
                return
            end if
        end do

        write(useUnit,'(a)') '[gate A] fpndsp == ndspev == bispev, all == splev separation-of-variables'

        1000 format('[gate A] ',a,' failed for ',a,': ',a)
        1100 format('[gate A] fpndsp /= bispev for ',a,'  (max |diff| = ',1pe12.3,')')
        1200 format('[gate A] ndspev /= bispev for ',a,'  (max |diff| = ',1pe12.3,')')
        1300 format('[gate A] fpndsp /= splev oracle for ',a,'  (max |diff| = ',1pe12.3,')')
    end function gate_fpndsp

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
        zin = reshape(daregr_z,[mx*my])    ! flat row-major (x slowest, y fastest) = regrid_nd's z contract

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

    !> @brief Build a clamped knot vector on [0,1] for order k with nk1 coefficients (uniform interior).
    subroutine clamped_unit_knots(k,nk1,t,n)
        integer(FP_SIZE), intent(in)  :: k,nk1
        real(FP_REAL),    intent(out) :: t(:)
        integer(FP_SIZE), intent(out) :: n
        integer(FP_SIZE) :: k1,nint,i
        k1   = k+1
        n    = nk1+k1
        nint = n-2*k1                  ! number of interior knots
        t(1:k1) = zero
        do i=1,nint
           t(k1+i) = real(i,FP_REAL)/real(nint+1,FP_REAL)
        end do
        t(n-k:n) = one
    end subroutine clamped_unit_knots

    !> @brief Gate E — fpndsp at dims=3 vs two independent references (the first genuine dims>2 path).
    !!
    !! Reference 1 (separable): coefficients c(i,j,l)=a(i)*b(j)*d(l) make the 3-D tensor spline
    !! factorize to sx(x)*sy(y)*sz(z); the reference grid is the outer product of three verified 1-D
    !! splev evaluations.  Reference 2 (non-separable): each z-coefficient slab is evaluated over the
    !! (x,y) plane with the gate-A-verified fpndsp(dims=2), then combined along z with an explicit
    !! B-spline basis (fpbspl) — catching axis-ordering / separation bugs the symmetric case can hide.
    !! Distinct per-axis orders and grid sizes, so an axis swap cannot pass unnoticed.
    logical function gate_fpndsp_3d(useUnit) result(ok)
        integer, intent(in) :: useUnit

        ! distinct per-axis orders, coefficient counts and eval-grid sizes (never equal: no axis hides)
        integer(FP_SIZE), parameter :: kx=3,   ky=2,   kz=4
        integer(FP_SIZE), parameter :: nk1x=6, nk1y=5, nk1z=7
        integer(FP_SIZE), parameter :: mx=7,   my=5,   mz=9
        integer(FP_SIZE), parameter :: maxn=nk1z+kz+1, maxm=mz, maxk1=kz+1
        integer(FP_SIZE), parameter :: nc3=nk1x*nk1y*nk1z, mout3=mx*my*mz, mxy=mx*my
        real(FP_REAL),    parameter :: tol=1.0e-12_FP_REAL

        real(FP_REAL)    :: t3(maxn,3),xg3(maxm,3),c3(nc3),z3(mout3),w3(maxm,maxk1,3)
        integer(FP_SIZE) :: lidx3(maxm,3),n3(3),k3(3),m3(3)
        ! separable reference (outer product of three 1-D splev evaluations)
        real(FP_REAL)    :: av(maxn),bv(maxn),dv(maxn),sx(mx),sy(my),sz(mz),ref
        ! non-separable reference (verified dims=2 plane evals, combined along z by splev)
        real(FP_REAL)    :: v2(mxy,nk1z),z2(mxy),c2(nk1x*nk1y),w2(mx,kx+1,2),vz(mz),cz(maxn)
        real(FP_REAL)    :: t2(maxn,2),xg2(mx,2),maxdiff
        integer(FP_SIZE) :: lidx2(mx,2),n2(2),k2(2),m2(2)
        integer(FP_SIZE) :: ix,iy,iz,i,nx,ny,nz,l,ic,iout
        integer(FP_FLAG) :: ier

        ok = .true.
        t3 = zero; xg3 = zero; xg2 = zero
        k3 = [kx,ky,kz]; m3 = [mx,my,mz]

        ! clamped knot vectors on [0,1]; n_a = nk1_a+k_a+1 (10,8,12)
        call clamped_unit_knots(kx,nk1x,t3(:,1),nx)
        call clamped_unit_knots(ky,nk1y,t3(:,2),ny)
        call clamped_unit_knots(kz,nk1z,t3(:,3),nz)
        n3 = [nx,ny,nz]

        ! eval grids: strictly interior, distinct sizes
        do i=1,mx; xg3(i,1) = (real(i,FP_REAL)-half)/real(mx,FP_REAL); end do
        do i=1,my; xg3(i,2) = (real(i,FP_REAL)-half)/real(my,FP_REAL); end do
        do i=1,mz; xg3(i,3) = (real(i,FP_REAL)-half)/real(mz,FP_REAL); end do

        ! 1-D coefficient vectors (moderate, deterministic, non-trivial)
        av = zero; bv = zero; dv = zero
        do i=1,nk1x; av(i) = sin(0.7_FP_REAL*i)+1.5_FP_REAL; end do
        do i=1,nk1y; bv(i) = cos(0.4_FP_REAL*i)+1.3_FP_REAL; end do
        do i=1,nk1z; dv(i) = 0.5_FP_REAL+0.13_FP_REAL*i;     end do

        ! coefficient/output flat indices are row-major: x slowest, z fastest
        ! c3((ix-1)*nk1y*nk1z+(iy-1)*nk1z+iz),  z3((ix-1)*my*mz+(iy-1)*mz+iz)

        ! ===== Reference 1: separable coefficients c(i,j,l)=a(i)*b(j)*d(l) =====
        do ix=1,nk1x
           do iy=1,nk1y
              do iz=1,nk1z
                 c3((ix-1)*nk1y*nk1z+(iy-1)*nk1z+iz) = av(ix)*bv(iy)*dv(iz)
              end do
           end do
        end do

        call fpndsp(3_FP_DIM,t3,n3,c3,k3,xg3,m3,z3,w3,lidx3)

        call splev(t3(1:nx,1),nx,av(1:nx),kx,xg3(1:mx,1),sx,mx,OUTSIDE_EXTRAPOLATE,ier)
        call splev(t3(1:ny,2),ny,bv(1:ny),ky,xg3(1:my,2),sy,my,OUTSIDE_EXTRAPOLATE,ier)
        call splev(t3(1:nz,3),nz,dv(1:nz),kz,xg3(1:mz,3),sz,mz,OUTSIDE_EXTRAPOLATE,ier)

        maxdiff = zero
        do ix=1,mx
           do iy=1,my
              do iz=1,mz
                 ref  = sx(ix)*sy(iy)*sz(iz)
                 iout = (ix-1)*my*mz+(iy-1)*mz+iz
                 maxdiff = max(maxdiff, abs(z3(iout)-ref))
              end do
           end do
        end do
        if (maxdiff>=tol) then
           ok = .false.; write(useUnit,3000) 'separable',maxdiff; return
        end if

        ! ===== Reference 2: non-separable coefficients =====
        do ix=1,nk1x
           do iy=1,nk1y
              do iz=1,nk1z
                 c3((ix-1)*nk1y*nk1z+(iy-1)*nk1z+iz) = sin(1.1_FP_REAL*ix+0.7_FP_REAL*iy+0.5_FP_REAL*iz)
              end do
           end do
        end do

        call fpndsp(3_FP_DIM,t3,n3,c3,k3,xg3,m3,z3,w3,lidx3)

        ! evaluate the (x,y) plane of each z-coefficient slab with the gate-A-verified dims=2 path
        n2 = [nx,ny]; k2 = [kx,ky]; m2 = [mx,my]
        t2(:,1) = t3(:,1); t2(:,2) = t3(:,2)
        xg2(1:mx,1) = xg3(1:mx,1); xg2(1:my,2) = xg3(1:my,2)
        do l=1,nk1z
           do ix=1,nk1x
              do iy=1,nk1y
                 c2((ix-1)*nk1y+iy) = c3((ix-1)*nk1y*nk1z+(iy-1)*nk1z+l)
              end do
           end do
           call fpndsp(2_FP_DIM,t2,n2,c2,k2,xg2,m2,z2,w2,lidx2)
           v2(:,l) = z2
        end do

        ! combine along z: for each (x,y) the nk1z slab values ARE the z-spline coefficients -> splev
        maxdiff = zero
        do ix=1,mx
           do iy=1,my
              ic = (ix-1)*my+iy
              cz = zero; cz(1:nk1z) = v2(ic,1:nk1z)
              call splev(t3(1:nz,3),nz,cz(1:nz),kz,xg3(1:mz,3),vz,mz,OUTSIDE_EXTRAPOLATE,ier)
              do iz=1,mz
                 iout = (ix-1)*my*mz+(iy-1)*mz+iz
                 maxdiff = max(maxdiff, abs(z3(iout)-vz(iz)))
              end do
           end do
        end do
        if (maxdiff>=tol) then
           ok = .false.; write(useUnit,3000) 'non-separable',maxdiff; return
        end if

        write(useUnit,'(a)') '[gate E] fpndsp(dims=3) == separable & nested-2D references'

        3000 format('[gate E] fpndsp(dims=3) /= ',a,' reference  (max |diff| = ',1pe12.3,')')
    end function gate_fpndsp_3d

    !> @brief Gate F — regrid_nd at dims=3 (least-squares fit on given knots) vs a closed-form oracle.
    !!
    !! The first genuine dims>2 FIT path: the generalized fpgrre_nd solve (alternating-direction
    !! reduction + ping-pong buffer + back-substitution + residual contraction). Data is a separable
    !! in-space polynomial f(x,y,w)=Px(x)*Py(y)*Pw(w) with deg Px/Py/Pw = kx/ky/kw, which lies EXACTLY
    !! in the tensor B-spline space for any clamped knots. Fit via regrid_nd(dims=3, iopt=-1) on
    !! hand-built clamped knots with nk1(d)<m(d) (genuinely over-determined, so the residual path runs);
    !! the binary knot-direction arbiter is provably untouched on the iopt<0 path. Two independent
    !! checks: the residual fp must vanish, and the fitted spline -- evaluated by the gate-A/E-verified
    !! fpndsp(dims=3) at strictly-interior probe points off the data grid -- must reproduce the
    !! closed-form f. Distinct per-axis orders and grid sizes, so an axis swap cannot pass unnoticed.
    logical function gate_regrid_nd_3d(useUnit) result(ok)
        integer, intent(in) :: useUnit

        integer(FP_SIZE), parameter :: kx=3, ky=2, kw=4
        integer(FP_SIZE), parameter :: nk1x=6, nk1y=5, nk1w=7
        integer(FP_SIZE), parameter :: mx=8, my=6, mw=10
        integer(FP_SIZE), parameter :: nxk=nk1x+kx+1, nyk=nk1y+ky+1, nwk=nk1w+kw+1   ! 10,8,12
        integer(FP_SIZE), parameter :: maxn=nwk, maxk1=kw+1, nc3=nk1x*nk1y*nk1w, mz3=mx*my*mw
        integer(FP_SIZE), parameter :: mpx=5, mpy=4, mpw=6, maxmp=mpw   ! interior probe grid
        real(FP_REAL),    parameter :: tol=1.0e-10_FP_REAL

        real(FP_REAL)    :: t3(maxn,3),xg3(mw,3),c3(nc3),z3(mz3),fpf,lo(3),hi(3),s
        real(FP_REAL)    :: wrk(LWRK)
        integer(FP_SIZE) :: iwrk(KWRK),n3(3),k3(3),m3(3),nest3(3)
        integer(FP_FLAG) :: ier
        ! probe / closed-form oracle
        real(FP_REAL)    :: xgp(maxmp,3),zp(mpx*mpy*mpw),wp(maxmp,maxk1,3)
        integer(FP_SIZE) :: lidxp(maxmp,3),mp3(3)
        real(FP_REAL)    :: x,y,w,maxdiff
        integer(FP_SIZE) :: ix,iy,iw,i,iout

        ok = .true.
        t3 = zero; xg3 = zero; xgp = zero
        k3 = [kx,ky,kw]; m3 = [mx,my,mw]; lo = zero; hi = one

        ! clamped knot vectors on [0,1]; nest = exact knot counts (least-squares on given knots)
        call clamped_unit_knots(kx,nk1x,t3(:,1),n3(1))
        call clamped_unit_knots(ky,nk1y,t3(:,2),n3(2))
        call clamped_unit_knots(kw,nk1w,t3(:,3),n3(3))
        nest3 = n3

        ! strictly-increasing data grids on [0,1] (endpoints included)
        do i=1,mx; xg3(i,1) = real(i-1,FP_REAL)/real(mx-1,FP_REAL); end do
        do i=1,my; xg3(i,2) = real(i-1,FP_REAL)/real(my-1,FP_REAL); end do
        do i=1,mw; xg3(i,3) = real(i-1,FP_REAL)/real(mw-1,FP_REAL); end do

        ! data: separable in-space polynomial, row-major (x slowest, w fastest)
        do ix=1,mx
           x = xg3(ix,1)
           do iy=1,my
              y = xg3(iy,2)
              do iw=1,mw
                 w = xg3(iw,3)
                 iout = (ix-1)*my*mw + (iy-1)*mw + iw
                 z3(iout) = poly_x(x)*poly_y(y)*poly_w(w)
              end do
           end do
        end do

        ! least-squares fit on the given clamped knots (iopt=-1): the arbiter-free dims=3 solve
        s = zero; wrk = zero; iwrk = 0
        call regrid_nd(-1_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,s,nest3, &
                       n3,t3,c3,fpf,wrk,LWRK,iwrk,KWRK,ier)
        if (ier/=FITPACK_OK) then
           ok = .false.; write(useUnit,'(a,i0)') '[gate F] regrid_nd(dims=3) fit failed, ier=',ier; return
        end if

        ! check 1: data exactly in space -> least-squares residual must vanish
        if (fpf>=1.0e-9_FP_REAL) then
           ok = .false.; write(useUnit,3000) 'nonzero residual fp', fpf; return
        end if

        ! check 2: reconstruction at strictly-interior probe points off the data grid vs closed-form f
        mp3 = [mpx,mpy,mpw]
        do i=1,mpx; xgp(i,1) = (real(i,FP_REAL)-half)/real(mpx,FP_REAL); end do
        do i=1,mpy; xgp(i,2) = (real(i,FP_REAL)-half)/real(mpy,FP_REAL); end do
        do i=1,mpw; xgp(i,3) = (real(i,FP_REAL)-half)/real(mpw,FP_REAL); end do

        call fpndsp(3_FP_DIM,t3,n3,c3,k3,xgp,mp3,zp,wp,lidxp)

        maxdiff = zero
        do ix=1,mpx
           x = xgp(ix,1)
           do iy=1,mpy
              y = xgp(iy,2)
              do iw=1,mpw
                 w = xgp(iw,3)
                 iout = (ix-1)*mpy*mpw + (iy-1)*mpw + iw
                 maxdiff = max(maxdiff, abs(zp(iout)-poly_x(x)*poly_y(y)*poly_w(w)))
              end do
           end do
        end do
        if (maxdiff>=tol) then
           ok = .false.; write(useUnit,3100) maxdiff; return
        end if

        write(useUnit,'(a)') '[gate F] regrid_nd(dims=3) least-squares fit reproduces in-space polynomial'

        3000 format('[gate F] regrid_nd(dims=3) ',a,' = ',1pe12.3)
        3100 format('[gate F] regrid_nd(dims=3) fit /= closed-form f  (max |diff| = ',1pe12.3,')')

    contains

        pure real(FP_REAL) function poly_x(x)
            real(FP_REAL), intent(in) :: x
            poly_x = 1.0_FP_REAL + 0.5_FP_REAL*x - 0.3_FP_REAL*x**2 + 0.2_FP_REAL*x**3
        end function poly_x
        pure real(FP_REAL) function poly_y(y)
            real(FP_REAL), intent(in) :: y
            poly_y = 0.8_FP_REAL - 0.4_FP_REAL*y + 0.6_FP_REAL*y**2
        end function poly_y
        pure real(FP_REAL) function poly_w(w)
            real(FP_REAL), intent(in) :: w
            poly_w = 1.2_FP_REAL + 0.3_FP_REAL*w - 0.5_FP_REAL*w**2 + 0.1_FP_REAL*w**3 + 0.15_FP_REAL*w**4
        end function poly_w

    end function gate_regrid_nd_3d

    !> @brief Gate G — regrid_nd at dims=3 SMOOTHING (s>0): the generalized knot-direction arbiter.
    !!
    !! The first dims>2 smoothing fit, so the first runtime exercise of the N-D argmax arbiter
    !! (new_knot_dimension_nd) AND of the p>0 discontinuity branch at dims>2. No legacy dims=3
    !! reference exists, so the oracle is independent: data is a separable Runge product (NOT in the
    !! spline space, so smoothing genuinely adds knots). After fitting: (1) success, (2) the smoothing
    !! identity |fp-s|<=tol*s, (3) knots distributed -- in particular axis 3 received knots, which the
    !! old literal-2 arbiter could never do, and (4) fp equals the unit-weight SSR recomputed through
    !! the gate-A/E/F-verified fpndsp.
    logical function gate_regrid_nd_3d_smoothing(useUnit) result(ok)
        integer, intent(in) :: useUnit

        integer(FP_SIZE), parameter :: kx=3, ky=3, kw=3
        integer(FP_SIZE), parameter :: mx=9, my=8, mw=7, mz=mx*my*mw
        integer(FP_SIZE), parameter :: nest=17, maxm=mx, maxk1=kx+1
        integer(FP_SIZE), parameter :: maxc=(nest-kx-1)**3, lwrk3=4000, kwrk3=200

        real(FP_REAL)    :: t3(nest,3),xg3(maxm,3),c3(maxc),z3(mz),zfit(mz),w3(maxm,maxk1,3)
        real(FP_REAL)    :: fp0cal,fp,s,lo(3),hi(3),ssr,wrk(lwrk3)
        integer(FP_SIZE) :: iwrk(kwrk3),n3(3),k3(3),m3(3),nest3(3),nmin3(3),lidx3(maxm,3)
        integer(FP_FLAG) :: ier
        integer(FP_SIZE) :: ix,iy,iw,i,iout

        ok = .true.
        t3 = zero; xg3 = zero
        k3 = [kx,ky,kw]; m3 = [mx,my,mw]; nest3 = nest; nmin3 = 2*(k3+1); n3 = nmin3
        lo = zero; hi = one

        ! strictly-increasing data grids on [0,1] (endpoints included)
        do i=1,mx; xg3(i,1) = real(i-1,FP_REAL)/real(mx-1,FP_REAL); end do
        do i=1,my; xg3(i,2) = real(i-1,FP_REAL)/real(my-1,FP_REAL); end do
        do i=1,mw; xg3(i,3) = real(i-1,FP_REAL)/real(mw-1,FP_REAL); end do

        ! data: separable Runge product (not in the spline space), row-major (x slowest, w fastest)
        do ix=1,mx
           do iy=1,my
              do iw=1,mw
                 iout = (ix-1)*my*mw + (iy-1)*mw + iw
                 z3(iout) = runge(xg3(ix,1))*runge(xg3(iy,2))*runge(xg3(iw,3))
              end do
           end do
        end do

        ! calibrate: a huge s returns the least-squares polynomial residual fp0
        wrk = zero; iwrk = 0
        call regrid_nd(0_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,1.0e30_FP_REAL,nest3, &
                       n3,t3,c3,fp0cal,wrk,lwrk3,iwrk,kwrk3,ier)
        if (.not.FITPACK_SUCCESS(ier) .or. fp0cal<=zero) then
           ok = .false.; write(useUnit,3000) 'calibration failed (ier/fp0)', fp0cal; return
        end if

        ! smoothing fit at s = 1% of the polynomial residual: forces knot addition on every axis
        s = 1.0e-2_FP_REAL*fp0cal
        t3 = zero; wrk = zero; iwrk = 0
        call regrid_nd(0_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,s,nest3, &
                       n3,t3,c3,fp,wrk,lwrk3,iwrk,kwrk3,ier)
        if (.not.FITPACK_SUCCESS(ier)) then
           ok = .false.; write(useUnit,'(a,i0)') '[gate G] regrid_nd(dims=3) smoothing failed, ier=',ier; return
        end if

        ! check 1: smoothing identity -- the solver converges to |fp-s| < tol*s (tol = 1e-3)
        if (abs(fp-s) > 1.0e-3_FP_REAL*s) then
           ok = .false.; write(useUnit,3000) '|fp-s| exceeds tol*s', abs(fp-s); return
        end if

        ! check 2: the new arbiter distributed knots -- axis 3 was refined (impossible for the old
        ! literal-2 arbiter), and at least two axes received interior knots
        if (n3(3)<=nmin3(3) .or. count(n3>nmin3)<2) then
           ok = .false.; write(useUnit,3100) n3(1),n3(2),n3(3); return
        end if

        ! check 3: fp equals the unit-weight SSR recomputed through the independent fpndsp
        call fpndsp(3_FP_DIM,t3,n3,c3,k3,xg3,m3,zfit,w3,lidx3)
        ssr = zero
        do ix=1,mx
           do iy=1,my
              do iw=1,mw
                 iout = (ix-1)*my*mw + (iy-1)*mw + iw
                 ssr  = ssr + (z3(iout)-zfit(iout))**2
              end do
           end do
        end do
        if (abs(ssr-fp) > 1.0e-9_FP_REAL*max(fp,one)) then
           ok = .false.; write(useUnit,3000) 'fp /= SSR via fpndsp', abs(ssr-fp); return
        end if

        write(useUnit,'(a)') '[gate G] regrid_nd(dims=3) smoothing converges, distributes knots, fp==SSR'

        3000 format('[gate G] regrid_nd(dims=3) smoothing ',a,' = ',1pe12.3)
        3100 format('[gate G] regrid_nd(dims=3) knots not distributed: n = ',i0,1x,i0,1x,i0)

    contains

        pure real(FP_REAL) function runge(t)
            real(FP_REAL), intent(in) :: t
            runge = one/(one + 20.0_FP_REAL*(t-half)**2)
        end function runge

    end function gate_regrid_nd_3d_smoothing

end module fitpack_grid_nd_tests
