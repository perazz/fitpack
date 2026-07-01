! **************************************************************************************************
!   FITPACK N-D gridded core — correctness gates for the dimension-generic regrid/fpregr/fpgrre/fpndsp
!
!   Gates A/D/D' pin the dims=2 behaviour of the (now sole) gridded routines against independent
!   oracles over the standard 11x11 `daregr` battery across all iopt modes and spline orders; Gates
!   E/F/G exercise the genuine dims>2 paths (eval contraction and fit) against independent references.
!
!   Gates:
!     A. fpndsp vs bispev vs a splev separation-of-variables oracle (gridded evaluation, dims=2)
!     D. regrid fit self-consistency (fp==SSR via bispev, s=0 interpolates, valid knots) across the
!                               cubic iopt-chain, quintic and iopt=-1 lsq — covers fpregr + fpgrre
!     D'. same fp==SSR oracle on a NON-SQUARE grid (mx /= my) — guards the z = (ny,nx) y-fast data
!                               convention so an x<->y axis/size swap cannot hide behind mx==my
!     E. fpndsp(dims=3) vs independent references (slice 2) — separable coeffs vs a product of 1-D
!                               splev, and non-separable coeffs vs the dims=2 path combined along z
!     F. regrid(dims=3) least-squares FIT on given knots (slice 3) — the generalized fpgrre solve
!                               recovers an in-space polynomial; checked vs the closed-form (fpndsp)
!
!   See todo/fitpack_nd_grids.md (slices 1–3) and Dierckx, Ch. 5 §5.4 (pp. 98-103).
! **************************************************************************************************
module fitpack_grid_nd_tests
    use fitpack_core
    use fitpack_gridded_splines, only: fitpack_gridded_spline
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

    !> @brief Correctness gates for the gridded core: dims=2 oracles (A/D/D') and dims>2 paths (E/F/G).
    logical function test_nd_grid_equivalence(iunit) result(success)
        integer, optional, intent(in) :: iunit
        integer :: useUnit

        useUnit = output_unit
        if (present(iunit)) useUnit = iunit

        success = .true.

        ! Gate D: regrid(dims=2) fit self-consistency (backbone; transitively covers fpregr + fpgrre)
        if (success) success = gate_regrid_nd(useUnit)

        ! Gate D': same, on a NON-SQUARE grid (mx /= my) so an x<->y data/size swap cannot hide
        if (success) success = gate_regrid_nd_nonsquare(useUnit)

        ! Gate A: fpndsp vs bispev
        if (success) success = gate_fpndsp(useUnit)

        ! Gate E: fpndsp at dims=3 vs two independent references (the first genuine dims>2 path)
        if (success) success = gate_fpndsp_3d(useUnit)

        ! Gate F: regrid at dims=3 — the first genuine dims>2 FIT (generalized fpgrre solve)
        if (success) success = gate_regrid_nd_3d(useUnit)

        ! Gate G: regrid at dims=3 SMOOTHING (s>0) — the generalized knot-direction arbiter
        if (success) success = gate_regrid_nd_3d_smoothing(useUnit)

        ! Gate H: the generic fitpack_gridded_spline class (dims=3/5 fit+eval, row_major flag, comm)
        if (success) success = gate_gridded_spline_class(useUnit)

        ! ---- Peripherals: each *_nd routine bit-for-bit vs its 2-D oracle at dims=2, then dims=3 ----
        ! Gate I: bispeu_nd (scattered evaluation)
        if (success) success = gate_bispeu_nd(useUnit)

        ! Gate J: pardtc_nd / parder_nd / pardeu_nd (partial derivatives)
        if (success) success = gate_parder_nd(useUnit)

        ! Gate K: dblint_nd (box/domain integral)
        if (success) success = gate_dblint_nd(useUnit)

        ! Gate L: profil_nd (cross-section: fix one axis)
        if (success) success = gate_profil_nd(useUnit)

        ! Gate M: the class peripheral methods (eval/dfdx/dfdx_ongrid/integral/cross_section/derivative_spline)
        if (success) success = gate_gridded_spline_peripherals(useUnit)

        if (success) write(useUnit,'(a)') '[test_nd_grid_equivalence] all N-D gridded gates OK'

    end function test_nd_grid_equivalence

    !> @brief The 2-D fitting battery (orders + smoothing) shared by all gates.
    pure function cases() result(cs)
        type(grid_case) :: cs(3)
        cs(1) = grid_case(3,3,0.22_FP_REAL,'cubic smoothing s=0.22')
        cs(2) = grid_case(3,3,0.0_FP_REAL ,'cubic interpolation s=0')
        cs(3) = grid_case(5,5,0.2_FP_REAL ,'quintic smoothing s=0.2')
    end function cases

    !> @brief Run a dims=2 gridded fit for one case, returning the fitted knots/coefficients in the
    !!        bivariate (nx,tx,ny,ty) form the evaluation gates expect.
    subroutine fit_regrid(gc,nx,tx,ny,ty,c,fp,ier)
        type(grid_case),  intent(in)  :: gc
        integer(FP_SIZE), intent(out) :: nx,ny
        real(FP_REAL),    intent(out) :: tx(NXEST),ty(NYEST),c(NCMAX),fp
        integer(FP_FLAG), intent(out) :: ier

        integer(FP_SIZE) :: mx,my,m2(2),n2(2),k2(2),nest2(2)
        real(FP_REAL)    :: lo(2),hi(2),wrk(LWRK),t2(NXEST,2),xg(size(daregr_x),2)
        integer(FP_SIZE) :: iwrk(KWRK)
        real(FP_REAL)    :: zin(size(daregr_x)*size(daregr_y))

        mx = size(daregr_x,kind=FP_SIZE)
        my = size(daregr_y,kind=FP_SIZE)
        zin = reshape(daregr_z,[mx*my])   ! flat row-major (x slowest, y fastest) = regrid's z contract
        m2 = [mx,my]; k2 = [gc%kx,gc%ky]; nest2 = [NXEST,NYEST]
        lo = [daregr_x(1),daregr_y(1)]; hi = [daregr_x(mx),daregr_y(my)]
        xg = zero; xg(1:mx,1) = daregr_x; xg(1:my,2) = daregr_y
        t2 = zero; wrk = zero; iwrk = 0

        call regrid(IOPT_NEW_SMOOTHING,2_FP_DIM,m2,xg,zin,lo,hi,k2,gc%s,nest2, &
                       n2,t2,c,fp,wrk,LWRK,iwrk,KWRK,ier)

        nx = n2(1); ny = n2(2)
        tx = zero; tx(1:n2(1)) = t2(1:n2(1),1)
        ty = zero; ty(1:n2(2)) = t2(1:n2(2),2)
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
            allocate(t(maxn,2),xg(maxm,2),w(maxk1,maxm,2),lidx(maxm,2))
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

    !> @brief Gate D — regrid(dims=2) fit is self-consistent across all iopt modes.
    !!
    !! With regrid now the sole gridded engine there is no legacy regrid to compare against, so
    !! this asserts the defining property of every fit instead: fp equals the unit-weight residual
    !! sum of squares recomputed independently on the data grid through bispev (which routes through
    !! fpndsp), an s=0 smoothing fit interpolates (fp collapses to zero), and the knots stay valid.
    !! It runs the cubic iopt=0->1 continuation chain, a quintic fit, and a least-squares fit on
    !! prescribed knots (iopt=-1), so it still exercises fpregr + fpgrre end-to-end.
    logical function gate_regrid_nd(useUnit) result(ok)
        integer, intent(in) :: useUnit

        integer(FP_SIZE) :: mx,my
        integer(FP_SIZE) :: ntst(2),kk(2),m2(2),nest2(2)
        real(FP_REAL)    :: ttst(NXEST,2),ctst(NCMAX),fptst,wrkt(LWRK),xgt(size(daregr_x),2)
        integer(FP_SIZE) :: iwrkt(KWRK)
        integer(FP_FLAG) :: ier_t
        real(FP_REAL)    :: zin(size(daregr_x)*size(daregr_y)),xb,xe,yb,ye,s,lo(2),hi(2)
        integer(FP_SIZE) :: kx,ky,iopt,icase,i
        character(len=40):: lbl

        ok = .true.
        mx = size(daregr_x,kind=FP_SIZE); my = size(daregr_y,kind=FP_SIZE)
        xb = daregr_x(1); xe = daregr_x(mx); yb = daregr_y(1); ye = daregr_y(my)
        zin = reshape(daregr_z,[mx*my])    ! flat row-major (x slowest, y fastest) = regrid's z contract

        m2 = [mx,my]; nest2 = [NXEST,NYEST]; lo = [xb,yb]; hi = [xe,ye]
        xgt = zero; xgt(1:mx,1) = daregr_x; xgt(1:my,2) = daregr_y
        wrkt = zero; iwrkt = 0

        ! ---- Phase A: cubic, iopt=0 then iopt=1 continuation chain (incl. s=0 interpolation) ----
        kx = 3; ky = 3; kk = [kx,ky]
        do icase=1,4
            select case (icase)
               case (1); iopt = 0; s = 10.0_FP_REAL
               case (2); iopt = 1; s = 0.22_FP_REAL
               case (3); iopt = 1; s = 0.1_FP_REAL
               case (4); iopt = 1; s = 0.0_FP_REAL
            end select
            call regrid(iopt,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                           ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
            write(lbl,'(a,i0)') 'cubic iopt-chain #',icase
            if (.not.fit_ssr_ok(useUnit,trim(lbl),kk,m2,xgt,zin,s,iopt,ntst,ttst,ctst,fptst,ier_t)) then
                ok = .false.; return
            end if
        end do

        ! ---- Phase B: quintic, fresh iopt=0 ----
        kx = 5; ky = 5; kk = [kx,ky]; iopt = 0; s = 0.2_FP_REAL
        call regrid(iopt,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                       ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
        if (.not.fit_ssr_ok(useUnit,'quintic fresh iopt=0',kk,m2,xgt,zin,s,iopt,ntst,ttst,ctst,fptst,ier_t)) then
            ok = .false.; return
        end if

        ! ---- Phase C: least-squares on given knots, iopt=-1 (s=0 here is lsq, NOT interpolation) ----
        kx = 3; ky = 3; kk = [kx,ky]; s = 0.0_FP_REAL
        ntst = [11,11]
        ttst = zero
        ttst(kx+2:kx+4,1) = [(half*(i-2),i=1,3)]
        ttst(ky+2:ky+4,2) = ttst(kx+2:kx+4,1)
        call regrid(-1_FP_FLAG,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                       ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
        if (.not.fit_ssr_ok(useUnit,'lsq given knots iopt=-1',kk,m2,xgt,zin,s,-1_FP_SIZE,ntst,ttst,ctst,fptst,ier_t)) then
            ok = .false.; return
        end if

        write(useUnit,'(a)') '[gate D] regrid(dims=2): fp==SSR, s=0 interpolates, valid knots (cubic chain, quintic, iopt=-1)'
    end function gate_regrid_nd

    !> @brief Gate D' — regrid(dims=2) on a NON-SQUARE grid (mx /= my).
    !!
    !! The square 11x11 daregr battery cannot expose an x<->y swap in the data tensor or in the
    !! per-axis sizing: with mx==my both indexings address the same elements. A non-square grid makes
    !! any such swap a hard failure, so this directly guards the z = (ny,nx) (y-fast) storage
    !! convention. The synthetic data is deliberately NON-symmetric under x<->y (different
    !! centres/weights) so even a pure transpose on a hypothetical square grid would be caught. The
    !! fp==SSR oracle would not hold under any axis confusion.
    logical function gate_regrid_nd_nonsquare(useUnit) result(ok)
        integer, intent(in) :: useUnit

        integer(FP_SIZE), parameter :: MX = 9, MY = 13
        real(FP_REAL)    :: xg1(MX),yg1(MY),zin(MX*MY)
        integer(FP_SIZE) :: ntst(2),kk(2),m2(2),nest2(2)
        real(FP_REAL)    :: ttst(NXEST,2),ctst(NCMAX),fptst,wrkt(LWRK),xgt(MY,2)
        integer(FP_SIZE) :: iwrkt(KWRK)
        integer(FP_FLAG) :: ier_t
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
        wrkt = zero; iwrkt = 0
        kk = [3_FP_SIZE,3_FP_SIZE]

        do icase=1,2
            select case (icase)
               case (1); iopt = 0; s = 0.05_FP_REAL   ! smoothing
               case (2); iopt = 0; s = 0.0_FP_REAL    ! interpolation
            end select
            call regrid(iopt,2_FP_DIM,m2,xgt,zin,lo,hi,kk,s,nest2, &
                           ntst,ttst,ctst,fptst,wrkt,LWRK,iwrkt,KWRK,ier_t)
            write(lbl,'(a,i0)') 'non-square 9x13 #',icase
            if (.not.fit_ssr_ok(useUnit,trim(lbl),kk,m2,xgt,zin,s,iopt,ntst,ttst,ctst,fptst,ier_t)) then
                ok = .false.; return
            end if
        end do

        write(useUnit,'(a)') '[gate D''] regrid(dims=2): fp==SSR, s=0 interpolates (non-square 9x13)'
    end function gate_regrid_nd_nonsquare

    !> @brief Oracle check for one dims=2 regrid fit: the residual sum of squares fp must match
    !!        the value recomputed independently on the data grid through bispev, an s=0 smoothing
    !!        fit must interpolate (fp≈0), and the knots must be non-decreasing on each axis.
    logical function fit_ssr_ok(useUnit,label,k,m,xg,zin,s,iopt,n,t,c,fp,ier) result(ok)
        integer,          intent(in) :: useUnit
        character(*),     intent(in) :: label
        integer(FP_SIZE), intent(in) :: k(2),m(2),iopt,n(2)
        real(FP_REAL),    intent(in) :: xg(:,:),zin(:),s,t(:,:),c(:),fp
        integer(FP_FLAG), intent(in) :: ier

        real(FP_REAL)    :: zeval(m(1)*m(2)),ewrk((k(1)+1)*m(1)+(k(2)+1)*m(2)),ssr
        integer(FP_SIZE) :: eiwrk(m(1)+m(2)),d
        integer(FP_FLAG) :: ier_e

        ok = .true.

        ! the fit itself must have succeeded
        if (.not.FITPACK_SUCCESS(ier)) then
            write(useUnit,2000) trim(label),'fit',ier; ok = .false.; return
        end if

        ! knots must be non-decreasing on each axis
        do d=1,2
           if (any(t(2:n(d),d)<t(1:n(d)-1,d))) then
              write(useUnit,2100) trim(label),d; ok = .false.; return
           end if
        end do

        ! fp must equal the unit-weight SSR recomputed independently through bispev (-> fpndsp)
        call bispev(t(1:n(1),1),n(1),t(1:n(2),2),n(2),c,k(1),k(2), &
                    xg(1:m(1),1),m(1),xg(1:m(2),2),m(2),zeval, &
                    ewrk,size(ewrk,kind=FP_SIZE),eiwrk,size(eiwrk,kind=FP_SIZE),ier_e)
        if (.not.FITPACK_SUCCESS(ier_e)) then
            write(useUnit,2000) trim(label),'bispev',ier_e; ok = .false.; return
        end if
        ssr = sum((zin(1:m(1)*m(2))-zeval)**2)
        if (abs(ssr-fp) > 1.0e-9_FP_REAL*max(fp,one)) then
            write(useUnit,2200) trim(label),abs(ssr-fp); ok = .false.; return
        end if

        ! a smoothing fit (iopt>=0) with s=0 is an interpolating spline: the residual collapses to zero
        if (iopt>=0 .and. s==zero .and. fp>1.0e-9_FP_REAL) then
            write(useUnit,2300) trim(label),fp; ok = .false.; return
        end if

        2000 format('[gate D] ',a,': ',a,' failed, ier=',i0)
        2100 format('[gate D] ',a,': knots not non-decreasing on axis ',i0)
        2200 format('[gate D] ',a,': fp /= SSR via bispev, |diff| = ',1pe12.3)
        2300 format('[gate D] ',a,': s=0 did not interpolate, fp = ',1pe12.3)
    end function fit_ssr_ok

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

        real(FP_REAL)    :: t3(maxn,3),xg3(maxm,3),c3(nc3),z3(mout3),w3(maxk1,maxm,3)
        integer(FP_SIZE) :: lidx3(maxm,3),n3(3),k3(3),m3(3)
        ! separable reference (outer product of three 1-D splev evaluations)
        real(FP_REAL)    :: av(maxn),bv(maxn),dv(maxn),sx(mx),sy(my),sz(mz),ref
        ! non-separable reference (verified dims=2 plane evals, combined along z by splev)
        real(FP_REAL)    :: v2(mxy,nk1z),z2(mxy),c2(nk1x*nk1y),w2(kx+1,mx,2),vz(mz),cz(maxn)
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

    !> @brief Gate F — regrid at dims=3 (least-squares fit on given knots) vs a closed-form oracle.
    !!
    !! The first genuine dims>2 FIT path: the generalized fpgrre solve (alternating-direction
    !! reduction + ping-pong buffer + back-substitution + residual contraction). Data is a separable
    !! in-space polynomial f(x,y,w)=Px(x)*Py(y)*Pw(w) with deg Px/Py/Pw = kx/ky/kw, which lies EXACTLY
    !! in the tensor B-spline space for any clamped knots. Fit via regrid(dims=3, iopt=-1) on
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
        real(FP_REAL)    :: xgp(maxmp,3),zp(mpx*mpy*mpw),wp(maxk1,maxmp,3)
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
        call regrid(-1_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,s,nest3, &
                       n3,t3,c3,fpf,wrk,LWRK,iwrk,KWRK,ier)
        if (ier/=FITPACK_OK) then
           ok = .false.; write(useUnit,'(a,i0)') '[gate F] regrid(dims=3) fit failed, ier=',ier; return
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

        write(useUnit,'(a)') '[gate F] regrid(dims=3) least-squares fit reproduces in-space polynomial'

        3000 format('[gate F] regrid(dims=3) ',a,' = ',1pe12.3)
        3100 format('[gate F] regrid(dims=3) fit /= closed-form f  (max |diff| = ',1pe12.3,')')

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

    !> @brief Gate G — regrid at dims=3 SMOOTHING (s>0): the generalized knot-direction arbiter.
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

        real(FP_REAL)    :: t3(nest,3),xg3(maxm,3),c3(maxc),z3(mz),zfit(mz),w3(maxk1,maxm,3)
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
        call regrid(0_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,1.0e30_FP_REAL,nest3, &
                       n3,t3,c3,fp0cal,wrk,lwrk3,iwrk,kwrk3,ier)
        if (.not.FITPACK_SUCCESS(ier) .or. fp0cal<=zero) then
           ok = .false.; write(useUnit,3000) 'calibration failed (ier/fp0)', fp0cal; return
        end if

        ! smoothing fit at s = 1% of the polynomial residual: forces knot addition on every axis
        s = 1.0e-2_FP_REAL*fp0cal
        t3 = zero; wrk = zero; iwrk = 0
        call regrid(0_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,s,nest3, &
                       n3,t3,c3,fp,wrk,lwrk3,iwrk,kwrk3,ier)
        if (.not.FITPACK_SUCCESS(ier)) then
           ok = .false.; write(useUnit,'(a,i0)') '[gate G] regrid(dims=3) smoothing failed, ier=',ier; return
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

        write(useUnit,'(a)') '[gate G] regrid(dims=3) smoothing converges, distributes knots, fp==SSR'

        3000 format('[gate G] regrid(dims=3) smoothing ',a,' = ',1pe12.3)
        3100 format('[gate G] regrid(dims=3) knots not distributed: n = ',i0,1x,i0,1x,i0)

    contains

        pure real(FP_REAL) function runge(t)
            real(FP_REAL), intent(in) :: t
            runge = one/(one + 20.0_FP_REAL*(t-half)**2)
        end function runge

    end function gate_regrid_nd_3d_smoothing

    !> @brief Gate H — the generic fitpack_gridded_spline class is pure marshalling over regrid/ndspev.
    !!
    !! Asserts that a dims=3 and a dims=5 fit+eval through the class reproduce a direct regrid+ndspev
    !! call bit-for-bit; that the row_major flag (a natural column-major array) yields the same fit as
    !! the flat row-major buffer; and that comm pack/expand round-trips.
    logical function gate_gridded_spline_class(useUnit) result(ok)
        integer, intent(in) :: useUnit
        ok = .true.
        if (ok) ok = class_matches_direct(useUnit,'dims=3',3_FP_DIM,[6_FP_SIZE,5_FP_SIZE,7_FP_SIZE],3_FP_SIZE)
        if (ok) ok = class_matches_direct(useUnit,'dims=5', &
                          5_FP_DIM,[4_FP_SIZE,3_FP_SIZE,4_FP_SIZE,3_FP_SIZE,4_FP_SIZE],2_FP_SIZE)
        if (ok) ok = class_row_major_flag(useUnit)
        if (ok) ok = class_comm_roundtrip(useUnit)
        if (ok) write(useUnit,'(a)') &
            '[gate H] fitpack_gridded_spline == direct regrid+ndspev; row_major flag; comm round-trip'
    end function gate_gridded_spline_class

    !> @brief Build a strictly-increasing coordinate table and a smooth flat (row-major) value grid.
    subroutine build_grid(dims,m,xg,zflat)
        integer(FP_DIM),  intent(in)  :: dims
        integer(FP_SIZE), intent(in)  :: m(:)
        real(FP_REAL),    intent(out) :: xg(:,:)
        real(FP_REAL),    intent(out) :: zflat(:)
        integer(FP_DIM)  :: d
        integer(FP_SIZE) :: p,idx(dims),k
        real(FP_REAL)    :: v
        xg = zero
        do d=1,dims
           do p=1,m(d)
              xg(p,d) = real(p-1,FP_REAL)/real(m(d)-1,FP_REAL) + 0.15_FP_REAL*real(d,FP_REAL)
           end do
        end do
        ! flat row-major (axis 1 slowest): iterate the multi-index with axis dims fastest
        idx = 1
        do p=1,product(m(1:dims))
           v = one
           do d=1,dims
              v = v + real(d,FP_REAL)*xg(idx(d),d) + xg(idx(d),d)**2      ! smooth, non-symmetric
           end do
           zflat(p) = v
           k = dims
           do
              idx(k) = idx(k)+1
              if (idx(k)<=m(k)) exit
              idx(k) = 1; k = k-1
              if (k<1) exit
           end do
        end do
    end subroutine build_grid

    !> @brief Fit+eval a grid through the class and compare bit-for-bit to a direct regrid+ndspev call.
    logical function class_matches_direct(useUnit,label,dims,m,korder) result(ok)
        integer,          intent(in) :: useUnit
        character(*),     intent(in) :: label
        integer(FP_DIM),  intent(in) :: dims
        integer(FP_SIZE), intent(in) :: m(:),korder

        type(fitpack_gridded_spline) :: obj
        integer(FP_SIZE) :: maxm,nz,k(dims),n2(dims)
        real(FP_REAL), allocatable :: xg(:,:),zflat(:),t2(:,:),c2(:),wrk(:),zc_dir(:),zc_cls(:)
        integer(FP_SIZE), allocatable :: iwrk(:)
        real(FP_REAL)    :: fp2,s
        integer(FP_FLAG) :: ier,ierc

        ok  = .true.
        maxm = maxval(m); nz = product(m(1:dims)); k = korder; s = 0.05_FP_REAL
        allocate(xg(maxm,dims),zflat(nz),zc_dir(nz),zc_cls(nz))
        call build_grid(dims,m,xg,zflat)

        ! ---- class path ----
        ierc = obj%new_fit(xg,zflat,m=m,order=korder,smoothing=s)
        if (.not.FITPACK_SUCCESS(ierc)) then
            write(useUnit,4000) trim(label),'class fit',ierc; ok = .false.; return
        end if
        zc_cls = obj%eval_ongrid(xg,m,ier)
        if (.not.FITPACK_SUCCESS(ier)) then
            write(useUnit,4000) trim(label),'class eval',ier; ok = .false.; return
        end if

        ! ---- direct regrid + ndspev with identical inputs (obj already sized the workspace) ----
        allocate(t2(maxval(obj%nest(1:dims)),dims),source=zero)
        allocate(c2(size(obj%c)),source=zero)
        allocate(wrk(obj%lwrk),source=zero)
        allocate(iwrk(obj%liwrk),source=0)
        n2 = 0
        call regrid(IOPT_NEW_SMOOTHING,dims,m,xg,zflat,obj%left(1:dims),obj%right(1:dims), &
                    k,s,obj%nest(1:dims),n2,t2,c2,fp2,wrk,obj%lwrk,iwrk,obj%liwrk,ier)
        if (.not.FITPACK_SUCCESS(ier)) then
            write(useUnit,4000) trim(label),'direct regrid',ier; ok = .false.; return
        end if

        if (any(n2/=obj%knots(1:dims)) .or. fp2/=obj%fp .or. &
            .not.all(c2(1:size(obj%c))==obj%c)) then
            write(useUnit,4100) trim(label); ok = .false.; return
        end if

        call ndspev(dims,t2,n2,c2,k,xg,m,zc_dir,wrk,obj%lwrk,iwrk,obj%liwrk,ier)
        if (.not.all(zc_dir==zc_cls)) then
            write(useUnit,4200) trim(label),maxval(abs(zc_dir-zc_cls)); ok = .false.; return
        end if

        4000 format('[gate H] ',a,': ',a,' failed, ier=',i0)
        4100 format('[gate H] ',a,': class fit /= direct regrid (knots/c/fp)')
        4200 format('[gate H] ',a,': class eval /= direct ndspev, max|diff| = ',1pe12.3)
    end function class_matches_direct

    !> @brief A natural column-major array (row_major=.false.) must give the same fit as the flat
    !!        row-major buffer of the same logical grid.
    logical function class_row_major_flag(useUnit) result(ok)
        integer, intent(in) :: useUnit

        integer(FP_DIM),  parameter :: dims = 3
        integer(FP_SIZE), parameter :: m1 = 6, m2 = 5, m3 = 7
        type(fitpack_gridded_spline) :: obj_flat,obj_nat
        real(FP_REAL) :: xg(max(m1,max(m2,m3)),dims),zflat(m1*m2*m3),anat(m1,m2,m3)
        integer(FP_SIZE) :: mm(dims),i1,i2,i3,p
        integer(FP_FLAG) :: e1,e2

        ok = .true.
        mm = [m1,m2,m3]
        call build_grid(dims,mm,xg,zflat)

        ! natural Fortran array A(i1,i2,i3) = f(x_i1,y_i2,w_i3), same values as the flat row-major grid
        p = 0
        do i1=1,m1
           do i2=1,m2
              do i3=1,m3
                 p = p+1                       ! p walks row-major (axis 1 slowest) == zflat order
                 anat(i1,i2,i3) = zflat(p)
              end do
           end do
        end do

        e1 = obj_flat%new_fit(xg,zflat,m=mm,smoothing=0.02_FP_REAL)                    ! flat row-major
        e2 = obj_nat %new_fit(xg,anat,row_major=.false._FP_BOOL,smoothing=0.02_FP_REAL) ! natural column-major

        if (.not.FITPACK_SUCCESS(e1) .or. .not.FITPACK_SUCCESS(e2)) then
            write(useUnit,'(a)') '[gate H] row_major: a fit failed'; ok = .false.; return
        end if
        if (any(obj_flat%knots(1:dims)/=obj_nat%knots(1:dims)) .or. obj_flat%fp/=obj_nat%fp .or. &
            .not.all(obj_flat%c==obj_nat%c)) then
            write(useUnit,'(a)') '[gate H] row_major=.false. fit /= flat row-major fit'; ok = .false.; return
        end if
    end function class_row_major_flag

    !> @brief comm pack/expand must round-trip a fitted class object.
    logical function class_comm_roundtrip(useUnit) result(ok)
        integer, intent(in) :: useUnit

        integer(FP_DIM),  parameter :: dims = 3
        integer(FP_SIZE), parameter :: mm(3) = [6,5,7]
        type(fitpack_gridded_spline) :: obj,obj2
        real(FP_REAL) :: xg(7,dims),zflat(6*5*7)
        real(FP_COMM), allocatable :: buffer(:)
        integer(FP_FLAG) :: e1

        ok = .true.
        call build_grid(dims,mm,xg,zflat)
        e1 = obj%new_fit(xg,zflat,m=mm,smoothing=0.03_FP_REAL)
        if (.not.FITPACK_SUCCESS(e1)) then
            write(useUnit,'(a)') '[gate H] comm: fit failed'; ok = .false.; return
        end if

        allocate(buffer(obj%comm_size()))
        call obj%comm_pack(buffer)
        call obj2%comm_expand(buffer)

        if (obj2%dims/=obj%dims .or. any(obj2%knots(1:dims)/=obj%knots(1:dims)) .or. &
            obj2%fp/=obj%fp .or. size(obj2%c)/=size(obj%c)) then
            write(useUnit,'(a)') '[gate H] comm: metadata round-trip mismatch'; ok = .false.; return
        end if
        if (.not.all(obj2%c==obj%c) .or. .not.all(obj2%z==obj%z)) then
            write(useUnit,'(a)') '[gate H] comm: bulk-array round-trip mismatch'; ok = .false.; return
        end if
    end function class_comm_roundtrip

    ! =================================================================================================
    ! PERIPHERAL GATES: scattered eval / derivatives / integral / cross-section
    ! =================================================================================================

    !> @brief Evaluate the nder-th derivative of the polynomial sum_i coef(i) x^(i-1).
    pure real(FP_REAL) function poly_deriv(coef,nder,x) result(v)
        real(FP_REAL),    intent(in) :: coef(:)
        integer(FP_SIZE), intent(in) :: nder
        real(FP_REAL),    intent(in) :: x
        integer(FP_SIZE) :: i,p,j
        real(FP_REAL)    :: term
        v = zero
        do i=1,size(coef,kind=FP_SIZE)
           p = i-1
           if (p<nder) cycle
           term = coef(i)
           do j=0,nder-1
              term = term*real(p-j,FP_REAL)
           end do
           v = v + term*x**(p-nder)
        end do
    end function poly_deriv

    !> @brief Definite integral over [a,b] of the polynomial sum_i coef(i) x^(i-1).
    pure real(FP_REAL) function poly_integral(coef,a,b) result(v)
        real(FP_REAL), intent(in) :: coef(:),a,b
        integer(FP_SIZE) :: i
        v = zero
        do i=1,size(coef,kind=FP_SIZE)
           v = v + coef(i)*(b**i - a**i)/real(i,FP_REAL)
        end do
    end function poly_integral

    !> @brief Build a separable dims=3 tensor spline on clamped [0,1] knots: coeffs c(i,j,l)=a(i)b(j)d(l),
    !!        so s(x,y,z) = sa(x)*sb(y)*sd(z) with sa/sb/sd the 1-D splines of a/b/d. Returns the tensor
    !!        knots/orders/coeffs and the per-axis 1-D coefficient vectors (for an independent splev oracle).
    subroutine build_separable_3d(kx,ky,kw,nk1x,nk1y,nk1w,t3,n3,k3,c3,ca,cb,cd)
        integer(FP_SIZE), intent(in)  :: kx,ky,kw,nk1x,nk1y,nk1w
        real(FP_REAL),    intent(out) :: t3(:,:)
        integer(FP_SIZE), intent(out) :: n3(3),k3(3)
        real(FP_REAL),    intent(out) :: c3(:),ca(:),cb(:),cd(:)
        integer(FP_SIZE) :: ix,iy,iw,iout
        k3 = [kx,ky,kw]
        t3 = zero
        call clamped_unit_knots(kx,nk1x,t3(:,1),n3(1))
        call clamped_unit_knots(ky,nk1y,t3(:,2),n3(2))
        call clamped_unit_knots(kw,nk1w,t3(:,3),n3(3))
        do ix=1,nk1x; ca(ix) = one + 0.30_FP_REAL*ix - 0.05_FP_REAL*ix*ix; end do
        do iy=1,nk1y; cb(iy) = 0.5_FP_REAL + 0.20_FP_REAL*iy;              end do
        do iw=1,nk1w; cd(iw) = 2.0_FP_REAL - 0.15_FP_REAL*iw + 0.02_FP_REAL*iw*iw; end do
        do ix=1,nk1x
           do iy=1,nk1y
              do iw=1,nk1w
                 iout = ((ix-1)*nk1y+(iy-1))*nk1w + iw
                 c3(iout) = ca(ix)*cb(iy)*cd(iw)
              end do
           end do
        end do
    end subroutine build_separable_3d

    !> @brief Gate I — bispeu_nd (scattered evaluation).
    !!
    !! dims=2: identical to bispeu bit-for-bit over the daregr fit battery. dims=3: a separable spline
    !! (coeffs a(i)b(j)d(l)) evaluated at scattered points must equal the product of the three 1-D splev
    !! evaluations. Distinct per-axis orders/sizes so an axis swap is a hard failure.
    logical function gate_bispeu_nd(useUnit) result(ok)
        integer, intent(in) :: useUnit
        integer(FP_SIZE), parameter :: NP = 7
        ! dims=2
        type(grid_case)  :: cs(3),gc
        integer(FP_SIZE) :: nx,ny,nc,kx,ky,icase,i,mx,my
        integer(FP_FLAG) :: ier,ier2
        real(FP_REAL)    :: tx(NXEST),ty(NYEST),c(NCMAX),fp
        real(FP_REAL)    :: px(NP),py(NP),zref(NP),ztest(NP),bwrk(64)
        real(FP_REAL)    :: t2(NXEST,2),xg2(2,NP)
        integer(FP_SIZE) :: n2(2),k2(2)
        ! dims=3
        integer(FP_SIZE), parameter :: kx3=3,ky3=2,kw3=4,nk1x=6,nk1y=5,nk1w=7
        integer(FP_SIZE), parameter :: nxk=nk1x+kx3+1,nyk=nk1y+ky3+1,nwk=nk1w+kw3+1,maxn=nwk
        real(FP_REAL)    :: t3(maxn,3),c3(nk1x*nk1y*nk1w),ca(nk1x),cb(nk1y),cd(nk1w)
        real(FP_REAL)    :: xg3(3,NP),znd(NP),fx(NP),fy(NP),fw(NP)
        real(FP_REAL)    :: sa(maxn),sb(maxn),sd(maxn)
        integer(FP_SIZE) :: n3(3),k3(3)
        real(FP_REAL)    :: maxdiff

        ok = .true.
        cs = cases()
        mx = size(daregr_x,kind=FP_SIZE); my = size(daregr_y,kind=FP_SIZE)

        ! ---- dims=2 : bit-for-bit vs bispeu ----
        do icase=1,size(cs)
            gc = cs(icase); kx = gc%kx; ky = gc%ky
            call fit_regrid(gc,nx,tx,ny,ty,c,fp,ier)
            if (.not.FITPACK_SUCCESS(ier)) then
                ok=.false.; write(useUnit,5000) 'fit',trim(gc%label),FITPACK_MESSAGE(ier); return
            end if
            nc = (nx-kx-1)*(ny-ky-1)
            do i=1,NP
               px(i) = daregr_x(1) + (daregr_x(mx)-daregr_x(1))*(real(i,FP_REAL)-half)/real(NP,FP_REAL)
               py(i) = daregr_y(1) + (daregr_y(my)-daregr_y(1))*(real(NP-i,FP_REAL)+half)/real(NP,FP_REAL)
            end do
            call bispeu(tx,nx,ty,ny,c,kx,ky,px,py,zref,NP,bwrk,size(bwrk,kind=FP_SIZE),ier)
            t2 = zero; t2(1:nx,1)=tx(1:nx); t2(1:ny,2)=ty(1:ny)
            n2 = [nx,ny]; k2 = [kx,ky]; xg2(1,:)=px; xg2(2,:)=py
            call bispeu_nd(2_FP_DIM,t2,n2,c(1:nc),k2,xg2,NP,ztest,ier2)
            if (.not.FITPACK_SUCCESS(ier) .or. .not.FITPACK_SUCCESS(ier2) .or. .not.all(zref==ztest)) then
                ok=.false.; write(useUnit,5100) trim(gc%label),maxval(abs(zref-ztest)); return
            end if
        end do

        ! ---- dims=3 : separable spline vs product of 1-D splev ----
        call build_separable_3d(kx3,ky3,kw3,nk1x,nk1y,nk1w,t3,n3,k3,c3,ca,cb,cd)
        do i=1,NP
           xg3(1,i) = (real(i,FP_REAL)-half)/real(NP,FP_REAL)
           xg3(2,i) = one - 0.7_FP_REAL*xg3(1,i)
           xg3(3,i) = 0.2_FP_REAL + 0.6_FP_REAL*xg3(1,i)
        end do
        call bispeu_nd(3_FP_DIM,t3,n3,c3,k3,xg3,NP,znd,ier)
        if (.not.FITPACK_SUCCESS(ier)) then
            ok=.false.; write(useUnit,'(a,i0)') '[gate I] bispeu_nd(dims=3) failed, ier=',ier; return
        end if
        sa=zero; sa(1:nk1x)=ca; sb=zero; sb(1:nk1y)=cb; sd=zero; sd(1:nk1w)=cd
        call splev(t3(1:n3(1),1),n3(1),sa(1:n3(1)),kx3,xg3(1,:),fx,NP,OUTSIDE_NEAREST_BND,ier)
        call splev(t3(1:n3(2),2),n3(2),sb(1:n3(2)),ky3,xg3(2,:),fy,NP,OUTSIDE_NEAREST_BND,ier)
        call splev(t3(1:n3(3),3),n3(3),sd(1:n3(3)),kw3,xg3(3,:),fw,NP,OUTSIDE_NEAREST_BND,ier)
        maxdiff = maxval(abs(znd - fx*fy*fw))
        if (maxdiff > 1.0e-12_FP_REAL) then
            ok=.false.; write(useUnit,5200) maxdiff; return
        end if

        write(useUnit,'(a)') '[gate I] bispeu_nd == bispeu (dims=2); dims=3 == separable splev product'
        5000 format('[gate I] ',a,' failed for ',a,': ',a)
        5100 format('[gate I] bispeu_nd /= bispeu for ',a,'  (max |diff| = ',1pe12.3,')')
        5200 format('[gate I] bispeu_nd(dims=3) /= splev product  (max |diff| = ',1pe12.3,')')
    end function gate_bispeu_nd

    !> @brief Gate J — pardtc_nd / parder_nd / pardeu_nd (partial derivatives).
    !!
    !! dims=2: all three bit-for-bit vs pardtc/parder/pardeu across derivative orders. dims=3: an
    !! in-space separable polynomial fitted with regrid(iopt=-1,s=0) is differentiated by parder_nd
    !! (grid) and pardeu_nd (scattered) and checked against the closed-form analytic derivative.
    logical function gate_parder_nd(useUnit) result(ok)
        integer, intent(in) :: useUnit
        integer(FP_SIZE), parameter :: NP = 6
        ! dims=2
        type(grid_case)  :: cs(3),gc
        integer(FP_SIZE) :: nx,ny,nc,ncnew,kx,ky,icase,i,mx,my,mz,ic
        integer(FP_FLAG) :: ier,ier2
        real(FP_REAL)    :: tx(NXEST),ty(NYEST),c(NCMAX),fp
        real(FP_REAL)    :: newc2(NCMAX),newcnd(NCMAX)
        real(FP_REAL)    :: zg2(size(daregr_x)*size(daregr_y)),zgnd(size(daregr_x)*size(daregr_y))
        real(FP_REAL)    :: px(NP),py(NP),zs2(NP),zsnd(NP)
        real(FP_REAL)    :: wrkp(4000),wrknd(4000)
        integer(FP_SIZE) :: iwrkp(400),iwrknd(400)
        real(FP_REAL)    :: t2(NXEST,2),xg2(size(daregr_x),2),xgp2(2,NP)
        integer(FP_SIZE) :: n2(2),k2(2),m2(2),nu2(2),dxy(2,3)
        ! dims=3
        integer(FP_SIZE), parameter :: k3v=3,nk1=6,mg=8,nkk=nk1+k3v+1
        integer(FP_SIZE), parameter :: mpx=5,mpy=4,mpw=6,maxmp=mpw,nc3=nk1**3
        real(FP_REAL),    parameter :: px_c(4)=[1.0_FP_REAL,0.5_FP_REAL,-0.3_FP_REAL,0.2_FP_REAL]
        real(FP_REAL),    parameter :: py_c(4)=[0.8_FP_REAL,-0.4_FP_REAL,0.6_FP_REAL,0.1_FP_REAL]
        real(FP_REAL),    parameter :: pw_c(4)=[1.2_FP_REAL,0.3_FP_REAL,-0.5_FP_REAL,0.15_FP_REAL]
        real(FP_REAL)    :: t3(nkk,3),xg3(mg,3),c3(nc3),z3(mg*mg*mg),fpf,lo(3),hi(3)
        real(FP_REAL)    :: xgp(maxmp,3),zp(mpx*mpy*mpw),ptsc(3,NP),zsc(NP)
        integer(FP_SIZE) :: n3(3),k3(3),m3(3),nest3(3),mp3(3),nu3(3)
        integer(FP_SIZE) :: ix,iy,iw,iout
        real(FP_REAL)    :: x,y,w,maxdiff

        ok = .true.
        cs = cases()
        mx = size(daregr_x,kind=FP_SIZE); my = size(daregr_y,kind=FP_SIZE); mz = mx*my
        dxy = reshape([1,0, 0,1, 2,1],[2,3])   ! (dx,dy) combos, all < 3

        ! ---- dims=2 : bit-for-bit vs pardtc / parder / pardeu ----
        do icase=1,size(cs)
            gc = cs(icase); kx = gc%kx; ky = gc%ky
            call fit_regrid(gc,nx,tx,ny,ty,c,fp,ier)
            if (.not.FITPACK_SUCCESS(ier)) then
                ok=.false.; write(useUnit,6000) 'fit',trim(gc%label),FITPACK_MESSAGE(ier); return
            end if
            nc = (nx-kx-1)*(ny-ky-1)
            t2 = zero; t2(1:nx,1)=tx(1:nx); t2(1:ny,2)=ty(1:ny)
            n2 = [nx,ny]; k2 = [kx,ky]; m2 = [mx,my]
            xg2 = zero; xg2(1:mx,1)=daregr_x; xg2(1:my,2)=daregr_y
            do i=1,NP
               px(i) = daregr_x(1) + (daregr_x(mx)-daregr_x(1))*(real(i,FP_REAL)-half)/real(NP,FP_REAL)
               py(i) = daregr_y(1) + (daregr_y(my)-daregr_y(1))*(real(NP-i,FP_REAL)+half)/real(NP,FP_REAL)
            end do
            xgp2(1,:)=px; xgp2(2,:)=py

            do ic=1,3
                nu2 = dxy(:,ic)
                if (nu2(1)>=kx .or. nu2(2)>=ky) cycle
                ncnew = (nx-kx-1-nu2(1))*(ny-ky-1-nu2(2))

                ! pardtc
                call pardtc(tx,nx,ty,ny,c,kx,ky,nu2(1),nu2(2),newc2,ier)
                call pardtc_nd(2_FP_DIM,t2,n2,c(1:nc),k2,nu2,newcnd(1:nc),ier2)
                if (.not.FITPACK_SUCCESS(ier) .or. .not.FITPACK_SUCCESS(ier2) .or. &
                    .not.all(newc2(1:ncnew)==newcnd(1:ncnew))) then
                    ok=.false.; write(useUnit,6100) 'pardtc',trim(gc%label),nu2(1),nu2(2); return
                end if

                ! parder (grid)
                call parder(tx,nx,ty,ny,c,kx,ky,nu2(1),nu2(2),daregr_x,mx,daregr_y,my,zg2, &
                            wrkp,size(wrkp,kind=FP_SIZE),iwrkp,size(iwrkp,kind=FP_SIZE),ier)
                call parder_nd(2_FP_DIM,t2,n2,c(1:nc),k2,nu2,xg2,m2,zgnd, &
                            wrknd,size(wrknd,kind=FP_SIZE),iwrknd,size(iwrknd,kind=FP_SIZE),ier2)
                if (.not.FITPACK_SUCCESS(ier) .or. .not.FITPACK_SUCCESS(ier2) .or. &
                    .not.all(zg2(1:mz)==zgnd(1:mz))) then
                    ok=.false.; write(useUnit,6100) 'parder',trim(gc%label),nu2(1),nu2(2); return
                end if

                ! pardeu (scattered)
                call pardeu(tx,nx,ty,ny,c,kx,ky,nu2(1),nu2(2),px,py,zs2,NP, &
                            wrkp,size(wrkp,kind=FP_SIZE),iwrkp,size(iwrkp,kind=FP_SIZE),ier)
                call pardeu_nd(2_FP_DIM,t2,n2,c(1:nc),k2,nu2,xgp2,NP,zsnd, &
                            wrknd,size(wrknd,kind=FP_SIZE),ier2)
                if (.not.FITPACK_SUCCESS(ier) .or. .not.FITPACK_SUCCESS(ier2) .or. .not.all(zs2==zsnd)) then
                    ok=.false.; write(useUnit,6100) 'pardeu',trim(gc%label),nu2(1),nu2(2); return
                end if
            end do
        end do

        ! ---- dims=3 : in-space polynomial, analytic derivative oracle ----
        t3 = zero; xg3 = zero; k3 = k3v; m3 = mg; lo = zero; hi = one
        call clamped_unit_knots(k3v,nk1,t3(:,1),n3(1))
        call clamped_unit_knots(k3v,nk1,t3(:,2),n3(2))
        call clamped_unit_knots(k3v,nk1,t3(:,3),n3(3))
        nest3 = n3
        do i=1,mg; xg3(i,1) = real(i-1,FP_REAL)/real(mg-1,FP_REAL); end do
        xg3(:,2) = xg3(:,1); xg3(:,3) = xg3(:,1)
        do ix=1,mg
           x = xg3(ix,1)
           do iy=1,mg
              y = xg3(iy,2)
              do iw=1,mg
                 w = xg3(iw,3)
                 iout = ((ix-1)*mg+(iy-1))*mg + iw
                 z3(iout) = poly_deriv(px_c,0_FP_SIZE,x)*poly_deriv(py_c,0_FP_SIZE,y)*poly_deriv(pw_c,0_FP_SIZE,w)
              end do
           end do
        end do
        wrknd = zero; iwrknd = 0
        call regrid(-1_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,zero,nest3,n3,t3,c3,fpf, &
                    wrknd,size(wrknd,kind=FP_SIZE),iwrknd,size(iwrknd,kind=FP_SIZE),ier)
        if (.not.FITPACK_SUCCESS(ier) .or. fpf>1.0e-9_FP_REAL) then
            ok=.false.; write(useUnit,6200) 'fit residual', fpf; return
        end if

        nu3 = [1_FP_SIZE,0_FP_SIZE,1_FP_SIZE]
        mp3 = [mpx,mpy,mpw]
        do i=1,mpx; xgp(i,1) = (real(i,FP_REAL)-half)/real(mpx,FP_REAL); end do
        do i=1,mpy; xgp(i,2) = (real(i,FP_REAL)-half)/real(mpy,FP_REAL); end do
        do i=1,mpw; xgp(i,3) = (real(i,FP_REAL)-half)/real(mpw,FP_REAL); end do
        call parder_nd(3_FP_DIM,t3,n3,c3,k3,nu3,xgp,mp3,zp, &
                       wrknd,size(wrknd,kind=FP_SIZE),iwrknd,size(iwrknd,kind=FP_SIZE),ier)
        if (.not.FITPACK_SUCCESS(ier)) then
            ok=.false.; write(useUnit,'(a,i0)') '[gate J] parder_nd(dims=3) failed, ier=',ier; return
        end if
        maxdiff = zero
        do ix=1,mpx
           x = xgp(ix,1)
           do iy=1,mpy
              y = xgp(iy,2)
              do iw=1,mpw
                 w = xgp(iw,3)
                 iout = ((ix-1)*mpy+(iy-1))*mpw + iw
                 maxdiff = max(maxdiff, abs(zp(iout) - &
                     poly_deriv(px_c,nu3(1),x)*poly_deriv(py_c,nu3(2),y)*poly_deriv(pw_c,nu3(3),w)))
              end do
           end do
        end do
        if (maxdiff > 1.0e-8_FP_REAL) then
            ok=.false.; write(useUnit,6300) 'parder_nd', maxdiff; return
        end if

        ! pardeu_nd at scattered interior points vs the same analytic derivative
        do i=1,NP
           ptsc(1,i) = (real(i,FP_REAL)-half)/real(NP,FP_REAL)
           ptsc(2,i) = one - 0.7_FP_REAL*ptsc(1,i)
           ptsc(3,i) = 0.2_FP_REAL + 0.6_FP_REAL*ptsc(1,i)
        end do
        call pardeu_nd(3_FP_DIM,t3,n3,c3,k3,nu3,ptsc,NP,zsc,wrknd,size(wrknd,kind=FP_SIZE),ier)
        if (.not.FITPACK_SUCCESS(ier)) then
            ok=.false.; write(useUnit,'(a,i0)') '[gate J] pardeu_nd(dims=3) failed, ier=',ier; return
        end if
        maxdiff = zero
        do i=1,NP
           maxdiff = max(maxdiff, abs(zsc(i) - &
               poly_deriv(px_c,nu3(1),ptsc(1,i))*poly_deriv(py_c,nu3(2),ptsc(2,i))*poly_deriv(pw_c,nu3(3),ptsc(3,i))))
        end do
        if (maxdiff > 1.0e-8_FP_REAL) then
            ok=.false.; write(useUnit,6300) 'pardeu_nd', maxdiff; return
        end if

        write(useUnit,'(a)') '[gate J] pardtc/parder/pardeu_nd == 2-D (dims=2); dims=3 == analytic derivative'
        6000 format('[gate J] ',a,' failed for ',a,': ',a)
        6100 format('[gate J] ',a,'_nd /= 2-D for ',a,' (nu=',i0,',',i0,')')
        6200 format('[gate J] dims=3 ',a,' = ',1pe12.3)
        6300 format('[gate J] ',a,'(dims=3) /= analytic derivative  (max |diff| = ',1pe12.3,')')
    end function gate_parder_nd

    !> @brief Gate K — dblint_nd (box/domain integral).
    !!
    !! dims=2: bit-for-bit vs dblint over the daregr fit battery (full and sub-box). dims=3: an in-space
    !! separable polynomial fit integrated over a box must equal the closed-form product of 1-D integrals.
    logical function gate_dblint_nd(useUnit) result(ok)
        integer, intent(in) :: useUnit
        ! dims=2
        type(grid_case)  :: cs(3),gc
        integer(FP_SIZE) :: nx,ny,nc,kx,ky,icase,i,mx,my,ib
        integer(FP_FLAG) :: ier
        real(FP_REAL)    :: tx(NXEST),ty(NYEST),c(NCMAX),fp
        real(FP_REAL)    :: iwrk2(NXEST+NYEST),r2d,rnd,xb,xe,yb,ye
        real(FP_REAL)    :: t2(NXEST,2)
        integer(FP_SIZE) :: n2(2),k2(2)
        real(FP_REAL)    :: box(4,2)
        ! dims=3
        integer(FP_SIZE), parameter :: k3v=3,nk1=6,mg=8,nkk=nk1+k3v+1,nc3=nk1**3
        real(FP_REAL),    parameter :: px_c(4)=[1.0_FP_REAL,0.5_FP_REAL,-0.3_FP_REAL,0.2_FP_REAL]
        real(FP_REAL),    parameter :: py_c(4)=[0.8_FP_REAL,-0.4_FP_REAL,0.6_FP_REAL,0.1_FP_REAL]
        real(FP_REAL),    parameter :: pw_c(4)=[1.2_FP_REAL,0.3_FP_REAL,-0.5_FP_REAL,0.15_FP_REAL]
        real(FP_REAL)    :: t3(nkk,3),xg3(mg,3),c3(nc3),z3(mg*mg*mg),fpf,lo(3),hi(3)
        real(FP_REAL)    :: xbb(3),xee(3),rint,rref,wrk3(4000)
        integer(FP_SIZE) :: iwrk3(400),n3(3),k3(3),m3(3),nest3(3),ix,iy,iw,iout
        real(FP_REAL)    :: x,y,w

        ok = .true.
        cs = cases()
        mx = size(daregr_x,kind=FP_SIZE); my = size(daregr_y,kind=FP_SIZE)

        ! ---- dims=2 : bit-for-bit vs dblint ----
        do icase=1,size(cs)
            gc = cs(icase); kx = gc%kx; ky = gc%ky
            call fit_regrid(gc,nx,tx,ny,ty,c,fp,ier)
            if (.not.FITPACK_SUCCESS(ier)) then
                ok=.false.; write(useUnit,7000) 'fit',trim(gc%label),FITPACK_MESSAGE(ier); return
            end if
            nc = (nx-kx-1)*(ny-ky-1)
            t2 = zero; t2(1:nx,1)=tx(1:nx); t2(1:ny,2)=ty(1:ny)
            n2 = [nx,ny]; k2 = [kx,ky]
            ! box 1 = full domain, box 2 = interior sub-box
            box(:,1) = [daregr_x(1),daregr_x(mx),daregr_y(1),daregr_y(my)]
            box(:,2) = [daregr_x(1)+0.2_FP_REAL*(daregr_x(mx)-daregr_x(1)), &
                        daregr_x(1)+0.8_FP_REAL*(daregr_x(mx)-daregr_x(1)), &
                        daregr_y(1)+0.3_FP_REAL*(daregr_y(my)-daregr_y(1)), &
                        daregr_y(1)+0.9_FP_REAL*(daregr_y(my)-daregr_y(1))]
            do ib=1,2
               xb=box(1,ib); xe=box(2,ib); yb=box(3,ib); ye=box(4,ib)
               r2d = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,iwrk2)
               rnd = dblint_nd(2_FP_DIM,t2,n2,c(1:nc),k2,[xb,yb],[xe,ye])
               if (r2d/=rnd) then
                   ok=.false.; write(useUnit,7100) trim(gc%label),ib,abs(r2d-rnd); return
               end if
            end do
        end do

        ! ---- dims=3 : in-space polynomial, closed-form integral ----
        t3 = zero; xg3 = zero; k3 = k3v; m3 = mg; lo = zero; hi = one
        call clamped_unit_knots(k3v,nk1,t3(:,1),n3(1))
        call clamped_unit_knots(k3v,nk1,t3(:,2),n3(2))
        call clamped_unit_knots(k3v,nk1,t3(:,3),n3(3))
        nest3 = n3
        do i=1,mg; xg3(i,1) = real(i-1,FP_REAL)/real(mg-1,FP_REAL); end do
        xg3(:,2) = xg3(:,1); xg3(:,3) = xg3(:,1)
        do ix=1,mg
           x = xg3(ix,1)
           do iy=1,mg
              y = xg3(iy,2)
              do iw=1,mg
                 w = xg3(iw,3)
                 iout = ((ix-1)*mg+(iy-1))*mg + iw
                 z3(iout) = poly_deriv(px_c,0_FP_SIZE,x)*poly_deriv(py_c,0_FP_SIZE,y)*poly_deriv(pw_c,0_FP_SIZE,w)
              end do
           end do
        end do
        wrk3 = zero; iwrk3 = 0
        call regrid(-1_FP_FLAG,3_FP_DIM,m3,xg3,z3,lo,hi,k3,zero,nest3,n3,t3,c3,fpf, &
                    wrk3,size(wrk3,kind=FP_SIZE),iwrk3,size(iwrk3,kind=FP_SIZE),ier)
        if (.not.FITPACK_SUCCESS(ier) .or. fpf>1.0e-9_FP_REAL) then
            ok=.false.; write(useUnit,7200) 'fit residual', fpf; return
        end if
        xbb = [0.2_FP_REAL,0.1_FP_REAL,0.3_FP_REAL]
        xee = [0.9_FP_REAL,0.8_FP_REAL,0.95_FP_REAL]
        rint = dblint_nd(3_FP_DIM,t3,n3,c3,k3,xbb,xee)
        rref = poly_integral(px_c,xbb(1),xee(1))*poly_integral(py_c,xbb(2),xee(2))*poly_integral(pw_c,xbb(3),xee(3))
        if (abs(rint-rref) > 1.0e-9_FP_REAL*max(one,abs(rref))) then
            ok=.false.; write(useUnit,7300) abs(rint-rref); return
        end if

        write(useUnit,'(a)') '[gate K] dblint_nd == dblint (dims=2); dims=3 == closed-form box integral'
        7000 format('[gate K] ',a,' failed for ',a,': ',a)
        7100 format('[gate K] dblint_nd /= dblint for ',a,' box ',i0,'  (|diff| = ',1pe12.3,')')
        7200 format('[gate K] dims=3 ',a,' = ',1pe12.3)
        7300 format('[gate K] dblint_nd(dims=3) /= closed-form  (|diff| = ',1pe12.3,')')
    end function gate_dblint_nd

    !> @brief Gate L — profil_nd (cross-section: fix one axis).
    !!
    !! dims=2: bit-for-bit vs profil for both iopt (fix x, fix y). dims=3: a separable spline sliced at a
    !! fixed value on axis 3 (and, separately, the middle axis 2) yields a dims=2 spline whose fpndsp
    !! evaluation must equal the separable product with that factor frozen (independent splev oracle).
    logical function gate_profil_nd(useUnit) result(ok)
        integer, intent(in) :: useUnit
        ! dims=2
        type(grid_case)  :: cs(3),gc
        integer(FP_SIZE) :: nx,ny,nc,kx,ky,icase,nkx1,nky1
        integer(FP_FLAG) :: ier,ier2
        real(FP_REAL)    :: tx(NXEST),ty(NYEST),c(NCMAX),fp,u
        real(FP_REAL)    :: cu2(NXEST),cund(NXEST)
        real(FP_REAL)    :: t2(NXEST,2)
        integer(FP_SIZE) :: n2(2),k2(2)
        ! dims=3
        integer(FP_SIZE), parameter :: kx3=3,ky3=2,kw3=4,nk1x=6,nk1y=5,nk1w=7
        integer(FP_SIZE), parameter :: nwk=nk1w+kw3+1,maxn=nwk
        integer(FP_SIZE), parameter :: mge=6
        real(FP_REAL)    :: t3(maxn,3),c3(nk1x*nk1y*nk1w),ca(nk1x),cb(nk1y),cd(nk1w)
        real(FP_REAL)    :: sa(maxn),sb(maxn),sd(maxn)
        integer(FP_SIZE) :: n3(3),k3(3)
        real(FP_REAL)    :: maxdiff

        ok = .true.
        cs = cases()

        ! ---- dims=2 : bit-for-bit vs profil (fix x -> iopt/ax to give f(y); fix y -> g(x)) ----
        do icase=1,size(cs)
            gc = cs(icase); kx = gc%kx; ky = gc%ky
            call fit_regrid(gc,nx,tx,ny,ty,c,fp,ier)
            if (.not.FITPACK_SUCCESS(ier)) then
                ok=.false.; write(useUnit,8000) 'fit',trim(gc%label),FITPACK_MESSAGE(ier); return
            end if
            nc = (nx-kx-1)*(ny-ky-1); nkx1 = nx-kx-1; nky1 = ny-ky-1
            t2 = zero; t2(1:nx,1)=tx(1:nx); t2(1:ny,2)=ty(1:ny)
            n2 = [nx,ny]; k2 = [kx,ky]

            ! fix x = u (profil iopt=0 -> f(y); profil_nd ax=1). u interior on axis x.
            u = half*(tx(kx+1)+tx(nx-kx))
            call profil(0_FP_SIZE,tx,nx,ty,ny,c,kx,ky,u,ny,cu2,ier)
            call profil_nd(1_FP_DIM,2_FP_DIM,t2,n2,c(1:nc),k2,u,cund(1:nky1),ier2)
            if (.not.FITPACK_SUCCESS(ier) .or. .not.FITPACK_SUCCESS(ier2) .or. &
                .not.all(cu2(1:nky1)==cund(1:nky1))) then
                ok=.false.; write(useUnit,8100) 'fix-x',trim(gc%label); return
            end if

            ! fix y = u (profil iopt=1 -> g(x); profil_nd ax=2). u interior on axis y.
            u = half*(ty(ky+1)+ty(ny-ky))
            call profil(1_FP_SIZE,tx,nx,ty,ny,c,kx,ky,u,nx,cu2,ier)
            call profil_nd(2_FP_DIM,2_FP_DIM,t2,n2,c(1:nc),k2,u,cund(1:nkx1),ier2)
            if (.not.FITPACK_SUCCESS(ier) .or. .not.FITPACK_SUCCESS(ier2) .or. &
                .not.all(cu2(1:nkx1)==cund(1:nkx1))) then
                ok=.false.; write(useUnit,8100) 'fix-y',trim(gc%label); return
            end if
        end do

        ! ---- dims=3 : separable spline, slice one axis, vs independent splev oracle ----
        call build_separable_3d(kx3,ky3,kw3,nk1x,nk1y,nk1w,t3,n3,k3,c3,ca,cb,cd)
        sa=zero; sa(1:nk1x)=ca; sb=zero; sb(1:nk1y)=cb; sd=zero; sd(1:nk1w)=cd

        ! slice fixing axis 3 (fast axis) at u -> dims=2 spline over (x,y); check vs sa(x)*sb(y)*sd(u)
        if (ok) ok = slice_ok(useUnit,'ax=3',3_FP_DIM,0.42_FP_REAL)
        ! slice fixing the MIDDLE axis 2 at u -> dims=2 spline over (x,w); check vs sa(x)*sb(u)*sd(w)
        if (ok) ok = slice_ok(useUnit,'ax=2',2_FP_DIM,0.37_FP_REAL)
        if (.not.ok) return

        write(useUnit,'(a)') '[gate L] profil_nd == profil (dims=2); dims=3 slice == separable oracle'
        8000 format('[gate L] ',a,' failed for ',a,': ',a)
        8100 format('[gate L] profil_nd /= profil (',a,') for ',a)

    contains

        !> Slice the separable 3-D spline at axis `ax`=u, evaluate the resulting 2-D spline on a grid,
        !! and compare to the separable product with the sliced factor frozen at u.
        logical function slice_ok(useUnit,label,ax,u) result(good)
            integer,          intent(in) :: useUnit
            character(*),     intent(in) :: label
            integer(FP_DIM),  intent(in) :: ax
            real(FP_REAL),    intent(in) :: u
            integer(FP_DIM)  :: da,db,d
            integer(FP_SIZE) :: nka,nkb,i,j,iout
            integer(FP_SIZE) :: n2s(2),k2s(2),m2s(2)
            real(FP_REAL)    :: cu(nk1x*nk1y*nk1w),t2s(maxn,2),xg2s(mge,2),zev(mge*mge)
            real(FP_REAL)    :: ewrk(4000),fac_u(1),fa(mge),fb(mge),refv
            integer(FP_SIZE) :: eiwrk(400)
            integer(FP_FLAG) :: e

            good = .true.
            ! surviving axes (in order), their knot columns and degrees
            da = 0; db = 0
            do d=1,3
               if (d==ax) cycle
               if (da==0) then; da = d; else; db = d; end if
            end do
            nka = n3(da)-k3(da)-1; nkb = n3(db)-k3(db)-1

            call profil_nd(ax,3_FP_DIM,t3,n3,c3,k3,u,cu(1:nka*nkb),e)
            if (.not.FITPACK_SUCCESS(e)) then
                good=.false.; write(useUnit,'(a,a)') '[gate L] profil_nd(dims=3) failed for ',label; return
            end if

            ! assemble the dims=2 cross-section spline and evaluate on an interior grid
            t2s = zero; t2s(1:n3(da),1)=t3(1:n3(da),da); t2s(1:n3(db),2)=t3(1:n3(db),db)
            n2s = [n3(da),n3(db)]; k2s = [k3(da),k3(db)]; m2s = [mge,mge]
            do i=1,mge; xg2s(i,1) = (real(i,FP_REAL)-half)/real(mge,FP_REAL); end do
            xg2s(:,2) = xg2s(:,1)
            call bispev(t2s(1:n2s(1),1),n2s(1),t2s(1:n2s(2),2),n2s(2),cu,k2s(1),k2s(2), &
                        xg2s(:,1),mge,xg2s(:,2),mge,zev,ewrk,size(ewrk,kind=FP_SIZE), &
                        eiwrk,size(eiwrk,kind=FP_SIZE),e)
            if (.not.FITPACK_SUCCESS(e)) then
                good=.false.; write(useUnit,'(a,a)') '[gate L] bispev(cross-section) failed for ',label; return
            end if

            ! oracle: the two surviving 1-D factors on the grid, times the frozen factor at u
            call one_d_factor(da,xg2s(:,1),mge,fa)
            call one_d_factor(db,xg2s(:,2),mge,fb)
            call one_d_factor(ax,[u],1_FP_SIZE,fac_u)

            maxdiff = zero
            do i=1,mge
               do j=1,mge
                  iout = (i-1)*mge + j
                  refv = fa(i)*fb(j)*fac_u(1)
                  maxdiff = max(maxdiff, abs(zev(iout)-refv))
               end do
            end do
            if (maxdiff > 1.0e-12_FP_REAL) then
                good=.false.; write(useUnit,8200) label,maxdiff
            end if
            8200 format('[gate L] dims=3 slice ',a,' /= separable oracle  (max |diff| = ',1pe12.3,')')
        end function slice_ok

        !> Evaluate the 1-D spline factor for axis `d` (a/b/d coefficients) at points xp(1:np).
        subroutine one_d_factor(d,xp,np,fv)
            integer(FP_DIM),  intent(in)  :: d
            integer(FP_SIZE), intent(in)  :: np
            real(FP_REAL),    intent(in)  :: xp(:)
            real(FP_REAL),    intent(out) :: fv(:)
            integer(FP_FLAG) :: e
            select case (d)
               case (1); call splev(t3(1:n3(1),1),n3(1),sa(1:n3(1)),kx3,xp,fv,np,OUTSIDE_NEAREST_BND,e)
               case (2); call splev(t3(1:n3(2),2),n3(2),sb(1:n3(2)),ky3,xp,fv,np,OUTSIDE_NEAREST_BND,e)
               case (3); call splev(t3(1:n3(3),3),n3(3),sd(1:n3(3)),kw3,xp,fv,np,OUTSIDE_NEAREST_BND,e)
            end select
        end subroutine one_d_factor

    end function gate_profil_nd

    !> @brief Gate M — the fitpack_gridded_spline peripheral methods are pure marshalling.
    !!
    !! A dims=3 class fit is exercised through every new method and compared to the direct core call:
    !! eval (scattered) vs bispeu_nd, dfdx_ongrid vs parder_nd, dfdx (scattered) vs pardeu_nd, and
    !! integral vs dblint_nd are all bit-for-bit; cross_section / derivative_spline reproduce
    !! profil_nd / pardtc_nd bit-for-bit and are functionally consistent (a sliced spline evaluated on a
    !! grid equals the full fit with that axis frozen; the derivative spline equals dfdx_ongrid).
    logical function gate_gridded_spline_peripherals(useUnit) result(ok)
        integer, intent(in) :: useUnit
        integer(FP_DIM),  parameter :: dims = 3
        integer(FP_SIZE), parameter :: m(3) = [7,6,8], korder = 3, NP = 5

        type(fitpack_gridded_spline) :: obj,sub,dsp
        integer(FP_SIZE) :: maxm,nz,nc,nc_sub,nc_new,i,lwrk,kwrk
        integer(FP_SIZE) :: k3(3),n3(3),nu(3)
        real(FP_REAL), allocatable :: xg(:,:),zflat(:),wrk(:),wrk2(:),cu_dir(:),newc_dir(:)
        integer(FP_SIZE), allocatable :: iwrk(:)
        real(FP_REAL) :: xp(dims,NP),fe_cls(NP),fe_dir(NP),fd_cls(NP),fd_dir(NP)
        real(FP_REAL) :: fg_cls(product(m)),fg_dir(product(m)),fg_dsp(product(m))
        real(FP_REAL) :: lo(dims),hi(dims),int_cls,int_dir,u,s
        real(FP_REAL) :: xgf(maxval(m),dims),zf_full(m(1)*m(2)),zf_sub(m(1)*m(2))
        integer(FP_FLAG) :: ier,ierc

        ok = .true.
        maxm = maxval(m); nz = product(m); k3 = korder; s = 0.05_FP_REAL
        allocate(xg(maxm,dims),zflat(nz))
        call build_grid(dims,m,xg,zflat)

        ierc = obj%new_fit(xg,zflat,m=m,order=korder,smoothing=s)
        if (.not.FITPACK_SUCCESS(ierc)) then
            write(useUnit,9000) 'class fit',ierc; ok=.false.; return
        end if
        n3 = obj%knots(1:dims); nc = product(n3-k3-1); nu = [1_FP_SIZE,0_FP_SIZE,1_FP_SIZE]

        ! scattered points strictly inside the fitted domain
        do i=1,NP
           xp(1,i) = obj%left(1) + (obj%right(1)-obj%left(1))*(real(i,FP_REAL)-half)/real(NP,FP_REAL)
           xp(2,i) = obj%left(2) + (obj%right(2)-obj%left(2))*(real(NP-i,FP_REAL)+half)/real(NP,FP_REAL)
           xp(3,i) = obj%left(3) + (obj%right(3)-obj%left(3))*half
        end do

        ! --- eval (scattered) vs direct bispeu_nd (+ single-point overload) ---
        fe_cls = obj%eval(xp,ier)
        call bispeu_nd(dims,obj%t,n3,obj%c,k3,xp,NP,fe_dir,ier)
        if (.not.all(fe_cls==fe_dir) .or. obj%eval(xp(:,1))/=fe_dir(1)) then
            write(useUnit,9100) 'eval'; ok=.false.; return
        end if

        ! --- dfdx_ongrid vs direct parder_nd ---
        lwrk = nc + maxm*(maxval(k3-nu)+1)*dims; kwrk = maxm*dims
        allocate(wrk(lwrk),source=zero); allocate(iwrk(kwrk),source=0_FP_SIZE)
        fg_cls = obj%dfdx_ongrid(xg,m,nu,ier)
        call parder_nd(dims,obj%t,n3,obj%c,k3,nu,xg,m,fg_dir,wrk,lwrk,iwrk,kwrk,ier)
        if (.not.all(fg_cls==fg_dir)) then
            write(useUnit,9100) 'dfdx_ongrid'; ok=.false.; return
        end if

        ! --- dfdx (scattered) vs direct pardeu_nd (+ single-point overload) ---
        allocate(wrk2(nc),source=zero)
        fd_cls = obj%dfdx(xp,nu,ier)
        call pardeu_nd(dims,obj%t,n3,obj%c,k3,nu,xp,NP,fd_dir,wrk2,nc,ier)
        if (.not.all(fd_cls==fd_dir) .or. obj%dfdx(xp(:,1),nu)/=fd_dir(1)) then
            write(useUnit,9100) 'dfdx'; ok=.false.; return
        end if

        ! --- integral vs direct dblint_nd ---
        lo = obj%left(1:dims)  + 0.2_FP_REAL*(obj%right(1:dims)-obj%left(1:dims))
        hi = obj%left(1:dims)  + 0.8_FP_REAL*(obj%right(1:dims)-obj%left(1:dims))
        int_cls = obj%integral(lo,hi)
        int_dir = dblint_nd(dims,obj%t,n3,obj%c,k3,lo,hi)
        if (int_cls/=int_dir) then
            write(useUnit,9100) 'integral'; ok=.false.; return
        end if

        ! --- cross_section vs direct profil_nd (bit-for-bit) + functional slice eval ---
        u = half*(obj%t(korder+1,3)+obj%t(n3(3)-korder,3))
        sub = obj%cross_section(3_FP_DIM,u,ier)
        nc_sub = (n3(1)-k3(1)-1)*(n3(2)-k3(2)-1)
        allocate(cu_dir(nc_sub))
        call profil_nd(3_FP_DIM,dims,obj%t,n3,obj%c,k3,u,cu_dir,ier)
        if (sub%dims/=2 .or. any(sub%knots(1:2)/=n3(1:2)) .or. any(sub%order(1:2)/=k3(1:2)) .or. &
            .not.all(sub%c==cu_dir)) then
            write(useUnit,9100) 'cross_section marshalling'; ok=.false.; return
        end if
        ! functional: sub on the (x,y) grid == full fit with axis 3 frozen at u
        zf_sub = sub%eval_ongrid(xg(:,1:2),m(1:2),ier)
        xgf = zero; xgf(1:m(1),1)=xg(1:m(1),1); xgf(1:m(2),2)=xg(1:m(2),2); xgf(1,3)=u
        zf_full = obj%eval_ongrid(xgf,[m(1),m(2),1_FP_SIZE],ier)
        if (maxval(abs(zf_sub-zf_full)) > 1.0e-12_FP_REAL) then
            write(useUnit,9200) 'cross_section', maxval(abs(zf_sub-zf_full)); ok=.false.; return
        end if

        ! --- derivative_spline vs direct pardtc_nd (bit-for-bit) + functional == dfdx_ongrid ---
        dsp = obj%derivative_spline(nu,ier)
        allocate(newc_dir(nc))
        call pardtc_nd(dims,obj%t,n3,obj%c,k3,nu,newc_dir,ier)
        nc_new = product((n3-2*nu)-(k3-nu)-1)
        if (any(dsp%order(1:dims)/=k3-nu) .or. any(dsp%knots(1:dims)/=n3-2*nu) .or. &
            .not.all(dsp%c==newc_dir(1:nc_new))) then
            write(useUnit,9100) 'derivative_spline marshalling'; ok=.false.; return
        end if
        fg_dsp = dsp%eval_ongrid(xg,m,ier)
        if (.not.all(fg_dsp==fg_dir)) then
            write(useUnit,9200) 'derivative_spline', maxval(abs(fg_dsp-fg_dir)); ok=.false.; return
        end if

        write(useUnit,'(a)') '[gate M] class peripherals == direct core calls (eval/dfdx/integral/cross_section/derivative_spline)'
        9000 format('[gate M] ',a,' failed, ier=',i0)
        9100 format('[gate M] class ',a,' /= direct core call')
        9200 format('[gate M] class ',a,' functional check failed  (max |diff| = ',1pe12.3,')')
    end function gate_gridded_spline_peripherals

end module fitpack_grid_nd_tests
