! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   (C) Federico Perini, 12/09/2023
!   Based on the netlib library by Paul Dierckx
!
! **************************************************************************************************
!> @brief Generic N-dimensional gridded smoothing-spline fitter.
!!
!! `fitpack_gridded_spline` fits a tensor-product B-spline \f$ s(x_1,\ldots,x_d) \f$ to data on a
!! rectangular grid of arbitrary domain dimension `dims` (the 3D+ analogue of the concrete 2-D
!! `fitpack_grid_surface`). One runtime-`dims` type holds a fit of any dimensionality up to
!! `MAX_IDIM`; the guts marshal onto the dimension-generic core routines `regrid` (fit) and
!! `ndspev` (grid evaluation).
!!
!! CANONICAL DATA LAYOUT. Internally the function values live in a flat, row-major buffer with
!! domain axis 1 the slowest-varying and axis `dims` the fastest — exactly the contract of `regrid`.
!! The rank-natural entry points accept an assumed-rank array and an optional `row_major` flag:
!!   - `row_major = .true.` (default): the array's Fortran (column-major) storage already equals the
!!     internal buffer, so it is copied as-is. For a rank-`d` array that means passing the axes in
!!     reverse, `z(m_d, ..., m_1)` — the direct generalization of the 2-D `z(iy,ix)` convention.
!!   - `row_major = .false.`: the array is a natural Fortran grid `z(m_1, ..., m_d)` (dimension `k`
!!     is domain axis `k`); it is permuted into the row-major buffer on ingest (a one-time copy).
!! A flat rank-1 `z(:)` plus explicit `m(:)` is always accepted (row-major, unambiguous).
!!
!! @see fitpack_grid_surface — concrete 2-D face; regrid, ndspev — N-D core routines
module fitpack_gridded_splines
    use fitpack_core, only: FP_REAL,FP_SIZE,FP_FLAG,FP_DIM,FP_BOOL,FP_COMM,MAX_IDIM,zero,one, &
                            IOPT_NEW_SMOOTHING,IOPT_OLD_FIT,IOPT_NEW_LEASTSQUARES, &
                            regrid,ndspev,FITPACK_SUCCESS,FITPACK_MESSAGE,fitpack_error_handling, &
                            get_smoothing,FITPACK_INPUT_ERROR, &
                            FP_COMM_SIZE,FP_COMM_PACK,FP_COMM_EXPAND
    use fitpack_fitters
    implicit none
    private

    public :: fitpack_gridded_spline

    !> @brief Tensor-product gridded smoothing-spline fitter for any domain dimension.
    type, extends(fitpack_fitter) :: fitpack_gridded_spline

        !> Domain (tensor) dimension, 1 <= dims <= MAX_IDIM. Metadata below is active over 1:dims.
        integer(FP_DIM)  :: dims = 0

        !> Per-axis spline degree
        integer(FP_SIZE) :: order(MAX_IDIM) = 3

        !> Per-axis point counts, estimated/actual knot counts
        integer(FP_SIZE) :: m(MAX_IDIM)     = 0
        integer(FP_SIZE) :: nest(MAX_IDIM)  = 0
        integer(FP_SIZE) :: knots(MAX_IDIM) = 0

        !> Per-axis domain boundaries
        real(FP_REAL) :: left(MAX_IDIM)  = zero
        real(FP_REAL) :: right(MAX_IDIM) = zero

        !> Bulk arrays (only these are allocatable). t(:,d)=knots on axis d; xg(:,d)=grid coords on
        !> axis d; z(:)=flat function values (row-major, axis 1 slowest). c/wrk/iwrk/fp are inherited.
        real(FP_REAL), allocatable :: t(:,:)
        real(FP_REAL), allocatable :: xg(:,:)
        real(FP_REAL), allocatable :: z(:)

        contains

            procedure :: destroy       => grid_destroy
            procedure :: new_points    => grid_new_points
            procedure :: new_fit       => grid_new_fit
            procedure :: fit           => grid_fit_automatic_knots
            procedure :: least_squares => grid_fit_least_squares
            procedure :: interpolate   => grid_fit_interpolating
            procedure :: eval_ongrid   => grid_eval_ongrid

            !> Parallel communication
            procedure :: comm_size     => grid_comm_size
            procedure :: comm_pack     => grid_comm_pack
            procedure :: comm_expand   => grid_comm_expand

    end type fitpack_gridded_spline

    contains

    !> @brief Release all allocated storage and reset the fit state.
    elemental subroutine grid_destroy(this)
        class(fitpack_gridded_spline), intent(inout) :: this
        integer :: ierr
        call this%destroy_base()
        this%dims  = 0
        this%order = 3
        this%m     = 0
        this%nest  = 0
        this%knots = 0
        this%left  = zero
        this%right = zero
        deallocate(this%t,  stat=ierr)
        deallocate(this%xg, stat=ierr)
        deallocate(this%z,  stat=ierr)
    end subroutine grid_destroy

    !> @brief Load new gridded data. `z` is assumed-rank: a rank-`d` array (2 <= d <= MAX_IDIM) whose
    !!        shape gives the per-axis point counts, or a flat rank-1 array with explicit `m`.
    !!
    !! @param[in]  xg         Grid coordinates, `xg(1:m(d),d)` strictly increasing on axis `d`.
    !! @param[in]  z          Function values (see the module header for the layout conventions).
    !! @param[in]  m          Per-axis point counts; REQUIRED when `z` is rank-1, ignored otherwise.
    !! @param[in]  row_major  Layout of a rank-`d` `z` (default .true.); see the module header.
    !! @param[in]  order      Optional per-axis spline degree (scalar, applied to all axes).
    subroutine grid_new_points(this,xg,z,m,row_major,order)
        class(fitpack_gridded_spline), intent(inout) :: this
        real(FP_REAL),    intent(in)           :: xg(:,:)
        real(FP_REAL),    intent(in)           :: z(..)
        integer(FP_SIZE), intent(in), optional :: m(:)
        logical(FP_BOOL), intent(in), optional :: row_major
        integer(FP_SIZE), intent(in), optional :: order

        integer(FP_DIM)  :: d
        integer(FP_SIZE) :: mm(MAX_IDIM),nz,maxm
        logical(FP_BOOL) :: rowmaj
        real(FP_REAL), allocatable :: zc(:)

        call this%destroy()

        rowmaj = .true.; if (present(row_major)) rowmaj = row_major

        !  ---- resolve dims, per-axis sizes m, and the flat (row-major) value buffer ----
        select rank (z)
            rank (1)
                ! flat entry: values already row-major; m is mandatory
                this%dims = size(m,kind=FP_DIM)
                mm(1:this%dims) = m
                allocate(zc(size(z)), source=reshape(z,[size(z)]))
            rank (2)
                call ranked_sizes(shape(z,kind=FP_SIZE),rowmaj,this%dims,mm)
                zc = flat_row_major(reshape(z,[size(z)]),shape(z,kind=FP_SIZE),rowmaj)
            rank (3)
                call ranked_sizes(shape(z,kind=FP_SIZE),rowmaj,this%dims,mm)
                zc = flat_row_major(reshape(z,[size(z)]),shape(z,kind=FP_SIZE),rowmaj)
            rank (4)
                call ranked_sizes(shape(z,kind=FP_SIZE),rowmaj,this%dims,mm)
                zc = flat_row_major(reshape(z,[size(z)]),shape(z,kind=FP_SIZE),rowmaj)
            rank (5)
                call ranked_sizes(shape(z,kind=FP_SIZE),rowmaj,this%dims,mm)
                zc = flat_row_major(reshape(z,[size(z)]),shape(z,kind=FP_SIZE),rowmaj)
            rank (6)
                call ranked_sizes(shape(z,kind=FP_SIZE),rowmaj,this%dims,mm)
                zc = flat_row_major(reshape(z,[size(z)]),shape(z,kind=FP_SIZE),rowmaj)
            rank default
                this%dims = 0
                return
        end select

        this%m(1:this%dims) = mm(1:this%dims)
        nz   = product(this%m(1:this%dims))
        maxm = maxval(this%m(1:this%dims))

        !  spline degree
        this%order = 3
        if (present(order)) this%order(1:this%dims) = order

        !  grid coordinates and domain bounds (per axis)
        allocate(this%xg(maxm,this%dims), source=zero)
        do d=1,this%dims
           this%xg(1:this%m(d),d) = xg(1:this%m(d),d)
           this%left(d)  = xg(1,d)
           this%right(d) = xg(this%m(d),d)
        end do

        !  flat function values (row-major)
        call move_alloc(zc, this%z)

        !  knot storage: overestimate 2*(order+1) => order+m+1 per axis
        do d=1,this%dims
           this%nest(d) = 2*(this%order(d) + this%m(d) + 1)
        end do
        allocate(this%t(maxval(this%nest(1:this%dims)),this%dims), source=zero)

        !  spline coefficients: product over axes of (nest-order-1)
        allocate(this%c(product(this%nest(1:this%dims)-this%order(1:this%dims)-1)), source=zero)

        this%fp   = zero
        this%iopt = IOPT_NEW_SMOOTHING

        call grid_prepare_workspace(this)

    end subroutine grid_new_points

    !> @brief Load new gridded data and perform a fit in one call.
    integer(FP_FLAG) function grid_new_fit(this,xg,z,m,row_major,order,smoothing) result(ierr)
        class(fitpack_gridded_spline), intent(inout) :: this
        real(FP_REAL),    intent(in)           :: xg(:,:)
        real(FP_REAL),    intent(in)           :: z(..)
        integer(FP_SIZE), intent(in), optional :: m(:)
        logical(FP_BOOL), intent(in), optional :: row_major
        integer(FP_SIZE), intent(in), optional :: order
        real(FP_REAL),    intent(in), optional :: smoothing

        call this%new_points(xg,z,m,row_major,order)
        ierr = this%fit(smoothing)
    end function grid_new_fit

    !> @brief Fit a smoothing spline with automatic knot placement (delegates to regrid).
    integer(FP_FLAG) function grid_fit_automatic_knots(this,smoothing,order) result(ierr)
        class(fitpack_gridded_spline), intent(inout) :: this
        real(FP_REAL),    intent(in), optional :: smoothing
        integer(FP_SIZE), intent(in), optional :: order

        integer(FP_SIZE) :: loop,nit
        real(FP_REAL)    :: smooth_now(3)

        call get_smoothing(this%smoothing,smoothing,nit,smooth_now)

        ! optional per-axis order change (scalar applied to all axes)
        if (present(order)) this%order(1:this%dims) = order

        ! the workspace need grows with the order; ensure it covers the current order
        call grid_prepare_workspace(this)

        do loop=1,nit
            this%smoothing = smooth_now(loop)
            call regrid(this%iopt,this%dims,this%m(1:this%dims),this%xg,this%z, &
                        this%left(1:this%dims),this%right(1:this%dims),this%order(1:this%dims), &
                        this%smoothing,this%nest(1:this%dims),this%knots(1:this%dims),this%t, &
                        this%c,this%fp,this%wrk,this%lwrk,this%iwrk,this%liwrk,ierr)
            if (FITPACK_SUCCESS(ierr)) this%iopt = IOPT_OLD_FIT
        end do
    end function grid_fit_automatic_knots

    !> @brief Least-squares fit on the current knots (optionally recomputing them first).
    integer(FP_FLAG) function grid_fit_least_squares(this,smoothing,reset_knots) result(ierr)
        class(fitpack_gridded_spline), intent(inout) :: this
        real(FP_REAL), intent(in), optional :: smoothing
        logical,       intent(in), optional :: reset_knots

        logical :: do_reset
        do_reset = .false.; if (present(reset_knots)) do_reset = reset_knots
        if (do_reset) then
            this%iopt = IOPT_NEW_SMOOTHING
            ierr = this%fit(smoothing)
            if (.not.FITPACK_SUCCESS(ierr)) return
        end if
        this%iopt = IOPT_NEW_LEASTSQUARES
        ierr = this%fit()
    end function grid_fit_least_squares

    !> @brief Interpolating fit (s = 0).
    integer(FP_FLAG) function grid_fit_interpolating(this) result(ierr)
        class(fitpack_gridded_spline), intent(inout) :: this
        this%iopt = IOPT_NEW_SMOOTHING
        ierr = this%fit(smoothing=zero)
    end function grid_fit_interpolating

    !> @brief Evaluate the fitted spline on a rectangular grid (delegates to ndspev).
    !!
    !! @param[in]  xg    Evaluation coordinates, `xg(1:m(d),d)` strictly increasing on axis `d`.
    !! @param[in]  m     Per-axis evaluation-grid sizes.
    !! @param[out] ierr  Optional error flag.
    !! @return     zeval Flat, row-major (axis 1 slowest) evaluated values, length `product(m)`.
    function grid_eval_ongrid(this,xg,m,ierr) result(zeval)
        class(fitpack_gridded_spline), intent(in) :: this
        real(FP_REAL),    intent(in)  :: xg(:,:)
        integer(FP_SIZE), intent(in)  :: m(:)
        integer(FP_FLAG), intent(out), optional :: ierr
        real(FP_REAL) :: zeval(product(m))

        integer(FP_SIZE) :: maxm,maxk1,lwrk,kwrk
        real(FP_REAL),    allocatable :: wrk(:)
        integer(FP_SIZE), allocatable :: iwrk(:)
        integer(FP_FLAG) :: ier

        maxm  = maxval(m)
        maxk1 = maxval(this%order(1:this%dims))+1
        lwrk  = maxm*maxk1*this%dims
        kwrk  = maxm*this%dims
        allocate(wrk(lwrk),source=zero); allocate(iwrk(kwrk),source=0_FP_SIZE)

        call ndspev(this%dims,this%t,this%knots(1:this%dims),this%c,this%order(1:this%dims), &
                    xg,m,zeval,wrk,lwrk,iwrk,kwrk,ier)

        call fitpack_error_handling(ier,ierr,'gridded spline eval_ongrid')
    end function grid_eval_ongrid

    ! =================================================================================================
    ! WORKSPACE + LAYOUT HELPERS
    ! =================================================================================================

    !> @brief (Re)size the fit workspace to regrid's requirement for the current dims/order/nest/grid.
    !!        Mirrors regrid's internal lwest/kwest so the type is a correctly-sized caller; grows the
    !!        arrays only when needed so an iopt=1 continuation keeps its state in wrk(1:2+dims).
    subroutine grid_prepare_workspace(this)
        class(fitpack_gridded_spline), intent(inout) :: this

        integer(FP_DIM)  :: dims,i
        integer(FP_SIZE) :: maxm,maxnest,maxk1,maxk2,mm,mynx,bufmax,lwrk,kwrk
        integer(FP_SIZE) :: nk1(MAX_IDIM)

        dims    = this%dims
        maxm    = maxval(this%m(1:dims))
        maxnest = maxval(this%nest(1:dims))
        maxk1   = maxval(this%order(1:dims))+1
        maxk2   = maxval(this%order(1:dims))+2
        nk1(1:dims) = this%nest(1:dims)-(this%order(1:dims)+1)

        bufmax = 0
        do i=1,dims-1
           bufmax = max(bufmax, product(nk1(1:i))*product(this%m(i+1:dims)))
        end do
        mm   = max(maxnest, maxm, product(this%m(2:dims)), &
                   product(nk1(1:dims-1)), maxval(nk1(1:dims)))
        mynx = 2*bufmax

        kwrk = (1+dims) + maxm*dims + maxnest*dims
        lwrk = (2+dims) + maxnest*dims + maxm*maxk1*dims + 2*maxnest*maxk2*dims + mm + mynx

        if (.not.allocated(this%iwrk)) then
            allocate(this%iwrk(kwrk),source=0)
        elseif (size(this%iwrk)<kwrk) then
            deallocate(this%iwrk); allocate(this%iwrk(kwrk),source=0)
        end if
        this%liwrk = size(this%iwrk,kind=FP_SIZE)

        if (.not.allocated(this%wrk)) then
            allocate(this%wrk(lwrk),source=zero)
        elseif (size(this%wrk)<lwrk) then
            deallocate(this%wrk); allocate(this%wrk(lwrk),source=zero)
        end if
        this%lwrk = size(this%wrk,kind=FP_SIZE)

    end subroutine grid_prepare_workspace

    !> @brief Map a rank-`d` array shape to per-axis domain sizes, honouring the row_major flag.
    !!        row_major=.true.: axes are reversed (z(m_d,...,m_1)); .false.: natural (z(m_1,...,m_d)).
    pure subroutine ranked_sizes(shp,row_major,dims,m)
        integer(FP_SIZE), intent(in)  :: shp(:)
        logical(FP_BOOL), intent(in)  :: row_major
        integer(FP_DIM),  intent(out) :: dims
        integer(FP_SIZE), intent(out) :: m(MAX_IDIM)
        integer(FP_DIM) :: d
        dims = size(shp,kind=FP_DIM)
        do d=1,dims
           if (row_major) then
              m(d) = shp(dims-d+1)     ! reversed: domain axis d = array dim (dims-d+1)
           else
              m(d) = shp(d)            ! natural : domain axis d = array dim d
           end if
        end do
    end subroutine ranked_sizes

    !> @brief Return the flat row-major (axis-1-slowest) value buffer from a column-major flatten.
    !!        row_major=.true.: the input already IS row-major, copy as-is. row_major=.false.: the
    !!        input is a natural Fortran grid z(m_1,...,m_d); permute (reverse axes) into row-major.
    pure function flat_row_major(zc,shp,row_major) result(zr)
        real(FP_REAL),    intent(in) :: zc(:)
        integer(FP_SIZE), intent(in) :: shp(:)
        logical(FP_BOOL), intent(in) :: row_major
        real(FP_REAL) :: zr(size(zc))

        integer(FP_DIM)  :: d,dims
        integer(FP_SIZE) :: p,rowidx,k,idx(MAX_IDIM)

        if (row_major) then
            zr = zc
            return
        end if

        dims = size(shp,kind=FP_DIM)
        idx(1:dims) = 1
        do p=1,size(zc)
           ! row-major destination index (axis 1 slowest) for the current multi-index idx
           rowidx = 0
           do d=1,dims
              rowidx = rowidx*shp(d) + (idx(d)-1)
           end do
           zr(rowidx+1) = zc(p)                 ! zc is column-major: p iterates axis 1 fastest
           ! increment idx with axis 1 fastest (matches zc's storage order)
           k = 1
           do
              idx(k) = idx(k)+1
              if (idx(k)<=shp(k)) exit
              idx(k) = 1; k = k+1
              if (k>dims) exit
           end do
        end do
    end function flat_row_major

    ! =================================================================================================
    ! PARALLEL COMMUNICATION
    ! =================================================================================================

    !> @brief Communication buffer size: base fields + fixed metadata (1 + 6*MAX_IDIM) + bulk arrays.
    elemental integer(FP_SIZE) function grid_comm_size(this)
        class(fitpack_gridded_spline), intent(in) :: this
        grid_comm_size = this%core_comm_size() &
                       + 1 + 6*MAX_IDIM &
                       + FP_COMM_SIZE(this%t)  &
                       + FP_COMM_SIZE(this%xg) &
                       + FP_COMM_SIZE(this%z)
    end function grid_comm_size

    !> @brief Pack the fitter into a communication buffer.
    pure subroutine grid_comm_pack(this,buffer)
        class(fitpack_gridded_spline), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_pack(buffer)
        pos = this%core_comm_size() + 1

        buffer(pos) = real(this%dims, FP_COMM);              pos = pos + 1
        buffer(pos:pos+MAX_IDIM-1) = real(this%order, FP_COMM); pos = pos + MAX_IDIM
        buffer(pos:pos+MAX_IDIM-1) = real(this%m,     FP_COMM); pos = pos + MAX_IDIM
        buffer(pos:pos+MAX_IDIM-1) = real(this%nest,  FP_COMM); pos = pos + MAX_IDIM
        buffer(pos:pos+MAX_IDIM-1) = real(this%knots, FP_COMM); pos = pos + MAX_IDIM
        buffer(pos:pos+MAX_IDIM-1) = this%left;                 pos = pos + MAX_IDIM
        buffer(pos:pos+MAX_IDIM-1) = this%right;                pos = pos + MAX_IDIM

        call FP_COMM_PACK(this%t,  buffer(pos:));   pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_PACK(this%xg, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%xg)
        call FP_COMM_PACK(this%z,  buffer(pos:))
    end subroutine grid_comm_pack

    !> @brief Expand the fitter from a communication buffer.
    pure subroutine grid_comm_expand(this,buffer)
        class(fitpack_gridded_spline), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)
        integer(FP_SIZE) :: pos

        call this%core_comm_expand(buffer)
        pos = this%core_comm_size() + 1

        this%dims  = nint(buffer(pos), FP_DIM);                        pos = pos + 1
        this%order = nint(buffer(pos:pos+MAX_IDIM-1), FP_SIZE);        pos = pos + MAX_IDIM
        this%m     = nint(buffer(pos:pos+MAX_IDIM-1), FP_SIZE);        pos = pos + MAX_IDIM
        this%nest  = nint(buffer(pos:pos+MAX_IDIM-1), FP_SIZE);        pos = pos + MAX_IDIM
        this%knots = nint(buffer(pos:pos+MAX_IDIM-1), FP_SIZE);        pos = pos + MAX_IDIM
        this%left  = buffer(pos:pos+MAX_IDIM-1);                       pos = pos + MAX_IDIM
        this%right = buffer(pos:pos+MAX_IDIM-1);                       pos = pos + MAX_IDIM

        call FP_COMM_EXPAND(this%t,  buffer(pos:));   pos = pos + FP_COMM_SIZE(this%t)
        call FP_COMM_EXPAND(this%xg, buffer(pos:));   pos = pos + FP_COMM_SIZE(this%xg)
        call FP_COMM_EXPAND(this%z,  buffer(pos:))
    end subroutine grid_comm_expand

end module fitpack_gridded_splines
